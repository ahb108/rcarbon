modelTest <- function(bins, ages, errors, yearRange, calCurves, nsim, edge=500, resOffsets=0 ,resErrors=0, runm=50, raw=FALSE, model="exponential", method="standard", normalised=FALSE, simnormalised=normalised, ncores=1){

    if (!model %in% c("exponential","uniform")){
        stop("Currently only exponential and uniform models are supported.")
    }
    print("Calibrating observed ages...")
    # Calibrate and bin observed ages
    observed <- rspd(bins=bins,ages=ages,errors=errors,yearRange=yearRange,calCurves=calCurves, method=method, normalised=normalised, ncores=ncores, verbose=FALSE)
    finalSPD <- observed[,"SPD"]
    # Simulation
    C14Interval <- range(ages)
    sim <- matrix(NA,nrow=length(finalSPD),ncol=nsim)
    print("Monte-Carlo test...")
    flush.console()
    pb <- txtProgressBar(min=1, max=nsim, style=3)
    for (s in 1:nsim){
        setTxtProgressBar(pb, s)
        if (model=="uniform"){
            randomDates <- round(runif(length(unique(bins)),rev(yearRange)[1]-edge,rev(yearRange)[2]+edge))
        } else if (model=="exponential"){
            plusoffset <- min(finalSPD[finalSPD!=0])/10000 
            finalSPD <- finalSPD+plusoffset #avoid log(0)
            fit <- lm(log(finalSPD)~observed[,1])
            time <- seq(min(observed[,1])-edge,max(observed[,1])+edge,1)
            est <-  exp(fit$coefficients[1]) * exp(time*fit$coefficients[2])
            pweights <- est/sum(est)
            randomDates <- round(sample(time,size=length(unique(bins)),prob=pweights))
        }
        randomSDs <- sample(size=length(randomDates), errors, replace=TRUE)
        simDates <- round(uncalibrate(randomDates,randomSDs)[,4:3])
        randomDates <- simDates[,1]
        randomSDs <- simDates[,2] 
        tmp <- calibrate(ages=randomDates,errors=randomSDs, resOffsets=0 ,resErrors=0, timeRange=yearRange, calCurves='intcal13', method=method, normalised=simnormalised, ncores=ncores, verbose=FALSE)
        tmp <- lapply(tmp, `[[`, 2)
        tmp <- lapply(tmp,`[`, 2)
        simDateMatrix <- do.call("cbind",tmp)
        sim[,s] <- apply(simDateMatrix,1,sum)
        sim[,s] <- sim[,s]/sum(sim[,s])
        sim[,s] <- sim[,s]+plusoffset
        tmp1 <- runMean(sim[,s],runm)
        sim[!is.na(tmp1),s] <- tmp1[!is.na(tmp1)]
    }
    close(pb)
    # Envelope, z-scores, global p-value
    lo <- apply(sim,1,quantile,prob=0.025)
    hi <- apply(sim,1,quantile,prob=0.975)
    Zsim <- t(apply(sim,1,scale))
    zLo <- apply(Zsim,1,quantile,prob=0.025,na.rm=TRUE)
    zHi <- apply(Zsim,1,quantile,prob=0.975,na.rm=TRUE)
    Zscore_empirical <- (finalSPD - apply(sim, 1, mean))/apply(sim, 1, sd)
    busts <- which(Zscore_empirical< zLo)
    booms <- which(Zscore_empirical> zHi)
    busts2 <- which(finalSPD< lo)
    booms2 <- which(finalSPD> hi)
    observedStatistic <- sum(c(zLo[busts] - Zscore_empirical[busts]),c(Zscore_empirical[booms]-zHi[booms]))
    expectedstatistic <- abs(apply(Zsim,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=zLo)) + apply(Zsim,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=zHi)
    pvalue <- 1 - c(length(expectedstatistic[expectedstatistic <= observedStatistic]))/c(length(expectedstatistic)+1)
    # Results
    result <- data.frame(calBP=observed[,1],SPD=finalSPD,lo=lo,hi=hi)
    if(raw==FALSE){
        res <- list(result=result, sim=NA, pval=pvalue, fit=fit)
    }
    if(raw==TRUE){
        res <- list(result=result, sim=sim, pval=pvalue, fit=fit)
    }
    class(res) <- "rspdModelTest"
    return(res)
}

regionTest <- function(x, bins, regions,  nsim, runm=NA, timeRange=NA, datenormalised=FALSE, raw=FALSE, verbose=TRUE){

    if (any(!is.na(timeRange))){
        if (nrow(x[[1]][["agegrid"]]) != length(seq(timeRange[1],timeRange[2],-1))){
            for (d in 1:length(x)){
                tmpag <- x[[d]][["agegrid"]]
                tmpag <- tmpag[tmpag$calBP <= timeRange[1] & tmpag$calBP <= timeRange[2], ]
                x[[d]][["agegrid"]] <- tmpag
            }
        }
    }
    ## Calculate SPDs per bin
    binNames <- unique(bins)
    binnedMatrix <- matrix(NA, nrow=nrow(x[[1]][["agegrid"]]), ncol=length(binNames))
    regionList <- numeric()
    if (verbose & length(binNames)>1){
        print("Binning by site/phase...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Binning by site/phase...")
    }
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        slist <- x[index]
        tmp <- lapply(lapply(slist, `[[`, 2),`[`,2)
        if (datenormalised){
            tmp <- lapply(tmp,FUN=function(x) x/sum(x))
        }
        if (length(binNames)>1){
            spd.tmp <- Reduce("+", tmp) / length(index)
        } else {
            spd.tmp <- Reduce("+", tmp)
        }
        binnedMatrix[,b] <- spd.tmp[,1]
        regionList[b] <- regions[index][1]
    }
    close(pb)
    ## Combine observed bins for focal region
    observedSPD <- vector("list",length=length(unique(regionList)))
    names(observedSPD) <- unique(regionList)
    for (d in 1:length(unique(regionList))){
        focus <- unique(regionList)[d]
        index <- which(regionList==focus)
        tmpSPD <- apply(binnedMatrix[,index], 1, sum)
        if (!is.na(runm)){
            tmp1 <- runMean(tmpSPD,runm)
            tmpSPD[!is.na(tmp1)] <- tmp1[!is.na(tmp1)]
        }
        tmpSPD <- tmpSPD / sum(tmpSPD)
        observedSPD[[d]] <- data.frame(calBP=x[[1]][["agegrid"]][,1], SPD=tmpSPD)
    }
    ## Simulate focal dataset but draw bins from all regions
    simulatedSPD <- vector("list",length=length(unique(regionList)))
    for (d in 1:length(unique(regionList))){
        simulatedSPD[[d]] <- matrix(NA, nrow=nrow(x[[1]][["agegrid"]]), ncol=nsim)
    }
    print("Permutation test...")
    flush.console()
    pb <- txtProgressBar(min=1, max=nsim, style=3)
    for (s in 1:nsim){
        setTxtProgressBar(pb, s)
        simRegionList <- sample(regionList)
        for (d in 1:length(unique(simRegionList))){
            focus <- unique(regionList)[d]
            index <- which(simRegionList==focus)
            tmpSPD <- apply(binnedMatrix[,index],1,sum)
            if (!is.na(runm)){
                tmp1 <- runMean(tmpSPD,runm)
                tmpSPD[!is.na(tmp1)] <- tmp1[!is.na(tmp1)]
            }           
            tmpSPD <- tmpSPD/sum(tmpSPD)
            simulatedSPD[[d]][,s] <- tmpSPD
        }
    }
    names(simulatedSPD) <- unique(regionList)
    close(pb)
    simulatedCIlist <- vector("list",length=length(unique(regionList)))
    for (d in 1:length(unique(regionList))){
        simulatedCIlist[[d]] <- cbind(apply(simulatedSPD[[d]],1,quantile,prob=c(0.025)), apply(simulatedSPD[[d]],1,quantile,prob=c(0.975)))
        names(simulatedCIlist) <- unique(regionList)
    }
    pValueList <- numeric(length=length(simulatedSPD))
    ## Further statistics
    for (a in 1:length(simulatedSPD)){
        zscoreMean <- apply(simulatedSPD[[a]],1,mean)
        zscoreSD <- apply(simulatedSPD[[a]],1,sd)
        tmp.sim <- t(apply(simulatedSPD[[a]],1,function(x){ return((x - mean(x))/sd(x)) }))
        tmp.sim[is.na(tmp.sim)] <- 0
        tmp.obs <- observedSPD[[a]]
        tmp.obs[,2] <- (tmp.obs[,2] - zscoreMean) / zscoreSD
        tmp.obs[is.na(tmp.obs[,2]),2] <- 0
        tmp.ci <- t(apply(tmp.sim,1, quantile, prob=c(0.025,0.975)))
        expectedstatistic <- abs(apply(tmp.sim,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=tmp.ci[,1])) + apply(tmp.sim,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=tmp.ci[,2])   
        lower <- tmp.obs[,2] - tmp.ci[,1]
        indexLow <- which(tmp.obs[,2] < tmp.ci[,1])
        higher <- tmp.obs[,2] - tmp.ci[,2]
        indexHi <- which(tmp.obs[,2] > tmp.ci[,2])
        observedStatistic <- sum(abs(lower[indexLow]))+sum(higher[indexHi])
        pValueList[[a]] <- 1
        if (observedStatistic>0){    
            pValueList[[a]] <- 1 - c(length(expectedstatistic[expectedstatistic <= observedStatistic]))/c(length(expectedstatistic)+1)
        }
        names(pValueList) <- unique(regionList)
    }        
    if(raw==FALSE){
        res <- list(observed=observedSPD,envelope=simulatedCIlist,pValueList=pValueList)
    } else {  
        res <- list(observed=observedSPD,envelope=simulatedCIlist,raw=simulatedSPD,pValueList=pValueList)
    }
    class(res) <- "rspdRegionTest"
    return(res)
}
