modelTest <- function(x, errors, bins, nsim, runm=NA, timeRange=NA, edge=500, raw=FALSE, model=c("exponential","custom"), predgrid=NA, method="standard", datenormalised=FALSE, ncores=1, fitonly=FALSE, verbose=TRUE){

    ## Bin observed dates
    if (verbose){ print("Aggregating observed dates...") }
    observed <- rspd(x=x, bins=bins, timeRange=timeRange, datenormalised=datenormalised, runm=runm, spdnormalised=FALSE, verbose=FALSE)
    finalSPD <- observed$grid$SPD
    ## Simulation
    sim <- matrix(NA,nrow=length(finalSPD),ncol=nsim)
    if (verbose & !fitonly){
        print("Monte-Carlo test...")
    flush.console()
        pb <- txtProgressBar(min=1, max=nsim, style=3)
    }
    coeffs <- NA
    time <- seq(min(observed$grid$calBP)-edge,max(observed$grid$calBP)+edge,1)
    if (model=="exponential"){
        plusoffset <- min(finalSPD[finalSPD!=0])/10000 
        finalSPD <- finalSPD+plusoffset #avoid log(0)
        fit <- lm(log(finalSPD)~observed$grid$calBP)
        coeffs <- fit$coefficients
        est <-  exp(coeffs[1]) * exp(time*coeffs[2])
        predgrid <- data.frame(calBP=time, PrDens=est)
    } else if (model=="custom"){
        if (length(predgrid)!=2){
            stop("If you choose a custom model, you must provide a proper predgrid argument (two-column data.frame of calBP and predicted densities).")
        }
        plusoffset <- 0
    }
    if (fitonly){
        print("Done (SPD and fitted model only).")
        res <- list(result=NA, sim=NA, pval=NA, osbSPD=observed, fit=predgrid, coefficients=coeffs)
        return(res)
    }
    cragrid <- pdUncal(predgrid, verbose=FALSE)
    obscras <- x$metadata$CRA
    cragrid$PrDens[cragrid$CRA > max(obscras) | cragrid$CRA < min(obscras)] <- 0
    for (s in 1:nsim){
        if (verbose){ setTxtProgressBar(pb, s) }
        randomDates <- sample(cragrid$CRA, replace=TRUE, size=length(unique(bins)), prob=cragrid$PrDens)
        randomSDs <- sample(size=length(randomDates), errors, replace=TRUE)
        tmp <- calibrate(ages=randomDates,errors=randomSDs, resOffsets=0 ,resErrors=0, timeRange=timeRange, calCurves='intcal13', method=method, normalised=datenormalised, compact=FALSE, ncores=ncores, verbose=FALSE)
        tmp <- lapply(tmp$grid,`[`,2)
        simDateMatrix <- do.call("cbind",tmp)
        sim[,s] <- apply(simDateMatrix,1,sum)
        sim[,s] <- sim[,s] + plusoffset
        if (!is.na(runm)){
            sim[,s] <- runMean(sim[,s], runm, edge="fill")
        }
    }
    if (verbose){ close(pb) }
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
    result <- data.frame(calBP=observed$grid$calBP,SPD=finalSPD,lo=lo,hi=hi)
    if(raw==FALSE){ sim <- NA }
    res <- list(result=result, sim=sim, pval=pvalue, fit=predgrid, coefficients=coeffs)
    class(res) <- "rspdModelTest"
    if (verbose){ print("Done.") }
    return(res)
}

rmarkTest <- function(x, bins, marks,  nsim, runm=NA, timeRange=NA, datenormalised=FALSE, raw=FALSE, verbose=TRUE){

    ## Calculate SPDs per bin
    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))
    regionList <- numeric()
    if (verbose & length(binNames)>1){
        print("Binning by site/phase...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Binning by site/phase...")
    }
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        slist <- x$grids[index]
        slist <- lapply(slist,FUN=function(x) merge(calyears,x, all.x=TRUE)) 
        slist <- rapply(slist, f=function(x) ifelse(is.na(x),0,x), how="replace")
        slist <- lapply(slist, FUN=function(x) x[with(x, order(-calBP)), ])
        tmp <- lapply(slist,`[`,"PrDens")
        if (datenormalised){   
            outofTR <- lapply(tmp,sum)==0 # date out of range
            tmpc <- tmp[!outofTR]
            if (length(tmpc)>0){
                tmp <- lapply(tmpc,FUN=function(x) x/sum(x))
            }
        }
        if (length(binNames)>1){
            spd.tmp <- Reduce("+", tmp) / length(index)
        } else {
            spd.tmp <- Reduce("+", tmp)
        }
        binnedMatrix[,b] <- spd.tmp[,1]
        regionList[b] <- marks[index][1]
    }
    if (verbose & length(binNames)>1){ close(pb) }
    ## Combine observed bins for focal region
    observedSPD <- vector("list",length=length(unique(regionList)))
    names(observedSPD) <- unique(regionList)
    for (d in 1:length(unique(regionList))){
        focus <- unique(regionList)[d]
        index <- which(regionList==focus)
        tmpSPD <- apply(binnedMatrix[,index], 1, sum)
        if (!is.na(runm)){
            tmpSPD <- runMean(tmpSPD, runm, edge="fill")
        }
        tmpSPD <- tmpSPD / sum(tmpSPD)
        observedSPD[[d]] <- data.frame(calBP=calyears, SPD=tmpSPD)
    }
    ## Simulate focal dataset but draw bins from all regions
    simulatedSPD <- vector("list",length=length(unique(regionList)))
    for (d in 1:length(unique(regionList))){
        simulatedSPD[[d]] <- matrix(NA, nrow=nrow(calyears), ncol=nsim)
    }
    if (verbose){
        print("Permutation test...")
        flush.console()
        pb <- txtProgressBar(min=1, max=nsim, style=3)
    }
    for (s in 1:nsim){
        if (verbose){ setTxtProgressBar(pb, s) }
        simRegionList <- sample(regionList)
        for (d in 1:length(unique(simRegionList))){
            focus <- unique(regionList)[d]
            index <- which(simRegionList==focus)
            tmpSPD <- apply(binnedMatrix[,index],1,sum)
            if (!is.na(runm)){
                tmpSPD <- runMean(tmpSPD, runm, edge="fill")
            }           
            tmpSPD <- tmpSPD/sum(tmpSPD)
            simulatedSPD[[d]][,s] <- tmpSPD
        }
    }
    names(simulatedSPD) <- unique(regionList)
    if (verbose){ close(pb) }
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
    class(res) <- "rspdMarkTest"
    if (verbose){ print("Done.") }
    return(res)
}



