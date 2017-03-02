modelTest <- function(x, errors, nsim, bins=NA, runm=NA, timeRange=NA, raw=FALSE, model=c("exponential","explog","custom"), predgrid=NA, method="standard", datenormalised=FALSE, ncores=1, fitonly=FALSE, verbose=TRUE){

    if (verbose){ print("Aggregating observed dates...") }
    if (is.na(bins[1])){
        samplesize <- length(x$grids)
    } else {
        samplesize <- length(unique(bins))
    }
    observed <- spd(x=x, bins=bins, timeRange=timeRange, datenormalised=datenormalised, runm=runm, spdnormalised=TRUE, verbose=FALSE)
    finalSPD <- observed$grid$PrDens
    ## Simulation
    sim <- matrix(NA,nrow=length(finalSPD),ncol=nsim)
    if (verbose & !fitonly){
        print("Monte-Carlo test...")
    flush.console()
        pb <- txtProgressBar(min=1, max=nsim, style=3)
    }
    coeffs <- NA
    time <- seq(max(observed$grid$calBP),min(observed$grid$calBP),-1)
    if (model=="exponential"){
        plusoffset <- 0
        fit <- nls(y ~ exp(a + b * x), data=data.frame(x=time, y=finalSPD), start=list(a=0, b=0))
        est <- predict(fit, list(x=time))
        predgrid <- data.frame(calBP=time, PrDens=est)
    } else if (model=="explog"){
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
    } else {
        stop("Specified model not one of current choices.")
    }
    if (fitonly){
        print("Done (SPD and fitted model only).")
        res <- list(result=NA, sim=NA, pval=NA, osbSPD=observed, fit=predgrid, coefficients=coeffs)
        return(res)
    }
    cragrid <- uncalibrate(as.CalGrid(predgrid), verbose=FALSE)
    obscras <- x$metadata$CRA
    cragrid$PrDens[cragrid$CRA > max(obscras) | cragrid$CRA < min(obscras)] <- 0
    for (s in 1:nsim){
        if (verbose){ setTxtProgressBar(pb, s) }
        randomDates <- sample(cragrid$CRA, replace=TRUE, size=samplesize, prob=cragrid$PrDens)
        randomSDs <- sample(size=length(randomDates), errors, replace=TRUE)
        tmp <- calibrate(ages=randomDates,errors=randomSDs, resOffsets=0 ,resErrors=0, timeRange=timeRange, calCurves='intcal13', method=method, normalised=datenormalised, ncores=ncores, verbose=FALSE, calMatrix=TRUE)
        simDateMatrix <- tmp$calmat
        sim[,s] <- apply(simDateMatrix,1,sum)
        sim[,s] <- sim[,s] + plusoffset
        sim[,s] <- (sim[,s]/sum(sim[,s])) * sum(predgrid$PrDens)
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
    result <- data.frame(calBP=observed$grid$calBP,PrDens=finalSPD,lo=lo,hi=hi)
    if(raw==FALSE){ sim <- NA }
    res <- list(result=result, sim=sim, pval=pvalue, fit=predgrid, coefficients=coeffs)
    class(res) <- "SpdModelTest"
    if (verbose){ print("Done.") }
    return(res)
}

permTest <- function(x, marks,  timeRange, nsim, propfc=NA, bins=NA, runm=NA, datenormalised=FALSE, raw=FALSE, verbose=TRUE){

    if (is.na(propfc)){ prfc <- 1 } else { prfc <- propfc }
    if (is.na(bins[1])){
        binNames <- bins <- "bin1"
    } else {
        binNames <- unique(bins)
    }
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    propdf <- data.frame(calBP=calyears, Denom=NA, ObsProp=NA, EnvHi=NA, EnvLo=NA)
    binnedMatrix <- matrix(nrow=nrow(calyears), ncol=length(binNames))
    GroupList <- numeric()
    if (verbose & length(binNames)>1){
        print("Summing observed groups...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3)
    }
    caldateTR <- as.numeric(x$metadata[1,c("StartBP","EndBP")])
    caldateyears <- seq(caldateTR[1],caldateTR[2],-1)
    check <- caldateTR[1] >= timeRange[1] & caldateTR[2] <= timeRange[2]
    ## Observed SPDs
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        if (length(x$calmatrix)>1){
            if (!check){
                stop("The time range of the calibrated dataset must be at least as large as the spd time range.")
            } else {
                tmp <- x$calmatrix[,index, drop=FALSE]
                if (datenormalised){
                    tmp <- apply(tmp,2,FUN=function(x) x/sum(x))
                }
                spdtmp <- rowSums(tmp)
                if (length(binNames)>1){
                    spdtmp <- spdtmp / length(index)
                }
                binnedMatrix[,b] <- spdtmp[caldateyears<=timeRange[1] & caldateyears>=timeRange[2]]
            }
        } else {
            slist <- x$grids[index]
            slist <- lapply(slist,FUN=function(x) merge(calyears,x, all.x=TRUE)) 
            slist <- rapply(slist, f=function(x) ifelse(is.na(x),0,x), how="replace")
            slist <- lapply(slist, FUN=function(x) x[with(x, order(-calBP)), ])
            tmp <- lapply(slist,`[`,2)
            if (datenormalised){   
                outofTR <- lapply(tmp,sum)==0 # date out of range
                tmpc <- tmp[!outofTR]
                if (length(tmpc)>0){
                    tmp <- lapply(tmpc,FUN=function(x) x/sum(x))
                }
            }
            if (length(binNames)>1){
                spdtmp <- Reduce("+", tmp) / length(index)
            } else {
                spdtmp <- Reduce("+", tmp)
            }
            if (length(binNames)>1){
                spdtmp <- spdtmp / length(index)
            }
            binnedMatrix[,b] <- spdtmp[,1]
        }
        GroupList[b] <- marks[index][1]
    }
    if (verbose & length(binNames)>1){ close(pb) }
    observedSPD <- vector("list",length=length(unique(GroupList)))
    names(observedSPD) <- unique(GroupList)
    for (d in 1:length(unique(GroupList))){
        focus <- unique(GroupList)[d]
        index <- which(GroupList==focus)
        tmpSPD <- apply(binnedMatrix[,index], 1, sum)
        if (!is.na(runm)){
            tmpSPD <- runMean(tmpSPD, runm, edge="fill")
        }
        if (focus==prfc){
            focd <- tmpSPD
        }
        if (d==1){
            dall <- tmpSPD
        } else {
            dall <- dall+tmpSPD
        }
        tmpSPD <- tmpSPD / sum(tmpSPD)
        observedSPD[[d]] <- data.frame(calBP=calyears, PrDens=tmpSPD)
    }
    propdf$ObsProp <- focd / dall
    propdf$Denom <- dall
    simprop <- matrix(nrow=nrow(calyears), ncol=nsim)
    ## Permutations
    simulatedSPD <- vector("list",length=length(unique(GroupList)))
    for (d in 1:length(unique(GroupList))){
        simulatedSPD[[d]] <- matrix(NA, nrow=nrow(calyears), ncol=nsim)
    }
    if (verbose){
        print("Permuting the groups...")
        flush.console()
        pb <- txtProgressBar(min=1, max=nsim, style=3)
    }
    for (s in 1:nsim){
        if (verbose){ setTxtProgressBar(pb, s) }
        simGroupList <- sample(GroupList)
        for (d in 1:length(unique(simGroupList))){
            focus <- unique(GroupList)[d]
            index <- which(simGroupList==focus)
            tmpSPD <- apply(binnedMatrix[,index],1,sum)
            if (!is.na(runm)){
                tmpSPD <- runMean(tmpSPD, runm, edge="fill")
            }
            if (focus==prfc){
                focd <- tmpSPD
            }
            if (d==1){
                focd <- dall <- tmpSPD
            } else {
                dall <- dall+tmpSPD
            }
            tmpSPD <- tmpSPD/sum(tmpSPD)
            simulatedSPD[[d]][,s] <- tmpSPD
        }
        simprop[,s] <- focd / dall
    }
    names(simulatedSPD) <- unique(GroupList)
    if (verbose){ close(pb) }
    tmp <- apply(simprop,1,FUN=function(x) quantile(x,c(0.025,0.975), na.rm=TRUE))
    ## Statistics
    propdf$EnvLo <- tmp[1,]
    propdf$EnvHi <- tmp[2,]
    simulatedCIlist <- vector("list",length=length(unique(GroupList)))
    for (d in 1:length(unique(GroupList))){
        simulatedCIlist[[d]] <- cbind(apply(simulatedSPD[[d]],1,quantile,prob=c(0.025)), apply(simulatedSPD[[d]],1,quantile,prob=c(0.975)))
        names(simulatedCIlist) <- unique(GroupList)
    }
    pValueList <- numeric(length=length(simulatedSPD))
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
        names(pValueList) <- unique(GroupList)
    }        
    res <- list(observed=observedSPD, envelope=simulatedCIlist, proportions=propdf, pValueList=pValueList)
    if (raw){ res$raw <- simulatedSPD }
    class(res) <- "SpdPermTest"
    if (verbose){ print("Done.") }
    return(res)
}


