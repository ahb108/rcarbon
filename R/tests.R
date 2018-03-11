 if(getRversion() >= "2.15.1")  utils::globalVariables(c("s"))



#' @title Monte-Carlo simulation test for SPDs 
#'
#' @description Comparison of an observed summed radiocarbon date distribution (aka SPD) with simulated outcomes from a theoretical model.
#'
#' @param x A vector of radiocarbon ages 
#' @param errors A vector of errors corresponding to each radiocarbon age
#' @param nsim Number of simulations
#' @param bins A vector indicating which bin each radiocarbon date is assigned to.
#' @param runm A number indicating the window size of the moving average to smooth both observed and simulated SPDs. If set to \code{NA} no moving average is applied.Default is \code{NA}. 
#' @param timeRange  A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param raw A logical variable indicating whether all permuted SPDs should be returned or not. Default is FALSE.
#' @param model A vector indicating the model to be fitted. Currently the acceptable options are \code{'uniform'}, \code{'linear'}, \code{'exponential'} and \code{'custom'}. 
#' @param predgrid A data.frame containing calendar years (column \code{calBP} and associated summed probabilities (column \code{PrDens}). Required when \code{model} is set to \code{'custom'}.
#' @param calCurves A vector of calibration curves (one between 'intcal13','shcal13' and 'marine13'; default is 'intcal13')
#' @param datenormalised If set to TRUE the total probability mass of each calibrated date will be made to sum to unity (the default in most radiocarbon calibration software). This argument will only have an effect if the dates in \code{x} were calibrated without normalisation (via normalised=FALSE in the \code{\link{calibrate}} function), in which case setting \code{datenormalised=TRUE} here will rescale each dates probability mass to sum to unity before aggregating the dates, while setting \code{datenormalised=FALSE} will ensure unnormalised dates are used for both observed and simulated SPDs. Default is FALSE.
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity for both observed and simulated data. 
#' @param ncores Number of cores used for for parallel execution. Default is 1.
#' @param fitonly A logical variable. If set to TRUE, only the the model fitting is executed and returned. Default is FALSE.
#' @param a Starter value for the exponential fit with the \code{\link{nls}} function using the formula \code{y ~ exp(a + b * x)} where \code{y} is the summed probability and \code{x} is the date. Default is 0.
#' @param b Starter value for the exponential fit with the \code{\link{nls}} function using the formula \code{y ~ exp(a + b * x)} where \code{y} is the summed probability and \code{x} is the date. Default is 0. 
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @details The function implements a modified version of Timpson et al (2014) Monte-Carlo test for comparing a theoretical or fitted statistical model to an observed summed radiocarbon date distribution (aka SPD). A variety of theoretical expectations can be compared to the observed distribution by setting the \code{model} argument, for example to fit basic \code{'uniform'} (the mean of the SPD), \code{'linear'} (fitted using the \code{\link{lm}} function) or \code{model='exponential'} models (fitted using the \code{\link{nls}} function). Models are fitted to the period spanned by \code{timeRange} although \code{x} can contain dates outside this range to mitigate possible edge effects (see also \code{bracket}). Alternatively, it is possible for the user to provide a model of their own by setting \code{model='custom'} and then supplying a two-column data.frame to \code{predgrid}. The chosen model is then 'uncalibrated' (see \code{\link{uncalibrate}}) and \emph{n} radiocarbon ages are randomly drawn, with \emph{n} equivalent to the number of dates or number of unique site/phase bins if the latter are supplied by the \code{bin} argument. The simulated dates are then calibrated and an SPD for each simulation. This process is repeated \code{nsim} times to produce a set of simulated expected probabilities densities per each calendar year. The probabilities are then z-transformed, and a 95\% critical envelope is computed. Local departures from the model are defined as instances where the observed SPD (which is also z-transformed) is outside such an envelope, while an estimate of the global significance of the observed SPD is also computed by comparing the total areas of observed and simulated SPDs that fall outside the simulation envelope. 
#'
#' @return An object of class \code{SpdModelTest} with the following elements
#' \itemize{
#' \item{\code{result}} {A four column data.frame containing the observed probability density (column \emph{PrDens}) and the lower and the upper values of the simulation envelope (columns \emph{lo} and \emph{hi}) for each calendar year column \emph{calBP}}
#' \item{\code{sim}} {A matrix containing the simulation results. Available only when \code{raw} is set to TRUE} 
#' \item{\code{pval}} {A numeric vector containing the p-value of the global significance test.}  
#' \item{\code{fit}} {A data.frame containing the probability densities of the fitted model for each calendar year within the time range of analysis}  
#' \item{\code{fitobject}} {Fitted model. Not available when \code{model} is \code{'custom'}}  
#' }
#'
#' @references 
#' Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., (2014). Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method. Journal of Archaeological Science 52, 549-557. doi:10.1016/j.jas.2014.08.011
#'
#'
#' @examples
#' ## Example with Younger Dryas period Near East, including site bins
#' \dontrun{
#' data(emedyd)
#' caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE, calMatrix=TRUE)
#' bins <- binPrep(sites=emedyd$SiteName, ages=emedyd$CRA, h=50)
#' nsim=5 #toy example
#' expnull <- modelTest(caldates, errors=emedyd$Error, bins=bins, nsim=nsim, runm=50,
#' timeRange=c(16000,9000), model="exponential", datenormalised=FALSE)
#' plot(expnull, xlim=c(16000,9000))
#' round(expnull$pval,4) #p-value
#' }
#' @import utils
#' @import stats
#' @import foreach
#' @import parallel
#' @import doParallel
#' @export

modelTest <- function(x, errors, nsim, bins=NA, runm=NA, timeRange=NA, raw=FALSE, model=c("exponential","uniform","linear","custom"), predgrid=NA, calCurves='intcal13', datenormalised=FALSE, spdnormalised=FALSE, ncores=1, fitonly=FALSE, a=0, b=0, verbose=TRUE){
    
    if (ncores>1&!requireNamespace("doParallel", quietly=TRUE)){	
	warning("the doParallel package is required for multi-core processing; ncores has been set to 1")
	ncores=1
    } else {
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      on.exit(stopCluster(cl))	
    }
    if (verbose){ print("Aggregating observed dates...") }
    if (is.na(bins[1])){
        samplesize <- nrow(x$metadata)
    } else {
        samplesize <- length(unique(bins))
    }
    observed <- spd(x=x, bins=bins, timeRange=timeRange, datenormalised=datenormalised, runm=runm, spdnormalised=spdnormalised, verbose=FALSE)
    finalSPD <- observed$grid$PrDens
    if (fitonly == TRUE) {nsim <- 1}
    ## Simulation
    sim <- matrix(NA,nrow=length(finalSPD),ncol=nsim)
    if (verbose & !fitonly){
        print("Monte-Carlo test...")
    flush.console()
        pb <- txtProgressBar(min=1, max=nsim, style=3)
    }
    time <- seq(timeRange[1],timeRange[2],-1)
    fit <- NA
    if (model=="exponential"){
        fit <- nls(y ~ exp(a + b * x), data=data.frame(x=time, y=finalSPD), start=list(a=a, b=b))
        est <- predict(fit, list(x=time))
        predgrid <- data.frame(calBP=time, PrDens=est)
    } else if (model=="uniform"){
        predgrid <- data.frame(calBP=time, PrDens=mean(finalSPD))
    } else if (model=="linear"){
        fit <- lm(y ~ x, data=data.frame(x=time, y=finalSPD))
        est <- predict(fit, list(x=time))
        predgrid <- data.frame(calBP=time, PrDens=est)
    } else if (model=="custom"){
        if (length(predgrid)!=2){
            stop("If you choose a custom model, you must provide a proper predgrid argument (two-column data.frame of calBP and predicted densities).")
        }
    } else {
        stop("Specified model not one of current choices.")
    }
    if (fitonly){
        print("Done (SPD and fitted model only).")
        res <- list(result=NA, sim=NA, pval=NA, osbSPD=observed, fit=predgrid, fitobject=fit)
        return(res)
    }
    cragrid <- uncalibrate(as.CalGrid(predgrid), calCurves=calCurves, compact=FALSE, verbose=FALSE)
    cragrid <- cragrid[cragrid$CRA <= max(x$metadata$CRA) & cragrid$CRA >= min(x$metadata$CRA),]

    if (ncores==1)
    {
    for (s in 1:nsim){ 
	if (verbose){ setTxtProgressBar(pb, s) } 
        randomDates <- sample(cragrid$CRA, replace=TRUE, size=samplesize, prob=cragrid$PrDens) 
        randomSDs <- sample(size=length(randomDates), errors, replace=TRUE) 
        tmp <- calibrate(x=randomDates,errors=randomSDs, timeRange=timeRange, calCurves=calCurves, normalised=datenormalised, ncores=1, verbose=FALSE, calMatrix=TRUE) 
        simDateMatrix <- tmp$calmatrix 
	sim[,s] <- apply(simDateMatrix,1,sum) 
        sim[,s] <- (sim[,s]/sum(sim[,s])) * sum(predgrid$PrDens[predgrid$calBP <= timeRange[1] & predgrid$calBP >= timeRange[2]]) 
        if (spdnormalised){ sim[,s] <- (sim[,s]/sum(sim[,s])) } 
        if (!is.na(runm)){ sim[,s] <- runMean(sim[,s], runm, edge="fill") }	
    	}
    }	    

    if (ncores>1)
    {	     
	print("Progress bar disabled for multi-core processing")
    	sim <- foreach (s = 1:nsim, .combine='cbind', .packages='rcarbon') %dopar% {
        # if (verbose){ setTxtProgressBar(pb, s) }
        randomDates <- sample(cragrid$CRA, replace=TRUE, size=samplesize, prob=cragrid$PrDens)
        randomSDs <- sample(size=length(randomDates), errors, replace=TRUE)
        tmp <- calibrate(x=randomDates,errors=randomSDs, timeRange=timeRange, calCurves=calCurves, normalised=datenormalised, ncores=1, verbose=FALSE, calMatrix=TRUE)
        simDateMatrix <- tmp$calmatrix
        aux <- apply(simDateMatrix,1,sum)
        aux <- (aux/sum(aux)) * sum(predgrid$PrDens[predgrid$calBP <= timeRange[1] & predgrid$calBP >= timeRange[2]])
        if (spdnormalised){ aux <- (aux/sum(aux)) }
        if (!is.na(runm)){
            aux <- runMean(aux, runm, edge="fill")
        }
	aux
   	 }
    	#stopCluster(cl)
     }

    if (verbose){ close(pb) }
    ## Envelope, z-scores, global p-value
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
    pvalue <- 1 - c(length(expectedstatistic[expectedstatistic <= observedStatistic])+1)/c(length(expectedstatistic)+1)
    # Results
    result <- data.frame(calBP=observed$grid$calBP,PrDens=finalSPD,lo=lo,hi=hi)
    if(raw==FALSE){ sim <- NA }
    res <- list(result=result, sim=sim, pval=pvalue, fit=predgrid, fitobject=fit)
    class(res) <- "SpdModelTest"
    if (verbose){ print("Done.") }
    return(res)
}


#' @title  Random mark permutation test for SPDs
#'
#' @description Global and local significance test for comparing shapes of multiple SPDs using random permutations. 
#'
#' @param x A \code{CalDates} class object containing the calibrated radiocarbon dates.
#' @param marks A numerical or character vector containing the marks associated to each radiocarbon date.
#' @param timeRange  A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param bins  A vector indicating which bin each radiocarbon date is assigned to.  
#' @param nsim Number of random permutations 
#' @param runm A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied.Default is NA.  
#' @param datenormalised If set to TRUE the total probability mass of each calibrated date will be made to sum to unity (the default in most radiocarbon calibration software). This argument will only have an effect if the dates in \code{x} were calibrated without normalisation (via normalised=FALSE in the \code{\link{calibrate}} function), in which case setting \code{datenormalised=TRUE} here will rescale each dates probability mass to sum to unity before aggregating the dates, while setting \code{datenormalised=FALSE} will ensure unnormalised dates are used for both observed and simulated SPDs. Default is FALSE.
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity. 
#' @param raw A logical variable indicating whether all permuted SPDs should be returned or not. Default is FALSE.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @details The function generates a distribution of expected SPDs by randomly shuffling the marks assigned to each \emph{bin} (see \code{\link{spd}} for details on binning). The resulting distribution of probabilities for each \emph{mark} (i.e. group of dates) for each calendar year is z-transformed, and a 95\% simulation envelope is computed. Local significant departures are defined as instances where the observed SPD (which is also z-transformed) is outside such envelope. A global significance is also computed by comparing the total "area" outside the simulation envelope in the observed and simulated data. 
#'
#' @return An object of class \code{SpdPermTest} with the following elements
#' \itemize{
#' \item{\code{observed}} {A list containing data.frames with the summed probability (column \emph{PrDens} for each calendar year (column \emph{calBP} for each mark/group}
#' \item{\code{envelope}} {A list containing matrices with the lower and upper bound values of the simulation envelope for each mark/group} 
#' \item{\code{pValueList}} {A list of p-value associated with each mark/group}  
#' }
#'
#' @references 
#' Crema, E.R., Habu, J., Kobayashi, K., Madella, M., (2016). Summed Probability Distribution of 14 C Dates Suggests Regional Divergences in the Population Dynamics of the Jomon Period in Eastern Japan. PLOS ONE 11, e0154809. doi:10.1371/journal.pone.0154809
#'
#' @examples
#' ## compare demographic trajectories in Netherlands and Denmark  
#' \dontrun{ 
#' data(euroevol)
#' nld.dnk = subset(euroevol,Country=="Netherlands"|Country=="Denmark")
#' bins = binPrep(nld.dnk$SiteID,nld.dnk$C14Age,h=200)
#' dates = calibrate(nld.dnk$C14Age,nld.dnk$C14SD,normalised=FALSE)
#' res = permTest(dates,marks=as.character(nld.dnk$Country),nsim=1000,
#' bins=bins,runm=200,timeRange=c(10000,4000))
#' round(res$pValueList,4) #extract p-values
#' par(mfrow=c(2,1))
#' plot(res,focalm="Netherlands",main="Netherlands")
#' plot(res,focalm="Denmark",main="Denmark")
#' }
#' @import utils
#' @import stats
#' @import foreach
#' @import parallel
#' @import doParallel
#' @export


permTest <- function(x, marks,  timeRange, nsim, bins=NA, runm=NA, datenormalised=FALSE, spdnormalised=FALSE, raw=FALSE, verbose=TRUE){

    if (is.na(bins[1])){ bins <- as.character(1:nrow(x$metadata)) }
    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    binnedMatrix <- matrix(nrow=nrow(calyears), ncol=length(binNames))
    GroupList <- vector()
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
        tmpSPD <- apply(binnedMatrix[,index,drop=FALSE], 1, sum)
        if (!is.na(runm)){
            tmpSPD <- runMean(tmpSPD, runm, edge="fill")
        }
        if (d==1){
            dall <- tmpSPD
        } else {
            dall <- dall+tmpSPD
        }
        if (spdnormalised){ tmpSPD <- tmpSPD / sum(tmpSPD) }
        observedSPD[[d]] <- data.frame(calBP=calyears, PrDens=tmpSPD)
    }
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
            tmpSPD <- apply(binnedMatrix[,index,drop=FALSE],1,sum)
            if (!is.na(runm)){
                tmpSPD <- runMean(tmpSPD, runm, edge="fill")
            }
            if (d==1){
                dall <- tmpSPD
            } else {
                dall <- dall+tmpSPD
            }
            if (spdnormalised){ tmpSPD <- tmpSPD/sum(tmpSPD) }
            simulatedSPD[[d]][,s] <- tmpSPD
        }
    }
    names(simulatedSPD) <- unique(GroupList)
    if (verbose){ close(pb) }
    ## Statistics
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
            pValueList[[a]] <- 1 - c(length(expectedstatistic[expectedstatistic <= observedStatistic])+1)/c(length(expectedstatistic)+1)
        }
        names(pValueList) <- unique(GroupList)
    }        
    res <- list(observed=observedSPD, envelope=simulatedCIlist, pValueList=pValueList)
    if (raw){ res$raw <- simulatedSPD }
    class(res) <- "SpdPermTest"
    if (verbose){ print("Done.") }
    return(res)
}



#' @title Spatial Permutation Test of summed probability distributions.
#'
#' @description This function carries out local spatial permutation tests of summed radiocarbon probability distributions in order to detect local deviations in growth rates (Crema et al 2017).
#' 
#' @param calDates  A \code{CalDates} class object.
#' @param timeRange A vector of length 2 indicating the start and end date of the analysis in cal BP
#' @param bins A vector indicating which bin each radiocarbon date is assigned to. Must have the same length as the number of radiocarbon dates. Can be created using the  \code{\link{binPrep}}) function. Bin names should follow the format "x_y", where x refers to a unique location (e.g. a site) and y is a integer value (e.g. "S023_1", "S023_2","S034_1", etc.).  
#' @param locations A \code{SpatialPoints} or a \code{SpatialPointsDataFrame} class object. Rownames of each point should much the first part of the bin names supplied (e.g. "S023","S034") 
#' @param breaks A vector of break points for defining the temporal slices.
#' @param spatialweights A \code{spatialweights} class object defining the spatial weights between the locations (cf. \code{\link{spweights}})
#' @param nsim The total number of simulations. Default is 1000.
#' @param runm The window size of the moving window average. Must be set to \code{NA} if the rates of change are calculated from the raw SPDs. 
#' @param permute Indicates whether the permutations should be based on the \code{"bins"} or the \code{"locations"}. Default is \code{"locations"}. 
#' @param ncores Number of cores used for for parallel execution. Default is 1.
#' @param datenormalised A logical variable indicating whether the probability mass of each date within \code{timeRange} is equal to 1. Default is FALSE. 
#' @param raw A logical variable indicating whether permuted sets of geometric growth rates for each location should be returned. Default is FALSE. 
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#'
#' @details The function consists of the following seven steps: 1) for each location (e.g. a site) generate a local SPD of radiocarbon dates, weighting the contribution of dates from neighbouring sites using a weight scheme provided by the \code{spatialweights} class object; 2) define temporal slices (using \code{breaks} as break values), then compute the total probability mass within each slice; 3) compute the rate of change between abutting temporal slices by using the formula: \eqn{(SPD_{t}/SPD_{t+1}^{1/\Delta t}-1)}; 4) randomise the location of individual bins or the entire sequence of bins associated with a given location and carry out steps 1 to 3; 5) repeat step 4 \code{nsim} times and generate, for each location, a distribution of growth rates under the null hypothesis (i.e. spatial independence); 6) compare, for each location, the observed growth rate with the distribution under the null hypothesis and compute the p-values; and 7) compute the false-discovery rate for each location.    
#'
#' @return A \code{spatialTest} class object
#'
#' @references
#' Crema, E.R., Bevan, A., Shennan, S. (2017). Spatio-temporal approaches to archaeological radiocarbon dates. Journal of Archaeological Science, 87, 1-9.
#' 
#' @seealso \code{\link{permTest}} for a non-spatial permutation test; \code{\link{plot.spatialTest}} for plotting; \code{\link{spweights}} for computing spatial weights; \code{\link{spd2gg}} for computing geometric growth rates.
#'
#' @examples
#' ## Reproduce Crema et al 2017 ##
#'\dontrun{
#' data(euroevol) #load data
#'
#' ## Subset only for 8000 to 5000 Cal BP (c7200-4200 C14BP)
#' edge=800
#' timeRange=c(8000,5000)
#' euroevol2=subset(euroevol,C14Age<=c(timeRange[1]-edge)&C14Age>=c(timeRange[2]-edge))
#'
#' ## define chronological breaks
#' breaks=seq(8000,5000,-500)
#'
#' ## Create a SpatialPoints class object 
#' library(sp)
#' sites = unique(data.frame(SiteID=euroevol2$SiteID,
#' Longitude=euroevol2$Longitude,Latitude=euroevol2$Latitude))
#' locations=data.frame(Longitude=sites$Longitude,Latitude=sites$Latitude)
#' rownames(locations)=sites$SiteID
#' coordinates(locations)<-c("Longitude","Latitude")
#' proj4string(locations)<- CRS("+proj=longlat +datum=WGS84")
#'
#' ## Compute Distance and Spatial Weights 
#' distSamples=spDists(locations,locations,longlat = TRUE)
#' spatialweights=spweights(distSamples,h=100) #using a kernal bandwidth of 100km
#'
#' ## Calibration and binning
#' bins=binPrep(sites=euroevol2$SiteID,ages=euroevol2$C14Age,h=200)  
#' calDates=calibrate(x=euroevol2$C14Age,errors=euroevol2$C14SD,
#' timeRange=timeRange,normalised=FALSE)
#'
#' ## Main Analysis (over 2 cores; requires doParallel package) 
#' ## NOTE: the number of simulations should be ideally larger 
#' ## to ensure a better resolution of the p/q-values.
#' res.locations=SPpermTest(calDates,timeRange=timeRange,bins=bins,locations=locations,
#' spatialweights=spatialweights,breaks=breaks,ncores=2,nsim=100,
#' permute="locations",datenormalised=FALSE)
#' 
#' ## Plot results
#' library(rworldmap)
#' base=getMap(resolution="low") #optionally add base map
#' #retrieve coordinate limits#
#' xrange=bbox(res.locations$locations)[1,] 
#' yrange=bbox(res.locations$locations)[2,]
#'
#' par(mfrow=c(2,2))  
#' par(mar=c(0.1,0.1,0,0.5))
#' plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange)
#' plot(res.locations,index=4,add=TRUE,option="raw",breakRange=c(-0.005,0.005)) 
#' plot(res.locations,option="rawlegend",breakRange=c(-0.005,0.005),rd=3)
#' par(mar=c(0.1,0.1,0,0.5))
#' plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange) 
#' plot(res.locations,index=4,add=TRUE,option="test")  
#' plot(res.locations,option="testlegend")
#' }
#' @import utils
#' @import stats
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import sp
#' @export
 

SPpermTest<-function(calDates, timeRange, bins, locations, breaks, spatialweights, nsim=1000, runm=NA,permute="locations",ncores=1,datenormalised=FALSE,verbose=TRUE,raw=FALSE)
{

###################################
#### Load Dependency Libraries ####
###################################
    if (ncores>1&!requireNamespace("doParallel", quietly=TRUE)){	
	warning("the doParallel package is required for multi-core processing; ncores has been set to 1")
	ncores=1
    }	

##################################
#### Initial warning messages ####
##################################

     if (!"CalDates" %in% class(calDates)){
        stop("calDates must be an object of class 'calDates'.")
    }
    if (length(bins)>1){
         if (any(is.na(bins))){
            stop("Cannot have NA values in bins.")
        }
        if (length(bins)!=length(calDates$grid)){
            stop("bins (if provided) must be the same length as x.")
        }
         } else {
        bins <- rep("0_0",length(calDates$grid))
    }

   if (!("SpatialPoints" %in% class(locations)[1]|"SpatialPointsDataFrame" %in% class(locations)[1])){
        stop("locations must be an object of class 'SpatialPoints' or 'SpatialPointsDataFrame'.")
    }

   locations.id=row.names(locations@coords)
    if (is.null(locations.id))
    {
        stop("locations must have rownames")
    }
    if (!all(range(timeRange)==range(breaks)))
    {
	stop("Range of breaks values must much match the temporal range defined by timeRange")
    }
   
    if (ncores>1&raw==TRUE)
    {
	warnings("raw==TRUE available only for ncores=1")
    	raw=FALSE
    }
#############################
#### Create binnedMatrix ####
#############################

    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))


    if (verbose & length(binNames)>1){
        print("Binning by site/phase...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Binning by site/phase...")
    }
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        slist <- calDates$grid[index]
        slist <- lapply(slist,FUN=function(x) merge(calyears,x, all.x=TRUE)) 
        slist <- rapply(slist, f=function(x) ifelse(is.na(x),0,x), how="replace")
        slist <- lapply(slist, FUN=function(x) x[with(x, order(-calBP)), ])
        tmp <- lapply(slist,`[`,2)
        if (datenormalised){
            tmp <- lapply(tmp,FUN=function(x) x/sum(x))
        }
        if (length(binNames)>1){
            spd.tmp <- Reduce("+", tmp) / length(index)
        } else {
            spd.tmp <- Reduce("+", tmp)
        }
	binnedMatrix[,b] <- spd.tmp[,1]
    }
    if (verbose & length(binNames)>1){ close(pb) }


################################
### Observed Data Subroutine ###
################################ 


## Aggregate by Locations ##
    origins=unlist(lapply(strsplit(binNames,"_"),function(x){x[[1]]}))

    if (!all(origins%in%locations.id))
     {
        stop("Missing bins or locations")
     }

    resMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nrow(binnedMatrix))
    
    for (x in 1:length(unique(locations.id))) 
        {
            index=which(origins==unique(locations.id)[x])
            if(length(index)>1)
                {resMatrix[x,]=apply(binnedMatrix[,index],1,sum)}
            if(length(index)==1)
                {resMatrix[x,]=binnedMatrix[,index]}
        }

 ## Aggregate by break s##

    nBreaks=length(breaks)-1
    obsMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nBreaks)
    timeSequence=timeRange[1]:timeRange[2]
    
    for (x in 1:nBreaks)
        {
            index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
            obsMatrix[,x]=apply(resMatrix[,index],1,sum)
        }

## Apply SpatialWeights ##

    obsGridVal=t(spatialweights$w)%*%obsMatrix

## Compute Rate of Change #3

    rocaObs=t(apply(obsGridVal,1,function(x,d){
		   L=length(x)
		   res=numeric(length=L-1)
		   for (i in 1:c(L-1))
			{
			res[i]=(x[i+1]/x[i])^(1/d)-1
			}
		   return(res)},
		   d=abs(breaks[2]-breaks[1])))

##############################
### Permutation Subroutine ###
############################## 

    if (ncores>1)
   	 {
          cl <- makeCluster(ncores)
          registerDoParallel(cl)
          print(paste("Running permutation test in parallel on ",getDoParWorkers()," workers...",sep=""))
	  sumcombine<-function(a,b)
		{
		list(a[[1]]+b[[1]],a[[2]]+b[[2]],a[[3]]+b[[3]])
		}
	  resultHiLoEq<-foreach (x=1:nsim,.combine= sumcombine) %dopar% {

            simGridVal<-matrix(NA,nrow=nrow(spatialweights$w),ncol=nBreaks)
            
	    ## Aggregate by Site ## 

            simResMatrix=matrix(0,nrow=length(unique(locations.id)),ncol=nrow(binnedMatrix))

            ## Randomly assigne bins to locations.id ##
           
	    if (permute=="bins")
	    {    
	    simOrigins=sample(origins)
            for (x in 1:length(unique(locations.id)))
                {                    
                    index=which(simOrigins==unique(locations.id)[x])
                    if(length(index)>1)
                        {simResMatrix[x,]=apply(binnedMatrix[,index],1,sum)}
                    if(length(index)==1)
                        {simResMatrix[x,]=binnedMatrix[,index]}
                }
            

            ## Aggregate by breaks ##
	    
            aggMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nBreaks)
            
            for (x in 1:nBreaks)
                {
                    index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
                    aggMatrix[,x]=apply(simResMatrix[,index],1,sum)
                }
		       

           ## Apply Weights ##

           simGridVal=t(spatialweights$w)%*%aggMatrix
	    }
	    if (permute=="locations")
	    {
	     simMatrix=obsMatrix[sample(nrow(obsMatrix)),]	
             simGridVal=t(spatialweights$w)%*%simMatrix
		
	    }


           ## Compute Rate of Change ##

           rocaSim=t(apply(simGridVal,1,function(x,d){
		   L=length(x)
		   res=numeric(length=L-1)
		   for (i in 1:c(L-1))
			{
			res[i]=(x[i+1]/x[i])^(1/d)-1
			}
		   return(res)},
		   d=abs(breaks[2]-breaks[1])))

	    lo=rocaObs<rocaSim	    
	    hi=rocaObs>rocaSim
	    eq=rocaObs==rocaSim

          return(list(hi,lo,eq))
	  }
        stopCluster(cl)

        lo=resultHiLoEq[[1]]
	hi=resultHiLoEq[[2]]
	eq=resultHiLoEq[[3]]
    
	} else {

    hi=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)
    lo=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)
    eq=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)

    print("Permutation test...")
    flush.console()
    
    if(raw)
    {
	rocaSimAll = array(NA,dim=c(nsim,nrow(spatialweights$w),nBreaks-1))
    }
    pb <- txtProgressBar(min = 1, max = nsim, style=3)

        for (s in 1:nsim)
        {
            setTxtProgressBar(pb, s)
	    simGridVal<-matrix(NA,nrow=nrow(spatialweights$w),ncol=nBreaks)
            ## Aggregate by Site ## 
            simResMatrix=matrix(0,nrow=length(unique(locations.id)),ncol=nrow(binnedMatrix))

            ## Randomly assign bins to locations
            if (permute=="bins")
	    {
	    simOrigins=sample(origins)
            



            for (x in 1:length(unique(locations.id)))
                {                    
                    index=which(simOrigins==unique(locations.id)[x])
                    if(length(index)>1)
                        {simResMatrix[x,]=apply(binnedMatrix[,index],1,sum)}
                    if(length(index)==1)
                        {simResMatrix[x,]=binnedMatrix[,index]}
                }
            

            ##Aggregate by breaks##
            aggMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nBreaks)
            
            for (x in 1:nBreaks)
                {
                    index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
                    aggMatrix[,x]=apply(simResMatrix[,index],1,sum)
                }
		       

           ##Apply Weights 
           simGridVal=t(spatialweights$w)%*%aggMatrix
	    }
           if (permute=="locations")
           {
	     simMatrix=obsMatrix[sample(nrow(obsMatrix)),]	
             simGridVal=t(spatialweights$w)%*%simMatrix
	   }



           ##Compute Rate of Change
           rocaSim=t(apply(simGridVal,1,function(x,d){
		   L=length(x)
		   res=numeric(length=L-1)
		   for (i in 1:c(L-1))
			{
			res[i]=(x[i+1]/x[i])^(1/d)-1
			}
		   return(res)},
		   d=abs(breaks[2]-breaks[1])))

	    hi=hi+(rocaObs>rocaSim)
	    lo=lo+(rocaObs<rocaSim)
	    eq=eq+(rocaObs==rocaSim)
	    rocaSimAll[s,,]=rocaSim
        }
    close(pb)
    }


############################
### Compute Significance ###
############################ 
    
    pvalHi=(lo+eq+1)/c(nsim+1)
    pvalLo=(hi+eq+1)/c(nsim+1)
    pval=pvalHi
    pval[which(pvalHi>pvalLo)]=pvalLo[which(pvalHi>pvalLo)]
    pval=pval*2
    if (max(pval)>1)
    {
    	 pval[which(pval>1)]=1
    }

    ## Compute False Discovery Rate ##
    
    qvalHi=apply(pvalHi,2,function(x){return(p.adjust(x,method="fdr"))})
    qvalLo=apply(pvalLo,2,function(x){return(p.adjust(x,method="fdr"))})
    qval=apply(pval,2,function(x){return(p.adjust(x,method="fdr"))})




    metadata=data.frame(npoints=length(unique(locations.id)),ndates=nrow(calDates$metadata),nbins=length(binNames),nsim=nsim,permutationType=permute,datenormalised=datenormalised,breaks=nBreaks,timeRange=paste(timeRange[1],"-",timeRange[2],sep=""),weights.h=spatialweights$h,weights.kernel=spatialweights$kernel)
  
    if(raw==FALSE){rocaSimAll=NA} 
    reslist=list(metadata=metadata,rocaSim=rocaSimAll,rocaObs=rocaObs,pval=pval,pvalHi=pvalHi,pvalLo=pvalLo,qval=qval,qvalLo=qvalLo,qvalHi=qvalHi,locations=locations)
    
    class(reslist) <- append(class(reslist),"spatialTest")
    return(reslist)
}



#' @title Point to point test of SPD 
#'
#' @description Test for evaluating the difference in the summed probability values associated with two points in time.
#'
#' @param x result of \code{\link{modelTest}} with raw=TRUE.
#' @param p1 calendar year (in BP) of start point.
#' @param p2 calendar year (in BP) of end point.
#' @param interactive if set to TRUE enables an interactive selection of p1 and p2 from a graphical display of the SPD. Disabled when \code{p1} and \code{p2} are defined.
#'
#' @details The function compares observed differences in the summed probability values associated with two points in time against a distribution of expected values under the null hypothesis defined with the \code{\link{modelTest}} function. The two points can be specified manually (assigning BP dates to the arguments \code{p1} and \code{p2}) or interactively (clicking on a SPD plot). Note that \code{\link{modelTest}} should be executed setting the argument \code{raw} to \code{TRUE} (default is \code{FALSE}.   
#'
#'
#' @return A list with: the BP dates for the two points and the p-value obtained from a two-sided test.
#'
#' @references 
#'  Edinborough, K., Porcic, M., Martindale, A., Brown, T.J., Supernant, K., Ames, K.M., (2017). Radiocarbon test for demographic events in written and oral history. PNAS 201713012. doi:10.1073/pnas.1713012114
#' @examples
#' ## Example with Younger Dryas period Near East, including site bins
#' \dontrun{
#' data(emedyd)
#' caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE, calMatrix=TRUE)
#' bins <- binPrep(sites=emedyd$SiteName, ages=emedyd$CRA, h=50)
#' nsim=10 #toy example
#' expnull <- modelTest(caldates, errors=emedyd$Error, bins=bins, nsim=nsim, runm=50,
#' timeRange=c(16000,9000), model="exponential", datenormalised=FALSE, raw=TRUE)
#' p2pTest(x=expnull,p1=13000,p2=12500) #non-interactive mode
#' p2pTest(x=expnull) #interactive mode
#' }
#' @seealso \code{\link{modelTest}}.
#' @import utils
#' @import stats
#' @export


p2pTest <- function(x,p1=NA,p2=NA,interactive=TRUE)
{

if (is.na(x$sim[1]))
{
 stop("x should be an SpdModelTest class object produced using modelTest() with raw=TRUE")
}

 if (!is.na(p1)&!is.na(p2))
         {
	 if (p2>p1){stop("the end point should be more recent than the start point")}
	 interactive=FALSE
	 }
 if (interactive)
 {
    plot(x)	
    print("select start point")
    p1=round(locator(n=1)$x[1])
    index1=match(p1,x$result$calBP)
    p1.y = x$result[index1,2]
    points(p1,p1.y,pch=20)

    print("select end point")
    p2=round(locator(n=1)$x[1])
    if (p2>p1)
    {
	print("the end point should be more recent than the start point")
    	while(p2>p1)
	{
         print("select end point")
         p2=round(locator(n=1)$x[1])
	}
    }	    
    index2=match(p2,x$result$calBP)
    p2.y = x$result[index2,2]
    points(p2,p2.y,pch=20)
    lines(x$result$calBP[index1:index2],x$result$PrDens[index1:index2],lwd=2)
 }

if (!interactive)
{
    index1=match(p1,x$result$calBP)
    p1.y = x$result[index1,2]
    index2=match(p2,x$result$calBP)
    p2.y = x$result[index2,2]
}

 if (!p1%in%x$result$calBP | !p2%in%x$result$calBP)
 {
	stop("p1 and p2 should be within the temporal range of the spd")
 }
 
 obs.diff = p1.y-p2.y
 sim.diff = x$sim[index1,]-x$sim[index2,] 
 nsim = ncol(x$sim) 
 lo = sum(obs.diff < sim.diff)
 hi = sum(obs.diff > sim.diff)
 eq = sum(obs.diff == sim.diff)
 pvalHi=(lo+eq+1)/c(nsim+1)
 pvalLo=(hi+eq+1)/c(nsim+1)
 pval=ifelse(pvalHi<pvalLo,pvalHi,pvalLo)*2
 
 return(list(p1=p1,p2=p2,pval=pval))	
}
