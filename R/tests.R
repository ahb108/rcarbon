 if(getRversion() >= "2.15.1")  utils::globalVariables(c("s","calBP"))



#' @title Monte-Carlo simulation test for SPDs 
#'
#' @description Comparison of an observed summed radiocarbon date distribution (aka SPD) with simulated outcomes from a theoretical model.
#'
#' @param x A \code{CalDates} object containing calibrated radiocarbon ages 
#' @param errors A vector of errors corresponding to each radiocarbon age
#' @param nsim Number of simulations
#' @param bins A vector indicating which bin each radiocarbon date is assigned to.
#' @param runm A number indicating the window size of the moving average to smooth both observed and simulated SPDs. If set to \code{NA} no moving average is applied.Default is \code{NA}. 
#' @param timeRange  A vector of length 2 indicating the start and end date of the analysis in cal BP. The fitting process is applied considering the SPD within the interval defined by this parameter. If no values are supplied the earliest and latest median calibrated dates of the observed data will be used.
#' @param backsight A single numeric value defining the distance in time between the focal year and the backsight year for computing the rate of change. Default is 10.
#' @param changexpr An expression for calculating the rate of change in SPD between the focal year and a backsight year. Available input options are t1 (the SPD for the focal year), t0 (the SPD for the backsight year), d (the distance between t0 and t1), and any other standard constants and mathematical operators. A sensible default is provided.
#' @param gridclip Whether the sampling of random dates is constrained to the observed range (TRUE) or not (FALSE). Default is TRUE.
#' @param raw A logical variable indicating whether all permuted SPDs should be returned or not. Default is FALSE.
#' @param model A vector indicating the model to be fitted. Currently the acceptable options are \code{'uniform'}, \code{'linear'}, \code{'exponential'} and \code{'custom'}. Default is \code{'exponential'}. 
#' @param method Method for the creation of random dates from the fitted model. Either \code{'uncalsample'} or \code{'calsample'}. Default is \code{'uncalsample'}. See below for details. 
#' @param predgrid A data.frame containing calendar years (column \code{calBP}) and associated summed probabilities (column \code{PrDens}). Required when \code{model} is set to \code{'custom'}.
#' @param normalised Whether the simulated dates should be normalised or not. Default based on whether x is normalised or not.  
#' @param datenormalised Argument kept for backward compatibility with previous versions.
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity for both observed and simulated data. 
#' @param edgeSize Controls edge effect by expanding the fitted model beyond the range defined by \code{timeRange}.
#' @param ncores Number of cores used for for parallel execution. Default is 1.
#' @param fitonly A logical variable. If set to TRUE, only the the model fitting is executed and returned. Default is FALSE.
#' @param a Starter value for the exponential fit with the \code{\link{nls}} function using the formula \code{y ~ exp(a + b * x)} where \code{y} is the summed probability and \code{x} is the date. Default is 0.
#' @param b Starter value for the exponential fit with the \code{\link{nls}} function using the formula \code{y ~ exp(a + b * x)} where \code{y} is the summed probability and \code{x} is the date. Default is 0. 
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @details The function implements a Monte-Carlo test for comparing a theoretical or fitted statistical model to an observed summed radiocarbon date distribution (aka SPD) and associated rates of changes. A variety of theoretical expectations can be compared to the observed distribution by setting the \code{model} argument, for example to fit basic \code{'uniform'} (the mean of the SPD), \code{'linear'} (fitted using the \code{\link{lm}} function) or \code{model='exponential'} models (fitted using the \code{\link{nls}} function). Models are fitted to the period spanned by \code{timeRange} although \code{x} can contain dates outside this range to mitigate possible edge effects (see also \code{bracket}). Alternatively, it is possible for the user to provide a model of their own by setting \code{model='custom'} and then supplying a two-column data.frame to \code{predgrid}. The function generates \code{nsim} theoretical SPDs from the fitted model via Monte-Carlo simulation, this is then used to define a 95\% critical envelope for each calendar year. The observed SPD is then compared against the simulation envelope; local departures from the model are defined as instances where the observed SPD is outside such an envelope, while an estimate of the global significance of the observed SPD is also computed by comparing the total areas of observed and simulated SPDs that fall outside the simulation envelope. The theoretical SPDs can be generated using two different sampling approaches defined by the parameter \code{method}. If \code{method} is set to \code{'uncalsample'} each date is drawn after the fitted model is backcalibrated as a whole and adjusted for a baseline expectation; if it is set to  \code{'calsample'} samples are drawn from the fitted model in calendar year then individually back calibrated and recalibrated (the approach of Timpson et al. 2014). For each simulation, both approaches produces \eqn{n} samples, with \eqn{n} equal to the number of bins or number of dates (when bins are not defined). Differences between these two approaches are particularly evident at dates coincident with steeper portions of the calibration curve. If more than one type of calibration curve is associated with the observed dates, at each Monte-Carlo iteration, the function randomly assigns each bin to one of the calibration curves with probability based on the proportion of dates within the bin associated to the specific curves. For example, if a bin is composed of four dates and three are calibrated with 'intcal20' the probability of that particular bin being assigned to 'intcal20' is 0.75.  
#' @note
#'\itemize{
#'\item {Windows users might receive a memory allocation error with larger time span of analysis (defined by the parameter \code{timeRange}). This can be avoided by increasing the memory limit with the \code{\link{memory.limit}} function.}
#'\item {Users experiencing a \code{Error: cannot allocate vector of size ...} error message can increase the memory size using the \code{Sys.setenv()}, for example: \code{Sys.setenv("R_MAX_VSIZE" = 16e9)}.} 
#'\item {The function currently supports only dates calibrated with 'intcal20',intcal13','intcal13nhpine16','shcal20','shcal13','shcal13shkauri16', 'marine20', and 'marine13'.}
#'}
#' 
#' @return An object of class \code{SpdModelTest} with the following elements
#' \itemize{
#' \item{\code{result}} {A four column data.frame containing the observed probability density (column \emph{PrDens}) and the lower and the upper values of the simulation envelope (columns \emph{lo} and \emph{hi}) for each calendar year column \emph{calBP}}
#' \item{\code{result.roc}} {A four column data.frame containing the observed rates of change (column \emph{roc}) and the lower and the upper values of the simulation envelope (columns \emph{lo.roc} and \emph{hi.roc}) for the mid point between two chronological blocks \emph{calBP}}
#' \item{\code{sim}} {A matrix containing the simulation results of the summed probabilities. Available only when \code{raw} is set to TRUE} 
#' \item{\code{sim.roc}} {A matrix containing the simulation results of the rate of change of summed probabilities. Available only when \code{raw} is set to TRUE} 
#' \item{\code{pval}} {A numeric vector containing the p-value of the global significance test for the summed probabilities}  
#' \item{\code{pval.roc}} {A numeric vector containing the p-value of the global significance test for the rates of change}  
#' \item{\code{fit}} {A data.frame containing the probability densities of the fitted model for each calendar year within the time range of analysis}  
#' \item{\code{fitobject}} {Fitted model. Not available when \code{model} is \code{'custom'}}  
#' \item{\code{n}} {Number of radiocarbon dates.}
#' \item{\code{nbins}}{Number of bins.}
#' \item{\code{nsim}}{Number of Monte-Carlo simulations.}
#' \item{\code{backsight}}{Backsight size.} 
#' }
#'
#' @references 
#' 
#' Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., (2014). Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method. Journal of Archaeological Science, 52, 549-557. doi:10.1016/j.jas.2014.08.011
#'
#'
#' @examples
#' ## Example with Younger Dryas period Near East, including site bins
#' \dontrun{
#' data(emedyd)
#' caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE)
#' bins <- binPrep(sites=emedyd$SiteName, ages=emedyd$CRA, h=50)
#' nsim=5 #toy example
#' expnull <- modelTest(caldates, errors=emedyd$Error, bins=bins, nsim=nsim, runm=50,
#' timeRange=c(16000,9000), model="exponential", datenormalised=FALSE)
#' plot(expnull, xlim=c(16000,9000))
#' round(expnull$pval,4) #p-value
#' summary(expnull)
#' }
#' @import utils
#' @import stats
#' @import doSNOW 
#' @import snow
#' @import foreach 
#' @import iterators 
#' @export

modelTest <- function(x, errors, nsim, bins=NA, runm=NA, timeRange=NA,backsight=50,changexpr=expression((t1/t0)^(1/d)-1),gridclip=TRUE, raw=FALSE, model=c("exponential"),method=c("uncalsample"),predgrid=NA, normalised=NA,datenormalised=NA, spdnormalised=FALSE, ncores=1, fitonly=FALSE, a=0, b=0, edgeSize=500,verbose=TRUE){
   
  caltimeRange =c(55000,0)
  if (any(x$metadata$CalCurve %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")))
  {
    caltimeRange =c(50000,0)
  }
    if (fitonly == TRUE) {nsim <- 1}
    if (ncores>1&!requireNamespace("doSNOW", quietly=TRUE)){	
	warning("the doSnow package is required for multi-core processing; ncores has been set to 1")
	ncores=1
    } else {
      cl <- snow::makeCluster(ncores)
      registerDoSNOW(cl)
      on.exit(snow::stopCluster(cl))	
    }
    if (!any(method%in%c("uncalsample","calsample")))
    {
	stop("The 'method' argument must be either 'uncalsample' or 'calsample'")
    }

    ccrange = c(max(medCal(x)),min(medCal(x)))

    if (anyNA(timeRange))
    {
	    timeRange=ccrange
    }

    if (is.na(normalised))
    {
	    normalised=FALSE    
	    if(x$metadata$Normalised[1]==TRUE)
	    {
		    normalised=TRUE
	    }
    }

    if (normalised!=x$metadata$Normalised[1])
    {
	warning("The normalisation setting of x and normalised are different")
    }

    if (!is.na(datenormalised))
    {
	    if (datenormalised!=normalised)
	    {
		    warning("'datenormalised' is not equal to 'normalised'. The datenormalised setting will be used for the normalisation setting of the calibration of simulated dates")
		    normalised=datenormalised
	    }
	    if (datenormalised!=x$metadata$Normalised[1])
	    {
		    if (x$metadata$Normalised[1])
		    {
			    warning("Input dates are normalised but datenormalised is set to FALSE. The datenormalised setting will be ignored")
			    normalised=datenormalised
		    }
		    if (!x$metadata$Normalised[1])
		    {

			    warning("Input dates are not normalised but datenormalised is set to TRUE. The datenormalised setting will be ignored")
			    normalised=datenormalised
		    }
	    }
    }



    calCurves = x$metadata$CalCurve
	if (!all(calCurves%in%c("intcal20","shcal20","marine20",'intcal13','intcal13nhpine16','shcal13','shcal13shkauri16','marine13')))
	{
	stop("modelTest() currently accepts only dates calibrated with the following calibration curves: 'intcal20','intcal13','intcal13nhpine16','shcal20','shcal13','shcal13shkauri16', 'marine20', and 'marine13'") 
	}

    unique.calCurves = as.character(sort(unique(calCurves)))
    ncc = length(unique.calCurves) #count number of unique calibration curves



    if (verbose){ print("Aggregating observed dates...") }

    #Generate matrix of sample sizes for each curve
    if (is.na(bins[1])){
        samplesize <- t(matrix(table(calCurves),nrow=ncc,ncol=nsim))
        colnames(samplesize) = names(table(calCurves))
    } else {
        samplesize <- curveSamples(bins=bins,calCurves=calCurves,nsim=nsim)
        if (ncc==1) {
          samplesize = matrix(samplesize,ncol=1,nrow=length(samplesize))
          colnames(samplesize) = names(table(calCurves))
        }
    }
    if (ncc>1) {samplesize=samplesize[,unique.calCurves]}

    
    # Create artificial bins in case bins are not supplied 
    if (is.na(bins[1])){ bins <- as.character(1:nrow(x$metadata)) }

    observed <- spd(x=x, bins=bins, timeRange=timeRange, runm=runm, spdnormalised=spdnormalised, verbose=FALSE, edgeSize=edgeSize)
    finalSPD <- observed$grid$PrDens
    
    ## Simulation
    sim <- matrix(NA,nrow=length(finalSPD),ncol=nsim)
    if (verbose & !fitonly){
	    print("Monte-Carlo test...")
	    flush.console()
    }
    fit.time <- seq(timeRange[1],timeRange[2],-1)
    pred.time <- fit.time
    if (gridclip) 
    {
	    st = max(ccrange[1],timeRange[1])+edgeSize
	    en = min(ccrange[2],timeRange[2])-edgeSize
	    pred.time <- seq(st,en,-1)
    }	

    fit <- NA
    if (model=="exponential"){
        fit <- nls(y ~ exp(a + b * x), data=data.frame(x=fit.time, y=finalSPD), start=list(a=a, b=b))
        est <- predict(fit, list(x=pred.time))
        predgrid <- data.frame(calBP=pred.time, PrDens=est)
    } else if (model=="uniform"){
        predgrid <- data.frame(calBP=pred.time, PrDens=mean(finalSPD))
    } else if (model=="linear"){
        fit <- lm(y ~ x, data=data.frame(x=fit.time, y=finalSPD))
        est <- predict(fit, list(x=pred.time))
        predgrid <- data.frame(calBP=pred.time, PrDens=est)
    } else if (model=="custom"){
	    if (length(predgrid)!=2){
		    stop("If you choose a custom model, you must provide a proper predgrid argument (two-column data.frame of calBP and predicted densities).")
	    }
	    if (!all(colnames(predgrid)%in%c("calBP","PrDens")))
	    {
		    stop("Column names in the predgrid argument should be 'calBP' and 'PrDens'")
	    }
    } else {
        stop("Specified model not one of current choices.")
    }
    if (fitonly){
	    print("Done (SPD and fitted model only).")
	    predgrid <- subset(predgrid,calBP<=timeRange[1]&calBP>=timeRange[2])
	    res <- list(result=NA, sim=NA, pval=NA, osbSPD=observed, fit=predgrid, fitobject=fit)
	    return(res)
    }


 # Add Extra Edges with PrDens=0   
    if (edgeSize>0)
    {
	predgrid = rbind.data.frame(data.frame(calBP=(max(predgrid$calBP)+edgeSize):c(predgrid$calBP[1]+1),PrDens=0),predgrid)
    	predgrid = rbind.data.frame(predgrid,data.frame(calBP=min(predgrid$calBP):(min(predgrid$calBP)-edgeSize),PrDens=0))
	if (any(predgrid$calBP<=0|predgrid$calBP>=caltimeRange[1]))
	    {
		    warning("edgeSize reduced")
		    predgrid = subset(predgrid, calBP<=caltimeRange[1]&calBP>=0)
	    }
    }

#  predgrid$PrDens = predgrid$PrDens/sum(predgrid$PrDens)

# Prepare Sampling Grid(s)
    cragrids = vector("list",length=ncc)
    for (i in 1:ncc)
    {		
	    tmp.grid <- uncalibrate(as.CalGrid(predgrid), calCurves=unique.calCurves[i], compact=FALSE, verbose=FALSE)
	    cragrids[[i]] <- tmp.grid

	    # Cllipping the uncalibrated grid
	    if (gridclip)
	    {
		    cragrids[[i]] <- tmp.grid[tmp.grid$CRA <= max(x$metadata$CRA) & tmp.grid$CRA >= min(x$metadata$CRA),]
	    }
    }


# Actual Method

    opts = NULL
    if (verbose)
    {
    if (ncores>1){ print(paste("Running in parallel on ",getDoParWorkers()," workers...",sep=""))}
    pb <- txtProgressBar(min=0, max=nsim, style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    }

    if (ncores==1)
    {
    for (s in 1:nsim){ 
	if (verbose){ setTxtProgressBar(pb, s) } 

    if (method=="uncalsample")
    {
	    randomDates <- vector("list",length=ncc)
	    ccurve.tmp <- numeric()
	    for (i in 1:ncc)
	    {       
		    randomDates[[i]] = sample(cragrids[[i]]$CRA,replace=TRUE,size=samplesize[s,i],prob=cragrids[[i]]$PrDens)
		    ccurve.tmp = c(ccurve.tmp,rep(unique.calCurves[i],samplesize[s,i]))
	    }

	    randomSDs <- sample(size=length(unlist(randomDates)), errors, replace=TRUE) 
    }


    if (method=="calsample")
    {
	    randomDates <- vector("list",length=ncc)
	    ccurve.tmp <- numeric()
	    for (i in 1:ncc)
	    {       
		    randomDates[[i]] = sample(cragrids[[i]]$CRA,replace=TRUE,size=samplesize[s,i],prob=cragrids[[i]]$Raw)
		    ccurve.tmp = c(ccurve.tmp,rep(unique.calCurves[i],samplesize[s,i]))
	    }

	    randomSDs <- sample(size=length(unlist(randomDates)), errors, replace=TRUE) 
    }


        tmp <- calibrate(x=unlist(randomDates),errors=randomSDs, timeRange=timeRange, calCurves=ccurve.tmp, normalised=normalised, ncores=1, verbose=FALSE, calMatrix=TRUE) 

        simDateMatrix <- tmp$calmatrix 
	sim[,s] <- apply(simDateMatrix,1,sum) 
        sim[,s] <- (sim[,s]/sum(sim[,s])) * sum(predgrid$PrDens[predgrid$calBP <= timeRange[1] & predgrid$calBP >= timeRange[2]]) 
        if (spdnormalised){ sim[,s] <- (sim[,s]/sum(sim[,s])) } 
        if (!is.na(runm)){ sim[,s] <- runMean(sim[,s], runm, edge="fill") }	
    	}
    }	    

    if (ncores>1)
    {	     
	    sim <- foreach (s = 1:nsim, .combine='cbind', .packages='rcarbon',.options.snow = opts) %dopar% {
		    
		    randomDates <- vector("list",length=ncc)
		    ccurve.tmp <- numeric()

		    if (method=="uncalsample")
		    {
			    for (i in 1:ncc)
			    {       
				    randomDates[[i]] = sample(cragrids[[i]]$CRA,replace=TRUE,size=samplesize[s,i],prob=cragrids[[i]]$PrDens)
				    ccurve.tmp = c(ccurve.tmp,rep(unique.calCurves[i],samplesize[s,i]))
			    }
		    }


		    if (method=="calsample")
		    {
			    for (i in 1:ncc)
			    {       
				    randomDates[[i]] = sample(cragrids[[i]]$CRA,replace=TRUE,size=samplesize[s,i],prob=cragrids[[i]]$Raw)
				    ccurve.tmp = c(ccurve.tmp,rep(unique.calCurves[i],samplesize[s,i]))
			    }
		    }

		    randomSDs <- sample(size=length(unlist(randomDates)), errors, replace=TRUE) 


		    tmp <- calibrate(x=unlist(randomDates),errors=randomSDs, timeRange=timeRange, calCurves=ccurve.tmp, normalised=normalised, ncores=1, verbose=FALSE, calMatrix=TRUE) 

		    simDateMatrix <- tmp$calmatrix
		    aux <- apply(simDateMatrix,1,sum)
		    aux <- (aux/sum(aux)) * sum(predgrid$PrDens[predgrid$calBP <= timeRange[1] & predgrid$calBP >= timeRange[2]])
		    if (spdnormalised){ aux <- (aux/sum(aux)) }
		    if (!is.na(runm)){
			    aux <- runMean(aux, runm, edge="fill")
		    }
		    aux
	    }
    }
    if (verbose){ close(pb) }



    ## rate of change subroutine
    timeSequence = timeRange[1]:timeRange[2]

    foo = function(spd,backsight,timeSequence,changexpr)
    {
	    obs=rep(NA,length(timeSequence))

	    for (i in 1:c(length(obs)-backsight))
	    {
		    d=backsight 	
		    t0 = spd[i]
		    t1 = spd[i+backsight]
		    obs[i+backsight] = eval(changexpr)
		    if (t1==0|t0==0){obs[i+backsight]=NA}
	    }
	    return(obs)
    }

    obs.roc = foo(finalSPD,backsight=backsight,timeSequence=timeSequence,changexpr=changexpr)
    sim.roc = apply(sim,2,foo,backsight=backsight,timeSequence=timeSequence,changexpr=changexpr)
    
    
    ## Envelope, z-scores, global p-value
    
    lo <- apply(sim,1,quantile,prob=0.025,na.rm=TRUE)
    hi <- apply(sim,1,quantile,prob=0.975,na.rm=TRUE)
    lo.roc = apply(sim.roc,1,quantile,prob=0.025,na.rm=TRUE)
    hi.roc = apply(sim.roc,1,quantile,prob=0.975,na.rm=TRUE)

    Zsim <- t(apply(sim,1,scale))
    zLo <- apply(Zsim,1,quantile,prob=0.025,na.rm=TRUE)
    zHi <- apply(Zsim,1,quantile,prob=0.975,na.rm=TRUE)
    Zsim.roc <- t(apply(sim.roc,1,scale))
    zLo.roc <- apply(Zsim.roc,1,quantile,prob=0.025,na.rm=TRUE)
    zHi.roc <- apply(Zsim.roc,1,quantile,prob=0.975,na.rm=TRUE)
    
    Zscore_empirical <- (finalSPD - apply(sim, 1, mean))/apply(sim, 1, sd)
    Zscore_empirical.roc <- (obs.roc - apply(sim.roc, 1, mean))/apply(sim.roc, 1, sd)

    busts <- which(Zscore_empirical< zLo)
    booms <- which(Zscore_empirical> zHi)
    busts2 <- which(finalSPD< lo)
    booms2 <- which(finalSPD> hi)
    busts.roc <- which(Zscore_empirical.roc < zLo.roc)
    booms.roc <- which(Zscore_empirical.roc > zHi.roc)
    busts2.roc <- which(obs.roc < lo.roc)
    booms2.roc <- which(obs.roc > hi.roc)

    observedStatistic <- sum(c(zLo[busts] - Zscore_empirical[busts]),c(Zscore_empirical[booms]-zHi[booms]))
    observedStatistic.roc <- sum(c(zLo.roc[busts.roc] - Zscore_empirical.roc[busts.roc]),c(Zscore_empirical.roc[booms.roc]-zHi.roc[booms.roc]))

    expectedstatistic <- abs(apply(Zsim,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=zLo)) + apply(Zsim,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=zHi)
    expectedstatistic.roc <- abs(apply(Zsim.roc,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=zLo.roc)) + apply(Zsim.roc,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=zHi.roc)

    pvalue <- c(length(expectedstatistic[expectedstatistic > observedStatistic])+1)/c(nsim+1)
    pvalue.roc <- c(length(expectedstatistic.roc[expectedstatistic.roc > observedStatistic.roc])+1)/c(nsim+1)


    # Results
    result <- data.frame(calBP=observed$grid$calBP,PrDens=finalSPD,lo=lo,hi=hi)
    result.roc <- data.frame(calBP=timeSequence,roc=obs.roc,lo.roc=lo.roc,hi.roc=hi.roc) #time is midpoint of transition


    predgrid <- subset(predgrid,calBP<=timeRange[1]&calBP>=timeRange[2])
    if(raw==FALSE){ sim <- NA; sim.roc<-NA }
    res <- list(result=result,result.roc=result.roc, sim=sim, sim.roc=sim.roc, pval=pvalue, pval.roc=pvalue.roc, fit=predgrid, fitobject=fit,nbins=length(unique(bins)),n=nrow(x$metadata),nsim=nsim,backsight=backsight)
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
#' @param backsight A single numeric value defining the distance in time between the focal year and the backsight year for computing the rate of change. Default is 10.
#' @param changexpr An expression for calculating the rate of change in SPD between the focal year and a backsight year. Available input options are t1 (the SPD for the focal year), t0 (the SPD for the backsight year), d (the distance between t0 and t1), and any other standard constants and mathematical operators. A sensible default is provided.
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
#' summary(res)
#' par(mfrow=c(2,1))
#' plot(res,focalm="Netherlands",main="Netherlands")
#' plot(res,focalm="Denmark",main="Denmark")
#' }
#' @import utils
#' @import stats
#' @export

permTest <- function(x, marks, timeRange, backsight=10,changexpr=expression((t1/t0)^(1/d)-1),nsim, bins=NA, runm=NA, datenormalised=FALSE, spdnormalised=FALSE, raw=FALSE, verbose=TRUE){

    if (is.na(bins[1])){ bins <- as.character(1:nrow(x$metadata)) }
    marks <- as.character(marks) 
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

    #compute rates of change
    timeSequence = timeRange[1]:timeRange[2]

    foo = function(spd,backsight,timeSequence,changexpr)
    {
	    obs=rep(NA,length(timeSequence))

	    for (i in 1:c(length(obs)-backsight))
	    {
		    d=backsight 	
		    t0 = spd[i]
		    t1 = spd[i+backsight]
		    obs[i+backsight] = eval(changexpr)
		    if (t1==0|t0==0){obs[i+backsight]=NA}
	    }
	    return(obs)
    }

  observedROC = simulatedROC = vector("list",length=length(observedSPD))

  for (i in 1:length(observedSPD))
  {
    tmp=foo(observedSPD[[i]][,2],backsight=backsight,timeSequence=timeSequence,changexpr=changexpr)
    observedROC[[i]] = data.frame(calBP=timeSequence,roc=tmp)
    names(observedROC)[[i]]=names(observedSPD)[[i]]
    simulatedROC[[i]] = apply(simulatedSPD[[i]],2,foo,backsight=backsight,timeSequence=timeSequence,changexpr=changexpr)
  }


    ## Simulation Envelope
    simulatedCIlist = simulatedCIlist.roc = vector("list",length=length(unique(GroupList)))
    for (d in 1:length(unique(GroupList))){
        simulatedCIlist[[d]] <- cbind(apply(simulatedSPD[[d]],1,quantile,prob=c(0.025),na.rm=TRUE), apply(simulatedSPD[[d]],1,quantile,prob=c(0.975),na.rm=TRUE))
        simulatedCIlist.roc[[d]] <- cbind(apply(simulatedROC[[d]],1,quantile,prob=c(0.025),na.rm=TRUE), apply(simulatedROC[[d]],1,quantile,prob=c(0.975),na.rm=TRUE))
        names(simulatedCIlist) <- unique(GroupList)
        names(simulatedCIlist.roc) <- unique(GroupList)
    }




    ## Compute Global P-value
    pValueList = pValueList.roc = numeric(length=length(simulatedSPD))

    for (a in 1:length(simulatedSPD)){
        zscoreMean <- apply(simulatedSPD[[a]],1,mean)
        zscoreSD <- apply(simulatedSPD[[a]],1,sd)
        zscoreMean.roc <- apply(simulatedROC[[a]],1,mean,na.rm=TRUE)
        zscoreSD.roc <- apply(simulatedROC[[a]],1,sd,na.rm=TRUE)

        tmp.sim <- t(apply(simulatedSPD[[a]],1,function(x){ return((x - mean(x))/sd(x)) }))
        tmp.sim.roc <- t(apply(simulatedROC[[a]],1,function(x){ return((x - mean(x))/sd(x)) }))
#         tmp.sim[is.na(tmp.sim)] <- 0

        tmp.obs <- observedSPD[[a]]
        tmp.obs[,2] <- (tmp.obs[,2] - zscoreMean) / zscoreSD
#         tmp.obs[is.na(tmp.obs[,2]),2] <- 0
        tmp.obs.roc <- observedROC[[a]]
        tmp.obs.roc[,2] <- (tmp.obs.roc[,2] - zscoreMean.roc) / zscoreSD.roc

        tmp.ci <- t(apply(tmp.sim,1, quantile, prob=c(0.025,0.975),na.rm=TRUE))
        tmp.ci.roc <- t(apply(tmp.sim.roc,1, quantile, prob=c(0.025,0.975),na.rm=TRUE))

        expectedstatistic <- abs(apply(tmp.sim,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=tmp.ci[,1])) + apply(tmp.sim,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=tmp.ci[,2])   
        expectedstatistic.roc <- abs(apply(tmp.sim.roc,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=tmp.ci.roc[,1])) + apply(tmp.sim.roc,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=tmp.ci.roc[,2])   
        lower <- tmp.obs[,2] - tmp.ci[,1]
        indexLow <- which(tmp.obs[,2] < tmp.ci[,1])
        higher <- tmp.obs[,2] - tmp.ci[,2]
        indexHi <- which(tmp.obs[,2] > tmp.ci[,2])
        lower.roc <- tmp.obs.roc[,2] - tmp.ci.roc[,1]
        indexLow.roc <- which(tmp.obs.roc[,2] < tmp.ci.roc[,1])
        higher.roc <- tmp.obs.roc[,2] - tmp.ci.roc[,2]
        indexHi.roc <- which(tmp.obs.roc[,2] > tmp.ci.roc[,2])

        observedStatistic <- sum(abs(lower[indexLow]))+sum(higher[indexHi])
        observedStatistic.roc <- sum(abs(lower.roc[indexLow.roc]))+sum(higher.roc[indexHi.roc])

        pValueList[[a]] <- 1
        pValueList.roc[[a]] <- 1

        if (observedStatistic>0){    
            pValueList[[a]] <-c(length(expectedstatistic[expectedstatistic > observedStatistic])+1)/c(nsim+1)
        }
        if (observedStatistic.roc>0){    
            pValueList.roc[[a]] <-c(length(expectedstatistic.roc[expectedstatistic.roc > observedStatistic.roc])+1)/c(nsim+1)
        }

        names(pValueList) <- unique(GroupList)
    } 

    res <- list(observed=observedSPD, envelope=simulatedCIlist, pValueList=pValueList)
    res$observed.roc=observedROC
    res$envelope.roc=simulatedCIlist.roc
    res$pValueList.roc=pValueList.roc

    metadata=vector("list",length=length(unique(GroupList)))
    names(metadata) = unique(GroupList)
	for (k in 1:length(unique(GroupList)))
	{
	        i=unique(GroupList)[k]
		tmp = c(sum(marks%in%i),length(unique(bins[which(marks%in%i)])))
		metadata[[k]]=tmp
	}
    res$nsim = nsim
    res$metadata = metadata    
    res$backsight = backsight
    if (raw)
    {
	    res$raw <- simulatedSPD
	    res$raw.row <- simulatedROC
    }
    class(res) <- "SpdPermTest"
    if (verbose){ print("Done.") }
    return(res)
}


#' @title Spatial Permutation Test of summed probability distributions.
#' @description This function is deprecated. Please use \code{\link{sptest}} instead.
#' @return A \code{spatialTest} class object
#' @name SPpermTest-deprecated
#' @usage SPpermTest(calDates, timeRange, bins, locations, breaks, 
#' spatialweights,rate=expression((t2/t1)^(1/d)-1), nsim=1000, runm=NA,permute="locations",
#' ncores=1,datenormalised=FALSE,verbose=TRUE,raw=FALSE)
#' @seealso \code{\link{rcarbon-deprecated}}
#' @keywords internal
NULL
#' @rdname rcarbon-deprecated
#' @section \code{SPpermTest}:
#' For \code{SPpermTest}, use \code{\link{sptest}}.
#' @export
SPpermTest<-function(calDates, timeRange, bins, locations, breaks, spatialweights, rate=expression((t2/t1)^(1/d)-1),nsim=1000, runm=NA,permute="locations",ncores=1,datenormalised=FALSE,verbose=TRUE,raw=FALSE)
{
.Deprecated("sptest")
 sptest(calDates=calDates, timeRange=timeRange, bins=bins, locations=locations, breaks=breaks, spatialweights=spatialweights,rate=rate, nsim=nsim, runm=runm,permute=permute,ncores=ncores,datenormalised=datenormalised,verbose=verbose,raw=raw)
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
#' @param rate An expression defining how the rate of change is calculated, where \code{t1} is the summed probability for a focal block, \code{t2} is the summed probability for next block, and \code{d} is the duration of the blocks. Default is a geometric growth rate (i.e \code{expression((t2/t1)^(1/d)-1)}).
#' @param nsim The total number of simulations. Default is 1000.
#' @param runm The window size of the moving window average. Must be set to \code{NA} if the rates of change are calculated from the raw SPDs. 
#' @param permute Indicates whether the permutations should be based on the \code{"bins"} or the \code{"locations"}. Default is \code{"locations"}. 
#' @param ncores Number of cores used for for parallel execution. Default is 1.
#' @param datenormalised A logical variable indicating whether the probability mass of each date within \code{timeRange} is equal to 1. Default is FALSE. 
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param raw A logical variable indicating whether permuted sets of geometric growth rates for each location should be returned. Default is FALSE. 
#'
#' @details The function consists of the following seven steps: 1) for each location (e.g. a site) generate a local SPD of radiocarbon dates, weighting the contribution of dates from neighbouring sites using a weight scheme provided by the \code{spatialweights} class object; 2) define temporal slices (using \code{breaks} as break values), then compute the total probability mass within each slice; 3) compute the rate of change between abutting temporal slices by using the formula: \eqn{(SPD_{t}/SPD_{t+1}^{1/\Delta t}-1)}; 4) randomise the location of individual bins or the entire sequence of bins associated with a given location and carry out steps 1 to 3; 5) repeat step 4 \code{nsim} times and generate, for each location, a distribution of growth rates under the null hypothesis (i.e. spatial independence); 6) compare, for each location, the observed growth rate with the distribution under the null hypothesis and compute the p-values; and 7) compute the false-discovery rate for each location.    
#'
#' @return A \code{spatialTest} class object
#'
#' @references
#' Crema, E.R., Bevan, A., Shennan, S. (2017). Spatio-temporal approaches to archaeological radiocarbon dates. Journal of Archaeological Science, 87, 1-9.
#' 
#' @seealso \code{\link{permTest}} for a non-spatial permutation test; \code{\link{plot.spatialTest}} for plotting; \code{\link{spweights}} for computing spatial weights; \code{\link{spd2rc}} for computing geometric growth rates.
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
#' ## Main Analysis (over 2 cores; requires doSnow package) 
#' ## NOTE: the number of simulations should be ideally larger 
#' ## to ensure a better resolution of the p/q-values.
#' res.locations=sptest(calDates,timeRange=timeRange,bins=bins,locations=locations,
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
#' plot(res.locations,index=4,add=TRUE,legend=TRUE,option="raw",breakRange=c(-0.005,0.005)) 
#' par(mar=c(0.1,0.1,0,0.5))
#' plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange) 
#' plot(res.locations,index=4,add=TRUE,legend=TRUE,option="test")  
#' }
#' @import utils
#' @import stats
#' @import doSNOW 
#' @import snow
#' @import foreach 
#' @import iterators 
#' @import sp
#' @export
 

sptest<-function(calDates, timeRange, bins, locations, breaks, spatialweights, rate=expression((t2/t1)^(1/d)-1),nsim=1000, runm=NA,permute="locations",ncores=1,datenormalised=FALSE,verbose=TRUE,raw=FALSE)
{

  #######################
  #### Load Packages ####
  #######################
  
  if (ncores>1&!requireNamespace('doSNOW',quietly = TRUE)){
    warning("the doSNOW package is required for multi-core processing; ncores has been set to 1")
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

	if (length(unique(abs(diff(breaks))))!=1)
	{
		stop("Unequal break intervals is not supported")
	}


	if (ncores>1&raw==TRUE)
	{
		warning("raw==TRUE available only for ncores=1")
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
			tmp <- lapply(tmp,FUN=function(x) {if(sum(x)!=0){return(x/sum(x))}else{return(x)}})
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

	rocaObs=t(apply(obsGridVal,1,function(x,d,rate){
				L=length(x)
				res=numeric(length=L-1)
				for (i in 1:c(L-1))
				{
					t2 = x[i+1]
					t1 = x[i]
					res[i] = eval(rate)
# 					res[i]=(x[i+1]/x[i])^(1/d)-1
					# If no spd for both period the rate is NA
					if (x[i+1]==0|x[i]==0)
					{
						res[i]=NA	
					}
				}
				return(res)},
				d=abs(breaks[2]-breaks[1]),
				rate=rate))
	#if single transition transpose matrix:
	if (nBreaks==2){rocaObs=t(rocaObs)} 

	##############################
	### Permutation Subroutine ###
	############################## 
	opts = NULL
	if (ncores>1)
	{
	  cl <- snow::makeCluster(ncores)
	  registerDoSNOW(cl)
	  if (verbose)
	  {
		print(paste("Running permutation test in parallel on ",getDoParWorkers()," workers...",sep=""))
		pb <- txtProgressBar(min=0, max=length(x), style=3)
		progress <- function(n) setTxtProgressBar(pb, n)
		opts <- list(progress = progress)
	  }
		sumcombine<-function(a,b)
		{
			list(a[[1]]+b[[1]],a[[2]]+b[[2]],a[[3]]+b[[3]])
		}
		resultHiLoEq<-foreach (x=1:nsim,.combine= sumcombine,.options.snow = opts) %dopar% {

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

			rocaSim=t(apply(simGridVal,1,function(x,d,rate){
						L=length(x)
						res=numeric(length=L-1)
						for (i in 1:c(L-1))
						{
							t2 = x[i+1]
							t1 = x[i]
							res[i] = eval(rate)
# 							res[i]=(x[i+1]/x[i])^(1/d)-1
							if (x[i+1]==0|x[i]==0)
							{
								res[i]=NA
							}
						}
						return(res)},
						d=abs(breaks[2]-breaks[1]),
						rate=rate))

			lo=rocaObs<rocaSim	    
			hi=rocaObs>rocaSim
			eq=rocaObs==rocaSim

			return(list(hi,lo,eq))
		}
		snow::stopCluster(cl)

		hi=resultHiLoEq[[1]]
		lo=resultHiLoEq[[2]]
		eq=resultHiLoEq[[3]]

	} else {

		hi=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)
		lo=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)
		eq=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)


		if(verbose){
			print("Permutation test...")
			flush.console()
		}

		if(raw)
		{
			rocaSimAll = array(NA,dim=c(nsim,nrow(spatialweights$w),nBreaks-1))
		}
		
		if(verbose){pb <- txtProgressBar(min = 1, max = nsim, style=3)}

		for (s in 1:nsim)
		{
			if(verbose){setTxtProgressBar(pb, s)}
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
			rocaSim=t(apply(simGridVal,1,function(x,d,rate){
						L=length(x)
						res=numeric(length=L-1)
						for (i in 1:c(L-1))
						{
							t2 = x[i+1]
							t1 = x[i]
							res[i] = eval(rate)
# 							res[i]=(x[i+1]/x[i])^(1/d)-1
							if (x[i+1]==0|x[i]==0)
							{
								res[i]=NA	
							}
						}
						return(res)},
						d=abs(breaks[2]-breaks[1]),
						rate=rate))

			# if single transition transpose matrix:
			if (nBreaks==2) {rocaSim=t(rocaSim)}

			hi=hi+(rocaObs>rocaSim)
			lo=lo+(rocaObs<rocaSim)
			eq=eq+(rocaObs==rocaSim)
			if(raw){rocaSimAll[s,,]=rocaSim}
		}
		if(verbose){close(pb)}
	}


	############################
	### Compute Significance ###
	############################ 
	
	pvalHi=(lo+eq+1)/c(nsim+1)
	pvalLo=(hi+eq+1)/c(nsim+1)
	pval=pvalHi
	pval[which(pvalHi>pvalLo)]=pvalLo[which(pvalHi>pvalLo)]
	pval=pval*2


	if (max(pval,na.rm=TRUE)>1)
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
#' @param plot if set to TRUE the function plots the location of p1 and p2 on the SPD. Default is FALSE.
#'
#' @details The function compares observed differences in the summed probability values associated with two points in time against a distribution of expected values under the null hypothesis defined with the \code{\link{modelTest}} function. The two points can be specified manually (assigning BP dates to the arguments \code{p1} and \code{p2}) or interactively (clicking on a SPD plot). Note that \code{\link{modelTest}} should be executed setting the argument \code{raw} to \code{TRUE} (default is \code{FALSE}).   
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
#' caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE)
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


p2pTest <- function(x,p1=NA,p2=NA,interactive=TRUE,plot=FALSE)
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
		if (plot)
		{	    
			plot(x)	    
			points(p1,p1.y,pch=20)    
			points(p2,p2.y,pch=20)
			lines(x$result$calBP[index1:index2],x$result$PrDens[index1:index2],lwd=2)
		}
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




#' @title Summarise a \code{SpdModelTest} class object
#'
#' @description \code{summary} method for class "\code{SpdModelTest}" 
#'
#' @param object A \code{SpdModelTest} class object produced using the \code{link{modelTest}} function.
#' @param type Specifies whether the summary should be based on SPD ('spd') or associated rates of change ('roc'). Default is 'spd'.
#' @param ... Ignored
#'
#' @details The summary function returns metadata (number of radiocarbon dates, bins, and simulations), the p-value of the global significance test, and the chronological interval of local positive and negative deviations from the simulation envelope.
#' @seealso \code{\link{modelTest}}.
#' @import utils
#' @import stats
#' @method summary SpdModelTest
#' @export

summary.SpdModelTest<-function(object,type='spd',...) {

	if (!type%in%c('spd','roc'))
	{
	 stop("The argument 'type' should be either 'spd' or 'roc'.")
	}

	cat("'modelTest()' function summary:\n")
	cat("\n") 
	cat(paste("Number of radiocarbon dates: ",object$n,"\n",sep=""))
	cat(paste("Number of bins: ",object$nbins,"\n",sep=""))
	if(type=='roc'){cat(paste("Backsight size: ",object$backsight,"\n",sep=""))}
	cat("\n")
	cat(paste("Statistical Significance computed using ",object$nsim," simulations. \n",sep=""))
	if(type=='spd'){cat(paste("Global p-value: ",round(object$pval,5),".\n",sep=""))}
	if(type=='roc'){cat(paste("Global p-value (rate of change): ",round(object$pval.roc,5),".\n",sep=""))}
	cat("\n")

	if (type=='spd')
	{
	obs <- object$result[,1:2]
	envelope <- object$result[,3:4]
	booms <- which(obs$PrDens>envelope[,2])
	busts <- which(obs$PrDens<envelope[,1])
	}

	if (type=='roc')
	{
	obs <- object$result.roc[,1:2]
	envelope <- object$result.roc[,3:4]
	booms <- which(obs$roc>envelope[,2])
	busts <- which(obs$roc<envelope[,1])
	}


	if (length(booms)>0)
	{
		cat(paste("Signficant positive local deviations at:\n"))
		i=1
		while (i < length(obs$calBP))
		{
			if(!is.na(obs[i,2]))
			{	
				if(obs[i,2]>envelope[i,2])
				{
					ss=obs$calBP[i]
					while(obs[i,2]>envelope[i,2])
					{
						ee=obs$calBP[i]
						i=i+1
						if (i>length(obs$calBP))
						{
							i = length(obs$calBP)
							ee=obs$calBP[i]
							break()
						}
					}
					if (ss!=ee)	
						cat(paste(ss,"~",ee," BP \n",sep=""))
					if (ss==ee)	
						cat(paste(ss," BP \n",sep=""))
				}
		}
			i = i+1
		}
	}


	if (length(booms)==0)
	{
		cat(paste("No significant positive local deviations"))
	}

	cat("\n")

	if (length(busts)>0)
	{
		cat(paste("Significant negative local deviations at:\n"))


		i=1
		while (i < length(obs$calBP))
		{
			if(!is.na(obs[i,2]))
			{	
				if(obs[i,2]<envelope[i,1])
				{
					ss=obs$calBP[i]
					while(obs[i,2]<envelope[i,1])
					{
						ee=obs$calBP[i]
						i=i+1
						if (i>length(obs$calBP))
						{
							i = length(obs$calBP)
							ee=obs$calBP[i]
							break()
						}
					}
					if (ss!=ee)	
						cat(paste(ss,"~",ee," BP \n",sep=""))
					if (ss==ee)	
						cat(paste(ss," BP \n",sep=""))
				}
			}
			i = i+1
		}
	}


	if (length(busts)==0)
	{
		cat(paste("No significant positive local deviations"))
	}
}


#' @title Summarise a \code{SpdPermTest} class object
#'
#' @description \code{summary} method for class "\code{SpdPermTest}" 
#'
#' @param object A \code{SpdPermTest} class object produced using the \code{link{permTest}} function.
#' @param type Specifies whether the summary should be based on SPD ('spd') or associated rates of change ('roc'). Default is 'spd'.
#' @param ... Ignored
#'
#' @details The summary function returns metadata (number of radiocarbon dates, bins, and simulations), the p-value of the global significance test, and the chronological interval of local positive and negative deviations from the simulation envelope for each set.
#' @seealso  \code{\link{permTest}}.
#' @import utils
#' @import stats
#' @method summary SpdPermTest
#' @export

summary.SpdPermTest<-function(object,type='spd',...) {

	if (!type%in%c('spd','roc'))
	{
	 stop("The argument 'type' should be either 'spd' or 'roc'.")
	}

	cat("'permTest()' function summary:\n")
	cat("\n") 
	cat(paste("Number of sets: ",length(object$observed),"\n",sep=""))
	if(type=='roc'){cat(paste("Backsight size: ",object$backsight,"\n",sep=""))}
	cat(paste("Statistical Significance computed using ",object$nsim," simulations. \n",sep=""))
	cat("\n")
	
	for (i in 1:length(object$observed))
	{

		if (type=='spd'){cat(paste("--- ",names(object$observed)[i]," ---\n",sep=""))}
		if (type=='roc'){cat(paste("--- ",names(object$observed.roc)[i]," ---\n",sep=""))}
		cat(paste("Number of radiocarbon dates:",object$metadata[[i]][1],"\n",sep=""))
		cat(paste("Number of bins:",object$metadata[[i]][2],"\n",sep=""))
		cat("\n")
		if(type=='spd'){
			cat(paste("Global p-value: ",round(object$pValueList[i],5),"\n",sep=""))
			obs <- object$observed[[i]]
			envelope <-object$envelope[[i]]
			booms <- which(obs[,2]>envelope[,2])
			busts <- which(obs[,2]<envelope[,1])
		}
		if(type=='roc'){
			cat(paste("Global p-value (rate of change): ",round(object$pValueList.roc[i],5),"\n",sep=""))
			obs <- object$observed.roc[[i]]
			envelope <-object$envelope.roc[[i]]
			booms <- which(obs[,2]>envelope[,2])
			busts <- which(obs[,2]<envelope[,1])
		}
		cat("\n")

		if (length(booms)>0)
		{
			cat(paste("Significant positive local deviations at:\n"))
			i=1
			while (i < length(obs[,1]))
			{	
				if(!is.na(obs[i,2]))
				{
				if(obs[i,2]>envelope[i,2])
				{
					ss=obs[i,1]
					while(obs[i,2]>envelope[i,2])
					{
						ee=obs[i,1]
						i=i+1
						if (i>length(obs[,1]))
						{
							i = length(obs[,1])
							ee=obs[i,1]
							break()
						}
					}
					if (ss!=ee)	
						cat(paste(ss,"~",ee," BP \n",sep=""))
					if (ss==ee)	
						cat(paste(ss," BP \n",sep=""))
				}
				}
				i = i+1
			}
		}


		if (length(booms)==0)
		{
			cat(paste("No significant positive local deviations"))
		}

		cat("\n")

		if (length(busts)>0)
		{
			cat(paste("Significant negative local deviations at:\n"))


			i=1
			while (i < length(obs[,1]))
			{
				if(!is.na(obs[i,2]))
				{
				if(obs[i,2]<envelope[i,1])
				{
					ss=obs[i,1]
					while(obs[i,2]<envelope[i,1])
					{
						ee=obs[i,1]
						i=i+1
						if (i>length(obs[,1]))
						{
							i = length(obs[,1])
							ee=obs[i,1]
							break()
						}
					}
					if (ss!=ee)	
						cat(paste(ss,"~",ee," BP \n",sep=""))
					if (ss==ee)	
						cat(paste(ss," BP \n",sep=""))
				}
				}
				i = i+1
			}
		}


		if (length(busts)==0)
		{
			cat(paste("No significant positive local deviations"))
		}

		cat("\n")
	}
}




