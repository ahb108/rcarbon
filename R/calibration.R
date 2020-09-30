#' @title Calibrate radiocarbon dates
#'
#' @description Function for calibrating radiocarbon dates or uncalibrated SPDs.
#'
#' @param x A vector of uncalibrated radiocarbon ages or a UncalGrid class object.
#' @param errors A vector of standard deviations corresponding to each estimated radiocarbon age.
#' @param ids An optional vector of IDs for each date.
#' @param dateDetails An optional vector of details for each date which will be returned in the output metadata. 
#' @param calCurves Either a character string naming a calibration curve already provided with the rcarbon package (currently 'intcal20','intcal13','intcal13nhpine16','shcal20','shcal13','shcal13shkauri16',''marine13','marine20' and 'normal'(i.e. no calibration) are possible; default is 'intcal20') or a custom calibration curve as three-column matrix or data.frame (calibrated year BP, uncalibrated age bp, standard deviation). Different existing curves can be specified per dated sample, but only one custom curve can be provided for all dates.
#' @param resOffsets A vector of offset values for any marine reservoir effect (default is no offset).
#' @param resErrors A vector of offset value errors for any marine reservoir effect (default is no offset).
#' @param timeRange Earliest and latest data to calibrate for, in calendar years. Posterior probabilities beyond this range will be excluded (the default is sensible in most cases).
#' @param normalised A logical variable indicating whether the calibration should be normalised or not. Default is TRUE.
#' @param F14C A logical variable indicating whether calibration should be carried out in F14C space or not. Default is FALSE.
#' @param eps Cut-off value for density calculation. Default is 1e-5.
#' @param calMatrix a logical variable indicating whether the age grid should be limited to probabilities higher than \code{eps}
#' @param ncores Number of cores/workers used for parallel execution. Default is 1 (>1 requires doSnow package).
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param ... ignored
#'
#' @details This function computes one or more calibrated radiocarbon ages using the method described in Bronk Ramsey 2008 (see also  Parnell 2017). It is possible to specify different calibration curves or reservoir offsets individually for each date, and control whether the resulting calibrated distribution is normalised to 1 under-the-curve or not. Calculations can also be executed in parallel to reduce computing time. The function was modified from the \code{BchronCalibrate} function in the \code{Bchron} package developed by A.Parnell (see references below).
#'
#' @return An object of class CalDates with the following elements:
#' \itemize{
#' \item{\code{metadata}} {A data.frame containing relevant information regarding each radiocarbon date and the parameter used in the calibration process.}
#' \item{\code{grids}} {A list of calGrid class objects, containing the posterior probabilities for each calendar year. The most memory-efficient way to store calibrated dates, as only years with non-zero probability are stored, but aggregation methods such as \code{spd()} may then take longer to extract and combine multiple dates. NA when the parameter calMatrix is set to TRUE.} 
#' \item{\code{calMatrix}} {A matrix of probability values, one row per calendar year in timeRange and one column per date. By storing all possible years, not just those with non-zero probability, this approach takes more memory, but speeds up spd() and is suggested whenever the latter is to be used. NA when the parameter calMatrix is set to FALSE.}  
#' }
#'
#' @references 
#' Bronk Ramsey, C. 2008. Radiocarbon dating: revolutions in understanding, \emph{Archaeometry} 50.2: 249-75. DOI: https://doi.org/10.1111/j.1475-4754.2008.00394.x \cr
#' Parnell, A. 2017. Bchron: Radiocarbon Dating, Age-Depth Modelling, Relative Sea Level Rate Estimation, and Non-Parametric Phase Modelling, R package: https://CRAN.R-project.org/package=Bchron
#' @examples
#' x1 <- calibrate(x=4000, errors=30)
#' plot(x1)
#' summary(x1)
#' # Example with a Marine Date, using a DeltaR of 300 and a DeltaR error of 30
#' x2 <- calibrate(x=4000, errors=30, calCurves='marine20', resOffsets=300, resErrors=30)
#' plot(x2)
#' @import stats 
#' @import utils 
#' @import doSNOW 
#' @import snow
#' @import foreach 
#' @import iterators 
#' @export

calibrate <- function (x, ...) {
  UseMethod("calibrate")
}

#' @rdname calibrate
#' @export

calibrate.default <- function(x, errors, ids=NA, dateDetails=NA, calCurves='intcal20', resOffsets=0 , resErrors=0, timeRange=c(55000,0), normalised=TRUE, F14C=FALSE, calMatrix=FALSE, eps=1e-5, ncores=1, verbose=TRUE, ...){

    if (ncores>1&!requireNamespace("doSNOW", quietly=TRUE)){	
	warning("the doSNOW package is required for multi-core processing; ncores has been set to 1")
	ncores=1
    }	
    	
    # age and error checks
    if (length(x) != length(errors)){
        stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
    }
    if (!is.na(ids[1]) & (length(x) != length(ids))){
        stop("Ages and errors (and ids/details/offsets if provided) must be the same length.")
    }
    if (any(is.na(x))|any(is.na(errors))){
        stop("Ages or errors contain NAs")
    }
    if (is.null(calCurves)|anyNA(calCurves)){
	stop("calCurves is NULL or contain NAs")
    }
  
    if (!any(class(calCurves) %in% c("matrix","data.frame"))& any(calCurves %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")) & (max(timeRange) > 50000 | min(timeRange)<0))
    {
      timeRange=c(50000,0)
      warning("timeRange value not supported with the selected curve(s); calibrating using timeRange=c(50000,0)")
    }
  
    if (any(class(calCurves) %in% c("matrix","data.frame")))
    {
      if (max(calCurves[,1])<timeRange[1] | min(calCurves[,1]) < timeRange[2] )
      {
        timeRange=c(max(calCurves[,1]),min(calCurves[,1]))
        warning(paste0("timeRange value not supported with the selected curve(s); calibrating using timeRange=c(",timeRange[1],",",timeRange[1],")"))
      }
    }

    if (F14C==TRUE&normalised==FALSE)
    {
      normalised=TRUE
      warning("normalised cannot be FALSE when F14C is set to TRUE, calibrating with normalised=TRUE")
    }
    # calCurve checks and set-up
    if (any(class(calCurves) %in% c("matrix","data.frame"))){
        cctmp <- as.matrix(calCurves)
        if (ncol(cctmp)!=3 | !all(sapply(cctmp,is.numeric))){
            stop("The custom calibration curve must have just three numeric columns.")
        } else {
            colnames(cctmp) <- c("CALBP","C14BP","Error")
            if (max(cctmp[,2]) < max(x) | min(cctmp[,2]) > min(x)){
                stop("The custom calibration curve does not cover the input age range.")
            }
            #cclist <- vector(mode="list", length=1)
            #cclist[[1]] <- cctmp
            #names(cclist) <- "custom"
            calCurves <- rep("custom",length(x))
            cclist2 <- vector(mode="list", length=1)
            names(cclist2) <- "custom"
            calBPrange = seq(max(cctmp[,1]),min(cctmp[,1]),-1)
            if (F14C)
            {
              F14 <- exp(cctmp[,2]/-8033) 
              F14Error <-  F14*cctmp[,3]/8033 
              calf14 <- approx(cctmp[,1], F14, xout=calBPrange)$y 
              calf14error <-  approx(cctmp[,1], F14Error, xout=calBPrange)$y 
              cclist2[[1]] <- list(calf14=calf14,calf14error=calf14error,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))
            } else {
              cclist2[[1]] = list(mu=stats::approx(cctmp[,1], cctmp[,2], xout = calBPrange)$y,tau2 = stats::approx(cctmp[,1], cctmp[,3], xout = calBPrange)$y^2,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))
            } 
        } 
    } else if (!all(calCurves %in% c("intcal13","intcal20","shcal13","shcal20","marine13","marine20","intcal13nhpine16","shcal13shkauri16","normal"))){
        stop("calCurves must be a character vector specifying one or more known curves or a custom three-column matrix/data.frame (see ?calibrate.default).")
  } else {
    tmp <- unique(calCurves)
    if (length(calCurves)==1){ calCurves <- rep(calCurves,length(x)) }
    #cclist <- vector(mode="list", length=length(tmp))
    #names(cclist) <- tmp
    cclist2 <- vector(mode="list", length=length(tmp))
    names(cclist2) <- tmp
    for (a in 1:length(tmp)){
      cctmp <- read_cal_curve_from_file(tmp[a])
      #cclist[[tmp[a]]] <- cctmp
      calBPrange = seq(max(cctmp[,1]),min(cctmp[,1]),-1)
      if (F14C)
      {
        F14 <- exp(cctmp[,2]/-8033) 
        F14Error <-  F14*cctmp[,3]/8033 
        calf14 <- approx(cctmp[,1], F14, xout=calBPrange)$y 
        calf14error <-  approx(cctmp[,1], F14Error, xout=calBPrange)$y 
        cclist2[[tmp[a]]] <- list(calf14=calf14,calf14error=calf14error,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))
      } else {
      cclist2[[tmp[a]]] = list(mu=stats::approx(cctmp[,1], cctmp[,2], xout = calBPrange)$y,tau2 = stats::approx(cctmp[,1], cctmp[,3], xout = calBPrange)$y^2,calBPrange=calBPrange,calBPindex=which(calBPrange<=timeRange[1]&calBPrange>=timeRange[2]))}
    }
  }
  ## container and reporting set-up
  reslist <- vector(mode="list", length=2)
  sublist <- vector(mode="list", length=length(x))
  if (calMatrix){
    calmBP <- seq(timeRange[1],timeRange[2],-1)
    calmat <- matrix(ncol=length(x), nrow=length(calmBP))
    rownames(calmat) <- calmBP
    calmat[] <- 0
  }
  if (is.na(ids[1])){
    ids <- as.character(1:length(x))
  } else {
    ids <- as.character(ids)
    if (any(duplicated(ids))){ stop("The values in the ids argument must be unique or left as defaults.") }
  }
  if (length(resOffsets)==1){ resOffsets <- rep(resOffsets,length(x)) }
  if (length(resErrors)==1){ resErrors <- rep(resErrors,length(x)) }
  names(sublist) <- ids
  names(reslist) <- c("metadata","grids")
  opts = NULL
  if (verbose){
    print("Calibrating radiocarbon ages...")
    flush.console()
    if (ncores>1){ print(paste("Running in parallel (standard calibration only) on ",getDoParWorkers()," workers...",sep=""))}
    pb <- txtProgressBar(min=0, max=length(x), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }
  # calibration
  `%dofun%` <- `%do%`
  if (ncores>1){
    # parallellised
    cl <- snow::makeCluster(ncores)
    registerDoSNOW(cl)
    `%dofun%` <- `%dopar%`
  }
  b <- NULL # Added to solve No Visible Binding for Global Variable NOTE
  age = x - resOffsets
  error = sqrt(errors^2 + resErrors^2)
  
  sublist <- foreach (b=1:length(x),.options.snow = opts) %dofun% {
    if (verbose & ncores==1) {setTxtProgressBar(pb, b)}
    #age <- x[b] - resOffsets[b]
    #error <- sqrt(errors[b]^2 + resErrors[b]^2)
    if (F14C==FALSE)
    {  
      dens <- BP14C_calibration(age[b], error[b], cclist2[[calCurves[b]]]$mu, cclist2[[calCurves[b]]]$tau2, eps)
    }
    if (F14C==TRUE)
    {
      dens <- F14C_calibration(age[b], error[b], cclist2[[calCurves[b]]]$calf14, cclist2[[calCurves[b]]]$calf14error, eps)
    }
    if (normalised){
      dens <- normalise_densities(dens, eps)
    }
    calBPrange = cclist2[[calCurves[b]]]$calBPrange
    calBP = calBPrange[cclist2[[calCurves[b]]]$calBPindex]
    PrDens = dens[cclist2[[calCurves[b]]]$calBPindex]
    if (anyNA(PrDens)){stop("One or more dates are outside the calibration range")}
    res = list(calBP = calBP[PrDens > 0],PrDens = PrDens[PrDens > 0])
    return(res)
  }
  if (ncores>1){
    snow::stopCluster(cl)
  }
  names(sublist) <- ids
  if (calMatrix){
    for (a in 1:length(sublist)){
      calmat[as.character(sublist[[a]]$calBP),a] <- sublist[[a]]$PrDens
    }
  }
  ## clean-up and results
  if (length(x)>1 & verbose){ close(pb) }
  df <- data.frame(DateID=ids, CRA=x, Error=errors, Details=dateDetails, CalCurve=calCurves,ResOffsets=resOffsets, ResErrors=resErrors, StartBP=timeRange[1], EndBP=timeRange[2], Normalised=normalised, F14C=F14C, CalEPS=eps, stringsAsFactors=FALSE)
  reslist[["metadata"]] <- df
  if (calMatrix){
    reslist[["grids"]] <- NA
    reslist[["calmatrix"]] <- calmat
  } else {
    reslist[["grids"]] <- lapply(sublist,function(x){tmp=data.frame(calBP=x[[1]],PrDens=x[[2]]);class(tmp)=append("calGrid",class(tmp));return(tmp)})
    reslist[["calmatrix"]] <- NA
  }
  class(reslist) <- c("CalDates",class(reslist))
  if (verbose){ print("Done.") }
  return(reslist)
}


#' @title Calibrate a UncalGrid class object
#'
#' @description Function for calibrating a SPD generated by summing uncalibrated dates.
#'
#' @param type Currently, two options are available. With "fast", a look-up is conducted for each calendar year that picks up the amplitude of the nearest CRA age. With "full", each CRA is calibrated individually and then summed (much slower).
#' @param datenormalised Controls for calibrated dates with probability mass outside the timerange of analysis. If set to TRUE the total probability mass within the time-span of analysis is normalised to sum to unity. Should be set to FALSE when the parameter \code{normalised} in \code{\link{calibrate}} is set to FALSE. Default is FALSE. 
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity. 
#' @rdname calibrate
#' @export

calibrate.UncalGrid <- function(x, errors=0, calCurves='intcal20', timeRange=c(55000,0), eps=1e-5, type="fast", datenormalised=FALSE, spdnormalised=FALSE, verbose=TRUE,...){

  if (length(errors)==1){
    errors <- rep(errors,length(x$CRA))
  }
  if (any(calCurves %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")) & (max(timeRange) > 50000 | min(timeRange)<0))
  {
    timeRange=c(50000,0)
    warning("timeRange value not supported with the selected curve(s); calibrating using timeRange=c(50000,0)")
  }
    if (class(calCurves) %in% c("matrix","data.frame")){
        calcurve <- as.matrix(calCurves)
        if (ncol(calcurve)!=3 | !all(sapply(calcurve,is.numeric))){
            stop("The custom calibration curve must have just three numeric columns.")
    } else {
      colnames(calcurve) <- c("CALBP","C14BP","Error")
    }
  } else if (any(class(calCurves)=="character")){
    calcurve <- read_cal_curve_from_file(calCurves)
  } else {
    stop("calCurves must be a character vector specifying a known curve or a custom three-column matrix/data.frame (see ?calibrate.default).")
  }
  if (type=="full"){
    caleach <- calibrate(x=x$CRA, errors=errors, normalised=datenormalised, calMatrix=TRUE, ...)
    for (d in 1:ncol(caleach$calmatrix)){
      caleach$calmatrix[,d] <- caleach$calmatrix[,d] *  x$PrDens[d]
    }
    res <- data.frame(calBP=timeRange[1]:timeRange[2], PrDens=rowSums(caleach$calmatrix))
  } else if (type=="fast"){
    if (datenormalised){
      warning('Cannot normalise dates using fast method, so leaving unnormalised.')
    }
    if (verbose){ print("Calibrating...") }
    CRAdates <- data.frame(approx(calcurve[,1:2], xout=seq(max(calcurve[,1]),min(calcurve[,1]),-1)))
    names(CRAdates) <- c("calBP","CRA")
    CRAdates$CRA <- round(CRAdates$CRA,0)
    res <- merge(CRAdates, x, by="CRA",all.x=TRUE, sort=FALSE)
    res <- res[with(res, order(-calBP)), c("calBP","PrDens")]
    res$PrDens[is.na(res$PrDens)] <- 0
  } else {
    stop("Type must be 'full' or 'fast'.")
  }
  res <- res[which(res$calBP<=timeRange[1] & res$calBP>=timeRange[2]),]
  if (spdnormalised){
    res[res$PrDens < eps,"PrDens"] <- 0
    res$PrDens <- res$PrDens/sum(res$PrDens)
  } else {
    res[res$PrDens < eps,"PrDens"] <- 0
  }
  res <- res[res$PrDens > 0,]
  class(res) <- c("CalGrid", class(res))   
  if (verbose){ print("Done.") }
  return(res)
}



#' @title Uncalibrate (back-calibrate) a calibrated radiocarbon date (or summed probability distribution).
#'
#' @description Function for uncalibrating one or more radiocarbon dates.
#'
#' @param x Either a vector of calibrated radiocarbon ages or an object of class CalGrid.
#' @param CRAerrors A vector of standard deviations corresponding to each estimated radiocarbon age (ignored if x is a CalGrid object).
#' @param roundyear An optional vector of IDs for each date (ignored if x is a CalGrid object).
#' @param  calCurves A string naming a calibration curve already provided with the rcarbon package (currently 'intcal20','intcal13','intcal13nhpine16','shcal20','shcal13','shcal13shkauri16',''marine13','marine20' and 'normal') or a custom curve provided as matrix/data.frame in three columns ("CALBP","C14BP","Error"). The default is the 'intcal20' curve and only one curve can currently be specified for all dates. 
#' @param  eps Cut-off value for density calculation (for CalGrid objects only).
#' @param  compact A logical variable indicating whether only uncalibrated ages with non-zero probabilities should be returned (for CalGrid objects only).
#' @param  verbose A logical variable indicating whether extra information on progress should be reported (for CalGrid objects only).
#' @param ... ignored

#' @details This function takes one or more calibrated calendars and looks-up the corresponding uncalibrated age, error of the stated calibration curve at that point. It also provides a randomised estimate of the uncalibrate age based on the curve error (and optionally also a hypothetical measurement error.
#'
#' @return A data.frame with specifying the original data, the uncalibrated age without the calibration curve error (ccCRA), the calibration curve error at this point in the curve (ccError), a randomised uncalibrated age (rCRA) given both the stated ccError and any further hypothesised instrumental error provided by the CRAerrors argument (rError). 
#'
#' @examples
#' # Uncalibrate two calendar dates
#' uncalibrate(c(3050,2950))
#' @import stats 
#' @import utils 
#' @export

uncalibrate <- function (x, ...) {
  UseMethod("uncalibrate", x)
}

#' @rdname uncalibrate
#' @export

uncalibrate.default <- function(x, CRAerrors=0, roundyear=TRUE, calCurves='intcal20', ...){
    
    if (length(CRAerrors)==1){ CRAerrors <- rep(CRAerrors,length(x)) } 
    ## calCurve checks and set-up
    if (class(calCurves) %in% c("matrix","data.frame")){
        calcurve <- as.matrix(calCurves)
        if (ncol(calcurve)!=3 | !all(sapply(calcurve,is.numeric))){
            stop("The custom calibration curve must have just three numeric columns.")
    } else {
      colnames(calcurve) <- c("CALBP","C14BP","Error")
      if (max(calcurve[,1]) < max(x) | min(calcurve[,1]) > min(x)){
        stop("The custom calibration curve does not cover the input age range.")
      }
    }
  } else if (!all(calCurves %in% c("intcal13","intcal20","shcal13","shcal20","marine13","marine20","intcal13nhpine16","shcal13shkauri16"))){
    stop("calCurves must be a character vector specifying one or more known curves or a custom three-column matrix/data.frame (see ?calibrate.default).")
  } else {
    calcurve <- read_cal_curve_from_file(calCurves)
  }
  dates <- data.frame(approx(calcurve, xout=x))
  colnames(dates) <- c("calBP", "ccCRA")
  calcurve.error <- approx(calcurve[,c(1,3)], xout=dates$calBP)$y
  dates$ccError <- calcurve.error
  #     dates$rCRA <- rnorm(nrow(dates), mean=dates$ccCRA, sd=dates$ccError)
  dates$rCRA <- rnorm(nrow(dates), mean=dates$ccCRA, sd=sqrt(dates$ccError^2+CRAerrors^2))
  dates$rError <- CRAerrors
  if (roundyear){ dates$rCRA <- round(dates$rCRA) }
  return(dates)
}

#' @rdname uncalibrate
#' @export

uncalibrate.CalGrid <- function(x, calCurves='intcal20', eps=1e-5, compact=TRUE, verbose=TRUE, ...){

    if (verbose){ print("Uncalibrating...") }
    names(x) <- c("calBP","PrDens")
    ## calCurve checks and set-up
    if (class(calCurves) %in% c("matrix","data.frame")){
        calcurve <- as.matrix(calCurves)
        if (ncol(calcurve)!=3 | !all(sapply(calcurve,is.numeric))){
            stop("The custom calibration curve must have just three numeric columns.")
    } else {
      colnames(calcurve) <- c("CALBP","C14BP","Error")
      if (max(calcurve[,1]) < max(x$calBP) | min(calcurve[,1]) > min(x$calBP)){
        stop("The custom calibration curve does not cover the input age range.")
      }
    }
  } else if (!all(calCurves %in% c("intcal13","intcal20","shcal13","shcal20","marine13","marine20","intcal13nhpine16","shcal13shkauri16"))){
    stop("calCurves must be a character vector specifying one or more known curves or a custom three-column matrix/data.frame (see ?calibrate.default).")
  } else {
    calcurve <- read_cal_curve_from_file(calCurves)
  }
  mycras <- uncalibrate(x$calBP,calCurves=calCurves)
  res <- data.frame(CRA=max(calcurve[,2]):min(calcurve[,2]), PrDens=0)
  
  h = x$PrDens/sum(x$PrDens)
  mu = mycras$ccCRA
  s = mycras$ccError
  k = res$CRA
  
  res$Raw=unlist(sapply(k,function(x,mu,s,h){return(sum(dnorm(x,mu,s)*h))},h=h,mu=mu,s=s))
  res$Base=unlist(sapply(k,function(x,mu,s){return(sum(dnorm(x,mu,s)))},mu=mu,s=s))
  
  res$Raw=res$Raw/sum(res$Raw)
  
  res$Raw[res$Raw < eps] <- 0
  res$PrDens[res$Base>0] <- res$Raw[res$Base>0] / res$Base[res$Base>0]
  if (compact){ res <- res[res$PrDens > 0,] }
  res$PrDens=res$PrDens/sum(res$PrDens)
  class(res) <- c("UncalGrid", class(res)) 
  if (verbose){ print("Done.") }
  return(res)
}


#' @title Convert data to class CalGrid. 
#'
#' @description Tries to coerce any two-column matrix or data.frame to a calibrated probability distribution (an object of class "CalGrid") for use by the rcarbon package. 
#' 
#' @param x A two-column \code{matrix} or \code{data.frame} class object.
#'
#' @return A CalGrid class object of probabilities or summed probabilities per calendar year BP.
#' @examples
#' df <- data.frame(calBP=5000:2000,PrDens=runif(length(5000:2000)))
#' mycalgrid <- as.CalGrid(df)
#' plot(mycalgrid)
#' @export
#' 
as.CalGrid <- function(x) {
  df <- as.data.frame(x)
  if (ncol(x) == 2){
    names(df) <- c("calBP", "PrDens")
  } else {
    stop("Input must be 2 columns.")
  }
  class(df) <- c("CalGrid", class(df)) 
  return(df)
}

#' @title Convert data to class UncalGrid. 
#'
#' @description Tries to coerce any two-column matrix or data.frame to a uncalibrated probability distribution (an object of class "UncalGrid") for use by the rcarbon package. 
#' 
#' @param x A two-column \code{matrix} or \code{data.frame} class object.
#'
#' @return A CalGrid class object of probabilities or summed probabilities per CRA.
#' @examples
#' df <- data.frame(CRA=5000:2000,PrDens=runif(length(5000:2000)))
#' mycalgrid <- as.UncalGrid(df)
#' @export

as.UncalGrid <- function(x) {
  df <- as.data.frame(x)
  if (ncol(x) == 2){
    names(df) <- c("CRA", "PrDens")
  } else {
    stop("Input must be 2 columns.")
  }
  class(df) <- c("UncalGrid", class(df)) 
  return(df)
}



#' @title Convert to a CalDates object 
#' @description Convert other calibrated date formats to an rcarbon CalDates object. 
#' @param x One or more calibrated dated to convert (currently only BchronCalibratedDates and oxcAARCalibratedDatesList objects are supported)
#' @return A CalDates object
#' @examples
#' \dontrun{ 
#' library(Bchron)
#' library(oxcAAR)
#' quickSetupOxcal()
#' dates <- data.frame(CRA=c(3200,2100,1900), Error=c(35,40,50))
#' bcaldates <- BchronCalibrate(ages=dates$CRA, ageSds=dates$Error, 
#' calCurves=rep("intcal13", nrow(dates)))
#' rcaldates <- rcarbon::calibrate(dates$CRA, dates$Error, calCurves=rep("intcal13"))
#' ocaldates <- oxcalCalibrate(c(3200,2100,1900),c(35,40,50),c("a","b","c"))
#' ## Convert to rcarbon format
#' caldates.b <- as.CalDates(bcaldates)
#' caldates.o <- as.CalDates(ocaldates)
#' ## Comparison plot
#' plot(rcaldates$grids[[2]]$calBP,rcaldates$grids[[2]]$PrDens, 
#' type="l", col="green", xlim=c(2300,1900))
#' lines(caldates.b$grids[[2]]$calBP,caldates.b$grids[[2]]$PrDens, col="red")
#' lines(caldates.o$grids[[2]]$calBP,caldates.o$grids[[2]]$PrDens, col="blue")
#' legend("topright", legend=c("rcarbon","Bchron","OxCal"), col=c("green","red","blue"), lwd=2)
#' }
#' @export
#' 
as.CalDates <- function(x){
  if (!any(class(x)%in%c("BchronCalibratedDates","oxcAARCalibratedDatesList"))){
    stop("Currently, x must be of class BchronCalibratedDates or oxcAARCalibratedDatesList")
  }
  if (any(class(x)=="BchronCalibratedDates")){	    
    methods <- "Bchron"
    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(x))
    names(sublist) <- names(x)
    names(reslist) <- c("metadata","grids")
    ## metadata
    df <- as.data.frame(matrix(ncol=11, nrow=length(x)), stringsAFactors=FALSE)
    names(df) <- c("DateID","CRA","Error","Details","CalCurve","ResOffsets","ResErrors","StartBP","EndBP","Normalised","CalEPS")
    df$DateID <- names(x)
    df$CRA <- as.numeric(unlist(lapply(X=x, FUN=`[[`, "ages")))
    df$Error <- as.numeric(unlist(lapply(X=x, FUN=`[[`, "ageSds")))
    df$CalCurve <- as.character(unlist(lapply(X=x, FUN=`[[`, "calCurves")))
    df$ResOffsets <- NA
    df$ResErrors <- NA
    df$StartBP <- NA
    df$EndBP <- NA
    df$Normalised <- TRUE
    reslist[["metadata"]] <- df
    ## grids
    for (i in 1:length(x)){
      tmp <- x[[i]]
      res <- data.frame(calBP=rev(tmp$ageGrid),PrDens=rev(tmp$densities))
      class(res) <- append(class(res),"calGrid")        
      sublist[[i]] <- res
    }
    reslist[["grids"]] <- sublist
    reslist[["calmatrix"]] <- NA
    class(reslist) <- c("CalDates",class(reslist))
    return(reslist)
  }
  if (any(class(x)=="oxcAARCalibratedDatesList")){
    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(x))
    names(sublist) <- names(x)
    names(reslist) <- c("metadata","grids")
    ## metadata
    df <- as.data.frame(matrix(ncol=11, nrow=length(x)), stringsAFactors=FALSE)
    names(df) <- c("DateID","CRA","Error","Details","CalCurve","ResOffsets","ResErrors","StartBP","EndBP","Normalised","CalEPS")
    df$DateID <- names(x)
    df$CRA <- as.numeric(unlist(lapply(X=x, FUN=`[[`, "bp")))
    df$Error <- as.numeric(unlist(lapply(X=x, FUN=`[[`, "std")))
    df$CalCurve=lapply(lapply(lapply(lapply(lapply(x,FUN=`[[`,"cal_curve"),FUN=`[[`,"name"),strsplit," "),unlist),FUN=`[[`,2)
    df$CalCurve=tolower(df$CalCurve)
    df$ResOffsets <- NA
    df$ResErrors <- NA
    df$StartBP <- NA
    df$EndBP <- NA
    df$Normalised <- TRUE
    df$CalEPS <- 0
    reslist[["metadata"]] <- df
    ## grids
    for (i in 1:length(x)){
      tmp <- x[[i]]$raw_probabilities  
      rr <- range(tmp$dates)
      res <- 	approx(x=tmp$dates,y=tmp$probabilities,xout=ceiling(rr[1]):floor(rr[2]))
      res$x <- abs(res$x-1950)
      res <- data.frame(calBP=res$x,PrDens=res$y)
      res$PrDens <- res$PrDens/sum(res$PrDens)
      class(res) <- append(class(res),"calGrid")        
      sublist[[i]] <- res
    }
    reslist[["grids"]] <- sublist
    reslist[["calmatrix"]] <- NA
    class(reslist) <- c("CalDates",class(reslist))
    return(reslist)
  }       
}


#' @export

"[.CalDates" <- function(x,i){
  
  if (nrow(x$metadata)==0){
    stop("No data to extract")
  }
  if(!missing(i)) {
    if (all(is.numeric(i)) | all(is.character(i)) | all(is.logical(i))){
      if (length(x$calmatrix)>1){
        res <- list(metadata=x$metadata[i,], grids=NA, calmatrix=x$calmatrix[,i])
      } else {
        res <- list(metadata=x$metadata[i,], grids=x$grids[i], calmatrix=NA)
      }
      class(res) <- c("CalDates", class(res))        
    } else {
      stop("i must be a numeric, character or logical vector of length(x)")
    }
    return(res)
  }           
}


#' @export

length.CalDates <- function(x,...)
{
  return(nrow(x$metadata))
}

#' @title Computes the highest posterior density interval (HPDI) for calibrated dates 
#' 
#' @param x A \code{CalDates} class object.
#' @param credMass Interval probability mass
#' 
#' @return A list of matrices with HPDI ranges and associated probabilities in cal BP. Note that the sum of the probability mass will be only approximately equal to the argument \code{credMass}.
#' @examples
#' x <- calibrate(c(3050,2950),c(20,20))
#' hpdi(x)
#' @seealso \code{\link{calibrate}}
#' @export


hpdi <- function(x, credMass=0.95){
  
  cl <- class(x)
  if (!"CalDates"%in%cl){
    stop("x must be of class CalDates")
  }
  n <- nrow(x$metadata)
  result <- vector("list",length=n)
  for (i in 1:n){
    if (length(x$calmatrix)>1){
      grd <- data.frame(calBP=as.numeric(row.names(x$calmatrix)),PrDens=x$calmatrix[,i])
      grd <- grd[grd$PrDens >0,]
    } else {
      grd <- x$grids[[i]]
    }
    sorted <- sort(grd$PrDens , decreasing=TRUE)
    heightIdx = min( which( cumsum( sorted) >= sum(grd$PrDens) * credMass ) )
    height = sorted[heightIdx]
    indices = which( grd$PrDens >= height )
    gaps <- which(diff(indices) > 1)
    starts <- indices[c(1, gaps + 1)]
    ends <- indices[c(gaps, length(indices))]
    p_interval <- as.integer(length(starts))
    for (j in 1:length(starts))
      {
        p_interval[j] <- sum(grd$PrDens[starts[j]:ends[j]])
      }
    result[[i]] <- cbind(startCalBP = grd$calBP[starts], endCalBP = grd$calBP[ends],prob = p_interval) 
  }  
  return(result)
}

#' @title Summarise a \code{CalDates} class object
#'
#' @description Returns summary statistics of calibrated dates.
#'
#' @param object A \code{CalDates} class object.
#' @param prob A vector containing probabilities for the higher posterior density interval. Default is \code{c(0.683,0.954)}, i.e. 1 and 2-Sigma range.
#' @param calendar Whether the summary statistics should be computed in cal BP (\code{"BP"}) or in BCAD (\code{"BCAD"}).
#' @param ... further arguments passed to or from other methods.
#' @return A \code{data.frame} class object containing the ID of each date, along with the median date and one and two sigma (or a user specified probability) higher posterior density ranges.
#' @method summary CalDates
#' @export 

summary.CalDates<-function(object,prob=NA,calendar="BP",...) {
  
  foo = function(x,i){if(nrow(x)>=i){return(x[i,1:2])}else{return(c(NA,NA))}}
  if (is.na(prob)) 
  {
    prob = c(0.683,0.954)
    pnames = c("OneSigma","TwoSigma")
  } else {
    pnames = paste("p",prob,sep="_")
  }	
  pnames=paste(pnames,calendar,sep="_")
  probMats = vector("list",length=length(prob))
  for (i in 1:length(prob))
  {
    cols = max(unlist(lapply(hpdi(object,prob[i]),nrow)))
    tmpMatrix=matrix(NA,ncol=cols,nrow=nrow(object$metadata))
    
    for (j in 1:cols)
    {
      tmp=t(sapply(hpdi(object,prob[i]),foo,i=j))
      if (calendar=="BCAD")
      {
        tmp = t(apply(tmp,1,BPtoBCAD))
      }
      tmpMatrix[,j]=apply(tmp,1,paste,collapse=" to ")
    }
    colnames(tmpMatrix)=paste(pnames[i],1:cols,sep="_")
    probMats[[i]]=tmpMatrix
  }
  
  med.dates=medCal(object)
  
  if (calendar=="BP")
  {res=data.frame(DateID=object$metadata$DateID,MedianBP=med.dates)}
  else
  {res=data.frame(DateID=object$metadata$DateID,BPtoBCAD(med.dates))
  colnames(res)[2]="MedianBC/AD"}
  for (k in 1:length(probMats))
  {
    res=cbind.data.frame(res,probMats[[k]])
  }
  return(res)
}


#' @title Computes the median date of calibrated dates 
#'
#' @description Function for generating a vector median calibrated dates from a \code{CalDates} class object.
#' 
#' @param x A \code{CalDates} class object.
#'
#' @return A vector of median dates in cal BP
#' @examples
#' x <- calibrate(c(3050,2950),c(20,20))
#' medCal(x)
#' @seealso \code{\link{calibrate}}, \code{\link{barCodes}}
#' @export
medCal <- function(x)
{
  ndates=nrow(x$metadata)
  meddates=numeric()
  if (is.na(x$calmatrix[1]))
  {
    for (i in 1:ndates)
    {
      tmp=x$grids[[i]]
      tmp$Cumul=cumsum(tmp$PrDens)	 
      meddates[i]=tmp[which.min(abs(tmp$Cumul-max(tmp$Cumul)/2)),1]
    }		
  } else
    
  {
    cumcal=apply(x$calmatrix,2,cumsum)
    for (i in 1:ndates)
    {
      index = which.min(abs(cumcal[,i]-max(cumcal[,i])/2))
      meddates[i]=as.numeric(rownames(cumcal)[index])
    }
  }
  return(meddates)
}

#' @title Creates mixed calibration curves.
#'
#' @description Function for generating mixed calibration curves (e.g. intcal20/marine20, intcal20/shcal20)  with user-defined proportions.
#' 
#' @param calCurve1 Name of the first curve, chosen from 'intcal20','intcal13','intcal13nhpine16','shcal20','shcal13','shcal13shkauri16',''marine13','marine20'. Default is 'intcal20'.
#' @param calCurve2 Name of the second curve, from the same list. Default is 'marine20'.
#' @param p The proportion of contribution of the first curve. Default is 1.
#' @param resOffsets Offset value for the marine reservoir effect to be applied to calCurve2. Default is 0.
#' @param resErrors Error of the marine reservoir effect offset to be applied to calCurve2. Default is 0.
#' @return A three-column matrix containing calibrated year BP, uncalibrated age bp, and standard deviation. To be used as custom calibration curve for the \code{\link{calibrate}} function.
#' @details The function is based on the \code{mix.calibrationcurves} function of the \code{clam} package.   
#' @references 
#' Blaauw, M. and Christen, J.A.. 2011. Flexible paleoclimate age-depth models using an autoregressive gamma process. \emph{Bayesian Analysis}, 6, 457-474.
#' Blaauw, M. 2018. clam: Classical Age-Depth Modelling of Cores from Deposits. R package version 2.3.1. https://CRAN.R-project.org/packacge=clam
#' Marsh, E.J., Bruno, M.C., Fritz, S.C., Baker, P., Capriles, J.M. and Hastorf, C.A., 2018. IntCal, SHCal, or a Mixed Curve? Choosing a 14 C Calibration Curve for Archaeological and Paleoenvironmental Records from Tropical South America. Radiocarbon, 60(3), pp.925-940.
#' @examples
#' myCurve <- mixCurves('intcal20','marine20',p=0.7,resOffsets=300,resErrors=20)
#' x <- calibrate(4000,30,calCurves=myCurve)
#' @seealso \code{\link{calibrate}}
#' @export



mixCurves <- function (calCurve1 = "intcal20", calCurve2 = "marine20", p = 0.5, resOffsets = 0, 
                       resErrors = 0)
{
  File1 <- paste(system.file("extdata", package = "rcarbon"),"/", calCurve1, ".14c", sep = "")
  File2 <- paste(system.file("extdata", package = "rcarbon"),"/", calCurve2, ".14c", sep = "")
  options(warn = -1)
  c1 <- readLines(File1, encoding = "UTF-8")
  c2 <- readLines(File2, encoding = "UTF-8")
  c1 <- c1[!grepl("[#]", c1)]
  c2 <- c2[!grepl("[#]", c2)]
  c1.con <- textConnection(c1)
  c2.con <- textConnection(c2)
  c1 <- as.matrix(read.csv(c1.con, header = FALSE,stringsAsFactors = FALSE))[, 1:3]
  c2 <- as.matrix(read.csv(c2.con, header = FALSE,stringsAsFactors = FALSE))[, 1:3]
  close(c1.con)
  close(c2.con)
  options(warn = 0)
  colnames(c2) <- c("CALBP", "C14BP", "Error")
  colnames(c1) <- c("CALBP", "C14BP","Error")
  c2.mu <- approx(c2[, 1], c2[, 2], c1[,1], rule = 2)$y + resOffsets
  c2.error <- approx(c2[, 1], c2[, 3], c1[,1], rule = 2)$y
  c2.error <- sqrt(c2.error^2 + resErrors^2)
  mu <- p * c1[, 2] + (1 - p) * c2.mu
  error <- sqrt(p * c1[, 3]^2 + (1 - p) * c2.error^2)
  res = cbind(c1[, 1], mu, error)
  colnames(res) = c("CALBP", "C14BP", "Error")
  return(res)
}
