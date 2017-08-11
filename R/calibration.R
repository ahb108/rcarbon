#' @title Calibrate radiocarbon dates
#'
#' @description Function for calibrating one or more radiocarbon dates.
#'
#' @param ages A vector of uncalibrated radiocarbon ages .
#' @param errors A vector of standard deviations corresponding to each estimated radiocarbon age.
#' @param ids An optional vector of IDs for each date.
#' @param dateDetails An optional vector of details for each date which will be returned in the output metadata. 
#' @param calCurves Either a string naming a calibration curve already provided with the rcarbon package (currently 'intcal13', 'shcal13' and 'marine13' are possible; default is 'intcal13') or a custom calibration curve with three columns (calibrated year BP, uncalibrated age bp, standard deviation).
#' @param resOffsets A vector of offset values for any marine reservoir effect (default is no offset).
#' @param resErrors A vector of offset value errors for any marine reservoir effect (default is no offset).
#' @param timeRange Earliest and latest data to calibrate for, in calendar years. Posterior probabilites beyond this range will be excluded (the default is sensible in most cases).
#' @param normalised A logical variable indicating whether the calibration should be normalised or not. Default is TRUE.
#' @param eps Cut-off value for density calculation. Default is 1e-5.
#' @param calMatrix a logical variable indicating whether the age grid should be limited to probabilities higher than \code{eps}
#' @param ncores Number of cores/workers used for parallel execution. Default is 1 (>1 requires doParallel package).
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @details This function computes one or more calibrated radiocarbon ages using the method described in Bronk Ramsey 2008 (albeit not in F14C space). It is possible to specify different calibration curves or reservoir offsets individually for each date, and control whether the resulting calibrated distribution is normalised to 1 under-the-curve or not. Calculations can also be executed in parallel to reduce computing time.
#'
#' @return An object of class CalDates with the following elements
#' \itemize{
#' \item{\code{metadata}} {A data.frame containing relevant information regarding each radiocarbon date and the parameter used in the calibration process.}
#' \item{\code{grids}} {A list of calGrid class objects, containing the posterior probabilities for each calendar year. The most memor-efficient way to store calibrated dates, as only years with non-zero probability are stored, but aggregation methods such as spd() may then take longer to extract and combine multiple dates. NA when the parameter calMatrix is set to TRUE.} 
#' \item{\code{calMatrix}} {A matrix of probability values, one row per calendar year in timeRange and one column per date. By storing all possible years, not just those with non-zero probabilty, this approach takes more memory, but speeds up spd() and is suggested whenever the latter is to be used. NA when the parameter calMatrix is set to FALSE.}  
#' }
#'
#' @references 
#' Bronk Ramsey, C. 2008. Radiocarbon dating: revolutions in understanding, \emph{Archaeometry} 50.2: 249–75. DOI: https://doi.org/10.1111/j.1475-4754.2008.00394.x 
#'
#' @examples
#' x1 <- calibrate(ages=4000, errors=30)
#' plot(x1)
#' # Example with a Marine Date, using a DeltaR of 300 and a DeltaR error of 30
#' x2 <- calibrate(ages=4000, errors=30, calCurves='marine13', resOffsets=300, resErrors=30)
#' plot(x2)
#' @export

calibrate <- function (x, ...) {
   UseMethod("calibrate")
}

#' @rdname calibrate
#' @export

calibrate.default <- function(ages, errors, ids=NA, dateDetails=NA, calCurves='intcal13', resOffsets=0 , resErrors=0, timeRange=c(50000,0), normalised=TRUE, calMatrix=FALSE, eps=1e-5, ncores=1, verbose=TRUE){

    # age and error checks
    if (length(ages) != length(errors)){
        stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
    }
    if (!is.na(ids[1]) & (length(ages) != length(ids))){
        stop("Ages and errors (and ids/details/offsets if provided) must be the same length.")
    }
    if (any(is.na(ages))|any(is.na(errors))){
        stop("Ages or errors contain NAs")
    }
    # calCurve checks and set-up
    if (!all(class(calCurves)=="character")){
        if (any(class(calCurves) %in% c("matrix","data.frame"))){
            cctmp <- as.matrix(calCurves)
            if (ncol(cctmp)!=3 | !all(sapply(mycc,is.numeric))){
                stop("The custom calibration curve must have just three numeric columns.")
            } else {
                colnames(cctmp) <- c("CALBP","C14BP","Error")
                if (max(cctmp[,2]) < max(ages) | min(cctmp[,2]) > min(ages)){
                    stop("The custom calibration curve does not cover the input age range.")
                }
                cclist <- vector(mode="list", length=1)
                cclist[[1]] <- cctmp
                names(cclist) <- "custom"
                calCurves <- rep("custom",length(ages))
            }
        } else {
            stop("calCurves must be a character vector specifying one or more known curves or a custom three-column matrix/data.frame (see ?calibrate.default).")
        }
    } else {
        if (!all(calCurves %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16"))){
            stop("calCurves must be a character vector specifying one or more known curves or a custom three-column matrix/data.frame (see ?calibrate.default).")
        } else {
            tmp <- unique(calCurves)
            if (length(calCurves)==1){ calCurves <- rep(calCurves,length(ages)) }
            cclist <- vector(mode="list", length=length(tmp))
            names(cclist) <- tmp
            for (a in 1:length(tmp)){
                calCurveFile <- paste(system.file("data", package="rcarbon"), "/", tmp[1],".14c", sep="")
                options(warn=-1)
                cctmp <- readLines(calCurveFile, encoding="UTF-8")
                cctmp <- cctmp[!grepl("[#]",cctmp)]
                cctmp <- as.matrix(read.csv(textConnection(cctmp), header=FALSE, stringsAsFactors=FALSE))[,1:3]
                options(warn=0)
                colnames(cctmp) <- c("CALBP","C14BP","Error")
                cclist[[tmp[a]]] <- cctmp
            }
        }
    }
    # container and reporting set-up
    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(ages))
    if (calMatrix){
        calmBP <- seq(timeRange[1],timeRange[2],-1)
        calmat <- matrix(ncol=length(ages), nrow=length(calmBP))
        rownames(calmat) <- calmBP
        calmat[] <- 0
    }
    if (is.na(ids[1])){
        ids <- as.character(1:length(ages))
    } else {
        ids <- as.character(ids)
    }
    if (length(resOffsets)==1){ resOffsets <- rep(resOffsets,length(ages)) }
    if (length(resErrors)==1){ resErrors <- rep(resErrors,length(ages)) }
    names(sublist) <- ids
    names(reslist) <- c("metadata","grids")
    if (length(ages)>1 & verbose){
        print("Calibrating radiocarbon ages...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(ages), style=3)
    }
    # calibration
    if (ncores>1){
        # parallellised
        require(doParallel)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        if (verbose){ print(paste("Running in parallel (standard calibration only) on ",getDoParWorkers()," workers...",sep=""))}
        sublist <- foreach (b=1:length(ages)) %dopar% {
            calcurve <- cclist[[calCurves[b]]]
            calBP <- seq(max(calcurve),min(calcurve),-1)
            age <- ages[b] - resOffsets[b]
            error <- errors[b] + resErrors[b]
            mu <- approx(calcurve[,1], calcurve[,2], xout=calBP)$y
            tau <- error^2 + approx(calcurve[,1], calcurve[,3], xout=calBP)$y^2
            dens <- dnorm(age, mean=mu, sd=sqrt(tau))
            dens[dens < eps] <- 0	
            if (normalised){
                dens <- dens/sum(dens)
                dens[dens < eps] <- 0
                dens <- dens/sum(dens)
            }
            res <- data.frame(calBP=calBP,PrDens=dens)
            res <- res[which(calBP<=timeRange[1]&calBP>=timeRange[2]),]
            res <- res[res$PrDens > 0,]
            class(res) <- append(class(res),"calGrid")
            return(res)
        }
        stopCluster(cl)
        names(sublist) <- ids
        if (calMatrix){
            for (a in 1:length(sublist)){
                calmat[as.character(sublist[[a]]$calBP),a] <- sublist[[a]]$PrDens
            }
        }
    } else {
        # single core
        for (b in 1:length(ages)){
            if (length(ages)>1 & verbose){ setTxtProgressBar(pb, b) }
            calcurve <- cclist[[calCurves[b]]]
            calBP <- seq(max(calcurve),min(calcurve),-1)
            age <- ages[b] - resOffsets[b]
            error <- errors[b] + resErrors[b]
            mu <- approx(calcurve[,1], calcurve[,2], xout=calBP)$y
            tau <- error^2 + approx(calcurve[,1], calcurve[,3], xout=calBP)$y^2
            dens <- dnorm(age, mean=mu, sd=sqrt(tau))
            dens[dens < eps] <- 0
            if (normalised){
                dens <- dens/sum(dens)
                dens[dens < eps] <- 0
                dens <- dens/sum(dens)
            }
            res <- data.frame(calBP=calBP,PrDens=dens)
            res <- res[which(calBP<=timeRange[1]&calBP>=timeRange[2]),]
            if (calMatrix){ calmat[,b] <- res$PrDens }
            res <- res[res$PrDens > 0,]
            class(res) <- append(class(res),"calGrid")
            sublist[[ids[b]]] <- res
        }
    }
    # clean-up and results
    if (length(ages)>1 & verbose){ close(pb) }
    df <- data.frame(DateID=ids, CRA=ages, Error=errors, Details=dateDetails, CalCurve=calCurves,ResOffsets=resOffsets, ResErrors=resErrors, StartBP=timeRange[1], EndBP=timeRange[2], Normalised=normalised, CalEPS=eps, stringsAsFactors=FALSE)
    reslist[["metadata"]] <- df
    if (calMatrix){
        reslist[["grids"]] <- NA
        reslist[["calmatrix"]] <- calmat
    } else {
        reslist[["grids"]] <- sublist
        reslist[["calmatrix"]] <- NA
    }
    class(reslist) <- c("CalDates",class(reslist))
    if (verbose){ print("Done.") }
    return(reslist)
}

#' @export

calibrate.UncalGrid <- function(x, errors=0, calCurves='intcal13', timeRange=c(50000,0), compact=TRUE, eps=1e-5, type="fast", datenormalised=FALSE, spdnormalised=FALSE, verbose=TRUE, ...){

    if (length(errors)==1){
        errors <- rep(errors,length(x$CRA))
    }
    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- calcurve[!grepl("[#]",calcurve)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")
    if (type=="full"){
        caleach <- calibrate(ages=x$CRA, errors=errors, method="standard", normalised=datenormalised, compact=FALSE,...)
        tmp <- lapply(caleach$grids,`[`,2)
        tmp <- lapply(1:length(tmp),FUN=function(i) tmp[[i]]*x$PrDens[i])
        tmp <- do.call("cbind",tmp)
        res <- data.frame(calBP=caleach$grids[[1]]$calBP, PrDens=apply(tmp,1,sum))
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
    if (compact){ res <- res[res$PrDens > 0,] }
    class(res) <- c("CalGrid", class(res))   
    if (verbose){ print("Done.") }
    return(res)
}

#' @title Uncalibrate (back-calibrate) a radiocarbon date.
#'
#' @description Function for uncalibrating one or more radiocarbon dates.
#'
#' @param x Either a vector of uncalibrated radiocarbon ages or an object of class CalGrid.
#' @param CRAerrors A vector of standard deviations corresponding to each estimated radiocarbon age (ignored if x is a CalGrid object).
#' @param roundyear An optional vector of IDs for each date (ignored if x is a CalGrid object).
#' @param  calCurves A string naming a calibration curve already provided with the rcarbon package (currently 'intcal13', 'shcal13' and 'marine13' are possible; default is 'intcal13' and only one can currently be specified for all dates). 
#' @param  eps Cut-off value for density calculation (for CalGrid objects only).
#' @param  compact A logical variable indicating whether only uncalibrated ages with non-zero probabilities should be returned (for CalGrid objects only).
#' @param  verbose A logical variable indicating whether extra information on progress should be reported (for CalGrid objects only).
#' @details This function takes one or more calibrated calendars and looks-up the corresponding uncalibrated age, error of the stated calibration curve at that point. It also provides a randomised estimate of the uncalibrate age based on the curve error (and optionally also a hypothetical measurement error.
#'
#' @return A data.frame with specifying the original data, the uncalibrated age without the calibration curve error (ccCRA), the calibration curve error at this point in the curve (ccError), a randomised uncalibrated age (rCRA) given both the stated ccError and any further hypothesised instrumental error provided by the CRAerrors argument (rError). 
#'
#' @examples
#' # Uncalibrate two calendar dates
#' uncalibrate(c(3050,2950))
#' @export

uncalibrate <- function (x, ...) {
   UseMethod("uncalibrate", x)
}


#' @export

uncalibrate.default <- function(x, CRAerrors=NA, roundyear=TRUE, calCurves='intcal13'){
    
    if (length(CRAerrors)==1){ CRAerrors <- rep(CRAerrors,length(x)) } 
    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- calcurve[!grepl("[#]",calcurve)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")
    dates <- data.frame(approx(calcurve, xout=x))
    colnames(dates) <- c("calBP", "ccCRA")
    calcurve.error <- approx(calcurve[,c(1,3)], xout=dates$calBP)$y
    dates$ccError <- calcurve.error
    dates$rCRA <- rnorm(nrow(dates), mean=dates$ccCRA, sd=dates$ccError)
    dates$rError <- CRAerrors
    if (roundyear){ dates$rCRA <- round(dates$rCRA) }
    return(dates)
}

#' @export

uncalibrate.CalGrid <- function(x, calCurves='intcal13', eps=1e-5, compact=TRUE, verbose=TRUE){

    if (verbose){ print("Uncalibrating...") }
    names(x) <- c("calBP","PrDens")
    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- calcurve[!grepl("[#]",calcurve)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")
    mycras <- uncalibrate(x$calBP)
    res <- data.frame(CRA=max(calcurve[,2]):min(calcurve[,2]), PrDens=0)
    tmp <- vector(mode="list",length=nrow(mycras))
    basetmp <- vector(mode="list",length=nrow(mycras))
    if (length(tmp)>1 & verbose){
        flush.console()
        pb <- txtProgressBar(min=1, max=length(tmp), style=3)
    }
    for (a in 1:length(tmp)){
        basetmp[[a]] <- dnorm(res$CRA, mean=mycras$ccCRA[a], sd=mycras$ccError[a])
        tmp[[a]] <- basetmp[[a]] * x$PrDens[a]
        if (verbose){ setTxtProgressBar(pb, a) }
    }
    if (verbose){ close(pb) }
    unscGauss <- do.call("cbind",tmp)
    res$Raw <- rowSums(unscGauss)
    res$Raw[res$Raw < eps] <- 0
    base <- do.call("cbind",basetmp)
    res$Base <- rowSums(base)
    res$Raw[res$Raw < eps] <- 0
    res$PrDens[res$Base>0] <- res$Raw[res$Base>0] / res$Base[res$Base>0]
    if (compact){ res <- res[res$PrDens > 0,] }
    class(res) <- c("UncalGrid", class(res)) 
    if (verbose){ print("Done.") }
    return(res)
}

#' @export

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

#' @export

as.CalDates <- function(x){
    cl <- class(x)
    if (cl!="BchronCalibratedDates"){
	    stop("x must be of class BchronCalibratedDates")
    }
    methods <- "Bchron"
    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(x))
    ids <- as.character(1:length(x))
    names(sublist) <- ids
    names(reslist) <- c("metadata","grids")
    ages <- unlist(lapply(x,function(x){return(x[[1]])}))
    errors <-  unlist(lapply(x,function(x){return(x[[2]])}))
    calCurves <- as.character(unlist(lapply(x,function(x){return(x[[3]])})))

    for (i in 1:length(x))
    {
	tmp <- x[[i]]
	res <- data.frame(calBP=rev(tmp$ageGrid),PrDens=rev(tmp$densities))
        class(res) <- append(class(res),"calGrid")        
	calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves[i],".14c", sep="")
        options(warn=-1)
        cctmp <- readLines(calCurveFile, encoding="UTF-8")
        cctmp <- cctmp[!grepl("[#]",cctmp)]
        cctmp <- as.matrix(read.csv(textConnection(cctmp), header=FALSE, stringsAsFactors=FALSE))[,1]
        options(warn=0)
        calBP <- seq(max(cctmp),min(cctmp),-1)
	rownames(res) <- match(res[,1],calBP)
	sublist[[ids[i]]] <- res
    }	    
	 
    df <- data.frame(DateID=ids, CRA=ages, Error=errors, Details=NA, CalCurve=calCurves,ResOffsets=NA, ResErrors=NA, StartBP=NA, EndBP=NA, CalMethod="Bchron", Normalised=TRUE, CalEPS=NA, stringsAsFactors=FALSE)
    reslist[["metadata"]] <- df
    reslist[["grids"]] <- sublist
    class(reslist) <- c("CalDates",class(reslist))
    return(reslist)
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
        result[[i]] <- cbind(startCalBP = grd$calBP[starts], endCalBP = grd$calBP[ends]) 
    }  
    return(result)
}


#' @title Summarise a CalDates object
#'
#' @export
summary.CalDates<-function(x,prob=NA,calendar="BP") {
	
	foo = function(x,i){if(nrow(x)>=i){return(x[i,])}else{return(c(NA,NA))}}
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
		cols = max(unlist(lapply(hpdi(x,prob[i]),nrow)))
		tmpMatrix=matrix(NA,ncol=cols,nrow=nrow(x$metadata))

		for (j in 1:cols)
		{
		tmp=t(sapply(hpdi(x,prob[i]),foo,i=j))
		if (calendar=="BC/AD")
		{
			tmp = t(apply(tmp,1,BPtoBCAD))
			
		}
                tmpMatrix[,j]=apply(tmp,1,paste,collapse=" to ")
		}
		colnames(tmpMatrix)=paste(pnames[i],1:cols,sep="_")
		probMats[[i]]=tmpMatrix
		}
      
        med.dates=medCal(x)

	if (calendar=="BP")
	{res=data.frame(DateID=x$metadata$DateID,MedianBP=med.dates)}
	else
	{res=data.frame(DateID=x$metadata$DateID,BPtoBCAD(med.dates))
	colnames(res)[2]="MedianBC/AD"}
	for (k in 1:length(probMats))
	{
	res=cbind.data.frame(res,probMats[[k]])
	}
return(res)
}



#' @export
medCal <- function(x)
{
	ndates=nrow(x$metadata)
	meddates=numeric()
	if (is.na(x$calmatrix))
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
