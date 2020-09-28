#' @title Rescale a numeric vector to a specified minimum and maximum 
#' @description Rescale a numeric vector to a specified minimum and maximum.  
#' @param x numeric vector to smooth.
#' @param type what kind of rescaling to perform. Current options are 'simple' (default) and 'normal' which produces a z-score and 'custom' for which the 'to' argument must be specified.
#' @param to numeric vector of length 2 specifying the minimum and maximum value to perform a linear rescale between (default is 0 and 1)
#' @param na.rm Set to TRUE,this removes NAs before rescaling.
#' @return A numeric vector of rescaled values.
#' @examples
#' reScale(15:200)
#' @import stats
#' @export

reScale <- function(x, type="simple", to=c(0,1), na.rm=TRUE){

    types <- c("simple","normal")
    if (!type %in% types){
        stop("The rescale type you have chosen is not currently an option.")
    }
    if (max(x)-min(x)==0){
        warning("All the values in x are the same, and will just be recentred on 0 if type='normal' or max(to) if type='simple'.")
        if (type=="normal"){ res <- rep(0,length(x)) } else { res <- rep(max(to), length(x)) }
        return(res)
    }
    if (na.rm){ x <- na.omit(x) }
    if (type=="normal"){
        res <- (x-mean(x))/sd(x)
    } else {
        xrange <- range(x)
        mfac <- (to[2] - to[1])/(xrange[2] - xrange[1])
        res <- to[1] + (x - xrange[1]) * mfac
    }
    return(res)
}

#' @title Calculate a running mean from a numeric vector. 
#' @description Calculate a running mean from a numeric vector.  
#' @param x numeric vector to smooth.
#' @param n the size of the window in which to smooth.
#' @param edge How to treat edge cases where a full window is unavailable. Current options are 'NA' to fill with NAs or 'fill' to fill with original values 
#' @return A numeric vector of smoothed values.
#' @examples
#' x <- rnorm(1000)
#' y <- c(1:1000)
#' plot(y,x, type="l")
#' lines(runMean(x,50), col="red")
#' @import stats
#' @export

runMean <- function(x, n, edge="NA"){
    res <- x
    tmp <- filter(res,rep(1/n,n), sides=2)
    if (edge == "fill"){
        res[!is.na(tmp)] <- tmp[!is.na(tmp)]
    } else {
        res <- tmp
    }
    return(res)
}

#' Smooth a numeric vector using a Gaussian window
#' 
#' @description Smooth a numeric vector using a Gaussian window
#' @param x numeric vector of values to smooth.
#' @param alpha numeric value controlling the size of the gaussian smoothing window. Proportional to the standard deviation of the Gaussian smoothing kernel where sd=(N-1)/(2*alpha) with N being the length of the input vector.
#' @param window a fraction between 0 and 1 representing the proportion of the input vector to include in the moving window.
#' @details Adapted from \code{smth.gaussian} in the \code{smoother} package. 
#' @references 
#' Hamilton, N. (2015). smoother: Functions Relating to the Smoothing of Numerical Data, R package version 1.1, https://CRAN.R-project.org/package=smoother
#' @examples
#' smoothGauss(runif(200),alpha=5)
#' @import stats
#' @export

smoothGauss <- function(x, alpha, window=0.1){

    windowLength <- as.integer(max(abs(window*length(x)),1))
    hw <- abs(windowLength / 2.0)
    w <- sapply(c(0:(windowLength-1)), function(x){
        n <- x - as.integer(hw)
        k <- -0.5 * (abs(alpha) * n / hw) ^2
        exp(1)^k
    })
    sizeW <- length(w)
    sizeD <- length(x)
    w <- w/sum(w)
    hkwL <- as.integer(sizeW/2) 
    hkwR <- sizeW - hkwL
    smthfun <- function(i){
        ix.d <- c((i-hkwL):(i+hkwR-1))
        ix.w <- which(ix.d %in% 1:sizeD)
        ix.d <- ix.d[ix.w]
        if (length(ix.w) != sizeW){
            W.nm <- w[ix.w] / sum(w[ix.w])
        } else {
            W.nm <- w
        }  
        D.nm <- x[ix.d]
        as.numeric(D.nm %*% W.nm)
    }
    res <- sapply(c(1:sizeD), FUN=smthfun)
    res[c(1:hkwL,(sizeD - hkwR + 1):sizeD)] <- NA # remove tails
    return(res)
}

#' @import utils
#' @keywords internal

rangecheck <- function(x, bins, timeRange, datenormalised=FALSE){
    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    caldateTR <- as.numeric(x$metadata[1,c("StartBP","EndBP")])
    caldateyears <- seq(caldateTR[1],caldateTR[2],-1)
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))
    for (b in 1:length(binNames)){
        index <- which(bins==binNames[b])
        if (length(x$calmatrix)>1){
                tmp <- x$calmatrix[,index, drop=FALSE]
                if (datenormalised){
                    tmp <- apply(tmp,2,FUN=function(x) x/sum(x))
                }
                spdtmp <- rowSums(tmp)
                if (length(binNames)>1){
                    spdtmp <- spdtmp / length(index)
                }
                binnedMatrix[,b] <- spdtmp[caldateyears<=timeRange[1] & caldateyears>=timeRange[2]]
            
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
    }
    return(sum(apply(binnedMatrix,2,sum)==0)/ncol(binnedMatrix)*100)
}

#' @title Convert BP dates to BC/AD format 
#' @description Converts calibrated BP dates to BC/AD dates, omitting `year 0' 
#' @param x A numerical vector (currently only basic checks that these numbers are in a sensible range). 
#' @return A vector with BC/BCE dates expressed as negative numbers and AD/CE dates as positive ones.
#' @examples
#' BPtoBCAD(4200)
#' @export

BPtoBCAD <- function(x){
    if (any(x < 0)){ stop("Post-bomb dates (<0 BP) are not currently supported.") }
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[x < 1950,2] <- 1950-res[x < 1950,1]
    res[x >= 1950,2] <- 1949-res[x >= 1950,1]
    return(res[,2])
}

#' @title Convert BC/AD dates to BP format
#' @description Converts BC/AD dates to BP format while handling the absence of 'year 0' 
#' @param x A numerical vector (currently only basic checks that these numbers are in a sensible range).
#' @return A vector with BP dates.
#' @examples
#' BCADtoBP(-1268)
#' @export

BCADtoBP <- function(x){
    if (any(x == 0)){ stop("0 BC/AD is not a valid year.") }
    if (any(x > 1950)){ stop("Post-bomb dates (> AD 1950) are not currently supported.") }
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[x > 0,2] <- abs(res[x > 0,1] - 1950)
    res[x < 0,2] <- abs(res[x < 0,1] - 1949)
    return(res[,2])
}


#' @title Computes the median date of each bin
#'
#' @description Function for generating a vector of median calibrated dates for each each bin.
#' 
#' @param x A \code{CalDates} class object.
#' @param bins vector containing the bin names associated with each radiocarbon date. Can be generated using \code{\link{binPrep}}.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @return A vector of median dates in cal BP
#' @examples
#' \dontrun{
#' #Load EUROEVOL Data
#' data(euroevol)
#' #Subset Danish Dates
#' denmark <- subset(euroevol,Country=="Denmark")
#' #Calibrate and Bin
#' denmarkDates <- calibrate(x=denmark$C14Age,errors=denmark$C14SD) 
#' denmarkBins <- binPrep(sites=denmark$SiteID,ages=denmark$C14Age,h=200) #200 years bin size
#' #Compute median date for each bin
#' binMed(x=denmarkDates,bins=denmarkBins)
#' }
#' @seealso \code{\link{binPrep}},\code{\link{barCodes}}
#' @import utils
#' @import stats
#' @export

binMed <- function(x,bins,verbose=TRUE){
	if (!"CalDates" %in% class(x)){
        stop("x must be an object of class 'CalDates'.")
    }
    if (length(bins)>1){
        nbins <- length(unique(bins))
        if (any(is.na(bins))){
            stop("Cannot have NA values in bins.")
        }
        if (length(bins)!=nrow(x$metadata)){
            stop("bins (if provided) must be the same length as x.")
        }
    } else {
        bins <- rep("0_0",nrow(x$metadata))
    }
    binNames <- unique(bins)
    caltimeRange =c(55000,0)
    if (any(x$metadata$CalCurve %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")))
    {
      caltimeRange =c(50000,0)
    }    
    calyears <- data.frame(calBP=seq(caltimeRange[1], 0,-1))
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))
    if (verbose){
        if (length(x$calmatrix)>1){
            print("Aggregating...")
        } else {
            print("Extracting and aggregating...")
        }
    }
    if (verbose & length(binNames)>1){
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3)
    }
    caldateTR <- as.numeric(x$metadata[1,c("StartBP","EndBP")])
    caldateyears <- seq(caldateTR[1],caldateTR[2],-1)
    check <- caldateTR[1] >= caltimeRange[1] & caldateTR[2] <= 0
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        if (length(x$calmatrix)>1){
            if (!check){
                stop("The time range of the calibrated dataset must be at least as large as the spd time range.")
            } else {
                tmp <- x$calmatrix[,index, drop=FALSE]
                spdtmp <- rowSums(tmp)
                if (length(binNames)>1){
                    spdtmp <- spdtmp / length(index)
                }
                binnedMatrix[,b] <- spdtmp[caldateyears<=caltimeRange[1] & caldateyears>=0]
            }
        } else {
            slist <- x$grids[index]
            slist <- lapply(slist,FUN=function(x) merge(calyears,x, all.x=TRUE)) 
            slist <- rapply(slist, f=function(x) ifelse(is.na(x),0,x), how="replace")
            slist <- lapply(slist, FUN=function(x) x[with(x, order(-calBP)), ])
            tmp <- lapply(slist,`[`,2)
            if (length(binNames)>1){
                spdtmp <- Reduce("+", tmp) / length(index)
            } else {
                spdtmp <- Reduce("+", tmp)
            }
            binnedMatrix[,b] <- spdtmp[,1]
        }
    }
close(pb)
print("Done")
cumcal=apply(binnedMatrix,2,cumsum)

medbins=numeric()
for (i in 1:nbins)
{
medbins[i] = calyears[which.min(abs(cumcal[,i]-max(cumcal[,i])/2)),1]
}
return(medbins)
}


#' @title Compute weights from distance matrix
#'
#' @description Function for computing a matrix of gaussian or fixed weights from distance matrix
#'
#' @param distmat a symmetric matrix of inter-site distances (in km). 
#' @param h parameter of the Gaussian distance decay function.
#' @param kernel indicates the type of weighting function, either 'fixed' or 'gaussian'. Default is 'gaussian'. 
#'
#' @details This function generates a weight matrix (required for the \code{\link{SPpermTest}}) function. When \code{kernel=="fixed"}, the weight \eqn{w_{ij}} between site \eqn{i} and \eqn{j} is equal to 1 when their interdistance \eqn{d_{ij}} is below \code{h}, and equal to 0 when  \eqn{d_{ij}>h}.When \code{kernel=="gaussian"}, the weight is calculated with formula exp(-d_{ij}^2/h^2).
#'
#' @return An object of class spatialweights
#'
#' @examples
#' lon <- c(11.3426,0.1278,0.1218)
#' lat <- c(44.4949,51.5074,52.2053)
#' library(sp)
#' d <- spDists(x=cbind(lon,lat),y=cbind(lon,lat))
#' spweights(d,h=100)
#' spweights(d,h=100,kernel="fixed")
#' @import stats
#' @export


spweights<-function(distmat,h=NULL,kernel="gaussian")
{
    w=matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
    kernels <- c("gaussian","fixed")
    if (!kernel %in% kernels){
                stop("The kernel you have chosen is not currently an option.")
           }
    if (is.null(h))
    {
                stop("Distance parameter h undefined")  
    }
    for (x in 1:nrow(distmat))
        {
            if (kernel=="gaussian")
                {w[x,]=exp(-distmat[x,]^2/h^2)}
            if (kernel=="fixed")
                {w[x,]=as.numeric(distmat[x,]<=h)}
        }
    res=list(w=w,h=h,kernel=kernel)
    class(res) <- append(class(res),"spatialweights")
    return(res)   
}

#' @title Computes rates of change from SPDs
#'
#' @description Function for computing rates of change between abutting user-defined time-blocks.
#'
#' @param spd Summed Probability Distribution obtained using the \code{\link{spd}} function. 
#' @param breaks A vector giving the breakpoints between the time-blocks.
#' @param backsight A single numeric value defining the distance in time between the focal year and the backsight year for computing the rate of change.
#' @param changexpr An expression defining how the rate of change is calculated, where \code{t1} is the summed probability for a focal block or year, \code{t0} is the summed probability for previous block or backsight year, and \code{d} is the duration of the block or the length of the backsight. Default is a geometric growth rate (i.e \code{expression((t1/t0)^(1/d)-1)}).
#' @details When the argument \code{breaks} is supplied the function aggregates the summed probability within each time-block and compared them across abutting blocks using the expression defined by \code{changexpr}. When the argument \code{backsight} is provided he expression is based on the comparison between the summed probability of each year and the associated backsight year.  
#'
#' @return An object of class \code{spdRC}.
#'
#' @examples
#' \dontrun{
#' data(emedyd)
#' caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE, calMatrix=TRUE)
#' bins <- binPrep(sites=emedyd$SiteName, ages=emedyd$CRA, h=50)
#' emedyd.spd <- spd(caldates,bins,timeRange=c(16000,9000),runm=100)
#' emedyd.gg <- spd2rc(emedyd.spd,breaks=seq(16000,9000,-1000))
#' emedyd.gg2 <- spd2rc(emedyd.spd,backsight=10)
#' plot(emedyd.gg)
#' plot(emedyd.gg2)
#' }

#' @import stats
#' @export

spd2rc <- function(spd,breaks=NULL,backsight=NULL,changexpr=expression((t1/t0)^(1/d)-1))
{
	if (is.null(breaks)&is.null(backsight))
	{
		stop('Either breaks or backsight should be provided')
	}

	if (!is.null(breaks)&!is.null(backsight))
	{
		stop('Both breaks and backsight cannot be provided')
	}

	if (is.null(backsight))
	{
		if (length(unique(round(abs(diff(breaks)))))!=1)
		{
			stop("Unequal break intervals is not supported")
		}
		nBreaks = length(breaks)-1
	}

	timeRange = eval(parse(text=spd$metadata[2]))
	timeSequence = timeRange[1]:timeRange[2]
	
	# if block based:
	if(!is.null(breaks))
	{
		type='blocks'
		obs=numeric(length=nBreaks)    
		for (x in 1:nBreaks)
		{
			index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
			obs[x]=sum(spd$grid[index,2])
		}

		res=numeric(length=nBreaks-1)
		for (i in 1:(nBreaks-1))
		{
			d=abs(breaks[i+1]-breaks[i]) 	
			t0 = obs[i]
			t1 = obs[i+1]
			res[i] = eval(changexpr)
			if (t1==0|t0==0){res[i]=NA}
			# 		res[i]=(obs[i+1]/obs[i])^(1/d)-1
		}
	}

	# if backsight based:
	if (!is.null(backsight))
	{
		type ='backsight'
		breaks = NA
		res=rep(NA,length(timeSequence))
		obs=NA

		for (i in 1:c(length(timeSequence)-backsight))
		{
			d=backsight 	
			t0 = spd$grid$PrDens[i]
			t1 = spd$grid$PrDens[i+backsight]
			res[i+backsight] = eval(changexpr)
			if (t1==0|t0==0){res[i+backsight]=NA}
		}
	}

	res=list(sumblock=obs,roca=res,breaks=breaks,timeSequence=timeSequence,type=type)
	class(res) <- append(class(res),"spdRC")
	return(res)   
}


#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils

curveSamples <- function(bins,calCurves,nsim)
{
	x = table(bins,calCurves)
	x = prop.table(x,1)
        x = replicate(nsim,table(factor(apply(x,1,function(x,y){sample(x=y,size=1,prob=x)},y=colnames(x)),levels=unique(calCurves))))
	return(t(x))
}


#' @title Sample random calendar dates
#'
#' @description Randomly samples calendar dates from each calibrated date or bin.
#' @param x A 'CalDates' class object. 
#' @param bins A vector containing the bin names associated with each radiocarbon date. If set to NA, binning is not carried out. 
#' @param nsim Number of sampling repetitions.
#' @param boot A logical value indicating whether bootstrapping is carried out (see details below). Default is FALSE.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @details The function randomly samples calendar dates based from calibrated probability distributions. When the \code{bins} argument is supplied a single calendar date is sampled from each bin. When \code{boot=TRUE}, dates (or bins) are randomly sampled with replacement before calendar dates are sampled. 
#' 
#' @return An object of class \code{simdates} with the following elements
#' \itemize{
#' \item{\code{sdates}} {A matrix containing the randomly sampled calendar dates, with rows containing each of the \code{nsim} repetitions.}
#' \item{\code{weight}} {A vector (or matrix in when \code{boot=TRUE}) containing the total area under the curve of each date, normalised to sum to unity. Notice this will be identical for all dates if the calibration is carried out with the argument \code{normalised} set to TRUE.}
#'}
#'
#' @import stats 
#' @import utils 
#' @export


sampleDates <- function(x,bins=NA,nsim,boot=FALSE,verbose=TRUE)
{
	# initial checks ####
	if (!"CalDates" %in% class(x)){
		stop("x must be an object of class 'CalDates'.")
	}
	if (length(bins)>1){
		if (any(is.na(bins))){
			stop("Cannot have NA values in bins.")
		}
		if (length(bins)!=nrow(x$metadata)){
			stop("bins (if provided) must be the same length as x.")
		}
	}
	else {
		bins <- 1:(length(x))
	}
	uni.bins = unique(bins)
	nbins = length(uni.bins)
	binList = vector('list',length=nbins)
	binLength=numeric()
	calmatrix=FALSE
	if (anyNA(x$grids)){calmatrix=TRUE}
	if (calmatrix)
	{
		warning("Processing time is slower when dates are calibrated using calMatrix=TRUE",immediate.=TRUE)
	}
	# Aggregate Bins ####
	if (verbose){
		flush.console()
		print("Aggregating...")
		pb <- txtProgressBar(min = 1, max = nbins,style = 3)
	}
	for (b in 1:nbins)
	{
		if (verbose){setTxtProgressBar(pb, b)}
		tmp.dates=x[which(bins==uni.bins[b])]
		if (calmatrix)
		{
			tmp.prob=tmp.dates$calmatrix
			binLength[b]=1
			if (length(tmp.dates)>1)
			{	
				tmp.prob=apply(tmp.dates$calmatrix,1,sum)
				binLength[b]=ncol(tmp.dates$calmatrix)
			}
			tmp.prob=tmp.prob[which(tmp.prob>0)]
			binList[[b]]=data.frame(calBP=as.numeric(names(tmp.prob)),PrDens=tmp.prob)
		} else {
			if (length(tmp.dates)==1)
			{
				binList[[b]]=tmp.dates$grids[[1]]
				binLength[b]=1
			} else {
				binLength[b]=length(tmp.dates)
				st.date = max(unlist(lapply(tmp.dates$grids,function(x){max(x$calBP)})))	
				end.date = min(unlist(lapply(tmp.dates$grids,function(x){min(x$calBP)}))) 
				tmp.prob = data.frame(calBP=st.date:end.date,PrDens=0)
				tmp.prob$PrDens=apply(sapply(tmp.dates$grids,function(x,r){
								     res = rep(0,length(r))
								     res[which(r%in%x$calBP)]=x$PrDens
								     return(res)},r=tmp.prob$calBP),1,sum)
				binList[[b]]=tmp.prob
			}
		}
	}
	if (verbose) {close(pb)}
	# Sample random dates ####
	if (!boot)
	{
		if (verbose){print("Simulating dates ...")}
		res=sapply(binList,function(x,nsim){sample(x$calBP,size=nsim,prob=x$PrDens,replace=TRUE)},nsim=nsim)
		weight=sapply(binList,function(x){return(sum(x$PrDens))})/binLength
		weight=weight/sum(weight)
	}
	if (boot)
	{
		res = matrix(NA,nrow=nsim,ncol=nbins)
		weight = matrix(NA,nrow=nsim,ncol=nbins)
		if (verbose){
			print("Bootstrapping...")
			pb <- txtProgressBar(min = 1, max = nsim,style = 3)
		}
		for (i in 1:nsim)
		{
			if (verbose){setTxtProgressBar(pb, i)}
			index=sample(nbins,replace=TRUE)
			res[i,]=sapply(binList[index],function(x,nsim){sample(x$calBP,size=1,prob=x$PrDens)})
			weight[i,]=sapply(binList[index],function(x){sum(x$PrDens)})/binLength[index]
			weight[i,]=weight[i,]/sum(weight[i,])
		}
		if (verbose) {close(pb)}
	}
	if (verbose){print("Done")}
	result=list(sdates=res,weight=weight)
	class(result) = c('simdates',class(result))
	return(result)
}

#' @title Gaussian weighting of dates relative to 
#' @description Rescale a numeric vector to a specified minimum and maximum.  
#' @param x A numeric vector or an object of class CalDates.
#' @param mean A single numeric value indicating the value to centre the Gaussian kernel on.
#' @param sd A single numeric value indicating the standard deviation of the Gaussian kernel to be used.
#' @param type The type of output to produce: currently either "weighted" (for a simple total weight value for each date) or "raw" (a list of reweighted calibrated radiocarbon probabilities for each calibrated date).
#' @return A numeric vector of weights (or optionally a list of reweighted calibrated radiocarbon probabilities).
#' @examples
#' ## Example weighting fo a set of dates versus a focal date of 5950 calBP
#' years <- seq(6500, 5500, -10)
#' plot(cbind(years, gaussW(years, 5950, 50)))
#' ## Example weighting of three calibrated dates  versus a focal date of 5950 calBP
#' dates <- calibrate(c(5280, 5180, 5080), c(30,30,30), normalised=FALSE)
#' gaussW(dates, 5950, 50)
#' ## Or the same with raw output
#' dateswt <- gaussW(dates, 5950, 50, type="raw")
#' head(dateswt[[1]])
#' @export
#' 
gaussW <- function(x, mean, sd, type="weights"){
    if (!class(x)[1] %in% c("numeric", "CalDates")){
        stop("Input must be a numeric vector of calibrated years BP or an object of class CalDates.")
    }
    if (!type %in% c("weights","raw")){ stop("type must be either weights or raw.") }
    base <- (dnorm(x=mean, mean=mean, sd=sd))
    if (class(x)[1]=="numeric"){
        res <- dnorm(x=x, mean=mean, sd=sd) / base
        return(res)
    } else {
        res <- vector(mode="list", length=length(x))
        wts <- vector(mode="numeric", length=length(x))
        for (a in 1: length(x)){
            wt <- dnorm(x=x$grids[[a]]$calBP, mean=mean, sd=sd) / base
            res[[a]] <- cbind(x$grids[[a]]$calBP, x$grids[[a]]$PrDens * wt)
            wts[a] <- sum(x$grids[[a]]$PrDens * wt)
        }
    }
    if (type=="raw"){
        return(res)
    } else if (type=="weights"){
        return(wts)
    }
}

rybcolourmap <- function(range, ...) {
    ## slightly modified from beachcolourmap() in the spatstat package
    col <- rybcolours(range, ...)
    z <- colourmap(col, range=range)
    return(z)
}

rybcolours <- function(range, sealevel=0, ncolours=100, nbeach=0){
    ## modified from beachcolours() in the spatstat package
    stopifnot(is.numeric(range) && length(range)==2)
    stopifnot(all(is.finite(range)))
    yr <- colorRampPalette(c("yellow","orangered","darkred"), space="rgb")
    cb <- colorRampPalette(c("blue","cyan","yellow"), space="rgb")
    depths <- range[1]
    peaks <- range[2]
    dv <- diff(range)/(ncolours - 1)
    epsilon <- nbeach * dv/2
    lowtide <- max(sealevel - epsilon, depths)
    hightide <-  min(sealevel + epsilon, peaks)
    countbetween <- function(a, b, delta) { max(0, round((b-a)/delta)) }
    nsea <- countbetween(depths, lowtide, dv)
    nbeach <- countbetween(lowtide,  hightide, dv)
    nland <- countbetween(hightide,  peaks, dv)
    colours <- character(0)
    if(nsea > 0)  colours <- cb(nsea) # cyan/blue
    if(nbeach > 0)  colours <- c(colours,rep("yellow",nbeach)) # yellow
    if(nland > 0)  colours <- c(colours, yr(nland)) # darkred/yellow
    return(colours)
}

#' @title Subsetting calibrated dates
#' @description Subsets calibrated dates (\code{CalDates} class object) based on Logical expressions of time intervals.  
#' @param x A CalDates class object
#' @param s Logical expression indicating dates to keep. The expression should include the term \code{BP} which refers to specific dates.
#' @param p Probability mass meeting the condition defined by \code{ss}.
#' @param ... Further arguments to be passed to or from other methods (ignored).
#' @details The function subsets \code{CalDates} class objects by identifying all dates that have a probability mass larger than \code{p} for a user defined logical expression of temporal interval containing the term \code{BP}, where \code{BP} refers to radiocarbon date. See examples for further detailes
#' @return A CalDates class object.
#' @examples
#' ## Generate some calibrated dates
#' x = calibrate(c(12100,5410,5320,3320),errors=c(20,20,30,30))
#' ## Subsets all dates that have a probability mass above 0.8 before 10000 BP
#' x2 = subset(x,BP>10000,p=0.8)
#' ## Subsets all dates that have a probability mass above 0.5 between 6000 and 6300 BP
#' x3 = subset(x,BP>6000&BP<6300,p=0.5)
#' @import stats
#' @export

subset.CalDates=function(x,s,p,...)
{
	index=rep(NA,length(x))
	s=substitute(s)
	calmat=anyNA(x$calmatrix)
	if (calmat)
	{
		for (i in 1:length(x))
		{
			tmp=x$grids[[i]]	 
			BP=tmp$calBP
			prdens=tmp$PrDens/sum(tmp$PrDens)   
			sumprob=sum(prdens[eval(s)])
			index[i]=sumprob>p
		}
	} else {
		BP=as.numeric(row.names(x$calmatrix))
		tmp=apply(x$calmatrix,2,function(x){x/sum(x)})
		index=apply(tmp[eval(s),],2,sum)>p 
	}
	return(x[which(index)])
}




#' @title Which indices for calibrated dates
#' @description Gives the TRUE indices of calibrated dates (\code{CalDates} class object) based on Logical expressions of time intervals.  
#' @param x A CalDates class object
#' @param s Logical expression indicating dates to keep. The expression should include the term \code{BP} which refers to specific dates.
#' @param p Probability mass meeting the condition defined by \code{ss}.
#' @details The function subsets \code{CalDates} class objects by identifying all dates that have a probability mass larger than \code{p} for a user defined logical expression of temporal interval containing the term \code{BP}, where \code{BP} refers to radiocarbon date. See examples for further detailes
#' @return A CalDates class object.
#' @examples
#' ## Generate some calibrated dates
#' x = calibrate(c(12100,5410,5320,3320),errors=c(20,20,30,30))
#' ## Subsets all dates that have a probability mass above 0.8 before 10000 BP
#' x2 = which.CalDates(x,BP>10000,p=0.8)
#' ## Subsets all dates that have a probability mass above 0.5 between 6000 and 6300 BP
#' x3 = which.CalDates(x,BP>6000&BP<6300,p=0.5)
#' @import stats
#' @export

which.CalDates=function(x,s,p)
{
	index=rep(NA,length(x))
	s=substitute(s)
	calmat=anyNA(x$calmatrix)
	if (calmat)
	{
		for (i in 1:length(x))
		{
			tmp=x$grids[[i]]	 
			BP=tmp$calBP
			prdens=tmp$PrDens/sum(tmp$PrDens)   
			sumprob=sum(prdens[eval(s)])
			index[i]=sumprob>p
		}
	} else {
		BP=as.numeric(row.names(x$calmatrix))
		tmp=apply(x$calmatrix,2,function(x){x/sum(x)})
		index=apply(tmp[eval(s),],2,sum)>p 
	}
	return(which(index))
}

#' @title Combine multiple CalDates Class Objects into one.
#'
#' @param ... \code{CalDates} class objects to be concatenated.
#' @param fixIDs logical. If set to TRUE, each date is given a new ID based on sequential integer. Default is FALSE
#' @return An object of class CalDates   
#' @examples  
#' x1 = calibrate(c(2000,3400),c(20,20),ids=1:2)
#' x2 = calibrate(c(4000,3000),c(30,30),calCurves=c('intcal20','marine20'),
#' resOffsets=c(0,30),resErrors=c(0,20),ids=3:4)
#' mcurve <- mixCurves('intcal20','marine20',p=0.7,resOffsets=300,resErrors=20)
#' x3 = calibrate(5300,20,calCurves=mcurve,ids=5)
#' x = combine(x1,x2,x3)
#' ## x$metadata
#' @seealso \code{\link{calibrate}}
#' @export
combine = function(...,fixIDs=FALSE)
{
	x=c(...)
	n = length(x)/3
	metadata.index=seq(from=1,length.out=n,by=3)
    grids.index = seq(from=2,length.out=n,by=3)
	calmatrices.index = seq(from=3,length.out=n,by=3)
	
	if (sum(is.na(x[calmatrices.index]))!=n&(sum(is.na(x[grids.index]))!=n))
	{
		stop('CalDates class object can be combined only if all elements were created with calMatrix=TRUE or calMatrix=FALSE')
	}
	
	if (all(!is.na(x[calmatrices.index])))
	{
		cmat=TRUE
	} else {cmat=FALSE}
	
	metadata=x[[metadata.index[1]]]
	grids=x[[grids.index[1]]]
	calmatrix=x[[calmatrices.index[1]]]
	
	for (i in 2:n)
	{
		metadata = rbind.data.frame(metadata,x[[metadata.index[i]]])
		if (cmat)
		{
			grids=NA
			calmatrix=cbind(calmatrix,x[[calmatrices.index[i]]])
		}
		
		if (!cmat)
		{
			calmatrix=NA
			grids = c(grids,x[[grids.index[i]]])
		}		
	}
	
	res=list(metadata=metadata,grids=grids,calmatrix=calmatrix)
	class(res)=c('CalDates','list')
	if (fixIDs){metadata$DateID=1:nrow(metadata)}
	if (any(duplicated(metadata$DateID)))
	{
	  stop('Date IDs must be unique. Consider setting ids explicitly when using calibrate() or set fixIDs to TRUE')
	}
	return(res)
}



#' @title Apply taphonomic corrections or other transformations to an SPD.
#'
#' @param x An object of class \code{CalSPD}, \code{compositeKDE} or \code{stackCalSPD}.
#' @param correction An expression for transforming the SPD. Available input terms include: CalBP, the vector of \code{calBP} year within the time range; and \code{PrDens}, a matching vector of summed probability. The default expression is the taphonomic correction formula proposed by Surovell et al 2009.
#' 
#' @return An object of the same class as x   
#' @examples  
#' \dontrun{
#'data(emedyd)
#'region1 = subset(emedyd,Region==1)
#'x = calibrate(x=region1$CRA, errors=region1$Error,normalised=FALSE)
#'bins = binPrep(sites=region1$SiteName, ages=region1$CRA,h=50)
#'region1.spd = spd(x=x,bins=bins,timeRange=c(16000,8000))
#'region1.spd.corrected = transformSPD(region1.spd)
#'}
#' @references 
#' Surovell, T.A., Finley, J.B., Smith, G.M., Brantingham, P.J., Kelly, R., 2009. Correcting temporal frequency distributions for taphonomic bias. Journal of Archaeological Science 36, 1715–1724.
#' @export

transformSPD = function(x,correction=expression(PrDens / (5.726442 * 10^6 * (CalBP+2176.4)^-1.3925309)))
{
  if (!any(class(x)%in%c('compositeKDE','CalSPD','stackCalSPD')))
  {
    stop("x must be of class 'CalSPD', 'compositeKDE', or 'stackCalSPD'")
  }
  
  if (any(class(x)%in%c('compositeKDE')))
  {
    x$res.matrix=apply(x$res.matrix,2,function(x,CalBP,expr){PrDens=x;return(eval(expr))},CalBP=x$timeRange[1]:x$timeRange[2],expr=correction)
  }
  
  if (any(class(x)%in%c('CalSPD')))
  {
    CalBP = x$grid$calBP
    PrDens=x$grid$PrDens
    x$grid$PrDens=eval(correction)
  }
  
  if (any(class(x)%in%c('stackCalSPD')))
  {
    for (i in 1:length(x$spds))
    {
    CalBP = x$spds[[i]]$grid$calBP
    PrDens= x$spds[[i]]$grid$PrDens
    x$spds[[i]]$grid$PrDens=eval(correction)
    }
  }
  
  return(x)
}






