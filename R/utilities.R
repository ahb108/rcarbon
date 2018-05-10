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


## quickMarks <- function(x, verbose=TRUE){

##     if (!"CalDates" %in% class(x)){
##         stop("Input must be of class \"CalDates\"")
##     }
##     df <- as.data.frame(matrix(ncol=8,nrow=nrow(x$metadata)), stringsasFactors=TRUE)
##     names(df) <- c("DateID","CRA","Error","qMed","q95s","q95e","q68s","q68e")
##     print("Extracting approximate values...")
##     if (nrow(x$metadata)>1 & verbose){
##         flush.console()
##         pb <- txtProgressBar(min=1, max=nrow(x$metadata), style=3)
##     }
##     for (a in 1:nrow(x$metadata)){
##         if (nrow(x$metadata)>1 & verbose){ setTxtProgressBar(pb, a) }
##         if (length(x$calmatrix)>1){
##             tmp <- data.frame(calBP=as.numeric(row.names(x$calmatrix)),PrDens=x$calmatrix[,a])
##             tmp <- tmp[tmp$PrDens >0,] 
##         } else {
##             tmp <- x$grids[[a]]
##         }
##         tmp <- tmp[tmp$PrDens>0,]
##         tmp <- tmp[with(tmp, order(-PrDens)), ]
##         tmp$Cumul <- cumsum(tmp$PrDens)
##         tmp$Cumul <- tmp$Cumul/max(tmp$Cumul)
##         tmp1 <- tmp[tmp$Cumul <= 0.95,]
##         df[a,"q95s"] <- min(tmp1$calBP)
##         df[a,"q95e"] <- max(tmp1$calBP)
##         wdth <- max(tmp1$calBP)-min(tmp1$calBP)
##         df[a,"q68s"] <- min(tmp1$calBP) + (wdth*0.25)
##         df[a,"q68e"] <- max(tmp1$calBP) - (wdth*0.25)
##         df[a,"qMed"] <- round(mean(c(df[a,"q95s"],df[a,"q95e"])),0)
##         df[a,"DateID"] <- as.character(x$metadata$DateID[a])
##         df[a,"CRA"] <- x$metadata$CRA[a]
##         df[a,"Error"] <- x$metadata$Error[a]
##     }
##     if (nrow(x$metadata)>1 & verbose){ close(pb) }
##     class(df) <- append(class(df),"quickMarks")
##     return(df)
## }



#' Smooth a numeric vector using a Gaussian window
#' 
#' @description Smooth a numeric vector using a Gaussian window
#' @param x numeric vector of values to smooth.
#' @param alpha numeric value controlling the size of the gaussian smoothing window. Proportional to the standard deviation of the Gaussian smoothing kernel where sd=(N-1)/(2*alpha) with N being the length of the input vector.
#' @param window a fraction between 0 and 1 representing the proportion of the input vector to include in the moving window.
#' @references
#' Adapted from \code{\link[smoother]{smth.gaussian}}
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
#' @param x A numerical vector (currently no checks that these numbers are in a sensible range). 
#' @return A vector with BC/BCE dates expressed as negative numbers and AD/CE dates as positive ones.
#' @examples
#' BPtoBCAD(4200)
#' @export

BPtoBCAD <- function(x){
    res <- matrix(c(x, rep(NA,length(x))), ncol=2)
    res[x < 1950,2] <- 1950-res[x < 1950,1]
    res[x >= 1950,2] <- 1949-res[x >= 1950,1]
    return(res[,2])
}

#' @title Convert BC/AD dates to BP format
#' @description Converts BC/AD dates to BP format while handling the absence of 'year 0' 
#' @param x A numerical vector (currently no checks that these numbers are in a sensible range).
#' @return A vector with BC/BCE dates expressed as negative numbers and AD/CE dates as positive ones.
#' @examples
#' BCADtoBP(-1268)
#' @export

BCADtoBP <- function(x){
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
    calyears <- data.frame(calBP=seq(50000, 0,-1))
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
    check <- caldateTR[1] >= 50000 & caldateTR[2] <= 0
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
                binnedMatrix[,b] <- spdtmp[caldateyears<=50000 & caldateyears>=0]
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

#' @title Compute geometric growth rates from SPDs
#'
#' @description Function for computing the geometric growth rates between abutting user-defined time-blocks.
#'
#' @param spd Summed Probability Distribution obtained using the \code{\link{spd}} function. 
#' @param breaks A vector giving the breakpoints between the time-blocks.
#'
#' @details The function computes the growth rate between abutting phases as \eqn{(X_{t}/X_{t+1})^{(1/d)}-1}, where \eqn{X_{t}} is the summed probability of radiocarbon dates in the block \eqn{t}, and \eqn{d} is the duration of the time-blocks. 
#'
#' @return An object of class \code{spdGG} containing the total summed probability for each time-block and the geometric growth rate between abutting blocks.
#'
#' @examples
#' \dontrun{
#' data(emedyd)
#' caldates <- calibrate(x=emedyd$CRA, errors=emedyd$Error, normalised=FALSE, calMatrix=TRUE)
#' bins <- binPrep(sites=emedyd$SiteName, ages=emedyd$CRA, h=50)
#' emedyd.spd <- spd(caldates,bins,timeRange=c(16000,9000))
#' emedyd.gg <- spd2gg(emedyd.spd,breaks=seq(16000,9000,-1000))
#' plot(emedyd.gg)
#' }

#' @import stats
#' @export

spd2gg <- function(spd,breaks)
{
	require(rcarbon)	
	if (length(unique(round(abs(diff(breaks)))))!=1)
	{
		stop("Unequal break intervals is not supported")
	}
	nBreaks = length(breaks)-1
	timeRange = eval(parse(text=spd$metadata[2]))
	timeSequence = timeRange[1]:timeRange[2]
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
		res[i]=(obs[i+1]/obs[i])^(1/d)-1
	}
	res=list(sumblock=obs,geomg=res,breaks=breaks)
	class(res) <- append(class(res),"spdGG")
	return(res)   
}



# 
# gaussW <- function(x, bw){
#     exp(-(x^2)/(2*(bw^2)))
# }
# 
# 
# 
# rybcolourmap <- function(range, ...) {
#   col <- rybcolours(range, ...)
#   z <- colourmap(col, range=range)
#   return(z)
# }
# 
# rybcolours <- function(range, sealevel=0, ncolours=100, nbeach=0){
#   stopifnot(is.numeric(range) && length(range)==2)
#   stopifnot(all(is.finite(range)))
#   yr <- colorRampPalette(c("yellow","orangered","darkred"), space="rgb")
#   cb <- colorRampPalette(c("blue","cyan","yellow"), space="rgb")
#   depths <- range[1]
#   peaks <- range[2]
#   dv <- diff(range)/(ncolours - 1)
#   epsilon <- nbeach * dv/2
#   lowtide <- max(sealevel - epsilon, depths)
#   hightide <-  min(sealevel + epsilon, peaks)
#   countbetween <- function(a, b, delta) { max(0, round((b-a)/delta)) }
#   nsea <- countbetween(depths, lowtide, dv)
#   nbeach <- countbetween(lowtide,  hightide, dv)
#   nland <- countbetween(hightide,  peaks, dv)
#   colours <- character(0)
#   if(nsea > 0)  colours <- cb(nsea) # cyan/blue
#   if(nbeach > 0)  colours <- c(colours,rep("yellow",nbeach)) # yellow
#   if(nland > 0)  colours <- c(colours, yr(nland)) # darkred/yellow
#   return(colours)
# }
# 
# spJitter <- function(pts, xamount, yamount=xamount){
#     proj <- NA
#     if (!is.na(proj4string(pts)) | proj4string(pts)!="NA"){
#         proj <- proj4string(pts)
#     }
#     if (class(pts) == "SpatialPointsDataFrame"){
#         df <- cbind(coordinates(pts),pts@data)
#         df[,1] <- jitter(df[,1],amount=xamount)
#         df[,2] <- jitter(df[,2],amount=yamount)
#         coordinates(df) <- df[,1:2]
#         proj4string(df) <- proj
#     } else if (class(pts) == "SpatialPoints"){
#         df <- coordinates(pts)
#         df[,1] <- jitter(df[,1],amount=xamount)
#         df[,2] <- jitter(df[,2],amount=yamount)
#         df <- SpatialPoints(df, proj4string=CRS(proj))
#     } else {
#         stop("Only works for SpatialPoints* at present.")
#     }
#     return(df)
# }
# 
# 
