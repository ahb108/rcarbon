#' @title Binning function of radiocarbon dates.  
#'
#' @description Prepare a set of bins for controlling the aggregation of radiocarbon dates
#' known to be from the same phase of same archaeological site (for use with \code{\link{spd}}). Used in cases where there is a concern that unusually high levels of sampling for radiocarbon at a given site or in a given site phase will impede comparison between sites or phases. 
#' 
#' @param sites a vector of character strings (or number to coerce to character) of all sites or site phases. If character strings are used these should not contain underscores (see also below)
#' @param ages a vector of uncalibrated conventional radiocarbon ages or a \code{CalDates} class object obtained using the \code{\link{calibrate}} function.
#' @param h a single numeric value passed to \code{\link{hclust}} control degree of grouping of similar ages in a phase site.
#'
#' @details If \code{ages} is a \code{CalDates} class object, median dates are used for the clustering.
#'
#' @return A vector of character strings with the same length of the object supplied for the argument \code{ages} identifying intra-site or intra-phase grouping, for use with \code{\link{spd}}.The character strings effectively provide a "name" for each "phase" within a "site", using sequential integers after an underscore. For example if a site named "S001" had four dates grouped into two bins with two dates each, the resulting vector would be "S001_1", "S001_1", "S001_2", and "S001_2".  

#' @seealso \code{\link{spd}} for generating SPD; \code{\link{binsense}} for sensitivity analysis pertaining the choice of the parameter \code{h}.
#' @import stats
#' @import utils
#' @export
binPrep <- function(sites, ages, h){
    
    if(any(grepl("_", sites)))
    {
        stop("sites should not contain elements with underscores")	
    }
    
    clusters <- rep(NA,length(sites))
    
    if (any(class(ages)%in%"CalDates"))
    {
	    ages = medCal(ages)
    }

    for  (x in 1:length(unique(sites))){
        index <- which(sites==unique(sites)[x])
        if (length(index)>1){
            clusters[index] <- paste(unique(sites)[x],cutree(hclust(dist(ages[index])),h=h),sep="_")
        }
        if (length(index)==1){
            clusters[index] <- paste(unique(sites)[x],"1",sep="_")
        }                
    }
    
    return(clusters)
}

#' @title Sampling function to select a maximum number of dates per site, bin or phase.
#'
#' @description Function to select a subset of uncalibrated radiocarbon dates up to a maximum sample size per site, bin or phase.
#' 
#' @param ages A vector of uncalibrated radiocarbon ages
#' @param errors A vector of uncalibrated radiocarbon errors (same length as ages)
#' @param bins A vector of labels corresponding to site names, ids, bins or phases (same length as ages)
#' @param size A single integer specifying the maximum number of desired dates for each label stated bin.
#' @param errorcap A single integer specifying the maximum allowable error (applied only in cases where a bin has more than the maximum number of allowable dates).
#' @param thresh A single numeric value between 0 and 1 specifying the approximate proportion (after rounding) of the resulting sample that will be chosen according to lowest date errors. At the extremes, O produces a simple random sample whereas 1 selectes the sample dates with the lowest errors. Ignored if method="random".
#' @param method The method to be applied where "random" simple selects a random sample, whereas "splitsample", picks some proportion (see thresh) of the sample to minimise errors, and randomly samples the rest. At present, these are the only two options.
#' @param seed Allows setting of a random seed to ensure reproducibility.
#'
#' @return A numeric vector of the row indices corresponding to those of the input data.
#' @examples
#'data(euroevol)
#'foursites <- euroevol[euroevol$SiteID %in% c("S2072","S4380","S6139","S9222"),]
#'table(as.character(foursites$SiteID))
#'## Thin so each site has 10 dates each max, with random selection
#'thinInds<- thinDates(ages=foursites$C14Age, errors=foursites$C14SD, bins=foursites$SiteID, size=10, method="random", seed=123)
#'tdates <- foursites[thinInds,]
#'tdates
#'## Same but choose the first 60% (i.e. 6 dates) from the lowest errors and then fill in the rest randomly.
#'thinInds<- thinDates(ages=foursites$C14Age, errors=foursites$C14SD, bins=foursites$SiteID, size=10, method="splitsample", thresh=0.6, seed=123)
#'tdates1 <- foursites[thinInds,]
#'tdates1
#' @seealso \code{\link{binPrep}}
#' @export
#'
thinDates <- function(ages, errors, bins, size, errorcap=NA, thresh=0.5, method="random", seed=NA){
    if (length(size)!=1){ stop("The size argument must be a single integer.") }
    if (thresh < 0 |  thresh > 1){ stop("The thresh argument must be between 0 and 1.") }
    df <- data.frame(ages=ages, errors=errors, bins=bins, RN=1:length(ages))
    allbins <- unique(as.character(bins))
    sitelist <- vector(mode="list", length=length(allbins))
    for (a in 1:length(allbins)){
        df1 <- df[df$bins==allbins[a],]
        if (nrow(df1)<=size){
            sitelist[[a]] <- df1$RN
        } else {
            if (method=="random"){
                if (!is.na(errorcap)){ df1 <- df1[df1$errors <= errorcap,] }
                if (!is.na(seed) & is.numeric(seed)){ set.seed(seed) }
                sitelist[[a]] <- sample(df1$RN,size)
            } else if (method=="splitsample"){
                if (!is.na(errorcap)){ df1 <- df1[df1$errors <= errorcap,] }
                rnks <- rank(df1$errors)
                check <- min(rnks[rnks>(size*thresh)])
                myinds <- which(rnks<=check)
                if (!is.na(seed) & is.numeric(seed)){ set.seed(seed) }
                res <- sample(df1$RN[myinds],size*thresh)
                if (length(res)<size){
                    inds <- df1$RN
                    if (!is.na(seed) & is.numeric(seed)){ set.seed(seed) }
                    therest <- sample(inds[!inds %in% res],size*(1-thresh))
                    sitelist[[a]] <- c(res,therest)
                } else {
                    sitelist[[a]] <- res
                }
            } else {
                stop("The method argument must currently be one of random or splitsample")
            }
        }
    }
    allres <- sort(unlist(sitelist))
    return(allres)
}

#' @title Summed probability distributions (SPD) of radiocarbon dates.  
#'
#' @description The function generates Summed probability distributions (SPD) of radiocarbon dates, with optional binning routine for controlling inter-site or inter-phase variation in sampling intensity.
#'
#' @param x A \code{CalDates} class object containing the calibrated radiocarbon dates.
#' @param timeRange A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param bins A vector containing the bin names associated with each radiocarbon date. If set to NA, binning is not carried out. 
#' @param datenormalised Controls for calibrated dates with probability mass outside the timerange of analysis. If set to TRUE the total probability mass within the time-span of analysis is normalised to sum to unity. Should be set to FALSE when the parameter \code{normalised} in \code{\link{calibrate}} is set to FALSE. Default is FALSE. 
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity. 
#' @param runm A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied. Default is NA  
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#'
#' @details The binning routine consists of computing summed probability distribution of all dates associated to a given bin, divided by the number of contributing dates. This controls for any striking differences in sampling intensity, and ensures that each site phase is equally contributing to the final SPD (see Timpson et al 2014 for details). Bins can be generated using the \code{\link{binPrep}}, whilst the sensitivity to parameter choice can be explored with the \code{\link{binsense}} function.
#'
#' @return An object of class \code{CalSPD} with the following elements
#' \itemize{
#' \item{\code{metadata}} {A data.frame containing relevant information regarding the parameters used to create the SPD as well as sample size and number of bins}
#' \item{\code{grid}} {A \code{CalGrid} class object containing the summed probability associated to each calendar year between \code{timeRange[1]} and \code{timeRange[2]}}
#'}
#'
#' @references 
#' Timpson, A., et al, (2014). Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method. Journal of Archaeological Science 52: 549-557. DOI:10.1016/j.jas.2014.08.011
#'
#' @seealso \code{\link{calibrate}} for calibrating radiocarbon dates; \code{\link{binPrep}} for preparing bins.
#' @import utils
#' @export

spd <- function(x, timeRange, bins=NA, datenormalised=FALSE, spdnormalised=FALSE, runm=NA, verbose=TRUE, edgeSize=500){
    
    defcall <- as.list(args(spd))
    defcall <- defcall[-length(defcall)]
    speccall <- as.list(match.call())
    speccall <- speccall[-1]
    i <- match(names(defcall), names(speccall))
    i <- is.na(i)
    if (any(i)){
        speccall[names(defcall)[which(i)]] <- defcall[which(i)]
    }
    speccall <- as.data.frame(lapply(speccall,deparse), stringsAsFactors=FALSE)
    speccall <- speccall[,names(defcall)]
    speccall$ndates <- nrow(x$metadata)
    speccall$nbins <- nrow(x$metadata)
    if (!"CalDates" %in% class(x)){
        stop("x must be an object of class 'CalDates'.")
    }
    if (length(bins)>1){
        speccall$nbins <- length(unique(bins))
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

    if(datenormalised)
    {
	    true.timeRange=timeRange
	    ccrange=range(medCal(x))
	    ccrange[2]=ccrange[2]+edgeSize
	    ccrange[1]=ccrange[1]-edgeSize
	    if (ccrange[1]<0|ccrange[2]>50000)
	    {
		stop("timeRange beyond calibration curve. Ensure that timeRange[1]+edgeSize is smaller than 50000 and timeRange[2]-edgeSize is larger than 0")
	    }
	    timeRange[1]=ccrange[2]
	    timeRange[2]=ccrange[1]
    }	
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
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
    check <- caldateTR[1] >= timeRange[1] & caldateTR[2] <= timeRange[2]
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
    }
    if (verbose & length(binNames)>1){ close(pb) }

    finalSPD <- apply(binnedMatrix,1,sum)

    if (datenormalised)
    {
	    timeRange=true.timeRange
    }	    

    
    if (!is.na(runm)){
        finalSPD <- runMean(finalSPD, runm, edge="fill")
    }
    res <- data.frame(calBP=calyears$calBP, PrDens=finalSPD)
    res <- res[res$calBP <= timeRange[1] & res$calBP >= timeRange[2],]
    if (spdnormalised){
        res$PrDens <- res$PrDens/sum(res$PrDens, na.rm=TRUE)
    }
    class(res) <- c("CalGrid", class(res))
    reslist <- vector("list",length=2)
    names(reslist) <- c("metadata","grid")
    reslist[["metadata"]] <- speccall
    reslist[["grid"]] <- res
    class(reslist) <- c("CalSPD",class(reslist))
    if (verbose){ print("Done.") }
    return(reslist)
}
