 if(getRversion() >= "2.15.1")  utils::globalVariables(c("focalyear"))


#' @title Binning function of radiocarbon dates.  
#'
#' @description Prepare a set of bins for controlling the aggregation of radiocarbon dates
#' known to be from the same phase of same archaeological site (for use with \code{\link{spd}}). Used in cases where there is a concern that unusually high levels of sampling for radiocarbon at a given site or in a given site phase will impede comparison between sites or phases. 
#' 
#' @param sites a vector of character strings (or number to coerce to character) of all sites or site phases. If character strings are used these should not contain underscores (see also below)
#' @param ages a vector of uncalibrated conventional radiocarbon ages or a \code{CalDates} class object obtained using the \code{\link{calibrate}} function.
#' @param h a single numeric value passed to \code{\link{cutree}} control degree of grouping of similar ages in a phase site.
#' @param method the agglomeration method to be used, passed on to \code{\link{hclust}}. Defaults to "complete" as in \code{\link{hclust}}.
#'
#' @details If \code{ages} is a \code{CalDates} class object, median dates are used for the clustering.
#'
#' @return A vector of character strings with the same length of the object supplied for the argument \code{ages} identifying intra-site or intra-phase grouping, for use with \code{\link{spd}}.The character strings effectively provide a "name" for each "phase" within a "site", using sequential integers after an underscore. For example if a site named "S001" had four dates grouped into two bins with two dates each, the resulting vector would be "S001_1", "S001_1", "S001_2", and "S001_2".  

#' @seealso \code{\link{spd}} for generating SPD; \code{\link{binsense}} for sensitivity analysis pertaining the choice of the parameter \code{h}.
#' @import stats
#' @import utils
#' @export
binPrep <- function(sites, ages, h, method = "complete"){
    
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
            clusters[index] <- paste(unique(sites)[x],cutree(hclust(dist(ages[index]), method = method),h=h),sep="_")
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
#' @param thresh A single numeric value between 0 and 1 specifying the approximate proportion (after rounding) of the resulting sample that will be chosen according to lowest date errors. At the extremes, O produces a simple random sample whereas 1 selects the sample dates with the lowest errors. Ignored if method="random".
#' @param method The method to be applied where "random" simple selects a random sample, whereas "splitsample", picks some proportion (see thresh) of the sample to minimise errors, and randomly samples the rest. At present, these are the only two options.
#' @param seed Allows setting of a random seed to ensure reproducibility.
#'
#' @return A numeric vector of the row indices corresponding to those of the input data.
#' @examples
#'data(euroevol)
#'foursites <- euroevol[euroevol$SiteID %in% c("S2072","S4380","S6139","S9222"),]
#'table(as.character(foursites$SiteID))
#'## Thin so each site has 10 dates each max, with random selection
#'thinInds<- thinDates(ages=foursites$C14Age, errors=foursites$C14SD, 
#' bins=foursites$SiteID, size=10, method="random", seed=123)
#'tdates <- foursites[thinInds,]
#'tdates
#'## Same but choose the first 60% (i.e. 6 dates) from the lowest errors 
#'## and then fill in the rest randomly.
#'thinInds<- thinDates(ages=foursites$C14Age, errors=foursites$C14SD, 
#'bins=foursites$SiteID, size=10, method="splitsample", thresh=0.6, seed=123)
#'tdates1 <- foursites[thinInds,]
#'tdates1
#' @seealso \code{\link{binPrep}}
#' @export
#'
thinDates <- function(ages, errors, bins, size, thresh=0.5, method="random", seed=NA){
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
                if (!is.na(seed) & is.numeric(seed)){ set.seed(seed) }
                sitelist[[a]] <- sample(df1$RN,size)
            } else if (method=="splitsample"){
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
#' @param edgeSize Extra margin in C14 Age time to handle edge effect when \code{datenormalise} is set to TRUE. Default is 500.
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

spd <- function(x,timeRange, bins=NA, datenormalised=FALSE, spdnormalised=FALSE, runm=NA, verbose=TRUE, edgeSize=500){
    
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
    caltimeRange =c(55000,0)
    if (any(x$metadata$CalCurve %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")))
    {
      caltimeRange =c(50000,0)
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
	    if (ccrange[1]<0|ccrange[2]>caltimeRange[1])
	    {
		stop(paste0("timeRange beyond calibration curve. Ensure that timeRange[1]+edgeSize is smaller than ", caltimeRange[1]," and timeRange[2]-edgeSize is larger than 0"))
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

#' @title Stacked Summed Probability Distribution
#' @description Generates and combines multiple SPDs based on a user defined grouping.
#' @param x A \code{CalDates} class object containing the calibrated radiocarbon dates.
#' @param timeRange A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param bins A vector containing the bin names associated with each radiocarbon date. If set to NA, binning is not carried out. 
#' @param group A vector containing the grouping variable.
#' @param datenormalised Controls for calibrated dates with probability mass outside the timerange of analysis. If set to TRUE the total probability mass within the time-span of analysis is normalised to sum to unity. Should be set to FALSE when the parameter \code{normalised} in \code{\link{calibrate}} is set to FALSE. Default is FALSE. 
#' @param runm A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied. Default is NA  
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param edgeSize Controls edge effect by expanding the fitted model beyond the range defined by \code{timeRange}.
#'
#' @return An object of class \code{stackCalSPD} 
#' @examples
#' \dontrun{
#'data(emedyd)
#'x = calibrate(x=emedyd$CRA, errors=emedyd$Error,normalised=FALSE)
#'bins = binPrep(sites=emedyd$SiteName, ages=emedyd$CRA,h=50)
#'res = stackspd(x=x,timeRange=c(16000,8000),bins=bins,group=emedyd$Region)
#'}
#' @import utils
#' @export

stackspd <- function(x, timeRange, bins=NA, group=NULL, datenormalised=FALSE, runm=NA, verbose=TRUE, edgeSize=500){

	if (is.null(group))
	{
		stop("The argument group must be provided")
	}

	stackLength = length(unique(group))


	## Main for Loop

	stackG = unique(group)
	stackL = length(stackG)
	stackedSPDs = vector("list",length=stackL)
	names(stackedSPDs) = stackG

	for (i in 1:stackL)
	{
		if (verbose)
		{
			print(paste0("Computing SPD for group ", i, " out of ", stackL))
		}
	k = which(group==stackG[i])	
	b = bins[k]
	if (anyNA(b)){b=NA}

	stackedSPDs[[i]]=spd(x[k], timeRange=timeRange, bins=b, datenormalised=datenormalised, spdnormalised=FALSE, runm=runm, verbose=verbose, edgeSize=edgeSize)
	}

	## Create General Metadata
	metadata = list(ndates=length(x),ngroups=stackL,timeRange=timeRange)
	
	## Final List
	result = list(metadata=metadata,spds=stackedSPDs)

	## Define Class
	class(result) <- c("stackCalSPD",class(result))
	return(result)
}

#' @title Composite Kernel Density Estimates of Radiocarbon Dates
#'
#' @description Computes a Composite Kernel Density Estimate (CKDE) from multiple sets of randomly sampled calendar dates.
#' @param x A \code{simdates} class object, generated using \code{\link{sampleDates}}.
#' @param timeRange A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param bw Kernel bandwidth to be used.
#' @param normalised A logical variable indicating whether the contribution of individual dates should be equal (TRUE), or weighted based on the area under the curve of non-normalised calibration (FALSE). Default is TRUE.
#' @details The function computes Kernel Density Estimates using randomly sampled calendar dates contained in a \code{simdates} class object (generated using the \code{simulate.dates()} function). The output contains \code{nsim} KDEs, where \code{nsim} is the argument used in \code{simulate.dates()}. The resulting object can be plotted to visualise a CKDE (cf Brown 2017), and if \code{boot} was set to \code{TRUE} in \code{sampleDates} its bootstrapped variant (cf McLaughlin 2018 for a similar analysis). The shape of the CKDE is comparable to an SPD generated from non-normalised dates when the argument \code{normalised} is set to FALSE.
#' @return An object of class \code{ckdeSPD} with the following elements
#' \itemize{
#' \item{\code{timeRange}} {The \code{timeRange} setting used.}
#' \item{\code{res.matrix}} {A matrix containing the KDE values with rows representing calendar dates.}
#'}
#'
#' @references 
#' Brown, W. A. 2017. The past and future of growth rate estimation in demographic temporal frequency analysis: Biodemographic interpretability and the ascendance of dynamic growth models. \emph{Journal of Archaeological Science}, 80: 96–108. DOI: https://doi.org/10.1016/j.jas.2017.02.003 \cr
#' McLaughlin, T. R. 2018. On Applications of Space–Time Modelling with Open-Source 14C Age Calibration. \emph{Journal of Archaeological Method and Theory}. DOI https://doi.org/10.1007/s10816-018-9381-3
#'
#' @examples
#' \dontrun{
#'data(emedyd)
#'x = calibrate(x=emedyd$CRA, errors=emedyd$Error,normalised=FALSE)
#'bins = binPrep(sites=emedyd$SiteName, ages=emedyd$CRA,h=50)
#'s = sampleDates(x,bins=bins,nsim=100,boot=FALSE)
#'ckdeNorm = ckde(s,timeRange=c(16000,9000),bw=100,normalised=TRUE)
#'plot(ckdeNorm,type='multiline',calendar='BCAD')
#' }
#'
#' @seealso \code{\link{sampleDates}}
#'
#' @import stats
#' @import utils
#' @export

ckde<- function(x,timeRange,bw,normalised=FALSE)
{

	if (!"simdates" %in% class(x)){
		stop("x must be an object of class 'simdates'.")
	}
	true.timeRange = c(max(x$sdates),min(x$sdates))
	nsim = nrow(x$sdates)
	boot = FALSE
	if (any(class(x$weight)%in%c('matrix','array'))){boot=TRUE}

	if (normalised)
	{
		raw.kde=apply(x$sdates,1,density,bw=bw,na.rm=TRUE,from=true.timeRange[1],to=true.timeRange[2])
	}
	if (!normalised)
	{
		if (length(unique(x$weight))==1)
		{
			warning("Simulated dates were generated from a dates calibrated with the argument `normalised` set to TRUE. The composite KDE will be executed with 'normalised' set to TRUE. To obtain a non-normalised composite KDE calibrate setting 'normalised` to FALSE")
		}
		if (!boot)
		{
			raw.kde=apply(x$sdates,1,density,bw=bw,na.rm=TRUE,from=true.timeRange[1],to=true.timeRange[2],weights=x$weight)
		} else {
			raw.kde = lapply(1:nrow(x$sdates),function(x,y,z,t1,t2,bw){
						 return(density(y[x,],bw=bw,na.rm=TRUE,weights=z[x,],from=t1,to=t2))},y=x$sdates,z=x$weight,t1=true.timeRange[1],t2=true.timeRange[2],bw=bw)
		}
	}
	res.matrix = matrix(NA,nrow=length(timeRange[1]:timeRange[2]),ncol=nsim)
	for (i in 1:nsim)
	{
		#check this line
		res.matrix[,i]=approx(x=raw.kde[[i]]$x,xout=timeRange[1]:timeRange[2],y=raw.kde[[i]]$y)$y
	}

	result = list(timeRange=timeRange,res.matrix=res.matrix)
	class(result) = c("compositeKDE",class(result))
	return(result)
}


#' @title Map the spatio-temporal intensity of a set of radiocarbon dates
#'
#' @description Function for mapping the spatio-temporal intensity of radiocarbon dates for a given geographical region in one or more tim.
#' @param x An object of class CalDates with calibrated radiocarbon ages.
#' @param coords A two column matrix of geographical coordinates from a a projected coordinate system (no checks are made for this) and with the same number of rows as length(x).
#' @param sbw A single numeric value for the spatial bandwidth to be applied around each raster cell, expressed as the standard deviation of a continuous Gaussian kernel (passed as the sigma argument to density.ppp()).
#' @param focalyears A vector of numeric values for focal years, in calBP, that will be timesteps at which date intensity maps will be produced.
#' @param tbw A single numeric value for the temporal bandwidth to be applied around each focal year, expressed as the standard deviation of a continuous Gaussian kernel.
#' @param win The bounding polygon for the mapping (must be an object of class 'owin', see the spatstat package)
#' @param cellres The cell or pixel resolution of the output raster maps.
#' @param outdir The output directory for timeslice maps and data that are saved to file.
#' @param bins  A vector of labels corresponding to site names, ids, bins or phases (same length as x)
#' @param backsight A single numeric value (which will be coerced to be positive) that specifies a comparison timestep in the past for a mapping of temporal change.
#' @param maskthresh A single numeric value for a lower-bound cut-off for all maps, based on a minimum required spatial intensity of all dates in x.
#' @param changexpr An expression for calculating the change in spatial intensity between the focal year and a backsight year (as defined via the backsight argument). Available input options are t1 (the spatial intensity for the focal year), t0 (the spatial intensity for the backsight year) and tk (the overall spatial intensity for all dates irrespective of year), plus any other standard constants and mathematical operators. A sensible default is provided.
#' @param spjitter Whether noise is applied to the spatial coordinates or not. Default is TRUE. 
#' @param amount Amount of jitter applied to the spatial coordinates when \code{spjitter=TRUE}. Default is d/5, where d is difference between the closest coordinates.    
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param ... ignored or passed to internal functions.
#'
#' @details This function computes one or more timeslice maps of the spatio-temporal kernel intensity of radiocarbon dates across a geographic region and for a specific focal year. The user specifies the arbitrary sizes of both the spatial and the temporal Gaussian kernels that will be used to summarise radiocarbon date intensity per grid cell per timestep. The results allow standardisation of colour ramps, etc. across timesteps and are amenable to plotting individually via plot.stKde and/or for output to png for animation.
#'
#' @return A list object of class stKde with the following elements:
#' \itemize{
#' \item {A series of list items storing some of the input parameters such as the focalyears, sbw, tbw, backsight, maskthresh} 
#' \item{\code{nonfocal}} {An im object mapping the basic spatial intensity of all dates, without reference to a focal year.} 
#' \item{\code{impaths}} {A character vector of the paths to the individual timeslices stored on file. Maps are not stored in memory (see spkde() for further details of what is stored).}
#' \item{\code{stats}} {A list of data.frames offering summary statistics on each of the different types of output surface across all timeslices. Used primarily to allow consistent colour ramps across time-slices.}
#' \item{\code{ppp}} {The ppp object for all dates and the observation window.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Example with a subset of English and Welsh dates from the Euroevol dataset
#' data(ewdates)
#' data(ewowin)
#' x <- calibrate(x=ewdates$C14Age, errors=ewdates$C14SD, normalised=FALSE)
#' ## Create centennial timeslices (also with site binning)
#' bins1 <- binPrep(sites=ewdates$SiteID, ages=ewdates$C14Age, h=50)
#' stkde1 <- stkde(x=x, coords=ewdates[,c("Eastings", "Northings")], win=ewowin, 
#' sbw=40000, cellres=2000, focalyears=seq(6500, 5000, -100), tbw=50, bins=bins1, 
#' backsight=200, outdir="im",amount=1)
#' ## Plot an example of all four basic outputs for 5900 calBP
#' dev.new(height=2.5, width=8)
#' par(mar=c(0.5, 0.5, 2.5, 2))
#' plot(stkde1, 5900, type="all")
#' }
#' 
#' @import spatstat 
#' @export
#' 
stkde <- function(x, coords, sbw, focalyears, tbw, win, cellres, outdir=".", bins=NA, backsight=NA, maskthresh=0, changexpr=expression((t1-t0)/tk),spjitter=TRUE,amount=NULL, verbose=TRUE, ...){
    if (!"CalDates" %in% class(x)){ stop("x must be of class CalDates.") }
    outpaths <- vector(mode="character", length=length(focalyears))
    forsumstats <- c("focal", "proportion", "change")
    statslist <- vector(mode="list", length=length(forsumstats))
    names(statslist) <- forsumstats
    df0 <- as.data.frame(matrix(ncol=10, nrow=length(focalyears)), stringsAsFactors=FALSE)
    names(df0) <- c("year", "min", "2.5", "25", "median", "75", "97.5", "max", "mean", "sd")
    for (g in 1:length(forsumstats)){ statslist[[g]] <- df0 }
    if (verbose){
        print("Mapping focal intensities and then writing to file...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(focalyears), style=3)
    }
    backsight <- abs(backsight)

    dpoints = FALSE
    if (spjitter==FALSE&anyDuplicated(coords)){dpoints=TRUE}
    if (spjitter|dpoints)
    {
	coords[,1]=jitter(coords[,1],amount=amount)
    	coords[,2]=jitter(coords[,2],amount=amount)
	if(dpoints){warning("Duplicated coordinates: spkde() has been executed with spjitter set to TRUE")}
    }	    

    if (any(inside.owin(coords[,1],coords[,2],win)==FALSE))
    {
	in.pts=inside.owin(coords[,1],coords[,2],win)    
        x=x[which(in.pts)]
	coords=coords[which(in.pts),]
	bins=bins[which(in.pts)]
	plural=" points were "
	if (sum(!in.pts)==1){plural=" point was "}
	warning(paste0(sum(!in.pts),plural,"rejected as lying outside the window of analysis. Consider using a smaller setting for the agrument 'amount'"))
    }
    
    for (a in 1:length(focalyears)){
        if (verbose){ setTxtProgressBar(pb, a) }
        focalyear <- spkde(x=x, coords=coords, win=win, sbw=sbw, cellres=cellres, focalyear=focalyears[a], tbw=tbw, bins=bins, backsight=backsight, verbose=FALSE, maskthresh=0,spjitter=FALSE)
        outpath <- paste(outdir, "/",focalyears[a],".rda", sep="")
        for (d in 1:length(names(statslist))){
            thisone <- names(statslist)[d]
            q <- quantile(focalyear[[thisone]], probs=c(0,0.025,0.25,0.5,0.75,0.975,1))
            statslist[[names(statslist)[d]]][a,] <- c(focalyears[a], q, mean(focalyear[[thisone]]), sd(focalyear[[thisone]]))
        }
        outpaths[a] <- outpath
        save(focalyear, file=outpath)
    }
    if (verbose){ close(pb) }
    resnames <- c("years", "sbw", "tbw", "backsight", "maskthresh", "nonfocal", "impaths", "stats")
    res <- vector(mode="list", length=length(resnames))
    names(res) <- resnames
    res[["focalyears"]] <- focalyears
    pppA <- ppp(coords[,1],coords[,2], window=win)
    res[["ppp"]] <- pppA
    res[["sbw"]] <- sbw
    res[["tbw"]] <- tbw
    res[["backsight"]] <- backsight
    res[["maskthresh"]] <- maskthresh
    ppAwts <- sapply(x$grids, function(x){ sum(x$PrDens) })
    alldens <- density(pppA, weights=ppAwts, eps=cellres, sigma=sbw)
    res[["nonfocal"]] <- alldens
    res[["impaths"]] <- outpaths
    names(res[["impaths"]]) <- focalyears
    res[["stats"]] <- statslist
    class(res) <- c("stKde", class(res)) 
    if (verbose){ print("Done.") }
    return(res)
}

#' @title Map the spatial intensity of a set of radiocarbon dates for a given focal year.
#'
#' @description Function for mapping the spatial intensity of radiocarbon dates for a given geographical region and focal year.
#'
#' @param x An object of class CalDates with calibrated radiocarbon ages.
#' @param coords A two column matrix of geographical coordinates from a a projected coordinate system (no checks are made for this) and with the same number of rows as length(x).
#' @param sbw A single numeric value for the spatial bandwidth to be applied around each raster cell, expressed as the standard deviation of a continuous Gaussian kernel (passed as the sigma argument to density.ppp()).
#' @param focalyear A single numeric value for the focal year for the intensity map.
#' @param tbw A single numeric value for the temporal bandwidth to be applied around each focal year, expressed as the standard deviation of a continuous Gaussian kernel.
#' @param win The bounding polygon for the mapping (must be an object of class 'owin', see the spatstat package)
#' @param cellres The cell or pixel resolution of the output raster maps.
#' @param bins  A vector of labels corresponding to site names, ids, bins or phases (same length as x)
#' @param backsight A single numeric value (which will be coerced to be positive) that specifies a comparison timestep in the past for a mapping of temporal change.
#' @param nsim How many bootstrap simulations to run (default is none).
#' @param maskthresh A single numeric value for a lower-bound cut-off for all maps, based on a minimum required spatial intensity of all dates in x.
#' @param changexpr An expression for calculating the change in spatial intensity between the focal year and a backsight year (as defined via the backsight argument). Available input options are t1 (the spatial intensity for the focal year), t0 (the spatial intensity for the backsight year) and tk (the overall spatial intensity for all dates irrespective of year), plus any other standard constants and mathematical operators. A sensible default is provided.
#' @param raw Whether to output the raw simulations (if nsim is set) or just the summaries (the latter is default).
#' @param spjitter Whether noise is applied to the spatial coordinates or not. Default is TRUE. 
#' @param amount Amount of jitter applied to the spatial coordinates when \code{spjitter=TRUE}. Default is d/5, where d is difference between the closest coordinates.   
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param ... ignored or passed to internal functions.
#'
#' @details This function is not really intended for general use, but rather as an internal function for stkde(). Most applications should use the latter, but spkde has been exported and made externally available, both because this function retains the result in memory (in contrast to stkde) and with a view to possible addition of bootstrap methods in the future. Some function arguments therefore remain experimental. The function computes timeslice maps of the spatio-temporal kernel intensity of radiocarbon dates across a geographic region for a specific focal year. The user specifies the arbitrary size of both the spatial and the temporal Gaussian kernels that will be used to summarise radiocarbon date intensity per grid cell.
#'
#' @return A list object of class spKde with the following elements:
#' \itemize{
#' \item {A series of list items storing some of the input parameters such as the focal year, sbw, tbw, backsight, maskthresh.} 
#' \item{\code{nonfocal}} {An im object mapping the basic spatial intensity of all dates, without reference to a focal year.} 
#' \item{\code{focal}} {An im object mapping the spatial intensity of dates for the focal year (i.e. weighted by how much each dates probability distribution overlaps with a Gaussian kernel centred on the focal year with a standard deviation of tbw).}
#' \item{\code{proportion}} {An im object mapping the proportional intensity of dates for the focal year (i.e. the focal surface divided by the nonfocal surface).}
#' \item{\code{change}} {An im object mapping the amount of change between the intensity of dates for the focal year and a chosen backsight year (i.e. as defined by changexpr).}
#' }
#'
#' @examples
#' \dontrun{
#' ## Example for the focal year 5600 calBP (also with site binning), 
#' ## using a subset of English and Welsh dates from the Euroevol dataset
#' data(ewdates)
#' data(ewowin)
#' x <- calibrate(x=ewdates$C14Age, errors=ewdates$C14SD, normalised=FALSE)
#' bins1 <- binPrep(sites=ewdates$SiteID, ages=ewdates$C14Age, h=50)
#' spkde1 <- spkde(x=x, coords=ewdates[,c("Eastings", "Northings")], win=ewowin, 
#' sbw=40000, cellres=2000, focalyear=5600, tbw=50, bins=bins1, backsight=200,amount=1)
#' plot(spkde1$focal)
#' plot(spkde1$proportion)
#' }
#' 
#' @import spatstat 
#' @export
#' 
spkde <- function(x, coords, sbw, focalyear, tbw, win, cellres, bins=NA, backsight=NA, nsim=NULL, maskthresh=0, changexpr=expression((t1-t0)/tk), raw=FALSE, spjitter=TRUE,amount=NULL,verbose=TRUE, ...){
    if (!"CalDates" %in% class(x)){ stop("x must be of class CalDates.") }
    basic <- c("nsim", "year", "backsight", "sbw", "tbw", "maskthresh")
    coreres <- c("nonfocal", "focal", "proportion", "change")
    extrares <- append(coreres,c("samplemean", "samplestdev", "sampleprop", "samplechange"))
    if (!is.null(nsim)){
        resnames <- append(basic,extrares)
    } else {
        resnames <- append(basic,coreres)
    }   
    if (is.null(nsim) & verbose){ print("Mapping spatial intensities...") }
    
    dpoints = FALSE
    if (spjitter==FALSE&anyDuplicated(coords)){dpoints=TRUE}
    if (spjitter|dpoints)
    {
	coords[,1]=jitter(coords[,1],amount=amount)
    	coords[,2]=jitter(coords[,2],amount=amount)
	if(dpoints){warning("Duplicated coordinates: spkde() has been executed with spjitter set to TRUE")}
    }	    
    
    if (any(inside.owin(coords[,1],coords[,2],win)==FALSE))
    {
	in.pts=inside.owin(coords[,1],coords[,2],win)    
        x=x[which(in.pts)]
	coords=coords[which(in.pts),]
	bins=bins[which(in.pts)]
	plural=" points were "
	if (sum(!in.pts)==1){plural=" point was "}
	warning(paste0(sum(!in.pts),plural,"rejected as lying outside the window of analysis. Consider using a smaller setting for the agrument 'amount'"))
    }

    pppA <- ppp(coords[,1],coords[,2], window=win)

    ppAwts <- sapply(x$grids, function(x){ sum(x$PrDens) })
    if (!is.na(bins)[1]){
        bindf0 <- data.frame(origrows=1:length(ppAwts), bins=bins, auc=ppAwts, stringsAsFactors=FALSE)
        bindf1 <- as.data.frame(table(bins), stringsAsFactors=FALSE)
        names(bindf1) <- c("bins","ndates")
        bindf <- merge(bindf0, bindf1, by="bins", all.x=TRUE)
        bindf$finalwt <- bindf$auc / bindf$ndates
        bindf <- bindf[with(bindf, order(origrows)), ]
        ppAwts <- bindf$finalwt
    }
    tk <- density(pppA, weights=ppAwts, eps=cellres, sigma=sbw)
    pppAwtsg <- gaussW(x, focalyear, tbw)
    if (!is.na(bins)[1]){ pppAwtsg <- pppAwtsg / bindf$ndates }
    t1 <- density(pppA, weights=pppAwtsg, eps=cellres, sigma=sbw)
    if (!is.null(nsim)){
        s1 <- sampleDates(x, nsim=nsim, bins=bins, verbose=verbose, ...)
        imlist <- vector(mode="list", length=nsim)
        win <- as.owin(win)
        coords <- as.matrix(coords)
        if (!is.na(backsight)){
            note <- "backsight must be NA or a single positive integer"
            if (length(backsight)!=1){ stop(note) }
            if (!is.numeric(backsight)){ stop(note) }
            if (backsight<=0){ stop(note) }
            bsBP <- focalyear + backsight
            bimlist <- vector(mode="list", length=nsim)
        }
        if (nsim>1 & verbose){
            print("Mapping spatial intensities...")
            flush.console()
            pb <- txtProgressBar(min=1, max=nsim, style=3)
        }
        for (a in 1:nsim){
            if (nsim > 1 & verbose){ setTxtProgressBar(pb, a) }
            if (is.matrix(s1$sdates)){
                fset <- cbind(s1$sdates[a,], s1$weight, coords[,1], coords[,2])
            } else {
                fset <- cbind(s1$sdates, s1$weight, coords[,1], coords[,2])
            }
            fWt <- gaussW(fset[,1], focalyear, tbw)
            cWt <- fset[,2] * fWt
            if (!is.na(bins)[1]){ cWt <- cWt / bindf$ndates }
            pppF <- ppp(fset[,3],fset[,4], window=win)
            imlist[[a]] <- density(pppF, weights=cWt, eps=cellres, sigma=sbw)
            if (!is.na(backsight)){   
                fWt <- gaussW(fset[,1], bsBP, tbw)
                cWt <- fset[,2] * fWt
                if (!is.na(bins)[1]){ cWt <- cWt / bindf$ndates }
                pppF <- ppp(fset[,3],fset[,4], window=win)
                bimlist[[a]] <- density(pppF, weights=cWt, eps=cellres, sigma=sbw)
            }
        }
        if (nsim > 1 & verbose){ close(pb) }
    }
    if (raw & !is.null(nsim)){
        if (verbose){ print("Done.") }
        return(imlist)
    } else {
        res <- vector(mode="list", length=length(resnames))
        names(res) <- resnames
        res[["year"]] <- focalyear
        res[["backsight"]] <- backsight
        res[["sbw"]] <- sbw
        res[["tbw"]] <- tbw
        res[["maskthresh"]] <- maskthresh
        res[["nonfocal"]] <- tk
        t1[as.matrix(tk) < maskthresh] <- NA
        res[["focal"]] <- t1
        tmp <- t1 / tk
        tmp[as.matrix(tk) < maskthresh] <- NA
        res[["proportion"]] <- tmp
        if (!is.null(nsim)){
            res[["nsim"]] <- nsim
            t1 <- im.apply(imlist, mean)
            t1[as.matrix(tk) < maskthresh] <- NA
            res[["samplemean"]] <- t1
            tmp <- im.apply(imlist, sd)
            tmp[as.matrix(tk) < maskthresh] <- NA
            res[["samplestdev"]] <- tmp
            tmp <- t1 / tk
            tmp[as.matrix(tk) < maskthresh] <- NA
            res[["sampleprop"]] <- tmp
        }
        if (!is.na(backsight)){
            bsBP <- focalyear + backsight
            bswt <- gaussW(x, bsBP, tbw)
            t0 <- density(pppA, weights=bswt, eps=cellres, sigma=sbw)
            tmp <- eval(changexpr)
            tmp[as.matrix(tk) < maskthresh] <- NA
            res[["change"]] <- tmp
            if (!is.null(nsim)){
                t0 <- im.apply(bimlist, mean)
                tmp <- eval(changexpr)
                tmp[as.matrix(tk) < maskthresh] <- NA
                res[["samplechange"]] <- tmp
            }
        }
        class(res) <- c("spKde", class(res))   
        if (verbose){ print("Done.") }
        return(res)
    }
}
