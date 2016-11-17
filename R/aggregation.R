#' Calibrate (or uncalibrate) a set of radiocarbon probability distributions that are already summed/aggregated.
#'
#' @param rspd object of class CalSPD (a summed distribution of calibrated calender date probabilities) or UncalSPD (a summed distribution of either uncalibrated conventional radiocarbon age probabilities)
#' @param calcurve character string indicating the radiocarbon calibration curve to be used or an object of class 'CalCurve')
#'
#' @return An object of class CalSPD or UncalSPD
#'
#' @examples
#' rspdCal()
rspdCal <- function(rspd, calcurve='intcal13'){ 

    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- cctmp[!grepl("[#]",cctmp)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")
    calBP.out <- seq(max(calcurve),min(calcurve),-1)
    CRAdates <- data.frame(approx(calcurve[,1:2], xout=calBP.out))
    names(CRAdates) <- c("calBP.out","CRA")
    CRAdates$CRA <- round(CRAdates$CRA,0)
    CRApdf <-  data.frame(CRA=ages,prob.out=dens)
    res <- merge(CRAdates,CRApdf,by="CRA",all.x=TRUE, sort=FALSE)
    res <- res[with(res, order(-calBP.out)), c("calBP.out","prob.out")]
    res$prob.out[is.na(res$prob.out)] <- 0
    return(res$prob.out)
}

#' Prepare a set of bins for controlling the aggregation of radiocarbon dates
#' known to be from the same phase of same archaeological site (for use with rspd)
#'
#' Used in cases where there is a concern that unusually high levels of sampling for radiocarbon at a given site or in agiven site phase will impede comparison between sites or phases. 
#' 
#' @param sites a vector of character strings (or number to coerce to character) of all sites or site phases
#' @param calcurve a vector uncalibrated conventional radiocarbon ages
#' @param h a single numeric value passed to hclust() to control degree of grouping of similar ages in a phase site.
#'
#' @return A vector of character strings of length(ages) that identifying intra-site or intra-phase grouping, for use with rspd()
#'
#' @examples
#' binPrep()
binPrep <- function(sites, ages, h=200){
    
    clusters <- rep(NA,length(sites))
    
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

rspd <- function(x, timeRange, bins=NA, datenormalised=FALSE, spdnormalised=TRUE, runm=NA, verbose=TRUE){

    defcall <- as.list(args(rspd))
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
    speccall$ndates <- length(x)
    speccall$nbins <- length(x)
    if (!"calDates" %in% class(x)){
        stop("x must be an object of class 'calDates'.")
    }
    if (length(bins)>1){
        speccall$nbins <- length(unique(bins))
        if (any(is.na(bins))){
            stop("Cannot have NA values in bins.")
        }
        if (length(bins)!=length(x)){
            stop("bins (if provided) must be the same lenght as x.")
        }
    } else {
        bins <- rep("0_0",length(x))
    }
    binNames <- unique(bins)
    binnedMatrix <- matrix(NA, nrow=nrow(x[[1]][["agegrid"]]), ncol=length(binNames))
    if (verbose & length(binNames)>1){
        print("Binning by site/phase...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Binning by site/phase...")
    }
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        slist <- x[index]
        tmp <- lapply(lapply(slist, `[[`, 2),`[`,2)
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
    if (verbose){ print("Aggregating...") }
    finalSPD <- apply(binnedMatrix,1,sum)
    if (!is.na(runm)){
        tmp1 <- runMean(finalSPD,runm)
        finalSPD[!is.na(tmp1)] <- tmp1[!is.na(tmp1)]
    }
    res <- data.frame(calBP=x[[1]][["agegrid"]][,1], SPD=finalSPD)
    if (spdnormalised){
        res$SPD <- res$SPD/sum(res$SPD, na.rm=TRUE)
    }
    res <- res[res$calBP <= timeRange[1] & res$calBP >= timeRange[2],]
    reslist <- vector("list",length=2)
    names(reslist) <- c("metadata","agegrid")
    reslist[["metadata"]] <- speccall
    reslist[["agegrid"]] <- res
    class(reslist) <- append(class(reslist),"CalSPD")
    if (verbose){ print("Done.") }
    return(reslist)
}

siteW <- function(dateID ,cra,error,siteID,timeRange){
    df <- data.frame(DateID=dateID,CRA=cra,Error=error,SiteID=siteID,Weights=1)
    allsiteids <- unique(siteID)
    if (any(is.na(allsiteids))){ stop("Cannot have missing site IDs") }
    for (a in 1:length(allsiteids)){
        cat(paste(a,";",sep=""))
        tmp <- df[df$SiteID==allsiteids[a] & !is.na(df$SiteID),]
        tmpc <- mcalibrate(tmp$CRA,tmp$Error, progress=FALSE, timeRange=timeRange)
        if (length(tmpc)>1){
            overlaps <- rep(NA,length(tmpc))
            for (i in 1:length(tmpc)){
                otherinds <- 1:length(tmpc)
                otherinds <- otherinds[-i]
                tmpo <- rep(NA,length(otherinds))
                for (j in 1:length(otherinds)){
                    tmpo[j] <- sum(abs(tmpc[[i]]$prob.out-tmpc[[otherinds[j]]]$prob.out))/2
                }
                overlaps[i] <- sum(tmpo)/length(tmpo)
            }
            df[df$SiteID==allsiteids[a],"Weights"] <- overlaps
        }
    }
    return(df)
}

