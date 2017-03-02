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
binPrep <- function(sites, ages, h){
    
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

    if (verbose){ print("Extracting...") }
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
    speccall$ndates <- length(x$grids)
    speccall$nbins <- length(x$grids)
    if (!"calDates" %in% class(x)){
        stop("x must be an object of class 'calDates'.")
    }
    if (length(bins)>1){
        speccall$nbins <- length(unique(bins))
        if (any(is.na(bins))){
            stop("Cannot have NA values in bins.")
        }
        if (length(bins)!=length(x$grids)){
            stop("bins (if provided) must be the same length as x.")
        }
    } else {
        bins <- rep("0_0",length(x$grids))
    }
    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))
    if (verbose & length(binNames)>1){
        print("Binning...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3)
    }
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
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
        finalSPD <- runMean(finalSPD, runm, edge="fill")
    }
    res <- data.frame(calBP=calyears$calBP, PrDens=finalSPD)
    if (spdnormalised){
        res$PrDens <- res$PrDens/sum(res$PrDens, na.rm=TRUE)
    }
    res <- res[res$calBP <= timeRange[1] & res$calBP >= timeRange[2],]
    class(res) <- append(class(res),"CalGrid")
    reslist <- vector("list",length=2)
    names(reslist) <- c("metadata","grid")
    reslist[["metadata"]] <- speccall
    reslist[["grid"]] <- res
    class(reslist) <- append(class(reslist),"CalSPD")
    if (verbose){ print("Done.") }
    return(reslist)
}

overlapW <- function(calDates, bins, verbose=TRUE){
    df <- data.frame(Weight=rep(1,length(myCalDates$grids)))
    binNames <- unique(bins)
    if (any(is.na(binNames))){ stop("Cannot have missing bin IDs") }
    if (verbose){
        print("Checking overlaps within each bin...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Checking overlaps within each bin...")
    }
    for (a in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, a) }
        index <- which(bins==binNames[a])
        binneddates <- calDates$metadata[index,"DateID"]
        tmpc <- calDates$grids[index]
        if (length(tmpc)>1){
            maxbp <- max(sapply(tmpc,FUN=function(x) max(x$calBP)))
            minbp <- min(sapply(tmpc,FUN=function(x) min(x$calBP)))
            calyears <- data.frame(calBP=seq(maxbp, minbp,-1))
            slist <- lapply(tmpc,FUN=function(x) merge(calyears,x, all.x=TRUE)) 
            slist <- rapply(slist, f=function(x) ifelse(is.na(x),0,x), how="replace")
            slist <- lapply(slist, FUN=function(x) x[with(x, order(-calBP)), ])
            for (i in 1:length(slist)){
                mydateid <- binneddates[i]
                focal <- slist[[i]]
                others <- slist[-i]
                overlaps <- rep(NA,length(others))
                for (j in 1:length(others)){             
                    tmp <- apply(cbind(focal$PrDens,others[[j]]$PrDens),1,min)
                    overlaps[j] <- sum(tmp)
                }
                wt <- sum(focal$PrDens) * (1 - (sum(overlaps) / (sum(focal$PrDens)*length(overlaps))))
                df$Weight[calDates$metadata$DateID==mydateid] <- wt
            }
        }
    }
    if (verbose & length(binNames)>1){ close(pb) }
    if (verbose){ print("Done.") }
    return(df)
}

spd <- function(x, timeRange, bins=NA, datenormalised=FALSE, spdnormalised=TRUE, runm=NA, verbose=TRUE){
    
    if (verbose){ print("Extracting...") }
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
    speccall$ndates <- length(x$grids)
    speccall$nbins <- length(x$grids)
    if (!"calDates" %in% class(x)){
        stop("x must be an object of class 'calDates'.")
    }
    if (length(bins)>1){
        speccall$nbins <- length(unique(bins))
        if (any(is.na(bins))){
            stop("Cannot have NA values in bins.")
        }
        if (length(bins)!=length(x$grids)){
            stop("bins (if provided) must be the same length as x.")
        }
    } else {
        bins <- rep("0_0",length(x$grids))
    }
    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))
    if (verbose & length(binNames)>1){
        print("Binning...")
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
    if (verbose){ print("Aggregating...") }
    finalSPD <- apply(binnedMatrix,1,sum)
    if (!is.na(runm)){
        finalSPD <- runMean(finalSPD, runm, edge="fill")
    }
    res <- data.frame(calBP=calyears$calBP, PrDens=finalSPD)
    if (spdnormalised){
        res$PrDens <- res$PrDens/sum(res$PrDens, na.rm=TRUE)
    }
    res <- res[res$calBP <= timeRange[1] & res$calBP >= timeRange[2],]
    class(res) <- append(class(res),"CalGrid")
    reslist <- vector("list",length=2)
    names(reslist) <- c("metadata","grid")
    reslist[["metadata"]] <- speccall
    reslist[["grid"]] <- res
    class(reslist) <- append(class(reslist),"CalSPD")
    if (verbose){ print("Done.") }
    return(reslist)
}
