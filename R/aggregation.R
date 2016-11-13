calspdf <- function(ages, dens, calCurves='intcal13'){

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

binPrep <- function(sites,dates,h=200){
    
    clusters <- rep(NA,length(sites))
    
    for  (x in 1:length(unique(sites))){
        index <- which(sites==unique(sites)[x])
        if (length(index)>1){
            clusters[index] <- paste(unique(sites)[x],cutree(hclust(dist(dates[index])),h=h),sep="_")
        }
        if (length(index)==1){
            clusters[index] <- paste(unique(sites)[x],"1",sep="_")
        }                
    }
    
    return(clusters)
}

rspd <- function(bins, ages, errors, resOffsets=0, resErrors=0, runm=50, yearRange,calCurves, eps=1e-5, method="standard", normalised=FALSE, ncores=1, verbose=TRUE){
    
    # Set-up
    if (any(is.na(bins))){ stop("Cannot have NA values in the bins.") }
    dummygrid <- calibrate(ages=100,errors=10,timeRange=yearRange) # dummy
    individualDatesMatrix <- matrix(NA,nrow=nrow(dummygrid[[1]][["agegrid"]]),ncol=length(ages))
    # Calibrate
    tmp <- calibrate(ages=ages,errors=errors,resOffsets=resOffsets,resErrors=resErrors,timeRange=yearRange,calCurves=calCurves, method=method, normalised=normalised, ncores=ncores, verbose=verbose)
    tmp <- lapply(tmp, `[[`, 2)
    tmp <- lapply(tmp,`[`, 2)
    individualDatesMatrix <- do.call("cbind",tmp)
    # Aggregate
    binNames <- unique(bins)
    binnedMatrix <- matrix(NA,nrow=nrow(individualDatesMatrix),ncol=length(binNames))
    if (verbose){
        print("Binning by site/phase...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Binning by site/phase...")
    }
    for (b in 1:length(binNames)){
            if (verbose){ setTxtProgressBar(pb, b) }
            index <- which(bins==binNames[b])
            if (length(index)>1){    
                spd.tmp <- apply(individualDatesMatrix[,index],1,sum)
                binnedMatrix[,b] <- spd.tmp/length(index)
            } else {
                binnedMatrix[,b] <- individualDatesMatrix[,index]
            }
        }
    if (verbose){ close(pb) }
    finalSPD <- apply(binnedMatrix,1,sum)
    finalSPD <- finalSPD/sum(finalSPD)
    tmp1 <- runMean(finalSPD,runm)
    finalSPD[!is.na(tmp1)] <- tmp1[!is.na(tmp1)]
    res <- data.frame(calBP=dummygrid[[1]][["agegrid"]][,1], SPD=finalSPD)
    class(res) <- append(class(res),"RSPD")
    return(res)    
}

siteW <- function(dateID,cra,error,siteID,timeRange){
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

