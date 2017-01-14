calibrate <- function(ages, errors, ids=NA, dateDetails=NA, calCurves='intcal13', resOffsets=0 , resErrors=0, timeRange=c(50000,0), method="standard", normalised=FALSE, compact=TRUE, dfs=100, eps=1e-5, ncores=1, verbose=TRUE){

    ## NB. add manualpage thanks to Bchron/Parnell
    if (length(ages) != length(errors)){
        stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
    }
    if (!is.na(ids[1]) & (length(ages) != length(ids))){
        stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
    }
    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(ages))
    if (is.na(ids[1])){
        ids <- as.character(1:length(ages))
    } else {
        ids <- as.character(ids)
    }
    if (length(dfs)==1){ dfs <- rep(dfs, length(ages)) }
    if (length(resOffsets)==1){ resOffsets <- rep(resOffsets,length(ages)) }
    if (length(resErrors)==1){ resErrors <- rep(resErrors,length(ages)) }
    names(sublist) <- ids
    names(reslist) <- c("metadata","grids")
    tmp <- unique(calCurves)
    if (length(calCurves)==1){ calCurves <- rep(calCurves,length(ages)) }
    cclist <- vector(mode="list", length=length(tmp))
    for (a in 1:length(tmp)){
        calCurveFile <- paste(system.file("data", package="rcarbon"), "/", tmp,".14c", sep="")
        options(warn=-1)
        cctmp <- readLines(calCurveFile, encoding="UTF-8")
        cctmp <- cctmp[!grepl("[#]",cctmp)]
        cctmp <- as.matrix(read.csv(textConnection(cctmp), header=FALSE, stringsAsFactors=FALSE))[,1:3]
        options(warn=0)
        colnames(cctmp) <- c("CALBP","C14BP","Error")
        cclist[[tmp[a]]] <- cctmp
    }
    if (length(ages)>1 & verbose){
            print("Calibrating radiocarbon ages...")
            flush.console()
            pb <- txtProgressBar(min=1, max=length(ages), style=3)
    }
    # Parallellised
    if (ncores>1){
        require(doParallel)
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        if (verbose){ print(paste("Running in parallel (standard calibration only) on ",getDoParWorkers()," workers...",sep=""))}
        sublist <- foreach (b=1:length(ages)) %dopar% {
            method <- "standard"
            calcurve <- cclist[[calCurves[b]]]
            calBP <- seq(max(calcurve),min(calcurve),-1)
            age <- ages[b] - resOffsets[b]
            error <- errors[b] + resErrors[b]
            mu <- approx(calcurve[,1], calcurve[,2], xout=calBP)$y
            tau <- error^2 + approx(calcurve[,1], calcurve[,3], xout=calBP)$y
            dens <- dnorm(age, mean=mu, sd=sqrt(tau))
            dens[dens < eps] <- 0
            if (normalised){ dens <- dens/sum(dens) }
            res <- data.frame(calBP=calBP,PrDens=dens)
            
            res <- res[which(calBP<=timeRange[1]&calBP>=timeRange[2]),]
            if (compact){ res <- res[res$PrDens > 0,] }
            class(res) <- append(class(res),"calGrid")
            return(res)
        }
        stopCluster(cl)
        names(sublist) <- ids
    } else {
        for (b in 1:length(ages)){
            if (length(ages)>1 & verbose){ setTxtProgressBar(pb, b) }
            calcurve <- cclist[[calCurves[b]]]
            calBP <- seq(max(calcurve),min(calcurve),-1)
            age <- ages[b] - resOffsets[b]
            error <- errors[b] + resErrors[b]
            methods <- c("standard","tDist","Bchron","CalPallike")
            if (!method %in% methods){
                stop("The method you have chosen is not currently an option.")
            }
            if (method=="standard"){
                mu <- approx(calcurve[,1], calcurve[,2], xout=calBP)$y
                tau <- error^2 + approx(calcurve[,1], calcurve[,3], xout=calBP)$y
                dens <- dnorm(age, mean=mu, sd=sqrt(tau))
                dens[dens < eps] <- 0
                if (normalised){
                    dens <- dens/sum(dens)
                    dens[dens < eps] <- 0
                    dens <- dens/sum(dens)
                }
                res <- data.frame(calBP=calBP,PrDens=dens)
            } else if (method=="tDist"){
                mu <- approx(calcurve[,1], calcurve[,2], xout=calBP)$y
                tau <- error^2 + approx(calcurve[,1], calcurve[,3], xout=calBP)$y
                dens <- dt((age - mu)/sqrt(tau), df=dfs)
                dens[dens < eps] <- 0
                if (normalised){
                    dens <- dens/sum(dens)
                    dens[dens < eps] <- 0
                    dens <- dens/sum(dens)
                }
                res <- data.frame(calBP=calBP,PrDens=dens)
            } else if (method=="Bchron"){
                tmp <- BchronCalibrate(ages=age,ageSds=error,calCurves=calCurves[b],eps=eps)
                calBPtmp <- rev(as.numeric(tmp[[1]][4][[1]]))
                prob <- rev(as.numeric(tmp[[1]][[5]]))
                dens <- rep(0,length=length(calBP))
                index <- which(calBP %in% calBPtmp)
                dens[index] <- prob
                dens[dens < eps] <- 0
                res <- data.frame(calBP=calBP,PrDens=dens)
            } else if (method=="CalPallike"){
                CRAdates <- data.frame(approx(calcurve, xout=calBP))
                names(CRAdates) <- c("calBP","CRA")
                CRAdates$CRA <- round(CRAdates$CRA,0)
                CRApdf <-  data.frame(CRA=max(CRAdates$CRA):min(CRAdates$CRA))
                CRApdf$PrDens <- dnorm(CRApdf$CRA, mean=age, sd=error)
                CRApdf$PrDens[CRApdf$PrDens < eps] <- 0
                res <- merge(CRAdates,CRApdf,by="CRA",all.x=TRUE, sort=FALSE)
                res <- res[with(res, order(-calBP)), c("calBP","PrDens")]
                if (normalised){
                    dens <- res$PrDens / sum(res$PrDens)
                    dens[dens < eps] <- 0
                    res$PrDens <- dens/sum(dens)                }
            }
            res <- res[which(calBP<=timeRange[1]&calBP>=timeRange[2]),]
            if (compact){ res <- res[res$PrDens > 0,] }
            class(res) <- append(class(res),"calGrid")
            sublist[[ids[b]]] <- res
        }
    }
    if (length(ages)>1 & verbose){ close(pb) }
    df <- data.frame(DateID=ids, CRA=ages, Error=errors, Details=dateDetails, CalCurve=calCurves,ResOffsets=resOffsets, ResErrors=resErrors, StartBP=timeRange[1], EndBP=timeRange[2], CalMethod=method, Normalised=normalised, CalEPS=eps, stringsAsFactors=FALSE)
    reslist[["metadata"]] <- df
    reslist[["grids"]] <- sublist
    class(reslist) <- append(class(reslist),"calDates")
    if (verbose){ print("Done.") }
    return(reslist)
}


uncalibrate <- function(calBP, CRAerrors=NA, roundyear=TRUE, calCurves='intcal13', method="standard", normalised=TRUE, eps=1e-5){ 

    if (length(CRAerrors)==1){ CRAerrors <- rep(CRAerrors,length(calBP)) } 
    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- calcurve[!grepl("[#]",calcurve)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")
    dates <- data.frame(approx(calcurve, xout=calBP))
    colnames(dates) <- c("calBP", "ccCRA")
    calcurve.error <- approx(calcurve[,c(1,3)], xout=dates$calBP)$y
    if (method == "standard"){
        dates$ccError <- calcurve.error
        dates$rCRA <- rnorm(nrow(dates), mean=dates$ccCRA, sd=dates$ccError)
        dates$rError <- CRAerrors
        if (roundyear){ dates$rCRA <- round(dates$rCRA) }
    } else if (method == "Cremaetal16"){
        if (is.na(CRAerrors[1])){
            stop("For this method you must provide a numeric error argument for the anticipated measurement error")
        }
        dates$ccError <- calcurve.error       
        dates$rCRA <- rnorm(nrow(dates), mean=dates$ccCRA, sd=CRAerrors)
        dates$rError <- sqrt(CRAerrors^2 + calcurve.error^2)
        if (roundyear){
            dates$rCRA <- round(dates$rCRA)
            dates$rError <- round(dates$rError)
        }
    } else {
        stop("Not one of the currently supplied methods." )   
    }
    return(dates)
}

reScale <- function(x, type="simple", crng=NULL, na.rm=TRUE){

    types <- c("simple","normal", "custom")
    if (!type %in% types){
        stop("The rescale type you have chosen is not currently an option.")
    }
    if (na.rm){ x <- na.omit(x) }
    if (type=="normal"){
        res <- (x-mean(x))/sd(x)
    } else if (type=="custom"){
        if (is.null(crng)){
            stop("For custom type you need to specify a crng.")
        } else {
            xrange <- range(x)
            mfac <- (crng[2] - crng[1])/(xrange[2] - xrange[1])
            res <- crng[1] + (x - xrange[1]) * mfac
        }
    } else {
       res <- (x-min(x))/(max(x) - min(x))
   }
    return(res)
}

gaussW <- function(x, bw){
    exp(-(x^2)/(2*(bw^2)))
}

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

quickMarks <- function(x, verbose=TRUE){

    if (!"calDates" %in% class(x)){
        stop("Input must be of class \"calDates\"")
    }
    df <- as.data.frame(matrix(ncol=8,nrow=length(x$grid)), stringsasFactors=TRUE)
    names(df) <- c("DateID","CRA","Error","qMed","q95s","q95e","q68s","q68e")
    print("Extracting approximate values...")
    if (length(x$grid)>1 & verbose){
        flush.console()
        pb <- txtProgressBar(min=1, max=length(x$grid), style=3)
    }
    for (a in 1:length(x$grid)){
        if (length(x$grid)>1 & verbose){ setTxtProgressBar(pb, a) }
        tmp <- x$grid[[a]]
        tmp <- tmp[tmp$PrDens>0,]
        tmp <- tmp[with(tmp, order(-PrDens)), ]
        tmp$Cumul <- cumsum(tmp$PrDens)
        tmp$Cumul <- tmp$Cumul/max(tmp$Cumul)
        tmp1 <- tmp[tmp$Cumul <= 0.95,]
        df[a,"q95s"] <- min(tmp1$calBP)
        df[a,"q95e"] <- max(tmp1$calBP)
        wdth <- max(tmp1$calBP)-min(tmp1$calBP)
        df[a,"q68s"] <- min(tmp1$calBP) + (wdth*0.25)
        df[a,"q68e"] <- max(tmp1$calBP) - (wdth*0.25)
        df[a,"qMed"] <- round(mean(c(df[a,"q95s"],df[a,"q95e"])),0)
        df[a,"DateID"] <- as.character(x$metadata$DateID[a])
        df[a,"CRA"] <- x$metadata$CRA[a]
        df[a,"Error"] <- x$metadata$Error[a]
    }
    if (length(x)>1 & verbose){ close(pb) }
    class(df) <- append(class(df),"quickMarks")
    return(df)
}

