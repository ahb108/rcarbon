calibrate <- function (x, ...) {
   UseMethod("calibrate")
}

calibrate.default <- function(ages, errors, ids=NA, dateDetails=NA, calCurves='intcal13', resOffsets=0 , resErrors=0, timeRange=c(50000,0), method="standard", normalised=FALSE, calMatrix=FALSE, dfs=100, oxpath=NULL, eps=1e-5, ncores=1, verbose=TRUE){

    if (length(ages) != length(errors)){
        stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
    }
    if (!is.na(ids[1]) & (length(ages) != length(ids))){
        stop("Ages and errors (and ids/details/offsets if provided) must be the same length.")
    }
    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(ages))
    if (calMatrix){
        calmBP <- seq(timeRange[1],timeRange[2],-1)
        calmat <- matrix(ncol=length(ages), nrow=length(calmBP))
        rownames(calmat) <- calmBP
    }
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
            if (!calMatrix){ res <- res[res$PrDens > 0,] }
            class(res) <- append(class(res),"calGrid")
            return(res)
        }
        stopCluster(cl)
        names(sublist) <- ids
        if (calMatrix){
            calmat <- sapply(grids, FUN=function(x) x$PrDens)
            rownames(calmat) <- calmBP
            sublist <- lapply(sublist, FUN=function(x) x[x$PrDens > 0, ])
        }
    } else {
        for (b in 1:length(ages)){
            if (length(ages)>1 & verbose){ setTxtProgressBar(pb, b) }
            calcurve <- cclist[[calCurves[b]]]
            calBP <- seq(max(calcurve),min(calcurve),-1)
            age <- ages[b] - resOffsets[b]
            error <- errors[b] + resErrors[b]
            methods <- c("standard","tDist","Bchron","OxCal","CalPallike")
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
            } else if (method=="OxCal"){
                if (is.null(oxpath)){ stop("You need to provide an oxpath argument.")
                } else {
                    if (b==1){
                        tmptxt <- capture.output(setOxcalExecutablePath(oxpath))
                    }
                    mydate <- oxcalCalibrate(age, error, ids[b])
                    years <- 1950-mydate[[1]]$raw_probabilities$dates
                    dens <- mydate[[1]]$raw_probabilities$probabilities
                    dens <- approx(years, dens, xout=calBP, rule=2)
                    res <- data.frame(calBP=calBP,PrDens=dens$y)
                }
            } else if (method=="Bchron"){
                tmp <- BchronCalibrate(ages=age,ageSds=error,calCurves=calCurves[b], eps=eps, dfs=dfs)
                calBPtmp <- rev(as.numeric(tmp[[1]][4][[1]]))
                prob <- rev(as.numeric(tmp[[1]][[5]]))
                dens <- rep(0,length=length(calBP))
                index <- which(calBP %in% calBPtmp)
                dens[index] <- prob
                dens[dens < eps] <- 0
                res <- data.frame(calBP=calBP,PrDens=dens)
            }
            res <- res[which(calBP<=timeRange[1]&calBP>=timeRange[2]),]
            if (calMatrix){ calmat[,b] <- res$PrDens }
            res <- res[res$PrDens > 0,]
            class(res) <- append(class(res),"calGrid")
            sublist[[ids[b]]] <- res
        }
    }
    if (length(ages)>1 & verbose){ close(pb) }
    df <- data.frame(DateID=ids, CRA=ages, Error=errors, Details=dateDetails, CalCurve=calCurves,ResOffsets=resOffsets, ResErrors=resErrors, StartBP=timeRange[1], EndBP=timeRange[2], CalMethod=method, Normalised=normalised, CalEPS=eps, stringsAsFactors=FALSE)
    reslist[["metadata"]] <- df
    reslist[["grids"]] <- sublist
    if (calMatrix){
        reslist[["calmatrix"]] <- calmat
    }
    class(reslist) <- append(class(reslist),"CalDates")
    if (verbose){ print("Done.") }
    return(reslist)
}

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

uncalibrate <- function (x, ...) {
   UseMethod("uncalibrate", x)
}

uncalibrate.default <- function(calBP, CRAerrors=NA, roundyear=TRUE, calCurves='intcal13', method="standard"){ 

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

uncalibrate.CalGrid <- function(calgrid, calCurves='intcal13', eps=1e-5, unifp="local", compact=TRUE, verbose=TRUE){

    if (verbose){ print("Uncalibrating...") }
    names(calgrid) <- c("calBP","PrDens")
    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- calcurve[!grepl("[#]",calcurve)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")
    mycras <- uncalibrate(calgrid$calBP)
    res <- data.frame(CRA=max(calcurve[,2]):min(calcurve[,2]), PrDens=0)
    tmp <- vector(mode="list",length=nrow(mycras))
    basetmp <- vector(mode="list",length=nrow(mycras))
    if (length(tmp)>1 & verbose){
        flush.console()
        pb <- txtProgressBar(min=1, max=length(tmp), style=3)
    }
    for (a in 1:length(tmp)){
        basetmp[[a]] <- dnorm(res$CRA, mean=mycras$ccCRA[a], sd=mycras$ccError[a])
        tmp[[a]] <- basetmp[[a]] * calgrid$PrDens[a]
        if (verbose){ setTxtProgressBar(pb, a) }
    }
    if (verbose){ close(pb) }
    unscGauss <- do.call("cbind",tmp)
    res$Raw <- rowSums(unscGauss)
    res$Raw[res$Raw < eps] <- 0
    if (unifp=="local"){
        base <- do.call("cbind",basetmp)
        res$Base <- rowSums(base)
        res$Raw[res$Raw < eps] <- 0
    } else if (unifp=="global"){
        data(UnifCalYears)
        res$Base <- UnifCalYears[UnifCalYears$CRA %in% res$CRA,"PrDens"]
    } else {
        stop("Options for unifp are 'local' or 'global'.")
    }
    res$PrDens[res$Base>0] <- res$Raw[res$Base>0] / res$Base[res$Base>0]
    if (compact){ res <- res[res$PrDens > 0,] }
    class(res) <- c("UncalGrid", class(res)) 
    if (verbose){ print("Done.") }
    return(res)
}

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

"[.calDates" <- function(x,i){
    
    if (nrow(x$metadata)==0){
        stop("No data to extract")
    }
    if(!missing(i)) {
        if (all(is.numeric(i)) | all(is.character(i)) | all(is.logical(i))){
            if (length(x$calmat>0)){
                res <- list(metadata=x$metadata[i,], grids=x$grids[i], calmatrix=x$calmatrix[,i])
            } else {
                res <- list(metadata=x$metadata[i,], grids=x$grids[i])
            }
            class(res) <- c("CalDates", class(res))        
        } else {
            stop("i must be a numeric, character or logical vector of length(x)")
        }
        return(res)
    }           
}
