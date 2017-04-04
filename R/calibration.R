calibrate <- function (x, ...) {
   UseMethod("calibrate")
}

calibrate.default <- function(ages, errors, ids=NA, dateDetails=NA, calCurves='intcal13', resOffsets=0 , resErrors=0, timeRange=c(50000,0), method="standard", normalised=FALSE, calMatrix=FALSE, dfs=100, oxpath=NULL, iter=50000,eps=1e-5, ncores=1, verbose=TRUE){

    if (length(ages) != length(errors)){
        stop("Ages and errors (and ids/date details/offsets if provided) must be the same length.")
    }
    if (!is.na(ids[1]) & (length(ages) != length(ids))){
        stop("Ages and errors (and ids/details/offsets if provided) must be the same length.")
    }
   if (any(is.na(ages))|any(is.na(errors))){
        stop("Ages or errors contain NAs")
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
            methods <- c("standard","tDist","OxCal","MCMC")
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
                }
           	    mydate <- oxcalSingleDate(ages=age,error=error,OxCalExecute=oxpath, calCurve=calCurves[b])
		    dens <- approx(mydate$years, mydate$dens, xout=calBP, rule=2)
                    res <- data.frame(calBP=calBP,PrDens=dens$y)
		    normalised <- TRUE
            } else if (method=="MCMC"){
           	    mydate <- jagsSingleCalibrate(age=age,error=error, calCurve=calCurves[b],iter=iter)
		    dens <- approx(mydate$calBP, mydate$PrDens, xout=calBP, rule=1)
                    res <- data.frame(calBP=calBP,PrDens=dens$y)
		    res$PrDens[is.na(res$PrDens)] <- 0
		    normalised <- TRUE
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
    class(reslist) <- c("CalDates",class(reslist))
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

uncalibrate.CalGrid <- function(calgrid, calCurves='intcal13', eps=1e-5, compact=TRUE, verbose=TRUE){

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
    base <- do.call("cbind",basetmp)
    res$Base <- rowSums(base)
    res$Raw[res$Raw < eps] <- 0
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



"[.CalDates" <- function(x,i){
    
    if (nrow(x$metadata)==0){
        stop("No data to extract")
    }
    if(!missing(i)) {
        if (all(is.numeric(i)) | all(is.character(i)) | all(is.logical(i))){
            if (length(x$calmatrix>0)){
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

oxcalSingleDate<-function(id="tmp01",ages,error,OxCalExecute,calCurve)
    {
        fn <- tempfile()
        fnRecieve <- paste(fn, ".js", sep = "")
        fn <- paste(fn, ".oxcal", sep = "")
	
        cat("Options(){};\n",file=fn,append=FALSE) #Start Sequence#
        cat("Plot(){\n",file=fn,append=TRUE) #Start Sequence#
        if (calCurve=="marine13")
            {
                cat('Curve("Marine13.14c");',file=fn,append=TRUE)
                cat(paste('Delta_R(',DeltaR,',',DeltaRsd,');\n',sep=""),file=fn,append=TRUE)
            }
        if (calCurve=="shcal13")
            {
                cat('Curve("ShCal13");',file=fn,append=TRUE)
            }

        if (calCurve=="intcal13")
            {
                cat('Curve("IntCal13");',file=fn,append=TRUE)
            }
        cat(paste('R_Date(','\"',id,'\",',ages,',',error,');\n',sep=""),file=fn,append=TRUE)
        cat('};\n',file=fn,append=TRUE)
        excecuter=paste(OxCalExecute,fn)
        system(excecuter)        
	result <- scan(fnRecieve, character(0), sep = "\n",quiet=T)


  probs <- as.double(na.omit(unlist(strsplit(stringr::str_match(result, "(ocd\\[\\d+\\].likelihood.prob=\\[)(.*)(\\];)")[, 3], ", "))))
  pstart <- as.double(na.omit(unlist(strsplit(stringr::str_match(result, "(ocd\\[\\d+\\].likelihood.start=)(.*)(;)")[, 3], ", "))))
  resolution <- as.double(na.omit(unlist(strsplit(stringr::str_match(result, "(ocd\\[\\d+\\].likelihood.resolution=)(.*)(;)")[, 3], ", "))))
  normaliser <- as.double(na.omit(unlist(strsplit(stringr::str_match(result, "(ocd\\[\\d+\\].likelihood.probNorm=)(.*)(;)")[, 3], ", "))))

  if(is.na(normaliser)) {normaliser <- 1}

  years <- seq(pstart, by = resolution, length.out = length(probs))
  res <- data.frame(years= 1950-years, dens = probs * normaliser)
	return(res)
    }


jagsSingleCalibrate<-function(age,error,calCurves='intcal13',init=NA,iter=50000)
{

 #   calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")

    calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calCurves,".14c", sep="")
    options(warn=-1)
    calcurve <- readLines(calCurveFile, encoding="UTF-8")
    calcurve <- calcurve[!grepl("[#]",calcurve)]
    calcurve <- as.matrix(read.csv(textConnection(calcurve), header=FALSE, stringsAsFactors=FALSE))[,1:3]
    options(warn=0)
    colnames(calcurve) <- c("CALBP","C14BP","Error")


    require(rjags) 
    calBP <- calcurve[,1]
    C14BP <- calcurve[,2]
    C14err <- calcurve[,3]
    dataList <- list(nDate=1, X=age, sigma=error, calBP=rev(calBP), C14BP=rev(C14BP),C14err=rev(C14err))

    if (is.na(init)){init=calBP[which(abs(C14BP-age)==min(abs(C14BP-age)))[1]]}

    ##Specify JAGS Model
    modelString="
      model{
      for (i in 1:nDate) {
      theta[i] ~ dunif(0,50000)
      mu[i] <- interp.lin(theta[i], calBP[], C14BP[])
      sigmaCurve[i] <- interp.lin(theta[i], calBP[], C14err[])
      tau[i] <- 1/(pow(sigma[i],2)+pow(sigmaCurve[i],2))
      X[i] ~ dnorm(mu[i],tau[i])
      twenty.year[i] <- 20*round(theta[i]/20)
      ten.year[i] <- 10*round(theta[i]/10)
      one.year[i] <- round(theta[i])
 	}
      
      }"

    initsList <- list(theta=init)
    jagsModel <- jags.model(file=textConnection(modelString),data=dataList,inits=initsList,n.chains=1,n.adapt=3000)
    update(jagsModel)
    codaSamples=coda.samples(jagsModel,variable.names=c("one.year"),n.iter=iter)
    
    BP <- as.numeric(names(table(codaSamples)))
    PrDens <- as.numeric(table(codaSamples)/iter)

    fullBP <- min(BP):max(BP)
    res <- data.frame(BP=fullBP,PrDens=0)	
    res <- merge(x=res,y=data.frame(BP=BP,PrDens=PrDens),by.x="BP",by.y="BP",all=TRUE)
    res$PrDens.y[which(is.na(res$PrDens.y))] <- 0
    finalRes <- data.frame(calBP=rev(res$BP),PrDens=rev(res$PrDens.y)) 	
    return(finalRes)
}



