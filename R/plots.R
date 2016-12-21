plot.calDate <- function(calDate, label=NA, calendar="BP", type="standard"){

    types <- c("standard","simple")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    yearsBP <- calDate[["grid"]]$calBP
    prob <- calDate[["grid"]]$PrDens
    cra <- calDate[["metadata"]]$CRA
    error <- calDate[["metadata"]]$Error
    calcurve <- calDate[["metadata"]]$CalCurve
    calendars <- c("BP","BCAD")
    if (!calendar %in% calendars){
        stop("The calendar you have chosen is not currently an option.")
    }        
    if (yearsBP[1] > yearsBP[length(yearsBP)]){
        yearsBP <- rev(yearsBP)
        prob <- rev(prob)
    }
    if (calendar=="BP"){
        plotyears <- yearsBP
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])*-1
        xlabel <- "Years cal BP"
    } else if (calendar=="BCAD"){
        plotyears <- 1950-yearsBP
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
        xlabel <- "Years BC/AD"
    } else {
        stop("Unknown calendar type")
    }
    yvals <- c(0,prob,0,0)
    xrng <- c(min(xvals[yvals>0.000001])-50,max(xvals[yvals>0.000001])+50)
    xticks <- 100*(xrng%/%100 + as.logical(xrng%%100))
    xticks <- seq(xticks[1]-100, xticks[2], 100)
    yrng <- c(min(yvals[yvals>0]),max(yvals[yvals>0])+(max(yvals[yvals>0])*2))
    par(mar = c(4,4,1,3)) #c(bottom, left, top, right)
    par(cex.lab=0.75)
    plot(xvals,yvals, type="n", xlab=xlabel, ylab="", ylim=yrng, xlim=xrng, xaxt='n', yaxt='n', cex.axis=0.75)
    axis(1, at=xticks, labels=abs(xticks),las=2, cex.axis=0.75)
    axis(4, cex.axis=0.75)
    polygon(xvals,yvals, col="grey50", border="grey50")
    if (type != "simple"){
        par(new=TRUE)
        cradf1 <- data.frame(CRA=50000:0,Prob=dnorm(50000:0, mean=cra, sd=error))
        cradf1 <- cradf1[cradf1$Prob>0.0001,]
        xlim <- xrng
        ylim <- c(cra-(12*error),cra+(8*error))    
        cradf1$RX <- reScale(cradf1$Prob, type="custom", crng=c(xlim[1],(xlim[1]+diff(xlim)*0.33)))
        yticks <- ylim[1]:ylim[2]
        yticks <- yticks[yticks %% 200 == 0]
        plot(cradf1$RX,cradf1$CRA,type="l", axes=FALSE, xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, col=rgb(144,238,144,120,maxColorValue=255))
        polygon(c(cradf1$RX,rev(cradf1$RX)),c(cradf1$CRA,rep(xlim[1],length(cradf1$CRA))), col=rgb(144,238,144,80,maxColorValue=255), border=NA)
        axis(side=2, at=yticks, labels=abs(yticks),las=2, cex.axis=0.75)
        mtext(side=2, line=3, "Radiocarbon Age", cex=0.75)   
        calCurveFile <- paste(system.file("data", package="rcarbon"), "/", calcurve,".14c", sep="")
        options(warn=-1)
        cc <- readLines(calCurveFile, encoding="UTF-8")
        cc <- cc[!grepl("[#]",cc)]
        cc <- read.csv(textConnection(cc), header=FALSE, stringsAsFactors=FALSE)
        options(warn=0)
        names(cc) <- c("BP","CRA","Error","D14C","Sigma")
        if (calendar=="BCAD"){
            tmp <- (xrng-1950)*-1
            cc$RX <- 1950-cc$BP
        } else {
            tmp <- xrng*-1
            cc$RX <- cc$BP*-1
        }
        cc <- cc[cc$BP >= min(tmp) & cc$BP < max(tmp),] 
        cc$Hi <- cc$CRA + 10
        cc$Lo <- cc$CRA - 10
        ccbox <- c((max(yrng)*0.2),(max(yrng)*0.9))
        polygon(c(cc$RX,rev(cc$RX)),c(cc$Hi,rev(cc$Lo)), col=rgb(255,140,0,120,maxColorValue=255), border=NA)
        lines(cc$RX,cc$CRA, col=rgb(255,140,0,60,maxColorValue=255))
    }
    if (!is.na(label)){
        legend("topright", label, bty="n", cex=0.75)
    }
}

plot.rspdModelTest <- function(rspdModelTest, yMax=NA, xlim=c(0,1), drawaxes=TRUE, ...){

    obs <- rspdModelTest$result[,1:2]
    resolution <- 1
    envelope <- rspdModelTest$result[,3:4]
    if (is.na(yMax)){
        yMax <- max(envelope,obs[,2])
    }
    booms <- which(obs[,2]>envelope[,2])
    busts <- which(obs[,2]<envelope[,1])
    baseline <- rep(0,nrow(obs))
    if (is.null(xlim)){
        xlim <- c(max(obs[,1]),min(obs[,1]))
    }
    if (drawaxes){
        plot(obs[,1],obs[,2],xlim=xlim,ylim=c(0,yMax), xlab="cal BP",ylab="Normalised Summed Probability",type="l",col=1,lwd=0.5,...)
    } else {
        plot(obs[,1],obs[,2],xlim=xlim,ylim=c(0,yMax), xlab="",ylab="",type="l", xaxt="n", yaxt="n",col=1,lwd=0.5,...)
    }
    boomPlot <- baseline
    boomPlot[booms] <- obs[booms,2]
    bustPlot <- baseline
    bustPlot[busts] <- obs[busts,2]
    boomBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (x in 1:length(boomPlot)){
        if (boomPlot[x]>0&state=="off"){
            counter <- counter+1
            boomBlocks <- c(boomBlocks,vector("list",1))
            boomBlocks[[counter]] <- vector("list",2)
            boomBlocks[[counter]][[1]] <- boomPlot[x]
            boomBlocks[[counter]][[2]] <- obs[x,1]
            state <- "on"
        }
        if (state=="on"){
            if (boomPlot[x]>0){
                boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[x])
                boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],obs[x,1])
            }
            if (boomPlot[x]==0){
                state <- "off"
            }
        }   
    }
    bustBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (x in 1:length(bustPlot)){
        if (bustPlot[x]>0&state=="off"){
            counter <- counter+1
            bustBlocks <- c(bustBlocks,vector("list",1))
            bustBlocks[[counter]] <- vector("list",2)
            bustBlocks[[counter]][[1]] <- bustPlot[x]
            bustBlocks[[counter]][[2]] <- obs[x,1]
            state <- "on"
        }
        if (state=="on"){
            if (bustPlot[x]>0){
                bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[x])
                bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],obs[x,1])
            }
            if (bustPlot[x]==0){
                state <- "off"
            }
        }   
    }
    if (length(booms)>0){
        for (x in 1:length(boomBlocks)){
            polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(+100,length(boomBlocks[[x]][[1]])),rep(-100,length(boomBlocks[[x]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
        }
    }  
    if (length(busts)>0){
        for (x in 1:length(bustBlocks)){
            polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(+100,length(bustBlocks[[x]][[1]])),rep(-100,length(bustBlocks[[x]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
        }
    }  
    polygon(x=c(obs[,1],rev(obs[,1])),y=c(envelope[,1],rev(envelope[,2])),col=rgb(0,0,0,0.2),border=NA)
    spdSmooth <- runMean(obs[,2],200/resolution)
    ## lines(obs[,1],spdSmooth,col=1,lwd=2.5,lty=1)
    if (drawaxes){
        axis(side=1,at=seq(max(obs[,1]),min(obs[,1]),-100),labels=NA,tck = -.01)
    }
}
    
plot.rspdRegionTest <- function(data, focalregion="1", xlim=NA, yMax=1, drawaxes=TRUE, ...){

    obs <- data$observed[[focalregion]]
    if (any(is.na(xlim))){ xlim <- c(max(obs[,1]),min(obs[,1])) }
    resolution <- 1
    envelope <- data$envelope[[focalregion]]
    if (is.na(yMax)){ yMax=max(as.numeric(envelope),obs[,2]) }    
    booms <- which(obs[,2]>envelope[,2])
    busts <- which(obs[,2]<envelope[,1])
    baseline <- rep(0,nrow(obs))
    if (drawaxes){
        plot(obs[,1],obs[,2],xlim=xlim, ylim=c(0,yMax), xlab="cal BP",ylab="Normalised Summed Probability",type="l",col=1,lwd=0.5,axes=FALSE,...)
        axis(side=1,padj=-1)
        axis(side=2,padj=1)
    } else {
        plot(obs[,1],obs[,2],xlim=xlim, ylim=c(0,yMax), xlab="",ylab="",type="l",col=1,lwd=0.5, axes=FALSE,...)
    }
    box()
    boomPlot <- baseline
    if (length(booms)>0){ boomPlot[booms]=obs[booms,2] }
    bustPlot <- baseline
    if (length(busts)>0){ bustPlot[busts]=obs[busts,2] }                 
    boomBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (x in 1:length(boomPlot)){
        if (boomPlot[x]>0&state=="off"){
            counter <- counter+1
            boomBlocks <- c(boomBlocks,vector("list",1))
            boomBlocks[[counter]] <- vector("list",2)
            boomBlocks[[counter]][[1]] <- boomPlot[x]
            boomBlocks[[counter]][[2]] <- obs[x,1]
            state <- "on"
        }
        if (state=="on"){
            if (boomPlot[x]>0){
                boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[x])
                boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],obs[x,1])
            }
            if (boomPlot[x]==0){
                state <- "off"
            }
        }    
    }
    bustBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (x in 1:length(bustPlot)){
        if (bustPlot[x]>0&state=="off"){
            counter <- counter+1
            bustBlocks <- c(bustBlocks,vector("list",1))
            bustBlocks[[counter]] <- vector("list",2)
            bustBlocks[[counter]][[1]] <- bustPlot[x]
            bustBlocks[[counter]][[2]] <- obs[x,1]
            state <- "on"
        }
        if (state=="on"){
            if (bustPlot[x]>0){
                bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[x])
                bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],obs[x,1])
            }
            if (bustPlot[x]==0){
                state <- "off"
            }
        }    
    } 
    if (length(booms)>0){
        for (x in 1:length(boomBlocks)){
            polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(+100,length(boomBlocks[[x]][[1]])),rep(-100,length(boomBlocks[[x]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
        }
    }
    if (length(busts)>0){
        for (x in 1:length(bustBlocks)){
            polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(+100,length(bustBlocks[[x]][[1]])),rep(-100,length(bustBlocks[[x]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
        }
    }  
    polygon(x=c(obs[,1], rev(obs[,1])), y=c(envelope[,1], rev(envelope[,2])), col=rgb(0,0,0,0.2), border=NA)
    spdSmooth <- runMean(obs[,2], 200/resolution)
    ## lines(obs[,1],spdSmooth,col=1,lwd=2.5,lty=1)
    if (drawaxes){
        axis(side=1, at=seq(max(obs[,1]), min(obs[,1]),-100), labels=NA, tck=-0.01)
    }
}

barCodes <- function(x, yrng=c(0,0.03), width=20, col=rgb(0,0,0,25,maxColorValue=255), border=NA, fixXorder=FALSE,...){

    if (!"quickMarks" %in% class(x)){
        stop("Input must be of class \"quickMarks\"")
    }
    barcodes <- x$qMed
    if (fixXorder){ barcodes <- barcodes*-1 }
    halfbw <- width/2
    for (a in 1:length(barcodes)){
        polygon(x=c(barcodes[a]-halfbw,barcodes[a]-halfbw,barcodes[a]+halfbw,barcodes[a]+halfbw,barcodes[a]-halfbw),y=c(yrng[1],yrng[2],yrng[2],yrng[1],yrng[1]), border=border, col=col, ...)
    }
}

crossHairs <- function(x, pch.pts=19, cex.pts=1, fixXorder=FALSE, rescaleY=FALSE,...){

    if (!"quickMarks" %in% class(x)){
        stop("Input must be of class \"quickMarks\"")
    }
    if (rescaleY){
        cra <- reScale(x$CRA)
        error <- x$Error / (max(x$CRA)-min(x$CRA))
    } else {
        cra <- x$CRA
        error <- x$Error
    }
    if (fixXorder){
        xstart <- x$q68s *-1
        xend <- x$q68e *-1
        xmed <- x$qMed *-1
    } else {
        xstart <- x$q68s
        xend <- x$q68e
        xmed <- x$qMed
    }
    for (a in 1:nrow(x)){
        lines(c(xstart[a],xend[a]), c(cra[a],cra[a]), ...)
        lines(c(xmed[a],xmed[a]), c(cra[a]-error[a],cra[a]+error[a]), ...)
        points(xmed[a],cra[a], pch=pch.pts, cex=cex.pts, ...)
    }
}
