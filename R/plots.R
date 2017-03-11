plot.CalDates <- function(calDates, ind=1, label=NA, calendar="BP", type="standard"){

    types <- c("standard", "simple", "auc")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    yearsBP <- calDates$grids[[ind]]$calBP
    prob <- calDates$grids[[ind]]$PrDens
    cra <- calDates$metadata$CRA[ind]
    error <- calDates$metadata$Error[ind]
    calcurve <- calDates$metadata$CalCurve[ind]
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
    if (type=="standard" | type=="auc"){
        if (type=="auc"){
            lines(xvals, yvals/sum(yvals), col="black", lty="dotted")
        }
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

plot.SpdModelTest <- function(test, calendar="BP", ylim=NA, xlim=NA, col.obs="black", lwd.obs=0.5, xaxs="i", yaxs="i", bbty="f", drawaxes=TRUE, ...){

    obs <- test$result[,1:2]
    if (calendar=="BP"){
        obs$Years <- obs$calBP
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
    } else if (calendar=="BCAD"){
        obs$Years <- 1950-obs$calBP
        xlabel <- "Years BC/AD"
        if (any(is.na(xlim))){ xlim <- c(min(obs$Years),max(obs$Years)) }
    } else {
        stop("Unknown calendar type")
    }    
    envelope <- test$result[,3:4]
    if (any(is.na(ylim))){ ylim <- c(0, max(envelope[,"hi"], obs$PrDens)*1.1) }
    booms <- which(obs$PrDens>envelope[,2])
    busts <- which(obs$PrDens<envelope[,1])
    baseline <- rep(0,nrow(obs))
    if (drawaxes & bbty != "n"){
        plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab=xlabel, ylab="Summed Probability", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, ...)
    } else if (bbty != "n"){
        plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab="", ylab="", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
    }
    box()
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
            boomBlocks[[counter]][[2]] <- obs[x,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (boomPlot[x]>0){
                boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[x])
                boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],obs[x,"Years"])
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
            bustBlocks[[counter]][[2]] <- obs[x,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (bustPlot[x]>0){
                bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[x])
                bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],obs[x,"Years"])
            }
            if (bustPlot[x]==0){
                state <- "off"
            }
        }   
    }
    if (length(booms)>0){
        for (x in 1:length(boomBlocks)){
            if (bbty=="f"){
                polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(+100,length(boomBlocks[[x]][[1]])),rep(-100,length(boomBlocks[[x]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }  
    if (length(busts)>0){
        for (x in 1:length(bustBlocks)){
            if (bbty=="f"){
                polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(+100,length(bustBlocks[[x]][[1]])),rep(-100,length(bustBlocks[[x]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }  
    polygon(x=c(obs[,"Years"],rev(obs[,"Years"])),y=c(envelope[,1],rev(envelope[,2])),col=rgb(0,0,0,0.2),border=NA)
    if (drawaxes & bbty != "n"){
        axis(side=1,at=seq(max(obs[,"Years"]),min(obs[,"Years"]),-100),labels=NA,tck = -.01)
    }
    bbp <- list(booms=boomBlocks, busts=bustBlocks)
    class(bbp) <- c("BBPolygons",class(bbp))
    if (bbty %in% c("n","b")){ return(bbp) }
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

plot.CalSPD <- function(spd, runm=NA, calendar="BP", type="standard", xlim=NA, ylim=NA, spdnormalised=FALSE, rescale=FALSE, fill.p="grey75", border.p=NA, xaxt='s', yaxt='s', ...){

    types <- c("standard","simple")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    spdvals <- spd$grid$PrDens
    if (!is.na(runm)){ spdvals <- runMean(spdvals, runm, edge="fill") }
    if (spdnormalised){ spdvals <- spdvals/sum(spdvals) }
    if (rescale){ spdvals <- reScale(spdvals) }
    if (any(is.na(ylim))){ ylim <- c(0,max(spdvals)*1.1) }
    if (calendar=="BP"){
        plotyears <- spd$grid$calBP
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(plotyears),min(plotyears)) }
    } else if (calendar=="BCAD"){
        plotyears <- 1950-spd$grid$calBP
        xlabel <- "Years BC/AD"
        if (any(is.na(xlim))){ xlim <- c(min(plotyears),max(plotyears)) }
    } else {
        stop("Unknown calendar type")
    }
    if (xaxt=='n'){ xlabel <- "" }
    if (type=="standard"){
        par(xaxs="i")
        par(yaxs="i")
        plot(plotyears, spdvals, xlim=xlim, ylim=ylim, type="l", col="white", ylab="", xlab=xlabel, xaxt=xaxt, yaxt=yaxt)
        polygon(c(plotyears,rev(plotyears)),c(spdvals,rep(0,length(spdvals))),border=border.p, col=fill.p)
    } else if (type=="simple"){
        plot(plotyears, spdvals, xlim=xlim, ylim=ylim, type="l", ylab="", xlab=xlabel, xaxt=xaxt, yaxt=yaxt,...)
    }
}

plot.CalGrid <- function(x, calendar="BP", fill.p="grey50", border.p=NA, xlim=NA, ylim=NA, cex.lab=0.75, cex.axis=cex.lab, mar=c(4,4,1,3),...){

    yearsBP <- x$calBP
    prob <- x$PrDens
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
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
        xlabel <- "Years cal BP"
    } else if (calendar=="BCAD"){
        plotyears <- 1950-yearsBP
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
        xlabel <- "Years BC/AD"
    } else {
        stop("Unknown calendar type")
    }
    yvals <- c(0,prob,0,0)
    if (calendar=="BP"){
        xrng <- c(max(plotyears)+50, min(plotyears)-50)
        xticks <- 100*(xrng%/%100 + as.logical(xrng%%100))
        xticks <- seq(xticks[1]-100, xticks[2], -100)
    } else {
        xrng <- c(min(plotyears)-50, max(plotyears)+50)
        xticks <- 100*(xrng%/%100 + as.logical(xrng%%100))
        xticks <- seq(xticks[1]-100, xticks[2], 100)
    }
    if (is.na(ylim[1])){ ylim <- c(0,max(yvals*1.1)) }
    if (is.na(xlim[1])){ xlim <- xrng }
    par(mar=mar) #c(bottom, left, top, right)
    par(cex.lab=cex.lab)
    plot(xvals,yvals, type="n", xlab=xlabel, ylab="", xlim=xlim, ylim=ylim, xaxt='n', yaxt='n', cex.axis=cex.axis,...)
    axis(1, at=xticks, labels=abs(xticks), las=2, cex.axis=cex.axis)
    axis(4, cex.axis=cex.axis)
    polygon(xvals,yvals, col=fill.p, border=border.p)
}

plot.UncalGrid <- function(x, type="adjusted", fill.p="grey50", border.p=NA, xlim=NA, ylim=NA, cex.lab=0.75, cex.axis=cex.lab, mar=c(4,4,1,3),...){

    years <- x$CRA
    if (type=="adjusted"){
        prob <- x$PrDens
    } else if (type=="raw"){
        prob <- x$Raw
    } else {
        stop("Currently, type must be 'raw' or 'adjusted'.")
    }
    if (years[1] > years[length(years)]){
        years <- rev(years)
        prob <- rev(prob)
    }
    plotyears <- years
    xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
    xlabel <- "C14 Years"
    yvals <- c(0,prob,0,0)
    xrng <- c(max(plotyears)+50, min(plotyears)-50)
    if (dist(xrng)>10000){
        xticks <- 1000*(xrng%/%1000 + as.logical(xrng%%1000))
        xticks <- seq(xticks[1]-1000, xticks[2], -1000)
    } else {
        xticks <- 100*(xrng%/%100 + as.logical(xrng%%100))
        xticks <- seq(xticks[1]-100, xticks[2], -100)
    }
    if (is.na(xlim[1])){ xlim <- xrng }
    if (is.na(ylim[1])){ ylim <- c(0,max(yvals*1.1)) }
    par(mar=mar) #c(bottom, left, top, right)
    par(cex.lab=cex.lab)
    plot(xvals,yvals, type="n", xlab=xlabel, ylab="", xlim=xlim, ylim=ylim, xaxt='n', yaxt='n', cex.axis=cex.axis,...)
    axis(1, at=xticks, labels=abs(xticks), las=2, cex.axis=cex.axis)
    axis(4, cex.axis=cex.axis)
    polygon(xvals,yvals, col=fill.p, border=border.p)
}

plot.SpdPermTest <- function(test, focalm="1", calendar="BP", xlim=NA, ylim=NA, col.obs="black", lwd.obs=0.5, xaxs="i", yaxs="i", bbty="f", drawaxes=TRUE, ...){

    obs <- test$observed[[focalm]]
    if (calendar=="BP"){
        obs$Years <- obs$calBP
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
    } else if (calendar=="BCAD"){
        obs$Years <- 1950-obs$calBP
        xlabel <- "Years BC/AD"
        if (any(is.na(xlim))){ xlim <- c(min(obs$Years),max(obs$Years)) }
    } else {
        stop("Unknown calendar type")
    }
    envelope <- test$envelope[[focalm]]
    if (any(is.na(ylim))){ ylim <- c(0, max(envelope[,2], obs$PrDens)*1.1) }
    booms <- which(obs$PrDens>envelope[,2])
    busts <- which(obs$PrDens<envelope[,1])
    baseline <- rep(0,nrow(obs))
    if (drawaxes & bbty != "n"){
        plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab=xlabel, ylab="Summed Probability", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
        axis(side=1,padj=-1)
        axis(side=2,padj=1)
    } else if (bbty != "n"){
        plot(obs$Years,obs$PrDens, xlim=xlim, ylim=ylim, xlab="", ylab="", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
    }
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
            boomBlocks[[counter]][[2]] <- obs[x,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (boomPlot[x]>0){
                boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[x])
                boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],obs[x,"Years"])
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
            bustBlocks[[counter]][[2]] <- obs[x,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (bustPlot[x]>0){
                bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[x])
                bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],obs[x,"Years"])
            }
            if (bustPlot[x]==0){
                state <- "off"
            }
        }    
    }
    if (length(booms)>0){
        for (x in 1:length(boomBlocks)){
            if (bbty=="f"){
                polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(+100,length(boomBlocks[[x]][[1]])),rep(-100,length(boomBlocks[[x]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }
    if (length(busts)>0){
        for (x in 1:length(bustBlocks)){
            if (bbty=="f"){
                polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(+100,length(bustBlocks[[x]][[1]])),rep(-100,length(bustBlocks[[x]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }  
    if (bbty != "n"){
        polygon(x=c(obs[,"Years"], rev(obs[,"Years"])), y=c(envelope[,1], rev(envelope[,2])), col=rgb(0,0,0,0.2), border=NA)
        box()
    }
    if (drawaxes & bbty != "n"){
        axis(side=1, at=seq(max(obs[,"Years"]), min(obs[,"Years"]),-100), labels=NA, tck=-0.01)
    }
    bbp <- list(booms=boomBlocks, busts=bustBlocks)
    class(bbp) <- c("BBPolygons",class(bbp))
    if (bbty %in% c("n","b")){ return(bbp) }
}

bbpolygons <- function(x, baseline=1, height=1, calendar="BP", border=NA, bg=NA, col.boom=rgb(0.7,0,0,0.2), col.bust=rgb(0,0,0.7,0.2), border.boom=NA, border.bust=NA){

    if (!calendar %in% c("BP","BCAD")){ stop("Unknown calendar type") }
    if (!grepl("BBPolygons",class(x)[1])){
        stop("Input must be of class BBPolygons.")
    } else {
        boomBlocks <- x$booms
        bustBlocks <- x$busts
        plotrng <- par("usr") #c(x1, x2, y1, y2)
        width <- (plotrng[4]-plotrng[3]) * height
        if (baseline=="top"){
            realbase <- plotrng[4]-width
        } else if (baseline=="bottom"){
            realbase <- plotrng[3]+width
        } else {
            realbase <- plotrng[3] + ((plotrng[4]-plotrng[3]) * baseline)
        }
        plotymin <- realbase-width
        plotymax <- realbase+width
        if (!is.na(bg)){
            polygon(c(plotrng[1], plotrng[1], plotrng[2], plotrng[2], plotrng[1]),c(plotymin, plotymax, plotymax, plotymin, plotymin), col=bg, border=border)
        }
        if (length(boomBlocks)>0){
            for (x in 1:length(boomBlocks)){
                if (calendar=="BP"){
                    polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(realbase+width,length(boomBlocks[[x]][[1]])),rep(realbase-width,length(boomBlocks[[x]][[1]]))),col=col.boom, border=border.boom)
                } else {
                    polygon(c((1950-boomBlocks[[x]][[2]]),rev((1950-boomBlocks[[x]][[2]]))),c(rep(realbase+width,length(boomBlocks[[x]][[1]])),rep(realbase-width,length(boomBlocks[[x]][[1]]))),col=col.boom, border=border.boom)
                }
            }
        }
        if (length(bustBlocks)>0){
            for (x in 1:length(bustBlocks)){
                if (calendar=="BP"){
                    polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(realbase+width,length(bustBlocks[[x]][[1]])),rep(realbase-width,length(bustBlocks[[x]][[1]]))),col=col.bust, border=border.bust)
                } else {
                    polygon(c((1950-bustBlocks[[x]][[2]]),rev((1950-bustBlocks[[x]][[2]]))),c(rep(realbase+width,length(bustBlocks[[x]][[1]])),rep(realbase-width,length(bustBlocks[[x]][[1]]))),col=col.bust, border=border.bust)
                }
            }
        }
        if(!is.na(border)){
            polygon(c(plotrng[1], plotrng[1], plotrng[2], plotrng[2], plotrng[1]),c(plotymin, plotymax, plotymax, plotymin, plotymin), border=border)
        }
    }
}
