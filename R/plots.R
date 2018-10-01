#' @title Plot calibrated dates
#'
#' @description Plot calibrated radiocarbon dates.
#' @param x \code{CalDates} class object containing calibrated radiocarbon dates.
#' @param ind Number indicating the index value of the calibrated radiocarbon date to be displayed. Default is 1.
#' @param label (optional) Character vector to be shown on the top-right corner of the display window.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param type Either \code{'standard'} or \code{'auc'}. If set to \code{'auc'}, displays both the normalised (dashed line) and unnormalised curves. Default is \code{'standard'}.
#' @param xlab (optional) Label for the x axis. If unspecified the default setting will be applied ("Year BP" or "Year BC/AD"). 
#' @param ylab (optional) Label for the y axis. If unspecified the default setting will be applied ("Radiocarbon Age"). 
#' @param axis4  Logical value indicating whether an axis of probabilities values should be displayed. Default is TRUE. 
#' @param HPD Logical value indicating whether intervals of higher posterior density should be displayed. Default is FALSE.
#' @param credMass A numerical value indicating the size of the higher posterior density interval. Default is 0.95 (i.e. 95\%).
#' @param customCalCurve A three column data.frame or matrix that allows you to pass and plot a custom calibration curve if you used one during calibration. You can currently only provide one such custom curve which is used for all dates.
#' @param ... Additional arguments affecting the plot. 
#'
#' @seealso \code{\link{calibrate}}
#'
#' @examples
#' x <- calibrate(x=c(3402,3490,4042),errors=c(20,20,30))
#' plot(x) #display the first date
#' plot(x,2) #displays the second date
#' plot(x,3, calendar="BCAD", HPD=TRUE) #display in BC/AD with higher posterior density interval
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export  

plot.CalDates <- function(x, ind=1, label=NA, calendar="BP", type="standard", xlab=NA, ylab=NA, axis4=TRUE, HPD=FALSE, credMass=0.95, customCalCurve=NA,...){

    types <- c("standard", "simple", "auc")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    if (length(x$calmatrix)>1){
        grd <- data.frame(calBP=as.numeric(row.names(x$calmatrix)),PrDens=x$calmatrix[,ind])
        grd <- grd[grd$PrDens >0,]
        yearsBP <- grd$calBP
        prob <- grd$PrDens
    } else {
        yearsBP <- x$grids[[ind]]$calBP
        prob <- x$grids[[ind]]$PrDens
    }
    cra <- x$metadata$CRA[ind]
    error <- x$metadata$Error[ind]
    calendars <- c("BP","BCAD")
    if (!calendar %in% calendars){
        stop("The calendar you have chosen is not currently an option.")
    }        
    if (yearsBP[1] > yearsBP[length(yearsBP)]){
        yearsBP <- rev(yearsBP)
        prob <- rev(prob)
    }
    yvals <- c(0,prob,0,0)
    if (calendar=="BP"){
        plotyears <- yearsBP
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
        if (is.na(xlab)){ xlabel <- "Years cal BP" } else { xlabel <- xlab } 
    } else if (calendar=="BCAD"){
        plotyears <- BPtoBCAD(yearsBP)
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
        if (is.na(xlab)){ 
		xlabel <- "Years BC/AD"
		if (all(range(plotyears)<0)) {xlabel <- "Years BC"}
		if (all(range(plotyears)>0)) {xlabel <- "Years AD"}
	} else { xlabel <- xlab }       
    } else {
        stop("Unknown calendar type")
    }
    xrng <- c(min(xvals[yvals>0.000001])-50,max(xvals[yvals>0.000001])+50)
    if (calendar=="BP"){ xlim <- rev(xrng) } else { xlim <- xrng }
    xticks <- 100*(xrng%/%100 + as.logical(xrng%%100))
    if (calendar=="BP"){
        xticks <- seq(xticks[1], xticks[2]-100, 100)
    } else {
        xticks <- seq(xticks[1]-100, xticks[2], 100)
    }
    yrng <- c(min(yvals[yvals>0]),max(yvals[yvals>0])+(max(yvals[yvals>0])*2))
    par(cex.lab=0.75)
    plot(xvals,yvals, type="n", xlab=xlabel, ylab="", ylim=yrng, xlim=xlim, xaxt='n', yaxt='n', cex.axis=0.75,...)
    
    xticksLab <- xticks
    if (calendar=="BCAD")
    {
      if (any(xticksLab==0)){xticksLab[which(xticksLab==0)]=1}
      xticks[which(xticks>1)]=xticks[which(xticks>1)]-1
    }
    axis(1, at=xticks, labels=abs(xticksLab), las=2, cex.axis=0.75)
    
    if (axis4){ axis(4, cex.axis=0.75) }
    if (!HPD){
    polygon(xvals,yvals, col="grey50", border="grey50")
    } else {
    polygon(xvals,yvals, col="grey82", border="grey82")
    hdres <- hpdi(x,credMass=credMass)[[ind]]
    if(calendar=="BCAD"){hdres=1950-hdres}
	for (i in 1:nrow(hdres))
	{
	 index <- which(xvals%in%hdres[i,1]:hdres[i,2])
         polygon(c(xvals[index],xvals[index[length(index)]],xvals[index[1]]),c(yvals[index],0,0), col="grey50", border="grey50")
	}
    }

    if (type=="standard" | type=="auc"){
        if (type=="auc"){
            lines(xvals, yvals/sum(yvals), col="black", lty="dotted")
        }
        par(new=TRUE)
        cradf1 <- data.frame(CRA=50000:0,Prob=dnorm(50000:0, mean=cra, sd=error))
        cradf1 <- cradf1[cradf1$Prob>0.0001,]
        ylim <- c(cra-(12*error),cra+(8*error))    
        cradf1$RX <- reScale(cradf1$Prob, to=c(xlim[1],(xlim[1]+diff(xlim)*0.33)))
        yticks <- ylim[1]:ylim[2]
        yticks <- yticks[yticks %% 200 == 0]
        plot(cradf1$RX,cradf1$CRA,type="l", axes=FALSE, xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, col=rgb(144,238,144,120,maxColorValue=255))
        polygon(c(cradf1$RX,rev(cradf1$RX)),c(cradf1$CRA,rep(xlim[1],length(cradf1$CRA))), col=rgb(144,238,144,80,maxColorValue=255), border=NA)
        axis(side=2, at=yticks, labels=abs(yticks),las=2, cex.axis=0.75)
        if (is.na(ylab)){
            mtext(side=2, line=3, "Radiocarbon Age", cex=0.75)
        } else {
            mtext(side=2, line=3, xlab, cex=0.75)
        }
        calcurvemetadata <- x$metadata$CalCurve[ind]
        calcurvecheck <- TRUE
        if (calcurvemetadata == "custom" & !class(customCalCurve) %in% c("data.frame","matrix")){
            calcurvecheck <- FALSE
        }
        if (calcurvecheck){
            if (calcurvemetadata == "custom"){
                cc <- as.data.frame(customCalCurve)[,1:3]
                names(cc) <- c("BP","CRA","Error")
            } else {
                calCurveFile <- paste(system.file("extdata", package="rcarbon"), "/", calcurvemetadata,".14c", sep="")
                options(warn=-1)
                cc <- readLines(calCurveFile, encoding="UTF-8")
                cc <- cc[!grepl("[#]",cc)]
		cc.con <- textConnection(cc)
                cc <- read.csv(cc.con, header=FALSE, stringsAsFactors=FALSE)
		close(cc.con)
                options(warn=0)
                names(cc) <- c("BP","CRA","Error","D14C","Sigma")
            }
            if (calendar=="BCAD"){
                tmp <- (xrng-1950)*-1
                cc$RX <- 1950-cc$BP
            } else {
                tmp <- xrng
                cc$RX <- cc$BP
            }
            cc <- cc[cc$BP >= min(tmp) & cc$BP < max(tmp),] 
            cc$Hi <- cc$CRA + 10
            cc$Lo <- cc$CRA - 10
            ccbox <- c((max(yrng)*0.2),(max(yrng)*0.9))
            polygon(c(cc$RX,rev(cc$RX)),c(cc$Hi,rev(cc$Lo)), col=rgb(255,140,0,120,maxColorValue=255), border=NA)
            lines(cc$RX,cc$CRA, col=rgb(255,140,0,60,maxColorValue=255))
        }
    }
    if (!is.na(label)){
        legend("topright", label, bty="n", cex=0.75)
    }
}



#' @title Plot result of Monte-Carlo simulation of observed versus modelled SPDs
#'
#' @description The function visualises the observed summed probability distribution of radiocarbon dates along with a simulation envelope for the null model and regions of positive and negative deviation.
#'
#' @param x A \code{SpdModelTest} class object generated using the \code{\link{modelTest}} function.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param ylim the y limits of the plot.
#' @param col.obs Line colour for the observed SPD
#' @param lwd.obs Line width for the observed SPD
#' @param xaxs The style of x-axis interval calculation (see \code{\link{par}})
#' @param yaxs The style of y-axis interval calculation (see \code{\link{par}})
#' @param bbty Display options; one between \code{'b'},\code{'n'},and \code{'f'}. See details below.
#' @param drawaxes A logical value determining whether the axes should be displayed or not. Default is TRUE.
#' @param ... Additional arguments affecting the plot

#' @details The argument \code{bbty} controls the display options of the Monte-Carlo Test. Default settings (\code{bbty='f'}) displays the observed SPD (solid black line), the simulation envelope of the fitted model (shaded grey polygon) and regions of significance positive (red semi-transparent rectangle) and negative (blue semi-transparent rectangle) deviation. The option \code{bbty='b'} removes the regions of positive/negative deviations, whilst the option \code{bbty='n'} displays the simulation envelope on existing plot. 

#' @seealso \code{\link{modelTest}}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export 

plot.SpdModelTest <- function(x, calendar="BP", ylim=NA, xlim=NA, col.obs="black", lwd.obs=0.5, xaxs="i", yaxs="i", bbty="f", drawaxes=TRUE, ...){

    obs <- x$result[,1:2]
    if (calendar=="BP"){
        obs$Years <- obs$calBP
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
    } else if (calendar=="BCAD"){
        obs$Years <- BPtoBCAD(obs$calBP)
	if (all(range(obs$Years)<0)){xlabel <- "Years BC"}
	if (all(range(obs$Years)>0)){xlabel <- "Years AD"}
        if (any(is.na(xlim))){xlim <- c(min(obs$Years),max(obs$Years)) }
    } else {
        stop("Unknown calendar type")
    }    
    envelope <- x$result[,3:4]
    if (any(is.na(ylim))){ ylim <- c(0, max(envelope[,"hi"], obs$PrDens)*1.1) }
    booms <- which(obs$PrDens>envelope[,2])
    busts <- which(obs$PrDens<envelope[,1])
    baseline <- rep(NA,nrow(obs))
    if (drawaxes & bbty != "n"){
        plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab=xlabel, ylab="Summed Probability", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
    } else if (bbty != "n"){
        plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab="", ylab="", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
    }
    if (drawaxes){
	box()
    axis(side=2)}
    boomPlot <- baseline
 	if (length(booms)>0){ boomPlot[booms]=obs[booms,2] }
    bustPlot <- baseline
    if (length(busts)>0){ bustPlot[busts]=obs[busts,2] }           
    boomBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (i in 1:length(boomPlot)){
        if (!is.na(boomPlot[i])&state=="off"){
            counter <- counter+1
            boomBlocks <- c(boomBlocks,vector("list",1))
            boomBlocks[[counter]] <- vector("list",2)
            boomBlocks[[counter]][[1]] <- boomPlot[i]
            boomBlocks[[counter]][[2]] <- obs[i,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (!is.na(boomPlot[i])){
                boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[i])
                boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],obs[i,"Years"])
            }
            if (is.na(boomPlot[i])){
                state <- "off"
            }
        }   
    }
    bustBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (i in 1:length(bustPlot)){
        if (!is.na(bustPlot[i])&state=="off"){
            counter <- counter+1
            bustBlocks <- c(bustBlocks,vector("list",1))
            bustBlocks[[counter]] <- vector("list",2)
            bustBlocks[[counter]][[1]] <- bustPlot[i]
            bustBlocks[[counter]][[2]] <- obs[i,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (!is.na(bustPlot[i])){
                bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[i])
                bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],obs[i,"Years"])
            }
            if (is.na(bustPlot[i])){
                state <- "off"
            }
        }   
    }
    if (length(booms)>0){
        for (i in 1:length(boomBlocks)){
            if (bbty=="f"){
                polygon(c(boomBlocks[[i]][[2]],rev(boomBlocks[[i]][[2]])),c(rep(+100,length(boomBlocks[[i]][[1]])),rep(-100,length(boomBlocks[[i]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }  
    if (length(busts)>0){
        for (i in 1:length(bustBlocks)){
            if (bbty=="f"){
                polygon(c(bustBlocks[[i]][[2]],rev(bustBlocks[[i]][[2]])),c(rep(+100,length(bustBlocks[[i]][[1]])),rep(-100,length(bustBlocks[[i]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }  
    polygon(x=c(obs[,"Years"],rev(obs[,"Years"])),y=c(envelope[,1],rev(envelope[,2])),col=rgb(0,0,0,0.2),border=NA)
    if (drawaxes & bbty != "n" & calendar=="BP"){
	rr <- range(pretty(obs[,"Years"]))    
        axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
        axis(side=1,at=pretty(obs[,"Years"]),labels=abs(pretty(obs[,"Years"])))
    } else if (drawaxes & bbty != "n" & calendar=="BCAD"){
	yy <-  obs[,"Years"]
       
	rr <- range(pretty(yy))    
        prettyTicks <- seq(rr[1],rr[2],+100)
	prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
        axis(side=1,at=prettyTicks, labels=NA,tck = -.01)

        py <- pretty(yy)
	pyShown <- py
	if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
	py[which(py>1)] <-  py[which(py>1)]-1
	axis(side=1,at=py,labels=abs(pyShown))
    }

    bbp <- list(booms=boomBlocks, busts=bustBlocks)
    class(bbp) <- c("BBPolygons",class(bbp))
    if (bbty %in% c("n","b")){ return(bbp) }
}

#' @title Plot the median values of calibrated radiocarbon dates or bins 
#'
#' @description Plot the median values of multiple calibrated radiocarbon dates or bins in a barcode-like strip.
#' @param x A vector containing median values obtained from \code{\link{medCal}} or \code{\link{binMed}}  
#' @param yrng y-axis range of the bars.
#' @param width width of the bars (optional) 
#' @param col color of the bars
#' @param border the color to draw the border. Use border = NA to omit borders.
#' @param ... Additional arguments affecting the plot
#'
#' 
#' @seealso \code{\link{medCal}}; \code{\link{binMed}}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @examples
#'\dontrun{
#' #Load EUROEVOL Data
#' data(euroevol)
#' 
#' #Subset Danish Dates
#' denmark <- subset(euroevol,Country=="Denmark")
#'
#' #Calibrate and Bin
#' denmarkDates <- calibrate(x=denmark$C14Age,errors=denmark$C14SD) 
#' denmarkBins <- binPrep(sites=denmark$SiteID,ages=denmark$C14Age,h=200) #200 years bin size
#'
#' #Compute median date for each bin
#' bm <- binMed(x=denmarkDates,bins=denmarkBins)
#'
#' #Compute median date for each date
#' dm <- medCal(denmarkDates)
#' 
#' #Compute SPD 
#' denmarkSPD <- spd(x=denmarkDates,bins=denmarkBins,timeRange=c(10000,4000))
#' 
#' #Plot SPD and barCodes of median dates
#' plot(denmarkSPD,runm=200)
#' barCodes(dm,yrng=c(0,0.01)) 
#'
#' #Plot SPD and barCodes of median bins in BC/AD
#' plot(denmarkSPD,runm=200,calendar="BCAD")
#' barCodes(BPtoBCAD(bm),yrng=c(0,0.01)) 
#'}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export

barCodes <- function(x, yrng=c(0,0.03), width=20, col=rgb(0,0,0,25,maxColorValue=255), border=NA, ...){
    barcodes <- x
    halfbw <- width/2
    for (a in 1:length(barcodes)){
        polygon(x=c(barcodes[a]-halfbw,barcodes[a]-halfbw,barcodes[a]+halfbw,barcodes[a]+halfbw,barcodes[a]-halfbw),y=c(yrng[1],yrng[2],yrng[2],yrng[1],yrng[1]), border=border, col=col, ...)
    }
}


# 
# crossHairs <- function(x, pch.pts=19, cex.pts=1, fixXorder=FALSE, rescaleY=FALSE,...){
# 
#     if (!"quickMarks" %in% class(x)){
#         stop("Input must be of class \"quickMarks\"")
#     }
#     if (rescaleY){
#         cra <- reScale(x$CRA)
#         error <- x$Error / (max(x$CRA)-min(x$CRA))
#     } else {
#         cra <- x$CRA
#         error <- x$Error
#     }
#     if (fixXorder){
#         xstart <- x$q68s *-1
#         xend <- x$q68e *-1
#         xmed <- x$qMed *-1
#     } else {
#         xstart <- x$q68s
#         xend <- x$q68e
#         xmed <- x$qMed
#     }
#     for (a in 1:nrow(x)){
#         lines(c(xstart[a],xend[a]), c(cra[a],cra[a]), ...)
#         lines(c(xmed[a],xmed[a]), c(cra[a]-error[a],cra[a]+error[a]), ...)
#         points(xmed[a],cra[a], pch=pch.pts, cex=cex.pts, ...)
#     }
# }
# 

#' @title Plot a summed probability distribution
#'
#' @description Plot a summed probability distribution (SPD) of radiocarbon dates 
#' @param x A \code{CalSPD} class object.
#' @param runm A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied. Default is NA  
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param type Either \code{'standard'} or \code{'simple'}. The former visualise the SPD as an area graph, while the latter as line chart. 
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param ylim the y limits of the plot.
#' @param ylab (optional) Label for the y axis. If unspecified the default setting will be applied ("Summed Probability") 
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPD is normalised to sum to unity. 
#' @param rescale  A logical variable indicating whether the SPD should be rescaled to range 0 to 1.
#' @param fill.p Fill colour for the SPD
#' @param border.p Border colour for the SPD
#' @param xaxt Whether the x-axis tick marks should be displayed (\code{xaxt='s'}, default) or not (\code{xaxt='n'}).
#' @param yaxt Whether the y-axis tick marks should be displayed (\code{xaxt='s'}, default) or not (\code{xaxt='n'}).
#' @param add Whether or not the new graphic should be added to an existing plot.
#' @param ... Additional arguments affecting the plot
#'
#'
#'
#' @seealso \code{\link{spd}}; \code{\link{plot.CalGrid}}
#' @examples
#' \dontrun{
#' data(emedyd)
#' levant <- emedyd[emedyd$Region=="1"|emedyd$Region=="2",]
#' bins <- binPrep(levant$SiteName, levant$CRA, h=50)
#' x <- calibrate(levant$CRA, levant$Error, normalised=FALSE)
#' spd.levant <- spd(x, bins=bins, timeRange=c(17000,8000))
#' spd.northernlevant <- spd(x[levant$Region=="2"], bins=bins[levant$Region=="2"],
#' timeRange=c(17000,8000))
#' plot(spd.levant, runm=50, xlim=c(16000,9000))
#' plot(spd.northernlevant, runm=50, add=TRUE, fill.p="black")
#' legend("topleft", legend=c("All Levant dates","Northern Levant only"), 
#' fill=c("grey75","black"), border=NA)
#' plot(spd.levant, runm=50, xlim=c(16000,9000), type="simple")
#' plot(spd.northernlevant, runm=50, col="red", type="simple", add=TRUE)
#'}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export 

plot.CalSPD <- function(x, runm=NA, calendar="BP", type="standard", xlim=NA, ylim=NA, ylab="Summed Probability", spdnormalised=FALSE, rescale=FALSE, fill.p="grey75", border.p=NA, xaxt='s', yaxt='s', add=FALSE,...){

    types <- c("standard","simple")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    spdvals <- x$grid$PrDens
    if (!is.na(runm)){ spdvals <- runMean(spdvals, runm, edge="fill") }
    if (spdnormalised){ spdvals <- spdvals/sum(spdvals) }
    if (rescale){ spdvals <- reScale(spdvals) }
    if (any(is.na(ylim))){ ylim <- c(0,max(spdvals)*1.1) }
    if (calendar=="BP"){
        plotyears <- x$grid$calBP
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(plotyears),min(plotyears)) }
    } else if (calendar=="BCAD"){
        plotyears <- BPtoBCAD(x$grid$calBP)
        xlabel <- "Years BC/AD"
	if (all(range(plotyears)<0)){xlabel <- "Years BC"}
	if (all(range(plotyears)>0)){xlabel <- "Years AD"}
        if (any(is.na(xlim))){ xlim <- c(min(plotyears),max(plotyears)) }
    } else {
        stop("Unknown calendar type")
    }
    if (xaxt=='n'){ xlabel <- "" }
    if (yaxt=='n'){ ylabel <- "" } else { ylabel <- ylab }
    if (type=="standard"){
        par(xaxs="i")
        par(yaxs="i")
        if (!add){
            plot(plotyears, spdvals, xlim=xlim, ylim=ylim, type="l", col="white", ylab=ylabel, xlab=xlabel, xaxt="n", yaxt=yaxt)
        }
        polygon(c(plotyears,rev(plotyears)),c(spdvals,rep(0,length(spdvals))),border=border.p, col=fill.p)
    } else if (type=="simple"){
        if (!add){
            plot(plotyears, spdvals, xlim=xlim, ylim=ylim, type="l", ylab="", xlab=xlabel, xaxt="n", yaxt=yaxt, ...)
        } else {
            lines(plotyears, spdvals, xlim=xlim, ylim=ylim, type="l", ylab="", xlab=xlabel, ...)
        }
    }
    box()
    if (calendar=="BP" & xaxt!="n"){
	rr <- range(pretty(plotyears))    
        axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
        axis(side=1,at=pretty(plotyears),labels=abs(pretty(plotyears)))
    } else if (calendar=="BCAD" & xaxt!="n"){
	yy <-  plotyears
        rr <- range(pretty(yy))    
        prettyTicks <- seq(rr[1],rr[2],+100)
	prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
        axis(side=1,at=prettyTicks, labels=NA,tck = -.01)
        py <- pretty(yy)
	pyShown <- py
	if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
	py[which(py>1)] <-  py[which(py>1)]-1
	axis(side=1,at=py,labels=abs(pyShown))
    }
}


#' @title Plot a summed probability distribution (from a CalGrid object)
#'
#' @description Plot a summed radiocarbon probability distribution. This is a basic function for plotting SPDs that have been constructed manually or by calibrating a summed or otherwise irregular CRA grid. In most instances, it is sensible to use \code{plot.CalSPD} instead.
#' 
#' @param x A "CalGrid" class object of summed probabilities per calendar year BP.
#' @param runm A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied. Default is NA  
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param fill.p Fill colour of the polygon depicting the summed probability distribution.
#' @param border.p Border colour of the polygon depicting the summed probability distribution.
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param ylim Adjust y axis limits (otherwise sensible default).
#' @param cex.lab Size of label text.
#' @param cex.axis Size of axis text.
#' @param mar Adjust margins around plot.
#' @param add Whether or not the new graphic should be added to an existing plot.
#' @param ... Additional arguments affecting the plot
#'
#' @seealso \code{\link{spd}}; \code{\link{plot.CalSPD}}
#' @examples
#' data(euroevol)
#' mycaldates <- calibrate(euroevol[1:10,"C14Age"], euroevol[1:10,"C14SD"], normalised=FALSE)
#' myspd <- spd(mycaldates, timeRange=c(8000,2000))
#' plot(myspd) #ordinary plot using \code{plot.CalSPD}
#' plot(myspd$grid) #working plot using the internal CalGrid object
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export 

plot.CalGrid <- function(x, runm=NA, calendar="BP", fill.p="grey50", border.p=NA, xlim=NA, ylim=NA, cex.lab=0.75, cex.axis=cex.lab, mar=c(4,4,1,3), add=FALSE,...){

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
        plotyears <- BPtoBCAD(yearsBP)
        xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
        xlabel <- "Years BC/AD"
	if (all(range(plotyears)<0)){xlabel<-"Years BC"}
	if (all(range(plotyears)>0)){xlabel<-"Years AD"}
    } else {
        stop("Unknown calendar type")
    }
    yvals <- c(0,prob,0,0)
    if (calendar=="BP"){
        xrng <- c(max(plotyears)+50, min(plotyears)-50)
        xticks <- 200*(xrng%/%100 + as.logical(xrng%%200))
        xticks <- seq(xticks[1]-200, xticks[2], -200)
    } else {
        xrng <- c(min(plotyears)-50, max(plotyears)+50)
        xticks <- 200*(xrng%/%100 + as.logical(xrng%%200))
        xticks <- seq(xticks[1]-200, xticks[2], 200)
    }
    if (is.na(ylim[1])){ ylim <- c(0,max(yvals*1.1)) }
    if (is.na(xlim[1])){ xlim <- xrng }
    if (!add){
        par(mar=mar) #c(bottom, left, top, right)
        par(cex.lab=cex.lab)
        plot(xvals,yvals, type="n", xlab=xlabel, ylab="", xlim=xlim, ylim=ylim, xaxt='n', yaxt='n',axes=F, cex.axis=cex.axis,...)
    }
    xticksLab <- xticks
    if (calendar=="BCAD"){
        if (any(xticksLab==0)){xticksLab[which(xticksLab==0)]=1}
        xticks[which(xticks>1)]=xticks[which(xticks>1)]-1
    }
    if (!add){
        axis(1, at=xticks, labels=abs(xticksLab), las=2, cex.axis=cex.axis)
        axis(2, cex.axis=cex.axis)
    }
    if (!is.na(runm)){ yvals <- runMean(yvals, runm, edge="fill") }
    polygon(xvals,yvals, col=fill.p, border=border.p)
    box()
}


#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils

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




#' @title Plot result of mark permutation test of SPDs
#'
#' @description Visualises the observed SPD along with the simulation envelope generated from \code{\link{permTest}}, with regions of positive and negative deviations highlighted in red and blue.
#'
#' @param x A \code{SpdPermTest} class object. Result of random mark permutation test (see \code{\link{permTest}})
#' @param focalm Value specifying the name of the focal mark (group) to be plotted. 
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param ylim the y limits of the plot.
#' @param col.obs Line colour for the observed SPD
#' @param col.env Colour for the simulation envelope
#' @param lwd.obs Line width for the observed SPD
#' @param xaxs The style of x-axis interval calculation (see \code{\link{par}})
#' @param yaxs The style of y-axis interval calculation (see \code{\link{par}})
#' @param bbty Display options; one between \code{'b'},\code{'n'},and \code{'f'}. See details in \code{\link{plot.SpdModelTest}}.
#' @param drawaxes A logical value determining whether the axes should be displayed or not. Default is TRUE.
#' @param ... Additional arguments affecting the plot
#'
#' @seealso \code{\link{permTest}}; \code{\link{plot.SpdModelTest}};
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export 

plot.SpdPermTest <- function(x, focalm="1", calendar="BP", xlim=NA, ylim=NA, col.obs="black", col.env=rgb(0,0,0,0.2), lwd.obs=0.5, xaxs="i", yaxs="i", bbty="f", drawaxes=TRUE, ...){

    obs <- x$observed[[focalm]]
    if (calendar=="BP"){
        obs$Years <- obs$calBP
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
    } else if (calendar=="BCAD"){
        obs$Years <- BPtoBCAD(obs$calBP)
	if (all(range(obs$Years)<0)){xlabel <- "Years BC"}
	if (all(range(obs$Years)>0)){xlabel <- "Years AD"}
        if (any(is.na(xlim))){ xlim <- c(min(obs$Years),max(obs$Years)) }
    } else {
        stop("Unknown calendar type")
    }
    envelope <- x$envelope[[focalm]]
    if (any(is.na(ylim))){ ylim <- c(0, max(envelope[,2], obs$PrDens)*1.1) }
    booms <- which(obs$PrDens>envelope[,2])
    busts <- which(obs$PrDens<envelope[,1])
    baseline <- rep(NA,nrow(obs))
    if (drawaxes & bbty != "n"){
        plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab=xlabel, ylab="Summed Probability", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
        #axis(side=1,padj=-1)
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
    for (i in 1:length(boomPlot)){
        if (!is.na(boomPlot[i])&state=="off"){
            counter <- counter+1
            boomBlocks <- c(boomBlocks,vector("list",1))
            boomBlocks[[counter]] <- vector("list",2)
            boomBlocks[[counter]][[1]] <- boomPlot[i]
            boomBlocks[[counter]][[2]] <- obs[i,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (!is.na(boomPlot[i])){
                boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[i])
                boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],obs[i,"Years"])
            }
            if (is.na(boomPlot[i])){
                state <- "off"
            }
        }    
    }
    bustBlocks <- vector("list")
    counter <- 0
    state <- "off"
    for (i in 1:length(bustPlot)){
        if (!is.na(bustPlot[i])&state=="off"){
            counter <- counter+1
            bustBlocks <- c(bustBlocks,vector("list",1))
            bustBlocks[[counter]] <- vector("list",2)
            bustBlocks[[counter]][[1]] <- bustPlot[i]
            bustBlocks[[counter]][[2]] <- obs[i,"Years"]
            state <- "on"
        }
        if (state=="on"){
            if (!is.na(bustPlot[i])){
                bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[i])
                bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],obs[i,"Years"])
            }
            if (is.na(bustPlot[i])){
                state <- "off"
            }
        }    
    }
    if (length(booms)>0){
        for (i in 1:length(boomBlocks)){
            if (bbty=="f"){
                polygon(c(boomBlocks[[i]][[2]],rev(boomBlocks[[i]][[2]])),c(rep(+100,length(boomBlocks[[i]][[1]])),rep(-100,length(boomBlocks[[i]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }
    if (length(busts)>0){
        for (i in 1:length(bustBlocks)){
            if (bbty=="f"){
                polygon(c(bustBlocks[[i]][[2]],rev(bustBlocks[[i]][[2]])),c(rep(+100,length(bustBlocks[[i]][[1]])),rep(-100,length(bustBlocks[[i]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
            } else if (bbty %in% c("s","b","n")){
            } else {
                stop("Incorrect bbty argument.")
            }
        }
    }  
    if (bbty != "n"){
        polygon(x=c(obs[,"Years"], rev(obs[,"Years"])), y=c(envelope[,1], rev(envelope[,2])), col=col.env, border=NA)
        box()
    }
    #if (drawaxes & bbty != "n"){
    #    axis(side=1, at=seq(max(obs[,"Years"]), min(obs[,"Years"]),-100), labels=NA, tck=-0.01)
    #}

    if (drawaxes & bbty != "n" & calendar=="BP"){
	rr <- range(pretty(obs[,"Years"]))    
        axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
        axis(side=1,at=pretty(obs[,"Years"]),labels=abs(pretty(obs[,"Years"])))
    } else if (drawaxes & bbty != "n" & calendar=="BCAD"){
	yy <-  obs[,"Years"]
	rr <- range(pretty(yy))    
        prettyTicks <- seq(rr[1],rr[2],+100)
	prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
        axis(side=1,at=prettyTicks, labels=NA,tck = -.01)

        py <- pretty(yy)
	pyShown <- py
	if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
	py[which(py>1)] <-  py[which(py>1)]-1
	axis(side=1,at=py,labels=abs(pyShown))
    }
    bbp <- list(booms=boomBlocks, busts=bustBlocks)
    class(bbp) <- c("BBPolygons",class(bbp))
    if (bbty %in% c("n","b")){ return(bbp) }
}


#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils


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


#' @title Bin sensitivity Plot
#' 
#' @description Visually explores how choosing different values for \code{h} in the \code{\link{binPrep}} function affects the shape of the SPD.
#' 
#' @param x A \code{CalDates} class object containing calibrated radiocarbon dates.
#' @param y A vector containing the locations ID (e.g. site ID) of each calibrated date to be used for the binning process. 
#' @param h A vector of numbers containing values for the \code{h} parameter to be used in the \code{\link{binPrep}} function. 
#' @param timeRange A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param binning  Either \code{'CRA'} or \code{'calibrated'}. Indicate whether the binning should be carried usinig the  14C age or using the median calibrated date. Default is \code{'CRA'}.
#' @param raw A logical variable indicating whether all  SPDs should be returned or not. Default is FALSE.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param legend A logical variable indicating whether the legend should be displayed. Default is TRUE
#' @param ... Additional arguments to be passed to the \code{\link{spd}} function. 
#' 
#' @seealso \code{\link{binPrep}}; \code{\link{spd}}
#' @examples
#' \dontrun{
#' data(euroevol)
#' #subset Danish dates
#' denmark=subset(euroevol,Country=="Denmark")
#' denmarkDates=calibrate(x=denmark$C14Age,errors=denmark$C14SD)
#' binsense(x=denmarkDates,y=denmark$SiteID,h=seq(0,200,20),timeRange=c(10000,4000),runm=200)
#' }

#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export

binsense <- function(x,y,h,timeRange,calendar="BP",binning='CRA',raw=F,verbose=T,legend=T,...)
{
  if (!calendar %in% c("BP","BCAD")){ stop("Unknown calendar type") }
  if (!binning %in% c("CRA","calibrated")) {stop("binning should be either 'CRA' or 'calibrated'")} 	
  years <- timeRange[1]:timeRange[2]
  xlab <- "Years cal BP"
  coln <- numeric(length=length(h))
  xr <- timeRange
  if (calendar=="BCAD")
  {
   years <- BPtoBCAD(years)
   xlab <- "Years BC/AD"
   xr <- range(years)
   if (all(xr<0)){xlab <- "Years BC"}
   if (all(xr>0)){xlab <- "Years AD"}
  }

  res <- matrix(NA,nrow=length(years),ncol=length(h))
  craAges <- x$metadata$CRA
  
  if (length(craAges)!=length(y))
  {
    stop("x and y have different lengths (each calibrated date in x should have a matching location in y)")
  }


  if (verbose)
	 {
         print("Computing SPDs...")
	 flush.console()
         pb <- txtProgressBar(min = 1, max =length(h), style=3)
	 }
  for (b in 1:length(h))
    {
    if (verbose){setTxtProgressBar(pb, b)}	    
    if (binning == "CRA"){bins <- binPrep(sites=y,ages=craAges,h=h[b])}
    if (binning == "calibrated"){ bins <- binPrep(sites=y,ages=x,h=h[b])}
    spdtmp <- spd(x,bins= bins,timeRange=timeRange,spdnormalised=T,verbose=F,...)
    res[,b] <- spdtmp$grid$PrDens
    coln[b] <- paste("h.",h[b],sep="")
    }

  if (verbose)
	{ 
	close(pb)
        print("Done.") 
	}	
  if (legend==TRUE){layout(matrix(c(1,1,2,2),2,2),widths=c(1,0.2))}

  plot(years,res[,1],xlim=xr,ylim=range(res),type="n",xlab=xlab,ylab="normalised SPD",axes=F)
  axis(side=2)
  if (calendar=="BP") {axis(1)}
  if (calendar=="BCAD")
  {
   xticksAt=pretty(years)
   xticksLab=xticksAt
   if (any(xticksLab==0)){xticksLab[which(xticksLab==0)]=1}
   if (any(xticksAt>1)){xticksAt[which(xticksAt>1)]=xticksAt[which(xticksAt>1)]-1}
   axis(side=1,at=xticksAt,labels=abs(xticksLab))
  }  


  cl=colorRampPalette(c("indianred","royalblue"),alpha=0.5)

  minRange=apply(res,1,min)
  maxRange=apply(res,1,max)

  polygon(c(years,rev(years)),c(minRange,rev(maxRange)),col="darkgrey",border="NA")

  for (x in 1:length(h))
   {
    lines(years,res[,x],col=cl(length(h))[x],lwd=0.5)
   }

  if (legend==TRUE)
  { 
  par(mar=c(6,2,6,2))	  
  z=matrix(1:100,nrow=1)
  x=1
  y=seq(h[1],h[length(h)],len=100) # supposing 3 and 2345 are the range of your data
  image(x,y,z,col=cl(100),axes=FALSE,xlab="",ylab="")
  axis(2,padj=1,cex.axis=0.7)
  mtext(side=2,"h",line=2,las=2,cex=0.7)
  }
  box()
if (raw) {
colnames(res)=coln	
res=cbind.data.frame(calBP=timeRange[1]:timeRange[2],res)
}
}



#' @title Plot results of the local spatial permutation test of summed probability distributions.
#'
#' @description Displays local growth rates, p-values, and q-values retrieved from a \code{spatialTest} class object.
#'
#' @param x A \code{spatialTest} class object
#' @param index A numerical value indicating which transition to display. Ignored when \code{option="rawlegend"} or  \code{option="testlegend"}. Default is 1.
#' @param option Indicates what to display. Either "\code{raw}" (the local growth rate) or "\code{test}" (the test results, i.e. q and p values). 
#' @param breakRange A vector of length 2 defining the minimum and maximum values of growth rate to be displayed in the legend. If set to NA its computed from data range (default).
#' @param breakLength A numerical vector defining the number of breaks for growth rates to be displayed in the legend.
#' @param rd Number of decimal places of the growth rate to be displayed in the Legend
#' @param baseSize Numerical value giving the amount by which points should be magnified relative to the default settings in R. Default is 0.5
#' @param plim Threshold value for the p-values. Default is 0.05.
#' @param qlim Threshold value for the q-values. Default is 0.05.
#' @param legend Logical values specifying whether the legend should be displayed or not. Default is FALSE. 
#' @param legSize Numerical value giving the amount by which points should be magnified relative to the default settings in R for the Legend. Default is 1.
#' @param location A single keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center" to specify the location of the Legend. Default is "bottomright".
#' @param ... Graphical parameters to be passed to methods.
#'
#' @details
#' The function displays a distribution map of local growth rates (when \code{option="raw"}), q- and p-values (when \code{option="test"}), and the associated legends (when \code{option="rawlegend"} or  \code{option="testlegend"}).
#'
#' @seealso \code{\link{SPpermTest}}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export 


plot.spatialTest<-function(x,index=1,option,breakRange=NA,breakLength=7,rd=5,baseSize=0.5,plim=0.05,qlim=0.05,legend=FALSE,legSize=1,location="bottomright",...)
{
	if (!any(class(x)%in%c("spatialTest")))
	{
        stop("x is not a spatialTest class object")
	}

        if (is.null(index)&option%in%c("raw","test"))
	{
        stop("index value missing")
	}

	if (!option%in%c("raw","test","rawlegend","testlegend"))
	{
        stop(paste("The option ",option," is not available",sep=""))
	}

        locations=x$locations

	if(all(is.na(breakRange))){breakRange=range((x$rocaObs[,index]))} 

	if (option=="raw")
	{
	breaks=seq(breakRange[1],breakRange[2],length.out=breakLength)
	outbreak=c(-Inf,breaks,Inf)
	classes=cut(x$rocaObs[,index], outbreak,labels=F)
	cols=colorRampPalette(c("blue","white","red"))(breakLength+1)
	classes=cols[classes]
	plot(locations,col=classes,pch=20,cex=baseSize,...)
	if (legend)
	{
	breaks=round(seq(breakRange[1],breakRange[2],length.out=breakLength),rd)
	cols=colorRampPalette(c("blue","white","red"))(breakLength+1)
	breaksLab=numeric(breakLength+1)
	breaksLab[1]= paste("<",breaks[1])
	for (j in 2:c(breakLength+1))
	{
	 breaksLab[j] = paste(breaks[j-1],"to", breaks[j])
	 if (j==c(breakLength+1)) {breaksLab[j] = paste(">",breaks[length(breaks)])}
	}
	legend(location,legend=breaksLab,col=cols,pch=20,bg="white",cex=legSize)
	}
	}


	if (option=="test")
	{
	nBreaks=ncol(x$rocaObs)
	plusPoints=locations[which(x$pvalHi[,index]>=0.5),]
	minusPoints=locations[which(x$pvalHi[,index]<0.5),]

	# Set Base
	plot(locations,col=NA,xlab="",ylab="",axes=FALSE,...)
	points(plusPoints,col="darkgrey",pch=20,cex=baseSize)
	points(minusPoints,col="darkgrey",pch=20,cex=baseSize)


	# Set Positive
	positive.index=which(x$pvalHi[,index]<=plim)

	if (length(positive.index)>0)
		{
		positive=locations[positive.index,]
		points(positive,pch=20,col="orange",cex=baseSize)
		qpositive.index=which(x$qvalHi[,index]<=qlim&x$pvalHi[,index]<=plim) #Originally based on qvalHi
		if (length(qpositive.index)>0)
			{
				qpositive=locations[qpositive.index,]
				points(qpositive,pch=20,col="red",cex=baseSize)

			}
		}
	negative.index=which(x$pvalLo[,index]<=plim)

	if (length(negative.index)>0)
		{
		negative=locations[negative.index,]
		points(negative,pch=20,col="cornflowerblue",cex=baseSize)
		qnegative.index=which(x$qvalLo[,index]<=qlim&x$pvalLo[,index]<=plim) #Originally based on qvalLo
		if (length(qnegative.index)>0)
			{
				qnegative=locations[qnegative.index,]
				points(qnegative,pch=20,col="darkblue",cex=baseSize)

			}

		}

	if (legend)
	{
	legend(location,legend=c(paste("positive deviation (p<",plim,")",sep=""),paste("positive deviation (q<",qlim,")",sep=""),paste("negative deviation (p<",plim,")",sep=""),paste("negative deviation (q<",qlim,")",sep="")),pch=20, col=c("orange","red","cornflowerblue","darkblue"),bg="white",cex=legSize)
	}

	}

}



#' @title Plot \code{spdGG} class objects 
#'
#' @description Plot calibrated geometric growth rates.
#' @param x \code{spdGG} class object containing geometric growth rates.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param ... Additional arguments affecting the plot. 
#'
#' @seealso \code{\link{spd2gg}}
#'
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export  

plot.spdGG<- function(x,calendar="BP",...)
{
	breaks=x$breaks
	obs=x$sumblock
	res=x$geomg
	par(mar=c(4,4,4,4))
	nn = paste(breaks[-length(breaks)],breaks[-1],sep=" to ")
	xxlab="Years cal BP"
	if (calendar=="BCAD")
	{
		bcad.breaks=BPtoBCAD(breaks)
		nn = paste(abs(bcad.breaks[-length(bcad.breaks)]),abs(bcad.breaks[-1]),sep="-")
		xxlab="Years BC/AD"
		if (all(range(bcad.breaks)<0)){xxlab="Years BC"}
		if (all(range(bcad.breaks)>0)){xxlab="Years AD"}

	}
	barplot(x$sumblock,names.arg=nn,ylab="Summed Probability",,space=0,col="bisque3",border=NA,xlab=xxlab,...)
	par(new=T)
	xx = 1:c(length(nn)-1)
	plot(0,0,xlim=c(0,length(nn)),ylim=range(res),axes=FALSE,xlab="",ylab="",type="n")
	lines(xx,res,lwd=2,col="darkgreen")
	points(xx,res,pch=20,col="darkgreen")
	axis(4,col="darkgreen", col.axis="darkgreen")
	mtext(side=4,"geometric growth rate",col="darkgreen",line=2)
	abline(h=0,lty=2,col="blue")
}

