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
#' @param credMass A numerical value indicating the size of the higher posterior density interval. Default is 0.95.
#' @param customCalCurve A three column data.frame or matrix that allows you to pass and plot a custom calibration curve if you used one during calibration. You can currently only provide one such custom curve which is used for all dates.
#' @param col The primary fill color for the calibrated date distribution. 
#' @param col2 The secondary colour fill color for the calibrated date distribution, used for regions outside the higher posterior interval. Ignored when \code{HPD=FALSE}.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex. Default is adjusted to 0.75.
#' @param cex.xylab The magnification to be used for x and y labels relative to the current setting of cex. Default is adjusted to 0.75.
#' @param cex.label The magnification to be used for the label on the top right corner defined by the argument \code{label}. Default is adjusted to 0.75.
#' @param add if set to \code{TRUE} the calibrated date is displayed over the existing plot. Default is \code{FALSE}. 
#' @param ... Additional arguments affecting the plot. 
#'
#' @seealso \code{\link{calibrate}}
#'
#' @examples
#' x <- calibrate(x=c(3402,3490,4042),errors=c(20,20,30))
#' plot(x) #display the first date
#' plot(x,2) #displays the second date
#' plot(x,3, calendar="BCAD", HPD=TRUE) #display in BC/AD with higher posterior density interval
#' @method plot CalDates
#' @export  


plot.CalDates <- function(x, ind=1, label=NA, calendar="BP", type="standard", xlab=NA, ylab=NA, axis4=TRUE, HPD=FALSE, credMass=0.95, customCalCurve=NA,add=FALSE,col='grey50',col2='grey82',cex.axis=0.75,cex.xylab=0.75,cex.label=0.75,...){

    types <- c("standard", "simple", "auc")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    caltimeRange =c(55000,0)
    if (any(x$metadata$CalCurve %in% c("intcal13","shcal13","marine13","intcal13nhpine16","shcal13shkauri16")))
    {
      caltimeRange =c(50000,0)
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

    if (!add)
    {
	    par(cex.lab=0.75)
	    plot(xvals,yvals, type="n", xlab=xlabel, ylab="", ylim=yrng, xlim=xlim, xaxt='n', yaxt='n', cex.lab=cex.xylab,cex.axis=cex.axis,...)
    }

    xticksLab <- xticks
    if (calendar=="BCAD")
    {
      if (any(xticksLab==0)){xticksLab[which(xticksLab==0)]=1}
      xticks[which(xticks>1)]=xticks[which(xticks>1)]-1
    }
    if(!add) {axis(1, at=xticks, labels=abs(xticksLab), las=2, cex.axis=cex.axis)}
    
    if (!add&axis4){ axis(4, cex.axis=cex.axis) }
    if (!HPD){
    polygon(xvals,yvals, col=col, border=col)
    } else {
    polygon(xvals,yvals, col=col2, border=col2)
    hdres <- hpdi(x,credMass=credMass)[[ind]]
    if(calendar=="BCAD"){hdres=1950-hdres}
	for (i in 1:nrow(hdres))
	{
	 index <- which(xvals%in%hdres[i,1]:hdres[i,2])
         polygon(c(xvals[index],xvals[index[length(index)]],xvals[index[1]]),c(yvals[index],0,0), col=col, border=col)
	}
    }

    if (type=="standard" | type=="auc"){
        if (type=="auc"){
            lines(xvals, yvals/sum(yvals), col="black", lty="dotted")
        }
    if(!add)
    {
	    par(new=TRUE)
	    cradf1 <- data.frame(CRA=caltimeRange[1]:0,Prob=dnorm(caltimeRange[1]:0, mean=cra, sd=error))
	    cradf1 <- cradf1[cradf1$Prob>0.0001,]
	    ylim <- c(cra-(12*error),cra+(8*error))    
	    cradf1$RX <- reScale(cradf1$Prob, to=c(xlim[1],(xlim[1]+diff(xlim)*0.33)))
	    yticks <- ylim[1]:ylim[2]
	    yticks <- yticks[yticks %% 200 == 0]
	    plot(cradf1$RX,cradf1$CRA,type="l", axes=FALSE, xlab=NA, ylab=NA, xlim=xlim, ylim=ylim, col=rgb(144,238,144,120,maxColorValue=255))
	    polygon(c(cradf1$RX,rev(cradf1$RX)),c(cradf1$CRA,rep(xlim[1],length(cradf1$CRA))), col=rgb(144,238,144,80,maxColorValue=255), border=NA)
	    axis(side=2, at=yticks, labels=abs(yticks),las=2, cex.axis=cex.axis)
	    if (is.na(ylab)){
		    title(line=3, ylab="Radiocarbon Age", cex.lab=cex.xylab)
	    } else {
		    title(line=3, ylab=ylab, cex=cex.xylab)
	    }
    }
        calcurvemetadata <- x$metadata$CalCurve[ind]
        calcurvecheck <- TRUE
        if (calcurvemetadata == "custom" & !any(class(customCalCurve) %in% c("data.frame","matrix","array"))){
            calcurvecheck <- FALSE
        }
        if (calcurvecheck&!add){
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
        legend("topright", label, bty="n", cex=cex.label)
    }
}

#' @title Plot multiple dates
#' 
#' @description Plot multiple radiocarbon dates.
#' @param x A CalDates class object with length > 1.
#' @param type Whether the calibrated dates are displayed as distributions (\code{'d'}) or as horizontal bars (\code{'b'}) spanning the HPD interval. Default is \code{'d'}.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param HPD Logical value indicating whether intervals of higher posterior density should be displayed. Default is FALSE.
#' @param credMass A numerical value indicating the size of the higher posterior density interval. Default is 0.95.
#' @param decreasing Whether dates should be plotted with decreasing order of median calibrated date (i.e. old to new; TRUE) or increasing order (i.e. new to old; FALSE). If set to NULL the dates plotted in the supplied order. Default is NULL 
#' @param label Whether the ID of each date should be displayed. Default is TRUE.
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param xlab (optional) Label for the x axis. If unspecified the default setting will be applied ("Year BP" or "Year BC/AD").
#' @param ylab (optional) Label for the y axis.
#' @param col.fill A vector of primary fill color for the calibrated date distribution. Default is 'grey50'.
#' @param col.fill2 A vector of secondary secondary colour fill color for the calibrated date distribution, used for regions outside the higher posterior interval. Ignored when \code{HPD=FALSE}. Default is 'grey82'.
#' @param col.line A vector of line color (ignored when \code{type} is set to \code{'d'}. Default is 1.
#' @param lwd Line width (ignored when \code{type} is set to \code{'d'}). Default is 1.
#' @param cex.lab The magnification to be used for x and y  labels relative to the current setting of cex. Default is adjusted to 1.
#' @param cex.id The magnification to be used the date labels relative to the current setting of cex. Default is adjusted to 1.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex. Default is adjusted to 1.
#' @param ydisp Whether the y axis should be displayed. Ignored when \code{type} is set to \code{'b'}. Default is FALSE
#' @param gapFactor Defines spacing between calibrated dates (when \code{type} is set to \code{'d'}) or the distance between the lines and the labels (when \code{type} is set to \code{'b'}) as proportion of individual y-axis ranges. Default is 0.2.
#' @param rescale Whether probability distributions should be rescaled (applicable only when \code{type} is set to \code{'d'}, default is FALSE).
#' @seealso \code{\link{calibrate}}
#'
#' @examples
#' data("emedyd")
#' tellaswad = subset(emedyd,SiteName=='Tell Aswad')
#' x = calibrate(tellaswad$CRA,tellaswad$Error,ids=tellaswad$LabID)
#' multiplot(x,HPD=TRUE,decreasing=TRUE,label=FALSE,gapFactor = 0.1)
#' multiplot(x,type='b',calendar='BCAD',cex.id = 0.5,lwd=2,gapFactor = 0.5)
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export  


multiplot<- function(x,type='d',calendar='BP',HPD=FALSE,credMass=0.95,decreasing=NULL,label=TRUE,xlim=NULL,xlab=NA,ylab=NA,col.fill='grey50',col.fill2='grey82',col.line='black',lwd=1,cex.id=1,cex.lab=1,cex.axis=1,ydisp=FALSE,gapFactor=0.2,rescale=FALSE)
{

	if(length(lwd)==1){lwd=rep(lwd,length(x))}
	if(length(col.line)==1){col.line=rep(col.line,length(x))}
	if(length(col.fill)==1){col.fill=rep(col.fill,length(x))}
	if(length(col.fill2)==1){col.fill2=rep(col.fill2,length(x))}


	if (!is.null(decreasing))
	{
	  ord = order(medCal(x),decreasing=decreasing)
		x = x[ord]
		col.line = col.line[ord]
		col.fill = col.fill[ord]
		col.fill2 = col.fill2[ord]
	}

	medDates = medCal(x)


	calendars <- c("BP","BCAD")
	if (!calendar %in% calendars){
		stop("The calendar you have chosen is not currently an option.")
	}        


	# Estimate general xlim
	if (is.null(xlim))
	{
		if(anyNA(x$grids))
		{
			tmp = apply(x$calmatrix,1,sum)
			st = as.numeric(names(tmp[which(tmp>0)[1]-1]))
			en = as.numeric(names(tmp[which(tmp>0)[length(which(tmp>0))]+1]))
		} else {
			st = max(unlist(lapply(x$grids,function(x){max(x$calBP)})))
			en = min(unlist(lapply(x$grids,function(x){min(x$calBP)})))
		}
		edge = 0.1*abs(st-en)
		xlim = c(st+edge,en-edge)
	}

	yearsBP = xlim[1]:xlim[2]



	if (calendar=="BP"){
		plotyears <- yearsBP
		xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
		if (is.na(xlab)){ xlabel <- "Years cal BP" } else { xlabel <- xlab } 
	} else if (calendar=="BCAD"){
		plotyears <- BPtoBCAD(yearsBP)
		xlim <- BPtoBCAD(xlim)
		xvals <- c(plotyears[1],plotyears,plotyears[length(plotyears)], plotyears[1])
		if (is.na(xlab)){ 
			xlabel <- "Years BC/AD"
			if (all(range(plotyears)<0)) {xlabel <- "Years BC"}
			if (all(range(plotyears)>0)) {xlabel <- "Years AD"}
		} else { xlabel <- xlab }       
	} else {
		stop("Unknown calendar type")
	}



	# Plot 
	if (type=='b')
	{
  		bse = hpdi(x,credMass=credMass)
		plot(0,0,xlim=xlim,ylim=c(0,length(bse)+1),axes=F,xlab=xlabel,ylab="",type='n')
		for (i in 1:length(bse))
		{
			tmp = matrix(bse[[i]][,-3],ncol=2)
			if (calendar=='BP')
			{
				apply(tmp,1,function(x,y,lwd,col){lines(c(x),c(y,y),lwd=lwd,col=col)},y=i,lwd=lwd[i],col=col.line[i])
				if(label){text(x=medDates[i],y=i+gapFactor,label=x$metadata$DateID[i],cex=cex.id)}
			}
			if (calendar=='BCAD')
			{
				apply(tmp,1,function(x,y,lwd,col){lines(BPtoBCAD(c(x)),c(y,y),lwd=lwd,col=col)},y=i,lwd=lwd[i],col=col.line[i])
				if(label){text(x=BPtoBCAD(medDates[i]),y=i+gapFactor,label=x$metadata$DateID[i],cex=cex.id)}
			}
		}
	}

	if (type=='d')
	{
		if(anyNA(x$grids))
		{
		  if (rescale)
		  {
		    x$calmatrix=apply(x$calmatrix,2,reScale)
		  }
			tmp = apply(x$calmatrix,1,max)
			ylim = as.numeric(c(0,tmp[which.max(tmp)]))
		} else {
		  if (rescale) 
		  {
		    for (i in 1:length(x))
		    {
		      x$grids[[i]]$PrDens = reScale(x$grids[[i]]$PrDens)
		    }
		  }
			ylim = c(0,max(unlist(lapply(x$grids,function(x){max(x$PrDens)}))))
		}

		# max ylim: combine ylim giving as a space 1/7 of the original distance
		gap = abs(diff(ylim))*gapFactor
		# generate ylim sequences:
		YLIMs = vector("list",length=length(x))
		YLIMs[[1]] = ylim
		for (i in 2:length(x))
		{
			YLIMs[[i]]=c(YLIMs[[i-1]][2]+gap,YLIMs[[i-1]][2]+gap+abs(diff(ylim)))
		}


		plot(0, 0, xlim=xlim, ylim=c(min(unlist(YLIMs)),max(unlist(YLIMs))+gap), type="n", ylab=ylab, xlab=xlabel, axes=F,cex.lab=cex.lab)



		for (i in 1:length(x))
		{
			tmpYlim= YLIMs[[i]]
			if (ydisp)
			{axis(2,at=c(tmpYlim[1],median(tmpYlim),max(tmpYlim)),labels=round(c(min(ylim),median(ylim),max(ylim)),2),las=2,cex.axis=cex.axis)}

			if(anyNA(x$grid))
			{
 				years = as.numeric(rownames(x$calmatrix))
				PrDens = x$calmatrix[,i]
				ii = which(PrDens>0)[1]-1
				jj = max(which(PrDens>0))+1
				years = years[ii:jj]
				PrDens = PrDens[ii:jj]
			} else {
				years=x$grid[[i]]$calBP
				PrDens = x$grid[[i]]$PrDens
			}

			if (calendar=='BCAD'){years=BPtoBCAD(years)}

			xvals = c(years,rev(years))
			yvals = c(PrDens+tmpYlim[1],rep(0,length(years))+tmpYlim[1])



			if (!HPD){
				polygon(xvals,yvals, col=col.fill[i], border=col.fill[i])
			} else {
				polygon(xvals,yvals, col=col.fill2[i], border=col.fill2[i])
				hdres <- hpdi(x,credMass=credMass)[[i]]
				for (j in 1:nrow(hdres))
				{
					if (calendar=='BCAD')
					{
						index <- which(xvals%in%BPtoBCAD(hdres[j,1]:hdres[j,2]))
					} else {
						index <- which(xvals%in%hdres[j,1]:hdres[j,2])
					}
					polygon(c(xvals[index],xvals[index[length(index)]],xvals[index[1]]),c(yvals[index],min(tmpYlim),min(tmpYlim)), col=col.fill[i], border=col.fill[i])
				}
			}


			if (label)
			{
				xx = medDates[i]
				if (calendar=='BCAD'){xx = BPtoBCAD(xx)}
				ylabel=ifelse(rescale,max(yvals)-0.5,max(yvals)+gap/2) 
				text(x=xx,y=ylabel,labels=x$metadata$DateID[i],cex=cex.id)
			}

		}

	}


	# draw x-axis
	if (calendar=="BP"){
		rr <- range(pretty(plotyears))    
		axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01,cex.axis=cex.axis)
		axis(side=1,at=pretty(plotyears),labels=abs(pretty(plotyears)),cex.axis=cex.axis)
	} else if (calendar=="BCAD"){
		yy <-  plotyears
		rr <- range(pretty(yy))    
		prettyTicks <- seq(rr[1],rr[2],100)
		prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
		axis(side=1,at=prettyTicks, labels=NA,tck = -.01,cex.axis=cex.axis)
		py <- pretty(yy)
		pyShown <- py
		if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
		py[which(py>1)] <-  py[which(py>1)]-1
		axis(side=1,at=py,labels=abs(pyShown),cex.axis=cex.axis)
	}








}














#' @title Plot result of Monte-Carlo simulation of observed versus modelled SPDs
#'
#' @description The function visualises the observed summed probability distribution of radiocarbon dates along with a simulation envelope for the null model and regions of positive and negative deviation.
#'
#' @param x A \code{SpdModelTest} class object generated using the \code{\link{modelTest}} function.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param type Whether to display SPDs ('spd') or rates of change ('roc'). Default is 'spd'.
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param ylim the y limits of the plot.
#' @param col.obs Line colour for the observed SPD.
#' @param lwd.obs Line width for the observed SPD.
#' @param xaxs The style of x-axis interval calculation (see \code{\link{par}})
#' @param yaxs The style of y-axis interval calculation (see \code{\link{par}})
#' @param bbty Display options; one between \code{'b'},\code{'n'},and \code{'f'}. See details below.
#' @param bbtyRet Whether details on the intervals of positive and negative deviations are returned. Default is TRUE.
#' @param drawaxes A logical value determining whether the axes should be displayed or not. Default is TRUE.
#' @param ... Additional arguments affecting the plot

#' @details The argument \code{bbty} controls the display options of the Monte-Carlo Test. Default settings (\code{bbty='f'}) displays the observed SPD (solid black line), the simulation envelope of the fitted model (shaded grey polygon) and regions of significance positive (red semi-transparent rectangle) and negative (blue semi-transparent rectangle) deviation. The option \code{bbty='b'} removes the regions of positive/negative deviations, whilst the option \code{bbty='n'} displays the simulation envelope on existing plot. 

#' @seealso \code{\link{modelTest}}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @method plot SpdModelTest
#' @export 

plot.SpdModelTest <- function(x, calendar="BP", type='spd', ylim=NA, xlim=NA, col.obs="black", lwd.obs=0.5, xaxs="i", yaxs="i", bbty="f", bbtyRet=TRUE, drawaxes=TRUE, ...){

	obs <- x$result[,1:2]
	envelope <- x$result[,3:4]
	
	if (!type%in%c('spd','roc'))
	{
	 stop("The argument 'type' should be either 'spd' or 'roc'.")
	}

	if (type=='roc')
	{
		obs <- x$result.roc[,1:2]
		envelope <- x$result.roc[,3:4]
		colnames(obs)<-c("calBP","PrDens")
		colnames(envelope)<-c("lo","hi")
	}

	if (calendar=="BP"){
		obs$Years <- obs$calBP
		xlabel <- "Years cal BP"
		if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
	} else if (calendar=="BCAD"){
		xlabel <- 'Years BC/AD'    
		obs$Years <- BPtoBCAD(obs$calBP)
		if (all(range(obs$Years)<0)){xlabel <- "Years BC"}
		if (all(range(obs$Years)>0)){xlabel <- "Years AD"}
		if (any(is.na(xlim))){xlim <- c(min(obs$Years),max(obs$Years)) }
	} else {
		stop("Unknown calendar type")
	}    


	if (any(is.na(ylim)))
		{
		       	ylim <- c(0, max(envelope[,"hi"], obs$PrDens,na.rm=TRUE)*1.1) 
			if (type=='roc') {ylim[1] <- min(c(envelope[,"lo"],obs$PrDens),na.rm=TRUE)}
		}
	
	booms <- which(obs$PrDens>envelope[,2])
	busts <- which(obs$PrDens<envelope[,1])
	baseline <- rep(NA,nrow(obs))


  	ylab = 'Summed Probability'
	if (type=='roc') 
	{
		ylab = 'Rate of Change'
		colpts = rep(col.obs,length(obs$PrDens))
		colpts[booms] = 'red'
		colpts[busts] = 'blue'
	}


	if (drawaxes & bbty != "n"){
			plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab=xlabel, ylab=ylab, type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE,...)
		if(type=='roc')
		{
			abline(h=0,lty=2,lwd=1)
		}

	} else if (bbty != "n")
	{
		plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab="", ylab="", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
		if (type=='roc')
		{
		abline(h=0,lty=2,lwd=1)
		}
	}




	if (drawaxes)
	{
		box()
		axis(side=2)
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
	if ((bbty %in% c("n","b")) & bbtyRet){ return(bbp) }
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
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex. Default is 1.
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
#' @method plot CalSPD
#' @export 

plot.CalSPD <- function(x, runm=NA, calendar="BP", type="standard", xlim=NA, ylim=NA, ylab="Summed Probability", spdnormalised=FALSE, rescale=FALSE, fill.p="grey75", border.p=NA, xaxt='s', yaxt='s', add=FALSE, cex.axis=1, ...){

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
			plot(plotyears, spdvals, xlim=xlim, ylim=ylim, type="l", col="white", ylab=ylabel, xlab=xlabel, xaxt="n", yaxt=yaxt, ...)
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
		axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01,cex.axis=cex.axis)
		axis(side=1,at=pretty(plotyears),labels=abs(pretty(plotyears)),cex.axis=cex.axis)
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
#' @method plot CalGrid
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



#' @title Plot a stack of SPDs
#'
#' @description Visualises multiple SPDs grouped as a \code{stackCalSPD} object.
#' 
#' @param x A \code{stackCalSPD} class object. Result of \code{\link{stackspd}} function.
#' @param type How to display the SPDs.Current options are \code{'stacked'},\code{'lines'}, '\code{'proportion'}. and \code{'multipanel'}. Default is \code{'stacked'}. 
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param spdnormalised A logical variable indicating whether the total probability mass of the SPDs are normalised to sum to unity. Default is FALSE. 
#' @param rescale  A logical variable indicating whether the SPDs should be rescaled to range 0 to 1. Default is FALSE.
#' @param runm A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied. Default is NA  
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param ylim the y limits of the plot.
#' @param xaxt Whether the x-axis tick marks should be displayed (\code{xaxt='s'}, default) or not (\code{xaxt='n'}).
#' @param yaxt Whether the y-axis tick marks should be displayed (\code{xaxt='s'}, default) or not (\code{xaxt='n'}).
#' @param gapFactor Defines spacing between SPDs as proportion of the y-axis range for multipanel plots. Default is 0.2.
#' @param col.fill Vector of fill color for the observed SPDs. The default color scheme is based on the Dark2 pallette of RColorBrewer package.
#' @param col.line Line colour for the observed SPDs.The default color scheme is based on the Dark2 palette of RColorBrewer package.
#' @param lwd.obs Line width for the observed SPDs. Default is 1.
#' @param lty.obs Line type for the observed SPDs. Default is 1.
#' @param cex.lab The magnification to be used for x and y  labels relative to the current setting of cex. Default is adjusted to 1.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex. Default is adjusted to 1.
#' @param legend Whether legend needs to be displayed. Item names will be retrieved from the values supplied in the argument \code{group} in \code{\link{stackspd}}. Default is TRUE.
#' @param legend.arg list of additional arguments to pass to \code{\link{legend}}; names of the list are used as argument names. Only used if \code{legend} is set to TRUE.
#' @param ylab a title for the y axis
#' @param ymargin multiplier for the maximum value on ylim range. Default is 1.1.
#' @param ... Additional arguments affecting the plot. 

#' @references 
#' Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. \url{https://CRAN.R-project.org/package=RColorBrewer}.
#' @examples
#' \dontrun{
#'data(emedyd)
#'emedyd = subset(emedyd,Region==1)
#'x = calibrate(x=emedyd$CRA, errors=emedyd$Error,normalised=FALSE)
#'bins = binPrep(sites=emedyd$SiteName, ages=emedyd$CRA,h=50)
#'res = stackspd(x=x,timeRange=c(16000,8000),bins=bins,group=emedyd$Region)
#'plot(res,type='stacked')
#'plot(res,type='lines')
#'plot(res,type='proportion')
#'plot(res,type='multipanel')
#'}
#' @method plot stackCalSPD
#' @export

plot.stackCalSPD <- function(x, type='stacked', calendar='BP', spdnormalised=FALSE,rescale=FALSE, runm=NA, xlim=NA, ylim=NA, xaxt='s', yaxt='s',gapFactor = 0.2, col.fill=NA, col.line=NA, lwd.obs=1, lty.obs=1, cex.lab=1, cex.axis=cex.lab,legend=TRUE,legend.arg=NULL, ylab=NA, ymargin=1.1, ...)
{
  # Based on Dark2 of RcolorBrewer
  colpal=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")
  

	#issue error messages
	if (!'stackCalSPD'%in%class(x))
	{
	 stop("The argument x should be a 'stackCalSPD' class object")
	}

	if (!type%in%c('multipanel','lines','stacked','proportion'))
	{
	 stop("The argument 'type' should be one between 'multipanel', 'lines', 'stacked'  or 'proportion'.")
	}
  
  # Calendar Setting
  if (!calendar %in% c("BP","BCAD")){
    stop("The calendar you have chosen is not currently an option.")
  }

  if (calendar=="BP"){
    plotyears <- x$metadata$timeRange[1]:x$metadata$timeRange[2]
    xlabel <- "Years cal BP"
    if (any(is.na(xlim))){ xlim <- c(max(plotyears),min(plotyears)) }
  } else if (calendar=="BCAD"){
    plotyears <- BPtoBCAD(x$metadata$timeRange[1]:x$metadata$timeRange[2])
    xlabel <- "Years BC/AD"
    if (all(range(plotyears)<0)){xlabel <- "Years BC"}
    if (all(range(plotyears)>0)){xlabel <- "Years AD"}
    if (any(is.na(xlim))){ xlim <- c(min(plotyears),max(plotyears)) }
  } else {
    stop("Unknown calendar type")
  }
  if (xaxt=='n'){ xlabel <- "" }
  if (yaxt=='n'){ ylabel <- "" } else { ylabel <- ylab }

  #pre-processing (create a single matrix for plotting)
  nsets = length(x$spds)
  if (nsets==1 & type=='multipanel') {type='stacked'}
  PrDens = matrix(NA,ncol=nsets,nrow=length(plotyears))
	
	for (i in 1:nsets)
	{
	  PrDens[,i]=x$spds[[i]][[2]][[2]]
	  if (!is.na(runm)){ PrDens[,i] <- runMean(PrDens[,i], runm, edge="fill") }
	}
	
	if (spdnormalised) {PrDens=apply(PrDens,2,function(x){x/sum(x)})}
	if (rescale){PrDens = apply(PrDens,2,reScale)}
	
	#Assign Colours and Line Types
	if (any(is.na(col.fill)))
	{
	  if (nsets<=8)
	  {
	    col.fill=colpal[1:nsets]
	  }
	  if (nsets > 8)
	  {
	    col.fill=sample(colors(),size=nsets,replace=TRUE)
	  warning("Color sequence randomised due to a large number of SPDs (>8). Consider selecting an appropriate color sequence manually")
	  }
	}
	
	if (any(is.na(col.line)))
	{
	  if (nsets<=8)
	  {
	    col.line=colpal[1:nsets]
	  }
	  if (nsets > 8)
	  {
	    col.line=sample(colors(),size=nsets,replace=TRUE)
	  warning("Color sequence randomised due to a large number of SPDs (>8). Consider selecting an appropriate color sequence manually")
	  }
	}
	
	if (any(is.na(lwd.obs))){lwd.obs=rep(1,nsets)}
	if (any(is.na(lty.obs))){lty.obs=rep(1,nsets)}

	#if only single values are provided
	if (length(lwd.obs)==1){lwd.obs=rep(lwd.obs,nsets)}
	if (length(lty.obs)==1){lty.obs=rep(lty.obs,nsets)}
	if (length(col.line)==1){col.line=rep(col.line,nsets)}
	if (length(col.fill)==1){col.fill=rep(col.fill,nsets)}
	




	#Plot

	#lines 
	if (type=='lines')
	{
	  if (any(is.na(ylim))){ ylim <- c(0,max(PrDens)*ymargin) }
	  if (is.na(ylab)){ylab='Summed Probability'}
	  plot(0, 0, xlim=xlim, ylim=ylim, type="l", ylab=ylab, xlab=xlabel, xaxt="n", yaxt=yaxt,cex.axis=cex.axis,cex.lab=cex.lab)
	  
	  for (i in 1:nsets)
	  {
	    lines(plotyears,PrDens[,i],col=col.line[i],lty=lty.obs[i],lwd=lwd.obs[i])
	  }
	  
	  if (legend)
	  {
	    if (is.null(legend.arg))
	    {
	      legend("topleft",legend=names(x$spds),col=col.line,lty=lty.obs,lwd=lwd.obs)
	    } else {
	      args.legend1 <- list("topleft", legend = names(x$spds),col=col.line,lty=lty.obs,lwd=lwd.obs)
	      args.legend1[names(legend.arg)] <- legend.arg
	      do.call("legend", args=args.legend1)
	    }
	  }
	}

	#stacked
	if (type=='stacked')
	{
		if (nsets>1)
		{	  
			PrDens = t(apply(PrDens,1,cumsum))
		}
		PrDens = cbind(0,PrDens)
	  
	  if (any(is.na(ylim))){ ylim <- c(0,max(PrDens)*ymargin) }
	  if (is.na(ylab)){ylab='Summed Probability'}
	  
	  plot(0, 0, xlim=xlim, ylim=ylim, type="l", ylab=ylab, xlab=xlabel, xaxt="n", yaxt=yaxt,cex.axis=cex.axis,cex.lab=cex.lab)
	  
	  for (i in 2:(nsets+1))
	  {
	    polygon(c(plotyears,rev(plotyears)),c(PrDens[,i],rev(PrDens[,i-1])),col=col.fill[i-1],lwd=lwd.obs[i-1],border=col.line[i-1],lty=lty.obs[i-1])
	  }
	  if (legend)
	  {
	    if (is.null(legend.arg))
	    {
	      legend("topleft",legend=names(x$spds),fill=col.fill)
	    } else {
	      
	      args.legend1 <- list("topleft", legend=names(x$spds),fill=col.fill)
	      args.legend1[names(legend.arg)] <- legend.arg
	      do.call("legend", args=args.legend1)
	    }
	  }
	}

	#proportion
	if (type=='proportion')
	{

		PrDens=prop.table(PrDens,1)
		if (nsets>1)
		{
			PrDens = t(apply(PrDens,1,cumsum))
		}
		PrDens = cbind(0,PrDens)
		if (is.na(ylab)){ylab='Relative Proportion'}
		plot(0, 0, xlim=xlim, ylim=c(0,1), type="l", ylab=ylab, xlab=xlabel, xaxt="n", yaxt=yaxt,cex.axis=cex.axis,cex.lab=cex.lab)

		for (i in 2:(nsets+1))
		{
			polygon(c(plotyears,rev(plotyears)),c(PrDens[,i],rev(PrDens[,i-1])),col=col.fill[i-1],lwd=lwd.obs[i-1],lty=lty.obs[i-1],border=col.fill[i-1])
		}
		if (legend)
		{
			if (is.null(legend.arg))
			{
				legend("topleft",legend=names(x$spds),fill=col.fill)
			} else {
				args.legend1 <- list("topleft", legend=names(x$spds),fill=col.fill)
				args.legend1[names(legend.arg)] <- legend.arg
				do.call("legend", args.legend1)
			}
		}
	}

	#multipanel
  if (type=='multipanel')
  {
    # individual ylim
    if (any(is.na(ylim))){ ylim <- c(0,max(PrDens)*ymargin) }
    if (is.na(ylab)){ylab='Summed Probability'}
    
    # max ylim: combine ylim giving as a space 1/7 of the original distance
    gap = abs(diff(ylim))*gapFactor
    # generate ylim sequences:
    YLIMs = vector("list",length=nsets)
    YLIMs[[1]] = ylim
    for (i in 2:nsets)
    {
      YLIMs[[i]]=c(YLIMs[[i-1]][2]+gap,YLIMs[[i-1]][2]+gap+abs(diff(ylim)))
    }
    
    plot(0, 0, xlim=xlim, ylim=c(min(unlist(YLIMs)),max(unlist(YLIMs))+gap), type="l", ylab=ylab, xlab=xlabel, axes=F,cex.lab=cex.lab)
    
    for (i in 1:nsets)
    {
      tmpYlim= YLIMs[[i]]
      axis(2,at=c(tmpYlim[1],median(tmpYlim),max(tmpYlim)),labels=round(c(min(ylim),median(ylim),max(ylim)),2),las=2,cex.axis=cex.axis)
      
      polygon(c(plotyears,rev(plotyears)),c(PrDens[,i]+tmpYlim[1],rep(0,length(plotyears))+tmpYlim[1]),col=col.fill[i],lwd=lwd.obs[i],border=col.line[i],lty=lty.obs[i])
    }
    
    if (legend)
    {
      if (is.null(legend.arg))
      {
        legend("topleft",legend=names(x$spds),fill=col.fill)
      } else {
        args.legend1 <- list("topleft", legend=names(x$spds),fill=col.fill)
        args.legend1[names(legend.arg)] <- legend.arg
        do.call("legend", args.legend1)
      }
    }
    
}


	# draw x-axis
	if (calendar=="BP" & xaxt!="n"){
	  rr <- range(pretty(plotyears))    
	  axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01,cex.axis=cex.axis)
	  axis(side=1,at=pretty(plotyears),labels=abs(pretty(plotyears)),cex.axis=cex.axis)
	} else if (calendar=="BCAD" & xaxt!="n"){
	  yy <-  plotyears
	  rr <- range(pretty(yy))    
	  prettyTicks <- seq(rr[1],rr[2],+100)
	  prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
	  axis(side=1,at=prettyTicks, labels=NA,tck = -.01,cex.axis=cex.axis)
	  py <- pretty(yy)
	  pyShown <- py
	  if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
	  py[which(py>1)] <-  py[which(py>1)]-1
	  axis(side=1,at=py,labels=abs(pyShown),cex.axis=cex.axis)
	}
}

#' @title Plot result of mark permutation test of SPDs
#'
#' @description Visualises the observed SPD along with the simulation envelope generated from \code{\link{permTest}}, with regions of positive and negative deviations highlighted in red and blue.
#'
#' @param x A \code{SpdPermTest} class object. Result of random mark permutation test (see \code{\link{permTest}})
#' @param focalm Value specifying the name of the focal mark (group) to be plotted. 
#' @param type Whether to display SPDs ('spd') or rates of change ('roc'). Default is 'spd'.
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
#' @method plot SpdPermTest
#' @export 

plot.SpdPermTest <- function(x, focalm=1, type='spd', calendar="BP", xlim=NA, ylim=NA, col.obs="black", col.env=rgb(0,0,0,0.2), lwd.obs=0.5, xaxs="i", yaxs="i", bbty="f", drawaxes=TRUE, ...){


	if (!type%in%c('spd','roc'))
	{
	 stop("The argument 'type' should be either 'spd' or 'roc'.")
	}

	obs <- x$observed[[focalm]]
	envelope <- x$envelope[[focalm]]
	ylab <- 'Summed Probability'

	if (type=='roc')
	{
		ylab <- 'Rate of Change'
		obs <- x$observed.roc[[focalm]]
		envelope <- x$envelope.roc[[focalm]]
		colnames(obs)<-c("calBP","PrDens")
		colnames(envelope)<-c("lo","hi")
	}

	if (calendar=="BP"){
		obs$Years <- obs$calBP
		xlabel <- "Years cal BP"
		if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
	} else if (calendar=="BCAD"){
		obs$Years <- BPtoBCAD(obs$calBP)
		xlabel <- "Years BC/AD"
		if (all(range(obs$Years)<0)){xlabel <- "Years BC"}
		if (all(range(obs$Years)>0)){xlabel <- "Years AD"}
		if (any(is.na(xlim))){ xlim <- c(min(obs$Years),max(obs$Years)) }
	} else {
		stop("Unknown calendar type")
	}


	if (any(is.na(ylim)))
		{
		       	ylim <- c(0, max(envelope[,2], obs$PrDens,na.rm=TRUE)) 
			if (type=='roc') {ylim[1] <- min(c(envelope[,1],obs$PrDens),na.rm=TRUE)}
		}

	booms <- which(obs$PrDens>envelope[,2])
	busts <- which(obs$PrDens<envelope[,1])
	baseline <- rep(NA,nrow(obs))
	if (drawaxes & bbty != "n"){
		plot(obs$Years, obs$PrDens, xlim=xlim, ylim=ylim, xlab=xlabel, ylab=ylab, type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
		axis(side=2,padj=1)
		if (type=='roc'){abline(h=0,lty=2)}
	} else if (bbty != "n"){
		plot(obs$Years,obs$PrDens, xlim=xlim, ylim=ylim, xlab="", ylab="", type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ...)
		if (type=='roc'){abline(h=0,lty=2)}
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
#' @param binning  Either \code{'CRA'} or \code{'calibrated'}. Indicate whether the binning should be carried using the 14C age or using the median calibrated date. Default is \code{'CRA'}.
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
#' @seealso \code{\link{sptest}}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @method plot spatialTest
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

	if (!option%in%c("raw","test"))
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



#' @title Plot \code{spdRC} class objects 
#'
#' @description Plot rates of change between time-blocks
#' @param x \code{spdRC} class object containing geometric growth rates.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param col.obs Line colour for the observed SPD
#' @param lwd.obs Line width for the observed SPD
#' @param xlim Limits for the x axis.
#' @param xaxs The style of x-axis interval calculation (see \code{\link{par}})
#' @param yaxs The style of y-axis interval calculation (see \code{\link{par}})
#' @param ... Additional arguments affecting the plot. 
#'
#' @seealso \code{\link{spd2rc}}
#'
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @method plot spdRC
#' @export  

plot.spdRC<- function(x,calendar="BP",col.obs="black", lwd.obs=0.5, xaxs="i", yaxs="i",xlim=NA,...)
{
	if (x$type=='blocks')
	{
		breaks=x$breaks
		obs=x$sumblock
		res=x$roca
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
		mtext(side=4,"rate of change",col="darkgreen",line=2)
		abline(h=0,lty=2,col="blue")
	}

	if (x$type=='backsight')
	{
		obs=data.frame(calBP=x$timeSequence, roc=x$roca)

		if (calendar=="BP"){
			obs$Years <- obs$calBP
			xlabel <- "Years cal BP"
			if (any(is.na(xlim))){ xlim <- c(max(obs$Years),min(obs$Years)) }
		} else if (calendar=="BCAD"){
			xlabel <- 'Years BC/AD'    
			obs$Years <- BPtoBCAD(obs$calBP)
			if (all(range(obs$Years)<0)){xlabel <- "Years BC"}
			if (all(range(obs$Years)>0)){xlabel <- "Years AD"}
			if (any(is.na(xlim))){xlim <- c(min(obs$Years),max(obs$Years)) }
		} else {
			stop("Unknown calendar type")
		}


			plot(obs$Years, obs$roc, xlim=xlim, xlab=xlabel, ylab='Rate of Change', type="l", col=col.obs, lwd=lwd.obs, xaxs=xaxs, yaxs=yaxs, axes=FALSE,...)

			abline(h=0,lty=2,lwd=1)
			box()
			axis(side=2)

			
			if (calendar=="BP"){
				rr <- range(pretty(obs[,"Years"]))    
				axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
				axis(side=1,at=pretty(obs[,"Years"]),labels=abs(pretty(obs[,"Years"])))
			} else if (calendar=="BCAD"){
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
	}
}


#' @title Plots a Composite Kernel Density Estimate of sampled radiocarbon dates.  
#' @param x A \code{compositeKDE} class object generated using the \code{\link{ckde}} function.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param type Either \code{envelope} or \code{multiline}. Default is \code{envelope}.
#' @param ylim the y limits of the plot.
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param fill.col Envelope color when \code{type='envelope'}. Default is 'lightgrey'.
#' @param interval Quantile interval for the envelope. Default is 0.95.
#' @param line.col Line color when \code{type='envelope'}. Default is 'black.
#' @param line.type Line type when \code{type='envelope'}. Default is 2.
#' @param multiline.alpha Alpha level for line transparency when \code{type='multiline'}. Default is 10/\code{nsim}, where \code{nsim} is the number of simulations. If \code{nsim} is smaller than 10, \code{multiline.alpha} will be set to 1.
#' @param multiline.col Line color when \code{type='multiline'}. Default is 'black'.
#' @param ... Additional arguments affecting the plot
#' @details Visualise a \code{compositeKDE} class object. If \code{type} is set \code{'envelope'} an envelope of the percentile interval defined by the parameter \code{interval} is shown along with the mean KDE. If \code{type} is set \code{'multiline'} all KDEs are shown. 
#' @seealso \code{\link{ckde}}; 
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @method plot compositeKDE
#' @export 
plot.compositeKDE <- function(x, calendar="BP", type='envelope', ylim=NA, xlim=NA, fill.col='lightgrey',interval=0.95,line.col='black',line.type=2, multiline.alpha=NA, multiline.col='black',...){

    types <- c("envelope","multiline")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    if (any(is.na(ylim))){ ylim <- c(0,max(x$res.matrix,na.rm=TRUE)*1.1)}
   
    plotyears = x$timeRange[1]:x$timeRange[2]
    
    if (calendar=="BP"){
        xlabel <- "Years cal BP"
        if (any(is.na(xlim))){ xlim <- c(max(plotyears),min(plotyears)) }
    } else if (calendar=="BCAD"){
        plotyears <- BPtoBCAD(plotyears)
        xlabel <- "Years BC/AD"
	if (all(range(plotyears)<0)){xlabel <- "Years BC"}
	if (all(range(plotyears)>0)){xlabel <- "Years AD"}
        if (any(is.na(xlim))){ xlim <- c(min(plotyears),max(plotyears)) }
    } else {
        stop("Unknown calendar type")
    }

    if (type=='multiline')
    {
	    plot(x=plotyears,y=x$avg.res,type='n',xlab=xlabel,ylab="Summed Probability",xlim=xlim,ylim=ylim,axes=FALSE,...)
	    if (is.na(multiline.alpha)){multiline.alpha=10/ncol(x$res.matrix)}
	    if (multiline.alpha>1){multiline.alpha=1}
	    mc = c(as.numeric(col2rgb(multiline.col)/255),multiline.alpha)
	    apply(x$res.matrix,2,lines,x=plotyears,col=rgb(mc[1],mc[2],mc[3],mc[4]))
    }

    if (type=='envelope')
    {
	avg=apply(x$res.matrix,1,mean,na.rm=TRUE)
	lo=apply(x$res.matrix,1,quantile,prob=(1-interval)/2,na.rm=TRUE)
   	hi=apply(x$res.matrix,1,quantile,prob=interval+(1-interval)/2,na.rm=TRUE)
	index = which(!is.na(hi))
	plot(x=plotyears,y=avg,type='n',xlab=xlabel,ylab="Summed Probability",xlim=xlim,ylim=ylim,axes=FALSE,...)
        polygon(x=c(plotyears[index],rev(plotyears[index])),y=c(lo[index],rev(hi[index])),border=NA,col=fill.col)
	lines(plotyears,avg,lty=line.type,col=line.col,lwd=2)
    }
    
    if(abs(diff(range(plotyears)))>10000){smallTickBreak=1000}
    if(abs(diff(range(plotyears)))<10000){smallTickBreak=100}
    if(abs(diff(range(plotyears)))<1000){smallTickBreak=10}
    
    if (calendar=="BP"){
	      rr <- range(pretty(plotyears))
	      rrSeq = seq(rr[2],rr[1],-smallTickBreak)
	      rrSeq = rrSeq[which(rrSeq>=min(plotyears)&rrSeq<=max(plotyears))]
        axis(side=1,at=rrSeq,labels=NA,tck = -.01)
        axis(side=1,at=axTicks(1),labels=axTicks(1))
    } else if (calendar=="BCAD"){
	      yy <-  plotyears
        rr <- range(pretty(yy))    
        prettyTicks <- seq(rr[1],rr[2],+smallTickBreak)
      	prettyTicks[which(prettyTicks>=0)] <-  prettyTicks[which(prettyTicks>=0)]-1
      	prettyTicks = prettyTicks[which(prettyTicks>=min(yy)&prettyTicks<=max(yy))]
        axis(side=1,at=prettyTicks, labels=NA,tck = -.01)
        py <- axTicks(1)
      	pyShown <- py
	if (any(pyShown==0)){pyShown[which(pyShown==0)]=1}
	py[which(py>1)] <-  py[which(py>1)]-1
	axis(side=1,at=py,labels=abs(pyShown))
    }
    axis(2)
    box()
}

#' @title Plot map(s) of the spatio-temporal intensity of radiocarbon dates.
#'
#' @description Plotting function for single or multiple maps of the spatio-temporal intensity of radiocarbon dates at given focal years.
#'
#' @param x An object of class stKde.
#' @param focalyears A vector of numeric values for focal years, in calBP, that will be timesteps at which date intensity maps will be plotted.
#' @param type A single character string stipulating which type of plot to create. Current options are "nonfocal", "mask", "focal", "proportion", "change" and "all".
#' @param plotdir Optional output directory for plots. If NULL, then only a single plot is made to the current device. If a valid output direcotry is provided then one or more timeslices maps are saved as png files (e.g. as source images for an animation).
#' @param imnames The format of the output files if output is as png files to a directory. The current two options are "basic" (labelling the images in basic sequence as preferred by animation software such as ffmpeg) or "byyear" (labelling the images by calBP year).
#' @param zlim  Numeric vector of length=2 which controlls the maximum or minimum of the colour ramp.
#' @param box  Logical. Plot a border around the map or not.
#' @param main Single character string specifying a main title. "auto" implies internal default titles are used.
#' @param col.mask The colour used to depict any areas that are being masked out.
#' @param ncolours The maximum number of colours to use in the colour ramp. 
#' @param ramptype What kind of treatment for the colour ramp. Current options are "raw" (do not try to standardise the ramps across timeslices),"std" (standardise each plot, typically by capturing the first 3sd in the main colour ramp and then outliers beyond that with the extreme colours of the ramp),"unl" (do not standardise but plot generalised high/low ramp labels) and "mmx" (scale to the minimum and maximum values across all timeslices).
#' @param withpts Plot with the original sample locations shown (current options are "y" and "n").
#' @param imdim Height and width of the plot to png in cm.
#' @param mars Margins of the plot to png.
#' @param ribbon Whether to plot the colour ramp legend or not.
#' @param ribargs Whether to plot the colour ramp legend or not.
#' @param cex.ribbon The size of the ribbon font.
#' @param cex.main The size of the title font.
#' @param pch.pts The symbols used for the plotted points.
#' @param col.pts The colours used for the plotted points.
#' @param cex.pts The size used for the plotted points.
#' @param tidydevs Logical for whether to clean up any open grpahics devices or not (default is TRUE).
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @param ... Additional arguments affecting the plot.
#'
#' @details This function plots to a screen device if a single focal year is stipulated and no output directory. Or if an output directory is stipulated, it plots one or more focal years to png files, with some basic formatting options and optional cross-year standardisation of the colour ramps (with a view to them being stills-for-video). For even more control of plotting, call this function one year at a time in a loop.
#' 
#' @examples
#' \dontrun{
#' ## Example with a subset of English and Welsh dates from the Euroevol dataset
#' data(ewdates)
#' data(ewowin)
#' dir.create(file.path("im"), showWarnings=FALSE)
#' x <- calibrate(x=ewdates$C14Age, errors=ewdates$C14SD, normalised=FALSE)
#' bins1 <- binPrep(sites=ewdates$SiteID, ages=ewdates$C14Age, h=50)
#' stkde1 <- stkde(x=x, coords=ewdates[,c("Eastings", "Northings")], 
#' win=ewowin, sbw=40000, cellres=2000, focalyears=seq(6500, 5000, -100),
#' tbw=50, bins=bins1, backsight=200, outdir="im")
#' 
#' ## Plot just the proportion surface for just 5900 calBP
#' plot(stkde1, 5900, type="proportion")
#'
#' ## Plot an example of all four basic outputs for just 5900 calBP
#' dev.new(height=2.5, width=8)
#' par(mar=c(0.5, 0.5, 2.5, 2))
#' plot(stkde1, 5900, type="all")
#'
#' ## Plot standardised change surfaces to a sub-directory called 
#' ## /png for all timeslices and save to file.
#' dir.create(file.path("png"), showWarnings=FALSE)
#' plot(stkde1, seq(6500,5000,-100), type="change", ramptype="std", withpts=TRUE, plotdir="png")
#'
#' ## Plot all four summary surfaces in one image, saving them to a sub-directory call 'pngall',
#' ## and with the output of the change map standardised to a common ramp 
#' ## (but leaving the focal and proportion maps unstandardised with simple ramp labelling)
#' dir.create(file.path("pngall"), showWarnings=FALSE)
#' plot(stkde1, seq(6500,5000,-100), type="all", ramptype=c("unl","unl","std"), imdim=cm(c(2.5,8)),
#' withpts=TRUE, plotdir="pngall")
#' }
#' 
#' @import spatstat 
#' @method plot stKde
#' @export
#' 
plot.stKde <- function(x, focalyears=NULL, type="focal", plotdir=NULL, imnames="byyear", zlim=NULL, box=FALSE, main="auto", col.mask="grey75", ncolours=256, ramptype="raw", imdim=c(10,10), mars=c(0.5, 0.5, 2.5, 2), ribbon=TRUE, ribargs=list(), cex.ribbon=0.5, withpts="n", cex.main=0.8, pch.pts=19, col.pts="grey50", cex.pts=0.1, tidydevs=TRUE, verbose=TRUE, ...){
    if (is.null(plotdir) & length(focalyears)>1){
        stop("With more that one focalyear, plots are written to file only, so plotdir must be specified.")
    }
    if (is.null(focalyears)){
        if (length(x$years)>1){
            warning("Plot object has more than one year of data but focalyears argument has not been specified...plotting first year by default.")
        }
        focalyears <- 1
    } else {
        focalyears <- as.character(focalyears)
    }
    maindef <- main
    types <- c("nonfocal", "mask", "focal", "proportion", "change", "all")
    if (!type %in% types){
        stop("The plot type you have chosen is not currently an option.")
    }
    if (x$maskthresh  > 0){
        mymask <- x$nonfocal < x$maskthresh
    }
    if (!is.null(plotdir)){
        if (imnames=="basic"){
            filenames <- paste(plotdir,"/", sprintf("img_%05d",1:length(focalyears)),".png",sep="")
        }  else if (imnames=="byyear"){
            filenames <- paste(plotdir,"/img_",focalyears,".png",sep="")
        } else {
            stop("imnames must by either basic or byyear.")
        }
    }
    if (verbose & length(focalyears)>1){
        print("Plotting to file...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(focalyears), style=3)
    }
    if (!all(ramptype %in% c("raw","std","unl","mmx"))){
        stop("Invalid value(s) for ramptype.")
    }
    if (length(ramptype)==3 & all(class(ramptype)=="character")){
        ramptypes <- ramptype
        names(ramptypes) <- c("focal","proportion","change")
    } else if (length(ramptype)==1 & all(class(ramptype)=="character")){ 
        ramptypes <- ramptype
        names(ramptypes) <- type
        if (type=="all"){
            ramptypes <- rep(ramptype,3)
            names(ramptypes) <- c("focal","proportion","change")
        }
    } else {
        stop("Incorrect number of values for ramptype argument.")
    }
    colfun <- spatstat.options("image.colfun")
    zlimcheck <- zlim
    for (a in 1:length(focalyears)){
        if (verbose & length(focalyears)>1){ setTxtProgressBar(pb, a) }
        load(x$impaths[focalyears[a]])
        if (main=="auto"){
            main.nonfocal <- paste("Overall Intensity\n(all years)",sep="")
            main.mask <- paste("Mask for Low Intensity Areas",sep="")
            main.focal <- paste("Focal Intensity\n(",focalyear$year, " calBP)",sep="")
            main.proportion <- paste("Focal Proportion\n(",focalyear$year, " calBP / all years)",sep="")
            main.change <- paste("Focal Change\n(",focalyear$year, " calBP, from ",focalyear$backsight, " yrs before)",sep="")
        }
        if (type=="nonfocal"){
            if (!is.null(plotdir)){
                png(filename=paste(plotdir,"/nonfocal.png",sep=""), res=300, height=imdim[1], width=imdim[2], units="cm")
                par(mar=mars)
            }
            plot(x$nonfocal, box=box, main="", ...)
            if (maindef=="auto"){
                title(main.nonfocal, cex.main=cex.main, ...)
            } else {
                title(maindef, cex.main=cex.main, ...)
            }
            if (withpts=="y"){ points(x$ppp, pch=pch.pts, cex=cex.pts, col=col.pts) }
            if (!is.null(plotdir)){ dev.off() }
        } else if (type=="mask"){
            if (x$maskthresh==0){
                stop("No mask present.")
            }   
            if (!is.null(plotdir)){
                png(filename=paste(plotdir,"/mask.png",sep=""), res=300, height=imdim[1], width=imdim[2], units="cm")
                par(mar=mars)
            }
            plot(mymask, box=box, main="", ...)
            if (maindef=="auto"){
                title(main.mask, cex.main=cex.main, ...)
            } else {
                title(maindef, cex.main=cex.main, ...)
            }
            if (!is.null(plotdir)){ dev.off() }
        } else if (type=="focal"){
            if (!is.null(plotdir)){
                png(filename=filenames[a], res=300, height=imdim[1], width=imdim[2], units="cm")
                par(mar=mars)
            }
            if (ramptypes["focal"]=="raw"){
                bc <- colfun(ncolours)
                plot(focalyear$focal, col=bc, box=box, main="", ...)
            } else if (ramptypes["focal"]=="unl"){
                bc <- colfun(ncolours)
                plot(focalyear$focal, col=bc, box=box, main="", ribargs=list(at=focalyear$focal$yrange, labels=c("low","high"), las=2), ...)
            } else if (ramptypes["focal"]=="mmx"){
                if (is.null(zlimcheck)){ zlim <- c(min(x$stats$focal[,'min']), max(x$stats$focal[,'max'])) }
                bc <- colourmap(colfun(ncolours), range=zlim)
                plot(focalyear$focal, col=bc, box=box, main="", ...)
            } else if (ramptypes["focal"]=="std"){
                m <- mean(x$stats$focal[,"mean"])
                stdev <- mean(x$stats$focal[,"sd"])
                nz <- 3
                zlim <- c((m-(nz*stdev)),(m+(nz*stdev)))
                xmax <- max(x$stats$focal[,'max'])
                xmin <- min(x$stats$focal[,'min'])
                if (xmax > zlim[2] & xmin < zlim[1]){
                    bc <- colourmap(colfun(ncolours-2), range=zlim)
                    bc <- colourmap(col=c("#000D7F",attr(bc,"stuff")$outputs,"#FEEF23"), breaks=c(xmin, attr(bc,"stuff")$breaks, xmax))
                } else if (xmax < zlim[2] & xmin > zlim[1]){
                    bc <- colourmap(colfun(ncolours), range=zlim)
                } else if (xmax < zlim[2] & xmin < zlim[1]){
                    bc <- colourmap(colfun(ncolours-1), range=zlim)
                    bc <- colourmap(col=c("#000D7F",attr(bc,"stuff")$outputs), breaks=c(xmin, attr(bc,"stuff")$breaks))
                } else if (xmax > zlim[2] & xmin > zlim[1]){
                    bc <- colourmap(colfun(ncolours-1), range=zlim)
                    bc <- colourmap(col=c(attr(bc,"stuff")$outputs,"#FEEF23"), breaks=c(attr(bc,"stuff")$breaks, xmax))
                }
                plot(focalyear$focal, col=bc, box=box, main="", ...)     
            }
            if (!is.na(col.mask) & (focalyear$maskthresh>0)){ plot(mymask, add=TRUE, col=c("white",col.mask)) }
            if (withpts=="y"){ points(x$ppp, pch=pch.pts, cex=cex.pts, col=col.pts) }
            if (maindef=="auto"){
                title(main.focal, cex.main=cex.main, ...)
            } else {
                title(maindef, cex.main=cex.main, ...)
            }
            if (!is.null(plotdir)){ dev.off() }
        } else if (type=="proportion"){
            if (!is.null(plotdir)){
                png(filename=filenames[a], res=300, height=imdim[1], width=imdim[2], units="cm")
                par(mar=mars)
            }
            if (ramptypes["proportion"]=="raw"){
                bc <- topo.colors(ncolours)
                plot(focalyear$proportion, col=bc, ramptype=ramptype, box=box, main="", ...)
            } else if (ramptypes["proportion"]=="mmx"){
                if (is.null(zlimcheck)){ zlim <- c(min(x$stats$proportion[,'min']), max(x$stats$proportion[,'max'])) }
                bc <- colourmap(topo.colors(ncolours), range=zlim)
                plot(focalyear$proportion, col=bc, ramptype=ramptype, box=box, main="", ...)
            } else if (ramptypes["proportion"]=="unl"){
                bc <- topo.colors(ncolours)
                plot(focalyear$proportion, col=bc, ramptype=ramptype, box=box, main="", ribargs=list(at=focalyear$proportion$yrange, labels=c("low","high"), las=2), ...)
            } else if (ramptypes["proportion"]=="std"){
                m <- mean(x$stats$proportion[,"mean"])
                stdev <- mean(x$stats$proportion[,"sd"])
                nz <- 3
                zlim <- c((m-(nz*stdev)),(m+(nz*stdev)))
                xmax <- max(x$stats$proportion[,'max'])
                xmin <- min(x$stats$proportion[,'min'])
                if (xmax > zlim[2] & xmin < zlim[1]){
                    bc <- colourmap(topo.colors(ncolours-2), range=zlim)
                    bc <- colourmap(col=c("#4C00FFFF",attr(bc,"stuff")$outputs,"#FFE0B3FF"), breaks=c(xmin, attr(bc,"stuff")$breaks, xmax))
                } else if (xmax < zlim[2] & xmin > zlim[1]){
                    bc <- colourmap(topo.colors(ncolours), range=zlim)
                } else if (xmax < zlim[2] & xmin < zlim[1]){
                    bc <- colourmap(topo.colors(ncolours-1), range=zlim)
                    bc <- colourmap(col=c("#4C00FFFF",attr(bc,"stuff")$outputs), breaks=c(xmin, attr(bc,"stuff")$breaks))
                } else if (xmax > zlim[2] & xmin > zlim[1]){
                    bc <- colourmap(topo.colors(ncolours-1), range=zlim)
                    bc <- colourmap(col=c(attr(bc,"stuff")$outputs,"#FFE0B3FF"), breaks=c(attr(bc,"stuff")$breaks, xmax))
                }
                plot(focalyear$proportion, col=bc, box=box, main="", ribbon=ribbon, ...)
            }
            if (!is.na(col.mask) & (focalyear$maskthresh>0)){ plot(mymask, add=TRUE, col=c("white",col.mask)) }
            if (withpts=="y"){ points(x$ppp, pch=pch.pts, cex=cex.pts, col=col.pts) }
            if (maindef=="auto"){
                title(main.proportion, cex.main=cex.main, ...)
            } else {
                title(maindef, cex.main=cex.main, ...)
            }
            if (!is.null(plotdir)){ dev.off() }
        } else if (type=="change"){
            if (!is.null(plotdir)){
                png(filename=filenames[a], res=300, height=imdim[1], width=imdim[2], units="cm")
                par(mar=mars)
            } 
            if (ramptypes["change"]=="raw"){
                if (is.null(zlimcheck)){ zlim <- c(min(focalyear$change),max(focalyear$change)) }
                bc <- rybcolourmap(zlim, ncolours=ncolours)
                plot(focalyear$change, col=bc, box=box, main="", ribbon=ribbon, ...)
            } else if (ramptypes["change"]=="unl"){
                if (is.null(zlim)){ zlim <- c(min(focalyear$change),max(focalyear$change)) }
                bc <- rybcolourmap(zlim, ncolours=ncolours)
                plot(focalyear$change, col=bc, box=box, main="", ribbon=ribbon, ribargs=list(at=focalyear$change$yrange, labels=c("low","high"), las=2), ...)
            } else if (ramptypes["change"]=="mmx"){
                if (is.null(zlimcheck)){ zlim <- c(min(x$stats$change[,'min']), max(x$stats$change[,'max'])) }
                bc <- colourmap(colfun(ncolours), range=zlim)
                plot(focalyear$change, col=bc, box=box, main="", ...)
            } else if (ramptypes["change"]=="std"){
                m <- mean(x$stats$change[,"mean"])
                stdev <- mean(x$stats$change[,"sd"])
                nz <- 3
                zlim <- c((m-(nz*stdev)),(m+(nz*stdev)))
                xmax <- max(x$stats$change[,'max'])
                xmin <- min(x$stats$change[,'min'])
                if (xmax > zlim[2] & xmin < zlim[1]){
                    bc <- rybcolourmap(zlim, ncolours=ncolours-2)
                    bc <- colourmap(col=c("blue",attr(bc,"stuff")$outputs,"darkred"), breaks=c(xmin, attr(bc,"stuff")$breaks, xmax))
                } else if (xmax < zlim[2] & xmin > zlim[1]){
                    bc <- rybcolourmap(zlim, ncolours=ncolours)
                } else if (xmax < zlim[2] & xmin < zlim[1]){
                    bc <- rybcolourmap(zlim, ncolours=ncolours-1)
                    bc <- colourmap(col=c("blue",attr(bc,"stuff")$outputs), breaks=c(xmin, attr(bc,"stuff")$breaks))
                } else if (xmax > zlim[2] & xmin > zlim[1]){
                    bc <- rybcolourmap(zlim, ncolours=ncolours-1)
                    bc <- colourmap(col=c(attr(bc,"stuff")$outputs,"darkred"), breaks=c(attr(bc,"stuff")$breaks, xmax))
                }
                plot(focalyear$change, col=bc, box=box, main="", ribbon=ribbon, ...)
            }        
            if (!is.na(col.mask) & x$maskthresh  > 0){ plot(mymask, add=TRUE, col=c("white",col.mask)) }
            if (withpts=="y"){ points(x$ppp, pch=pch.pts, cex=cex.pts, col=col.pts) }
            if (maindef=="auto"){
                title(main.change, cex.main=cex.main, ...)
            } else {
                title(maindef, cex.main=cex.main, ...)
            }
            if (!is.null(plotdir)){ dev.off() }
        } else if (type=="all"){
            if (!is.null(plotdir)){
                png(filename=filenames[a], res=300, height=imdim[1], width=imdim[2], units="cm")
                par(mar=mars)
            } 
            par(mfrow=c(1,4))
            plot(x, focalyear$year, type="focal", ramptype=ramptype, ...)
            plot(x, focalyear$year, type="nonfocal", ramptype=ramptype, ...)
            plot(x, focalyear$year, type="proportion", ramptype=ramptype, ...)
            plot(x, focalyear$year, type="change", ramptype=ramptype, ...)
            if (!is.null(plotdir)){ dev.off() }
        }
    }
    if (verbose & length(focalyears)>1){
        close(pb)
        if (!is.null(plotdir) & tidydevs){ while (!is.null(dev.list())){ dev.off() } }
        print("Done.")
    }
}
