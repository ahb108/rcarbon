#' @title Sample random calendar dates. 
#'
#' @description Randomly samples calendar dates from each calibrated date or bin.
#' @param x A 'CalDates' class object. 
#' @param bins A vector containing the bin names associated with each radiocarbon date. If set to NA, binning is not carried out. 
#' @param nsim Number of sampling repetitions.
#' @param boot A logical value indicating whether bootstrapping is carried out (see details below). Default is FALSE.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' @details The function randomly samples calendar dates based from calibrated probability distributions. When the \code{bins} argument is supplied a single calendar date is sampled from each bin. When \code{boot=TRUE}, dates (or bins) are randomly sampled with replacement before calendar dates are sampled. 
#' 
#' @return An object of class \code{simdates} with the following elements
#' \itemize{
#' \item{\code{sdates}} {A matrix containing the randomly sampled calendar dates, with rows containing each of the \code{nsim} repetitions.}
#' \item{\code{weight}} {A vector (or matrix in when \code{boot=TRUE}) containing the total area under the curve of each date, normalised to sum to unity. Notice this will be identical for all dates if the calibration is carried out with the argument \code{normalised} set to TRUE.}
#'}
#'
#' @import stats 
#' @import utils 
#' @export


sampleDates <- function(x,bins,nsim,boot=FALSE,verbose=TRUE)
{
	if (!"CalDates" %in% class(x)){
		stop("x must be an object of class 'CalDates'.")
	}
	if (length(bins)>1){
		if (any(is.na(bins))){
			stop("Cannot have NA values in bins.")
		}
		if (length(bins)!=nrow(x$metadata)){
			stop("bins (if provided) must be the same length as x.")
		}
	}
	else {
		bins <- 1:(length(x))
	}

	uni.bins = unique(bins)
	nbins = length(uni.bins)
	binList = vector('list',length=nbins)
	binLength=numeric()

	calmatrix=FALSE
	if (anyNA(x$grids)){calmatrix=TRUE}

	if (calmatrix)
	{
	warning("Processing time is slower when dates are calibrated using calMatrix=TRUE",immediate.=TRUE)
	}

	if (verbose){
		flush.console()
		print("Aggregating...")
		pb <- txtProgressBar(min = 1, max = nbins,style = 3)
	}

	for (b in 1:nbins)
	{
		if (verbose){setTxtProgressBar(pb, b)}
		tmp.dates=x[which(bins==uni.bins[b])]
		if (calmatrix)
		{
			tmp.prob=tmp.dates$calmatrix
			binLength[b]=1
			if (length(tmp.dates)>1)
			{	
			tmp.prob=apply(tmp.dates$calmatrix,1,sum)
			binLength[b]=ncol(tmp.dates$calmatrix)
			}
			tmp.prob=tmp.prob[which(tmp.prob>0)]
			binList[[b]]=data.frame(calBP=as.numeric(names(tmp.prob)),PrDens=tmp.prob)
		} else {
			if (length(tmp.dates)==1)
			{
				binList[[b]]=tmp.dates$grids[[1]]
				binLength[b]=1
			} else {
				binLength[b]=length(tmp.dates)
				st.date = max(unlist(lapply(tmp.dates$grids,function(x){max(x$calBP)})))	
				end.date = min(unlist(lapply(tmp.dates$grids,function(x){min(x$calBP)}))) 
				tmp.prob = data.frame(calBP=st.date:end.date,PrDens=0)
				tmp.prob$PrDens=apply(sapply(tmp.dates$grids,function(x,r){
								     res = rep(0,length(r))
								     res[which(r%in%x$calBP)]=x$PrDens
								     return(res)},r=tmp.prob$calBP),1,sum)
				binList[[b]]=tmp.prob

			}

		}
	}

	if (verbose) {close(pb)}

	if (!boot)
	{
		if (verbose){print("Simulating dates ...")}
		res=sapply(binList,function(x,nsim){sample(x$calBP,size=nsim,prob=x$PrDens,replace=TRUE)},nsim=nsim)
		weight=sapply(binList,function(x){return(sum(x$PrDens))})/binLength
		weight=weight/sum(weight)
	}

	if (boot)
	{
		res = matrix(NA,nrow=nsim,ncol=nbins)
		weight = matrix(NA,nrow=nsim,ncol=nbins)
		if (verbose){
			print("Bootstrapping...")
			pb <- txtProgressBar(min = 1, max = nsim,style = 3)
		}

		for (i in 1:nsim)
		{
			
			if (verbose){setTxtProgressBar(pb, i)}
			index=sample(nbins,replace=TRUE)
			res[i,]=sapply(binList[index],function(x,nsim){sample(x$calBP,size=1,prob=x$PrDens)})
			weight[i,]=sapply(binList[index],function(x){sum(x$PrDens)})/binLength[index]
			weight[i,]=weight[i,]/sum(weight[i,])
		}
		if (verbose) {close(pb)}
	}
	if (verbose){print("done")}
	result=list(sdates=res,weight=weight)
	class(result) = c('simdates',class(result))
	return(result)
}



#' @title Composite Kernel Density Estimate SPD 
#'
#' @description Computes a Composite Kernel Density Estimate (CKDE) from multiple sets of randomly sampled calendar dates.
#' @param x A \code{simdates} class object, generated using \code{\link{sampleDates}}.
#' @param timeRange A vector of length 2 indicating the start and end date of the analysis in cal BP.
#' @param bw Kernel bandwith to be used.
#' @param normalised A logical variable indicating whether the contribution of individual dates should be equal (TRUE), or weighted based on the area under the curve of non-normalised calibration (FALSE). Default is TRUE.
#' @details The function computes Kernel Density Estimates using randomly sampled calendar dates contained in a \code{simdates} class object (generated using the \code{simulate.dates()} function). The output contains \code{nsim} KDEs, where \code{nsim} is the argument used in \code{simulate.dates()}. The resulting object can be plotted to visualise a CKDE (see Brown 2017), and if \code{boot} was set to \code{TRUE} in \code{simylate.dates()} its bootstraped variant (cf McLaughlin 2018 for a similar analysis). The shape of the CKDE is comparable to an SPD generated from non-normalised dates when the argument \code{normalised} is set to FALSE.
#' @return An object of class \code{ckdeSPD} with the following elements
#' \itemize{
#' \item{\code{timeRange}} {The \code{timeRange} setting used.}
#' \item{\code{res.matrix}} {A matrix containing the KDE values with rows representing calendar dates.}
#'}
#'
#' @references 
#' Brown, W. A. 2017. The past and future of growth rate estimation in demographic temporal frequency analysis: Biodemographic interpretability and the ascendance of dynamic growth models. \emph{Journal of Archaeological Science}, 80: 96–108. DOI: https://doi.org/10.1016/j.jas.2017.02.003 \cr
#' McLaughlin, T. R. 2018. On Applications of Space–Time Modelling with Open-Source 14C Age Calibration. \emph{Journal of Archaeological Method and Theory}. DOI https://doi.org/10.1007/s10816-018-9381-3
#'
#' @examples
#'data(emedyd)
#'x = calibrate(x=emedyd$CRA, errors=emedyd$Error,normalised=FALSE)
#'bins = binPrep(sites=emedyd$SiteName, ages=emedyd$CRA,h=50)
#'s = sampleDates(x,bins=bins,nsim=1000,boot=FALSE)
#'ckdeNorm = ckde(s1,timeRange=c(16000,9000),bw=100,normalised=TRUE)
#'ckdeNNorm = ckde(s1,timeRange=c(16000,9000),bw=100,normalised=FALSE)
#'plot(ckdeNorm)
#'plot(ckdeNNorm)
#'
#' @seealso \code{\link{sampleDates}}
#'
#' @import stats
#' @import utils
#' @export


ckde<- function(x,timeRange,bw,normalised=FALSE)
{

	if (!"simdates" %in% class(x)){
		stop("x must be an object of class 'simdates'.")
	}
	true.timeRange = c(max(x$sdates),min(x$sdates))
	nsim = nrow(x$sdates)
	boot = FALSE
	if (class(x$weight)=='matrix'){boot=TRUE}

	if (normalised)
	{
		raw.kde=apply(x$sdates,1,density,bw=bw,na.rm=TRUE,from=true.timeRange[1],to=true.timeRange[2])
		}
	if (!normalised)
	{
		if (length(unique(x$weight))==1)
		{
		warning("Simulated dates were generated from a dates calibrated with the argument `normalised` set to TRUE. The composite KDE will be executed with 'normalised' set to TRUE. To obtain a non-normalised composite KDE calibrate setting 'normalised` to FALSE")
		}
		if (!boot)
		{
			raw.kde=apply(x$sdates,1,density,bw=bw,na.rm=TRUE,from=true.timeRange[1],to=true.timeRange[2],weights=x$weight)
		} else {
			raw.kde = lapply(1:nrow(x$sdates),function(x,y,z,t1,t2,bw){
						 return(density(y[x,],bw=bw,na.rm=TRUE,weights=z[x,],from=t1,to=t2))},y=x$sdates,z=x$weight,t1=true.timeRange[1],t2=true.timeRange[2],bw=bw)
		}
	}
  res.matrix = matrix(NA,nrow=length(timeRange[1]:timeRange[2]),ncol=nsim)
  for (i in 1:nsim)
  {
	  #check this line
    res.matrix[,i]=approx(x=raw.kde[[i]]$x,xout=timeRange[1]:timeRange[2],y=raw.kde[[i]]$y)$y
  }

  result = list(timeRange=timeRange,res.matrix=res.matrix)
  class(result) = c("compositeKDE",class(result))
  return(result)
}

 
#' @description Plot a composite Kernel Density Estimate of sample dates.  
#' @param x A \code{compositeKDE} class object generated using the \code{ckde()} function.
#' @param calendar Either \code{'BP'} or \code{'BCAD'}. Indicate whether the calibrated date should be displayed in BP or BC/AD. Default is  \code{'BP'}.
#' @param type Either \code{envelope} or \code{multiline}. Default is \code{envelope}.
#' @param interval Percentile interval of the simulation envelope. Default is 0.95 (i.e. 95%).
#' @param xlim the x limits of the plot. In BP or in BC/AD depending on the choice of the parameter \code{calender}. Notice that if BC/AD is selected BC ages should have a minus sign (e.g. \code{c(-5000,200)} for 5000 BC to 200 AD).
#' @param alpha Alpha level for line transparency when \code{type='multiline'}. Default is 10/\code{nsim}, where \code{nsim} is the number of simulations. If \code{nsim} is smaller than 10, \code{alpha} will be set to 1.
#' @param fill.col Envelope color when \code{type='envelope'}. Default is 'lightgrey'.
#' @param line.col Line color when \code{type='envelope'}. Default is 'black.
#' @param line.type Line type when \code{type='envelope'}. Default is 2.
#' @param multiline.col Line color when \code{type='multiline'}. Default is 'black'.
#' @param ... Additional arguments affecting the plot
#' @details Visualise a \code{compositeKDE} class object. If \code{type} is set \code{'envelope'} an envelope of the percentile iterval defined by the parameter \code{interval} is shown along with the mean KDE. If \code{type} is set \code{'multiline'} all KDEs are shown. 
#' @seealso \code{\link{ckde}}; 
#' @examples
#' \dontrun{
#' data(emedyd)
#' levant <- emedyd[emedyd$Region=="1"|emedyd$Region=="2",]
#' bins <- binPrep(levant$SiteName, levant$CRA, h=50)
#' x <- calibrate(levant$CRA, levant$Error, normalised=FALSE)
#'}
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
#' @export 

plot.compositeKDE <- function(x, calendar="BP", type='envelope', ylim=NA, xlim=NA, alpha=NA, fill.col='lightgrey',line.col='black',line.type=2, interval=0.95,  multiline.col='black',...){

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
	    if (is.na(alpha)){alpha=10/ncol(x$res.matrix)}
	    if (alpha>1){alpha=1}
	    mc = c(as.numeric(col2rgb(multiline.col)/255),alpha)
	    apply(x$res.matrix,2,lines,x=plotyears,col=rgb(mc[1],mc[2],mc[3],mc[4]))
    }

    if (type=='envelope')
    {
	avg=apply(x$res.matrix,1,mean)   
	lo=apply(x$res.matrix,1,quantile,prob=(1-interval)/2)
   	hi=apply(x$res.matrix,1,quantile,prob=interval+(1-interval)/2)
	index = which(!is.na(hi))
	plot(x=plotyears,y=avg,type='n',xlab=xlabel,ylab="Summed Probability",xlim=xlim,ylim=ylim,axes=FALSE,...)
        polygon(x=c(plotyears[index],rev(plotyears[index])),y=c(lo[index],rev(hi[index])),border=NA,col=fill.col)
	lines(plotyears,avg,lty=line.type,col=line.col,lwd=2)
    }
    if (calendar=="BP"){
	rr <- range(pretty(plotyears))    
        axis(side=1,at=seq(rr[2],rr[1],-100),labels=NA,tck = -.01)
        axis(side=1,at=pretty(plotyears),labels=abs(pretty(plotyears)))
    } else if (calendar=="BCAD"){
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
    axis(2)
    box()
}




