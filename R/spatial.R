## Compute weights from distance matrix
defineNeighbour<-function(distmat,h=NULL,kernel="gaussian")
{
    w=matrix(NA,nrow=nrow(distmat),ncol=ncol(distmat))
    kernels <- c("gaussian","fixed")
    if (!kernel %in% kernels){
                stop("The kernel you have chosen is not currently an option.")
           }
    if (is.null(h))
    {
                stop("Distance parameter h undefined")  
    }
    for (x in 1:nrow(distmat))
        {
            if (kernel=="gaussian")
                {w[x,]=exp(-distmat[x,]^2/h^2)}
            if (kernel=="fixed")
                {w[x,]=as.numeric(distmat[x,]<=h)}
        }
    res=list(w=w,h=h,kernel=kernel)
    class(res) <- append(class(res),"spatialweights")
    return(res)   
}



## Compute Great Arc Distances Matrix 

greatArcDist<-function(Latitude,Longitude,verbose=FALSE)
    {

	if (length(Latitude)!=length(Longitude))
	{
        stop("Latitude and Longitude must have the same length.")
	}

        n=length(Latitude)
        res<-matrix(0,nrow=n,ncol=n)

	 if (verbose)
	 {
         print("Calculating great-arc distances...")
         flush.console()
         pb <- txtProgressBar(min = 1, max =n, style=3)
	 }
        for (i in 1:n)
            {
                 if (verbose){setTxtProgressBar(pb, i)}
                for (j in 1:n)
                    {
			if ((i>j)&((Latitude[i]!=Latitude[j])|(Longitude[i]!=Longitude[j])))
			{   	
                                res[i,j]<- c(60 * (180/pi) * acos(sin((pi/180)*Latitude[i]) * sin((pi/180)*Latitude[j])
								  + cos((pi/180)*Latitude[i]) * cos((pi/180)*Latitude[j]) *
									  cos((pi/180)*(Longitude[j] - Longitude[i])))) * 1852/1000

                        }
		    }
            }
        if (verbose)
	{ 
	close(pb)
        print("Done.") 
	}	
       return(as.matrix(as.dist(res)))
    }


##Core spatialSPD function##

spSPDpermtest<-function(calDates, timeRange, bins, locations, breaks, spatialweights, nsim=1000, runm=NA, verbose=TRUE,permute="locations",ncores=1,datenormalised=FALSE)
{

###################################
#### Load Dependency Libraries ####
###################################
    require(sp)	
    require(fdrtool)
    require(rcarbon)
    if (ncores>1) {require(doParallel)}

##################################
#### Initial warning messages ####
##################################

     if (!"CalDates" %in% class(calDates)){
        stop("calDates must be an object of class 'calDates'.")
    }
    if (length(bins)>1){
         if (any(is.na(bins))){
            stop("Cannot have NA values in bins.")
        }
        if (length(bins)!=length(calDates$grid)){
            stop("bins (if provided) must be the same length as x.")
        }
         } else {
        bins <- rep("0_0",length(calDates$grid))
    }

   if (!("SpatialPoints" %in% class(locations)[1]|"SpatialPointsDataFrame" %in% class(locations)[1])){
        stop("locations must be an object of class 'SpatialPoints' or 'SpatialPointsDataFrame'.")
    }

   locations.id=row.names(locations@coords)
    if (is.null(locations.id))
    {
        stop("locations must have rownames")
    }

   

#############################
#### Create binnedMatrix ####
#############################

    binNames <- unique(bins)
    calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
    binnedMatrix <- matrix(NA, nrow=nrow(calyears), ncol=length(binNames))


    if (verbose & length(binNames)>1){
        print("Binning by site/phase...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(binNames), style=3, title="Binning by site/phase...")
    }
    for (b in 1:length(binNames)){
        if (verbose & length(binNames)>1){ setTxtProgressBar(pb, b) }
        index <- which(bins==binNames[b])
        slist <- calDates$grid[index]
        slist <- lapply(slist,FUN=function(x) merge(calyears,x, all.x=TRUE)) 
        slist <- rapply(slist, f=function(x) ifelse(is.na(x),0,x), how="replace")
        slist <- lapply(slist, FUN=function(x) x[with(x, order(-calBP)), ])
        tmp <- lapply(slist,`[`,2)
        if (datenormalised){
            tmp <- lapply(tmp,FUN=function(x) x/sum(x))
        }
        if (length(binNames)>1){
            spd.tmp <- Reduce("+", tmp) / length(index)
        } else {
            spd.tmp <- Reduce("+", tmp)
        }
	binnedMatrix[,b] <- spd.tmp[,1]
    }
    if (verbose & length(binNames)>1){ close(pb) }


################################
### Observed Data Subroutine ###
################################ 


## Aggregate by Locations ##
    origins=unlist(lapply(strsplit(binNames,"_"),function(x){x[[1]]}))

    if (!all(origins%in%locations.id))
     {
        stop("Missing bins or locations")
     }

    resMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nrow(binnedMatrix))
    
    for (x in 1:length(unique(locations.id))) 
        {
            index=which(origins==unique(locations.id)[x])
            if(length(index)>1)
                {resMatrix[x,]=apply(binnedMatrix[,index],1,sum)}
            if(length(index)==1)
                {resMatrix[x,]=binnedMatrix[,index]}
        }

 ## Aggregate by break s##

    nBreaks=length(breaks)-1
    obsMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nBreaks)
    timeSequence=timeRange[1]:timeRange[2]
    
    for (x in 1:nBreaks)
        {
            index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
            obsMatrix[,x]=apply(resMatrix[,index],1,sum)
        }

## Apply SpatialWeights ##

    obsGridVal=t(spatialweights$w)%*%obsMatrix

## Compute Rate of Change #3

    rocaObs=t(apply(obsGridVal,1,function(x,d){
		   L=length(x)
		   res=numeric(length=L-1)
		   for (i in 1:c(L-1))
			{
			res[i]=(x[i+1]/x[i])^(1/d)-1
			}
		   return(res)},
		   d=abs(breaks[2]-breaks[1])))

##############################
### Permutation Subroutine ###
############################## 

    if (ncores>1)
   	 {	
          cl <- makeCluster(ncores)
          registerDoParallel(cl)
          print(paste("Running permutation test in parallel on ",getDoParWorkers()," workers...",sep=""))
	  sumcombine<-function(a,b)
		{
		list(a[[1]]+b[[1]],a[[2]]+b[[2]],a[[3]]+b[[3]])
		}
	  resultHiLoEq<-foreach (x=1:nsim,.combine= sumcombine) %dopar% {

            simGridVal<-matrix(NA,nrow=nrow(spatialweights$w),ncol=nBreaks)
            
	    ## Aggregate by Site ## 

            simResMatrix=matrix(0,nrow=length(unique(locations.id)),ncol=nrow(binnedMatrix))

            ## Randomly assigne bins to locations.id ##
           
	    if (permute=="bins")
	    {    
	    simOrigins=sample(origins)
            for (x in 1:length(unique(locations.id)))
                {                    
                    index=which(simOrigins==unique(locations.id)[x])
                    if(length(index)>1)
                        {simResMatrix[x,]=apply(binnedMatrix[,index],1,sum)}
                    if(length(index)==1)
                        {simResMatrix[x,]=binnedMatrix[,index]}
                }
            

            ## Aggregate by breaks ##
	    
            aggMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nBreaks)
            
            for (x in 1:nBreaks)
                {
                    index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
                    aggMatrix[,x]=apply(simResMatrix[,index],1,sum)
                }
		       

           ## Apply Weights ##

           simGridVal=t(spatialweights$w)%*%aggMatrix
	    }
	    if (permute=="locations")
	    {
	     simMatrix=obsMatrix[sample(nrow(obsMatrix)),]	
             simGridVal=t(spatialweights$w)%*%simMatrix
		
	    }


           ## Compute Rate of Change ##

           rocaSim=t(apply(simGridVal,1,function(x,d){
		   L=length(x)
		   res=numeric(length=L-1)
		   for (i in 1:c(L-1))
			{
			res[i]=(x[i+1]/x[i])^(1/d)-1
			}
		   return(res)},
		   d=abs(breaks[2]-breaks[1])))

	    lo=rocaObs<rocaSim	    
	    hi=rocaObs>rocaSim
	    eq=rocaObs==rocaSim

          return(list(hi,lo,eq))
	  }
        stopCluster(cl)

        lo=resultHiLoEq[[1]]
	hi=resultHiLoEq[[2]]
	eq=resultHiLoEq[[3]]
    
	} else {

    hi=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)
    lo=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)
    eq=matrix(0,nrow=nrow(spatialweights$w),ncol=nBreaks-1)

    print("Permutation test...")
    flush.console()

    pb <- txtProgressBar(min = 1, max = nsim, style=3)

        for (s in 1:nsim)
        {
            setTxtProgressBar(pb, s)
	    simGridVal<-matrix(NA,nrow=nrow(spatialweights$w),ncol=nBreaks)
            ## Aggregate by Site ## 
            simResMatrix=matrix(0,nrow=length(unique(locations.id)),ncol=nrow(binnedMatrix))

            ## Randomly assign bins to locations
            if (permute=="bins")
	    {
	    simOrigins=sample(origins)
            



            for (x in 1:length(unique(locations.id)))
                {                    
                    index=which(simOrigins==unique(locations.id)[x])
                    if(length(index)>1)
                        {simResMatrix[x,]=apply(binnedMatrix[,index],1,sum)}
                    if(length(index)==1)
                        {simResMatrix[x,]=binnedMatrix[,index]}
                }
            

            ##Aggregate by breaks##
            aggMatrix=matrix(NA,nrow=length(unique(locations.id)),ncol=nBreaks)
            
            for (x in 1:nBreaks)
                {
                    index=which(timeSequence<=breaks[x]&timeSequence>breaks[x+1])
                    aggMatrix[,x]=apply(simResMatrix[,index],1,sum)
                }
		       

           ##Apply Weights 
           simGridVal=t(spatialweights$w)%*%aggMatrix
	    }
           if (permute=="locations")
           {
	     simMatrix=obsMatrix[sample(nrow(obsMatrix)),]	
             simGridVal=t(spatialweights$w)%*%simMatrix
	   }



           ##Compute Rate of Change
           rocaSim=t(apply(simGridVal,1,function(x,d){
		   L=length(x)
		   res=numeric(length=L-1)
		   for (i in 1:c(L-1))
			{
			res[i]=(x[i+1]/x[i])^(1/d)-1
			}
		   return(res)},
		   d=abs(breaks[2]-breaks[1])))

	    hi=hi+(rocaObs>rocaSim)
	    lo=lo+(rocaObs<rocaSim)
	    eq=eq+(rocaObs==rocaSim)
	    
        }
    close(pb)
    }


############################
### Compute Significance ###
############################ 
    
    pvalHi=(lo+eq+1)/c(nsim+1)
    pvalLo=(hi+eq+1)/c(nsim+1)
    pval=pvalHi
    pval[which(pvalHi>pvalLo)]=pvalLo[which(pvalHi>pvalLo)]
    pval=pval*2
    if (max(pval)>1)
    {
    	 pval[which(pval>1)]=1
    }

    ## Compute False Discovery Rate (q-value) ##
    
    qvalHi=apply(pvalHi,2,function(x){return(fdrtool(x,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval)})
    qvalLo=apply(pvalLo,2,function(x){return(fdrtool(x,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval)})
    qval=apply(pval,2,function(x){return(fdrtool(x,statistic="pvalue",plot=FALSE,verbose=FALSE)$qval)})

    metadata=data.frame(npoints=length(unique(locations.id)),ndates=nrow(calDates$metadata),nbins=length(binNames),nsim=nsim,permutationType=permute,datenormalised=datenormalised,breaks=nBreaks,timeRange=paste(timeRange[1],"-",timeRange[2],sep=""),weights.h=spatialweights$h,weights.kernel=spatialweights$kernel)
   
    reslist=list(metadata=metadata,rocaObs=rocaObs,pval=pval,pvalHi=pvalHi,pvalLo=pvalLo,qval=qval,qvalLo=qvalLo,qvalHi=qvalHi,locations=locations)
    
    class(reslist) <- append(class(reslist),"spatialTest")
    return(reslist)
}


print.spSPD<-function(x)
{
print(x$metadata)
}

##Plot function for spatial SPD#


plot.spSPD<-function(x,index=NULL,basemap=TRUE,baseSize=0.5)
{
	if (!any(class(x)%in%c("spatialTest")))
	{
        stop("x is not a spatialTest class object")
	}

        if (is.null(index))
	{
        stop("index value missing")
	}

        require(sp)	
        locations=x$locations
	projection=strsplit(proj4string(locations),split=" ")[[1]]

	if (basemap) 
         {
	library(rworldmap)
	library(raster)
	library(maptools)
	if (!"+proj=longlat"%in%projection&!"+datum=WGS84"%in%projection)
	{
        stop("basemap available only for LatLong data with WGS84")
	}
	base=getMap(resolution="low")
         }

		
		nBreaks=ncol(x$rocaObs)


	plusPoints=locations[which(x$pvalHi[,index]>0.5),]
	minusPoints=locations[which(x$pvalHi[,index]<0.5),]

	# Set Base 
	par(mar=c(0.1,0.1,0,0.5))
	plot(locations,col=NA,xlab="",ylab="",axes=FALSE)

	if(basemap) {plot(base,col="grey75",border="beige",add=TRUE)}
	points(plusPoints,col="gold",pch=20,cex=0.5)
	points(minusPoints,col="lightblue",pch=20,cex=0.5)


	# Set Positive
	positive.index=which(x$pvalLo[,index]<=0.05)
	
	if (length(positive.index)>0)
		{
		positive=locations[positive.index,]
		points(positive,pch=20,col="orange",cex=0.5)
		qpositive.index=which(x$qvalLo[,index]<=0.05&x$pvalLo[,index]<=0.05) #Originally based on qvalHi
		if (length(qpositive.index)>0)
			{
				qpositive=locations[qpositive.index,]
				points(qpositive,pch=20,col="red",cex=0.5)

			}
		}
	negative.index=which(x$pvalHi[,index]<=0.05)
	
	if (length(negative.index)>0)
		{
		negative=locations[negative.index,]
		points(negative,pch=20,col="cornflowerblue",cex=0.5)
		qnegative.index=which(x$qvalHi[,index]<=0.05&x$pvalHi[,index]<=0.05) #Originally based on qvalLo
		if (length(qnegative.index)>0)
			{
				qnegative=locations[qnegative.index,]
				points(qnegative,pch=20,col="darkblue",cex=0.5)

			}

		}

}

rybcolourmap <- function(range, ...) {
  col <- rybcolours(range, ...)
  z <- colourmap(col, range=range)
  return(z)
}

rybcolours <- function(range, sealevel=0, ncolours=100,nbeach=0){
  stopifnot(is.numeric(range) && length(range)==2)
  stopifnot(all(is.finite(range)))
  yr <- colorRampPalette(c("yellow","orangered","darkred"), space="rgb")
  cb <- colorRampPalette(c("blue","cyan","yellow"), space="rgb")
  depths <- range[1]
  peaks <- range[2]
  dv <- diff(range)/(ncolours - 1)
  epsilon <- nbeach * dv/2
  lowtide <- max(sealevel - epsilon, depths)
  hightide <-  min(sealevel + epsilon, peaks)
  countbetween <- function(a, b, delta) { max(0, round((b-a)/delta)) }
  nsea <- countbetween(depths, lowtide, dv)
  nbeach <- countbetween(lowtide,  hightide, dv)
  nland <- countbetween(hightide,  peaks, dv)
  colours <- character(0)
  if(nsea > 0)  colours <- cb(nsea) # cyan/blue
  if(nbeach > 0)  colours <- c(colours,rep("yellow",nbeach)) # yellow
  if(nland > 0)  colours <- c(colours, yr(nland)) # darkred/yellow
  return(colours)
}

spJitter <- function(pts, xamount, yamount=xamount){
    proj <- NA
    if (!is.na(proj4string(pts)) | proj4string(pts)!="NA"){
        proj <- proj4string(pts)
    }
    if (class(pts) == "SpatialPointsDataFrame"){
        df <- cbind(coordinates(pts),pts@data)
        df[,1] <- jitter(df[,1],amount=xamount)
        df[,2] <- jitter(df[,2],amount=yamount)
        coordinates(df) <- df[,1:2]
        proj4string(df) <- proj
    } else if (class(pts) == "SpatialPoints"){
        df <- coordinates(pts)
        df[,1] <- jitter(df[,1],amount=xamount)
        df[,2] <- jitter(df[,2],amount=yamount)
        df <- SpatialPoints(df, proj4string=CRS(proj))
    } else {
        stop("Only works for SpatialPoints* at present.")
    }
    return(df)
}




