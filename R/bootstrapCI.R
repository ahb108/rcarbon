# bootstrapCI <- function(x,sites,h,boot=100,timeRange=c(7000,3000),type="sites",verbose=TRUE,ncores=1,runm=NA,...)
# {
# if (ncores>1&!requireNamespace("doParallel", quietly=TRUE)){	
#         warning("the doParallel package is required for multi-core processing; ncores has been set to 1")
#         ncores=1
#     }	
# 
#         cralist <- x$metadata$CRA
# 
##obs
# obs.bins <- binPrep(ages=cralist,h=h,sites=sites)
# obs.spd <- spd(x=x,bins=obs.bins,timeRange=timeRange,spdnormalised=TRUE,verbose=F,runm=runm,...)
# 
##bootstrap
# resmatrix <- matrix(NA,nrow=nrow(obs.spd$grid),ncol=boot)
# 
# if (verbose){
#         print("Bootstrapping..")
#         flush.console()
#         pb <- txtProgressBar(min=1, max=boot, style=3, title="Checking overlaps within each bin...")
#     }
# 
# if (type == "spd")
# {
#         errors <- x$metadata$Error
#         datenormalised <- x$metadata$Normalised[1]
#         samplesize <- length(unique(obs.bins))
#           predgrid <- obs.spd$grid
#         cragrid <- uncalibrate(as.CalGrid(predgrid), verbose=FALSE)
#         obscras <- x$metadata$CRA
#         cragrid$PrDens[cragrid$CRA > max(obscras) | cragrid$CRA < min(obscras)] <- 0
# }
# 
# 
# for (b in 1:boot)
# {
#   if (verbose){ setTxtProgressBar(pb, b) }
#         
#   if (type=="sites")
#     {
#      boot.sites<-sample(unique(sites),replace=TRUE)
#      osites <- numeric()
#      newsites<-numeric()
#      index <-numeric()
#      k=0
#         for (j in 1:length(boot.sites))
#         {
#           k=k+1
#           i<-which(sites==boot.sites[j])
#           index=c(index,i)
#           osites=c(osites,sites[i])
#           newsites=c(newsites,rep(k,length(i)))
#         }
#      boot.caldates=x[index]
#      boot.cralist<-boot.caldates$metadata$CRA
#      boot.bins <- binPrep(ages=boot.cralist,h=h,sites=newsites)
#      boot.spd <- spd(x=boot.caldates,bins=boot.bins,timeRange=timeRange,spdnormalised=TRUE,verbose=F,runm=runm,...)
#      resmatrix[,b]=boot.spd$grid$PrDens
#     }
# 
#   if (type=="bins")
#     {
#      boot.bins<-sample(unique(obs.bins),replace=TRUE)
#      k=0
#      obins <- numeric()
#      newbins <- numeric()
#      index <- numeric()
#      for (j in 1:length(boot.bins))
#         {
#           k=k+1
#           i<-which(obs.bins==boot.bins[j])
#           index=c(index,i)
#           obins=c(obins,obs.bins[i])
#           newbins=c(newbins,rep(k,length(i)))
#         }
#      boot.caldates=x[index]
#      boot.spd <- spd(x=boot.caldates,bins=newbins,timeRange=timeRange,spdnormalised=TRUE,verbose=F,runm=runm,...)
#      resmatrix[,b]=boot.spd$grid$PrDens
#     }
#   if (type=="spd")
#   {
#         randomDates <- sample(cragrid$CRA, replace=TRUE, size=samplesize, prob=cragrid$PrDens)
#         randomSDs <- sample(size=length(randomDates), errors, replace=TRUE)
#         tmp <- calibrate(x=randomDates,errors=randomSDs, resOffsets=0 ,resErrors=0, timeRange=timeRange, calCurves='intcal13', normalised=datenormalised, ncores=ncores, verbose=FALSE, calMatrix=TRUE)
#         simDateMatrix <- tmp$calmatrix
#         resmatrix[,b] <- apply(simDateMatrix,1,sum)
#         resmatrix[,b] <- (resmatrix[,b]/sum(resmatrix[,b])) 
#         if (!is.na(runm)){
#             resmatrix[,b] <- runMean(resmatrix[,b], runm, edge="fill")
#         }
#     }
#   }
#    
# if (verbose){ close(pb)
# print("Done.") }
# return(list(obs.spd=obs.spd,bootmatrix=resmatrix))
# }
# 
# bootPlot <- function(x,ci=0.95,calendar="BP")
# {
# if (!calendar %in% c("BP","BCAD")){ stop("Unknown calendar type") }
# timeRange=rev(range(x$obs.spd$grid$calBP))
#   
#   years <- timeRange[1]:timeRange[2]
#   xlab <- "Years BP"
#   xr <- timeRange
#   if (calendar=="BCAD")
#   {
#    years <- 1950 - years
#    xlab <- "Years BC/AD"
#    xr <- range(years)
#   }
#         
# 
# yrange=range(c(x$bootmatrix,x$obs.spd$grids$PrDens))
# plot(years,x$obs.spd$grids$PrDens,ylim=yrange,xlim=xr,xlab=xlab,type="n",axes=F)
# 
# lo<-apply(x$bootmatrix,1,quantile,prob=0+(1-ci)/2)
# hi<-apply(x$bootmatrix,1,quantile,prob=ci+(1-ci)/2)
# 
# polygon(c(years,rev(years)),c(lo,rev(hi)),border="lightgrey",col="lightgrey")
# lines(years,x$obs.spd$grid$PrDens,lwd=2)
# 
# axis(side=2)
# 
# if (calendar=="BP")
# {
# axis(side=1)
# }
# 
# if (calendar=="BCAD")
#   {
#    xticksAt=pretty(years)
#    xticksLab=xticksAt
#    if (any(xticksLab==0)){xticksLab[which(xticksLab==0)]=1}
#    if (any(xticksAt>1)){xticksAt[which(xticksAt>1)]=xticksAt[which(xticksAt>1)]-1}
#    axis(side=1,at=xticksAt,labels=xticksLab)
#   }  
# box()
# 
# }
# 
# 
# 
