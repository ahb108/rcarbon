#' Compute Great Arc Distances Matrix 
#'
#' This function computes and returns a Great Arc distance matrix from Latitude and Longitude data.
#'
#' @param Latitude A vector of Latitudes in decimal degrees
#' @param Longitude A vector of Longitudes in decimal degrees
#' @param progress If set to TRUE shows the progress bar is shown (default is FALSE)

#' @details Distances are calculated following the standard trigonometry formula
#'
#' @return An object of class matrix with distances in km

greatArcDist<-function(Latitude,Longitude,progress=FALSE)
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


