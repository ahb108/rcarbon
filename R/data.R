#' @title Radiocarbon dates from the EUROEVOL database
#' @description Radiocarbon dates (n=14,053) and site coordinates (n=4,213) from the EUROEVOL project database. Sites without radiocarbon dates (n=544), phase-codes, and other data have been omitted.

#' @format A data.frame with the following variables:
#' \describe{
#' \item{\code{C14ID}}{ID of each radiocarbon date}
#' \item{\code{C14Age}}{Radiocarbon age in 14C years BP}
#' \item{\code{C14SD}}{Radiocarbon age error}
#' \item{\code{LabCode}}{Labcode of the radiocarbon date}
#' \item{\code{Material}}{Material of the dated sample}
#' \item{\code{SiteID}}{ID of the site from which the sample has been recovered}
#' \item{\code{Latitude}}{Latitude of the sampling site in decimal degrees}
#' \item{\code{Longitude}}{Longitude of the sampling site in decimal degrees}
#' \item{\code{Country}}{Country where the sampling site is located}
#'}
#' @source Manning, K., Colledge, S., Crema, E., Shennan, S., Timpson, A., 2016. The Cultural Evolution of Neolithic Europe. EUROEVOL Dataset 1: Sites, Phases and Radiocarbon Data. Journal of Open Archaeology Data 5. doi:10.5334/joad.40
#' @references
#' Shennan, S., Downey, S.S., Timpson, A., Edinborough, K., Colledge, S., Kerig, T., Manning, K., Thomas, M.G., 2013. Regional population collapse followed initial agriculture booms in mid-Holocene Europe. Nature Communications 4, ncomms3486. doi:10.1038/ncomms3486
#'
#' Timpson, A., Colledge, S., Crema, E., Edinborough, K., Kerig, T., Manning, K., Thomas, M.G., Shennan, S., 2014. Reconstructing regional population fluctuations in the European Neolithic using radiocarbon dates: a new case-study using an improved method. Journal of Archaeological Science 52, 549-557. doi:10.1016/j.jas.2014.08.011
#'
#' @examples
#' \dontrun{
#' data(euroevol)
#' Ireland <- subset(euroevol,Country=="Ireland")
#' bins <- binPrep(Ireland$SiteID,Ireland$C14Age,h=200)
#' x <- calibrate(Ireland$C14Age,Ireland$C14SD)
#' spd.ireland <- spd(x,bins=bins,runm=200,timeRange=c(8000,4000))
#' plot(spd.ireland)
#'}
"euroevol"

#' @title Radiocarbon dates for the Eastern Mediterranean around the Younger Dryas
#' @description Radiocarbon dates (n=1915) and site coordinates (n=201) from a paper considering the relationship between human activity in the eastern Mediterranean/Middle East and early Holocene climate change, including the Younger Dryas.   

#' @format A data.frame with the following variables:
#' \describe{
#' \item{\code{LabID}}{Laboratory ID assigned to each radiocarbon date (where known)}
#' \item{\code{CRA}}{Radiocarbon age in 14C years BP}
#' \item{\code{Error}}{Radiocarbon age error}
#' \item{\code{Material}}{Material of the dated sample}
#' \item{\code{Species}}{Species of the dated sample (where identified)}
#' \item{\code{SiteName}}{Name of the site from which the sample has been recovered}
#' \item{\code{Country}}{Country where the sampling site is located}
#' \item{\code{Longitude}}{Longitude of the sampling site in decimal degrees}
#' \item{\code{Latitude}}{Latitude of the sampling site in decimal degrees}
#' \item{\code{Region}}{One of three analytical regions (1=southern Levant, 2=Northern Levant, 3= South-central Anatolia}
#'}
#' @source Palmisano, A., Bevan, A. and S. Shennan 2017. Data and code for demographic trends in the paper "Human responses and non-responses to climatic variations during the Last Glacial-Interglacial transition in the eastern Mediterranean", UCL Discovery Archive 1570274. doi:10.14324/000.ds.1570274.
#' @references
#' Roberts, N., Woodbridge, J., Bevan, A., Palmisano, A., Shennan, S. and E. Asouti 2017. Human responses and non-responses to climatic variations during the Last Glacial-Interglacial transition in the eastern Mediterranean. Quaternary Science Reviews, 184, 47-67. doi:10.1016/j.quascirev.2017.09.011.
#' 
#' @examples
#' \dontrun{
#' data(emedyd)
#' northernlevant <- emedyd[emedyd$Region=="2",]
#' bins <- binPrep(northernlevant$SiteName, northernlevant$CRA, h=50)
#' x <- calibrate(northernlevant$CRA, northernlevant$Error, normalised=FALSE)
#' spd.northernlevant <- spd(x, bins=bins, runm=50, timeRange=c(17000,8000))
#' plot(spd.northernlevant)
#'}
"emedyd"


#' @title Subset of EUROEVOL radiocarbon dates from Great Britain
#' @description Radiocarbon dates (n=2,324) and site coordinates (n=652) from England and Wales collected from the EUROEVOL project database. See \cite{\link{euroevol}} for more details regarding the source data.
#'
#' @format A data.frame with the following variables:
#' \describe{
#' \item{\code{C14ID}}{ID of each radiocarbon date}
#' \item{\code{C14Age}}{Radiocarbon age in 14C years BP}
#' \item{\code{C14SD}}{Radiocarbon age error}
#' \item{\code{LabCode}}{Labcode of the radiocarbon date}
#' \item{\code{Material}}{Material of the dated sample}
#' \item{\code{SiteID}}{ID of the site from which the sample has been recovered}
#' \item{\code{Eastings}}{Easting coordinates of the sampling site in meters (OSGB 1936 epsg:27700) }
#' \item{\code{Northings}}{Northing coordinates of the sampling site in meters (OSGB 1936 epsg:27700)}
#'}
#' @source Manning, K., Colledge, S., Crema, E., Shennan, S., Timpson, A., 2016. The Cultural Evolution of Neolithic Europe. EUROEVOL Dataset 1: Sites, Phases and Radiocarbon Data. Journal of Open Archaeology Data 5. doi:10.5334/joad.40
#'
"ewdates"


#' @title Polygonal window of England and Wales 
#' @description An \code{\link{owin}} class polygonal window of England and Wales.
#' @format An \code{owin} class object. 
#' @examples
#' \dontrun{
#' data(ewowin)
#' # Obtained from rworldmap:
#' library(maptools)
#' library(rgeos)
#' library(rworldmap)
#' bng <- CRS("+init=epsg:27700")
#' sbrit <- getMap(resolution="high")
#' sbrit <- spTransform(sbrit[sbrit$GEOUNIT =="United Kingdom",], bng)
#' tmp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(c(130000,130000,310000,
#' 425000,700000,700000,130000), c(0,300000,560000,598000,300000,0,0)))), 
#' "1")), 1:1, proj4string=bng)
#' sbrit <- gIntersection(sbrit,tmp) 
#' ewowin <- as.owin(sbrit)
#' }
#' @source South,A.2011.rworldmap: A New R package for Mapping Global Data. The R Journal Vol. 3/1 : 35-43.
"ewowin"

