#' @title rcarbon:Calibration and analysis of radiocarbon dates
#' 
#' @description The \pkg{rcarbon} package handles the calibration and analysis of radiocarbon, often but not exclusively for the purposes of archaeological research. It includes functions not only for basic calibration, uncalibration and plotting of one or more dates, but also a statistical framework for building demographic and related longitudinal inferences from aggregate radiocarbon date lists. 

#' @details
#'
#' Core functions in the \pkg{rcarbon} package can be grouped as follows:
#' \describe{
#'   \item{\strong{Calibration Functions}}{\code{\link{calibrate}} and \code{\link{uncalibrate}} enable the calibration and back-calibration for a variety of curves.}
#'   \item{\strong{Aggregation Functions}}{\code{\link{spd}} generates a summed probability distribution (SPD) of radiocarbon dates; \code{\link{binPrep}} can be used to define clusters of radiocarbon dates associated with the same context/phase}
#'     \item{\strong{Statistical Test Functions}}{\code{\link{modelTest}} compares the observed SPD against a variety of theoretical models (most typically an exponential curve) using the Monte-Carlo approach; \code{\link{p2pTest}} compares observed differences in SPD between two user-specified points in time against differences expected from a theoretical model; \code{\link{permTest}} compares two or more SPDs and test for the null hypothesis that all sets are derived from the same population; \code{\link{SPpermTest}} identifies, for defined intervals, locations with significantly higher or lower growth rate in the SPD compared to the pan-regional trend in the data}
#' }
#' @note
#' Up-to-date development version, bug-reports, and further information concerning  the \pkg{rcarbon} package can be found on GitHub (\url{https://github.com/ahb108/rcarbon}).
#' To see the preferred citation for the package, type citation("rcarbon").
#'
#'@references See individual functions for references.
#'
#' @author The \pkg{rcarbon} is developed and maintained by Andrew Bevan and Enrico Crema
#' @docType package
#' @name rcarbon
NULL
