% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregation.R
\name{stackspd}
\alias{stackspd}
\title{Stacked Summed Probability Distribution}
\usage{
stackspd(
  x,
  timeRange,
  bins = NA,
  group = NULL,
  datenormalised = FALSE,
  runm = NA,
  verbose = TRUE,
  edgeSize = 500
)
}
\arguments{
\item{x}{A \code{CalDates} class object containing the calibrated radiocarbon dates.}

\item{timeRange}{A vector of length 2 indicating the start and end date of the analysis in cal BP.}

\item{bins}{A vector containing the bin names associated with each radiocarbon date. If set to NA, binning is not carried out.}

\item{group}{A character or factor vector containing the grouping variable.}

\item{datenormalised}{Controls for calibrated dates with probability mass outside the timerange of analysis. If set to TRUE the total probability mass within the time-span of analysis is normalised to sum to unity. Should be set to FALSE when the parameter \code{normalised} in \code{\link{calibrate}} is set to FALSE. Default is FALSE.}

\item{runm}{A number indicating the window size of the moving average to smooth the SPD. If set to \code{NA} no moving average is applied. Default is NA}

\item{verbose}{A logical variable indicating whether extra information on progress should be reported. Default is TRUE.}

\item{edgeSize}{Controls edge effect by expanding the fitted model beyond the range defined by \code{timeRange}.}
}
\value{
An object of class \code{stackCalSPD}
}
\description{
Generates and combines multiple SPDs based on a user defined grouping.
}
\examples{
\dontrun{
data(emedyd)
x = calibrate(x=emedyd$CRA, errors=emedyd$Error,normalised=FALSE)
bins = binPrep(sites=emedyd$SiteName, ages=emedyd$CRA,h=50)
res = stackspd(x=x,timeRange=c(16000,8000),bins=bins,group=emedyd$Region)
}
}
