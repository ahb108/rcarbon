% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.spatialTest}
\alias{plot.spatialTest}
\title{Plot results of the local spatial permutation test of summed probability distributions.}
\usage{
\method{plot}{spatialTest}(
  x,
  index = 1,
  option,
  breakRange,
  breakLength = 7,
  rd = 5,
  baseSize = 0.5,
  plim = 0.05,
  qlim = 0.05,
  legend = FALSE,
  legSize = 1,
  location = "bottomright",
  ...
)
}
\arguments{
\item{x}{A \code{spatialTest} class object}

\item{index}{A numerical value indicating which transition to display. Ignored when \code{option="rawlegend"} or  \code{option="testlegend"}. Default is 1.}

\item{option}{Either "\code{raw}" to display local growth rates, "\code{test}" to display the test results (i.e. q and p values), or "\code{return}" to return a \code{sf} class object containing all relevant information for the given \code{index} value.}

\item{breakRange}{A vector of length 2 defining the minimum and maximum values of growth rate to be displayed in the legend. If set to NA its computed from data range (default).}

\item{breakLength}{A numerical vector defining the number of breaks for growth rates to be displayed in the legend.}

\item{rd}{Number of decimal places of the growth rate to be displayed in the Legend}

\item{baseSize}{Numerical value giving the amount by which points should be magnified relative to the default settings in R. Default is 0.5}

\item{plim}{Threshold value for the p-values. Default is 0.05.}

\item{qlim}{Threshold value for the q-values. Default is 0.05.}

\item{legend}{Logical values specifying whether the legend should be displayed or not. Default is FALSE.}

\item{legSize}{Numerical value giving the amount by which points should be magnified relative to the default settings in R for the Legend. Default is 1.}

\item{location}{A single keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center" to specify the location of the Legend. Default is "bottomright".}

\item{...}{Graphical parameters to be passed to methods.}
}
\description{
Displays local growth rates, p-values, and q-values retrieved from a \code{spatialTest} class object.
}
\details{
The function displays a distribution map of local growth rates (when \code{option="raw"}), q- and p-values (when \code{option="test"}), and the associated legends (when \code{option="rawlegend"} or  \code{option="testlegend"}).
}
\seealso{
\code{\link{sptest}}
}
