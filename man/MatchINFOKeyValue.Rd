% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/firevat_config.R
\name{MatchINFOKeyValue}
\alias{MatchINFOKeyValue}
\title{An internal function for parsing INFO column.}
\usage{
MatchINFOKeyValue(INFO.list)
}
\arguments{
\item{INFO.list}{A list of INFO values}
}
\value{
A vector of matching INFO values
}
\description{
Check whether an INFO entry is given with single key or key-value pair.
Return values accordingly.
}
\keyword{internal}
