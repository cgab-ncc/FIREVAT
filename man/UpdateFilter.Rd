% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/firevat_filter.R
\name{UpdateFilter}
\alias{UpdateFilter}
\title{UpdateFilter}
\usage{
UpdateFilter(vcf.filter, param.values)
}
\arguments{
\item{vcf.filter}{A list from MakeFilterFromConfig}

\item{param.values}{A numeric vector contains filtering value
(same length with length(vcf.config.filter))}
}
\value{
Updated vcf.filter (list)
}
\description{
Update filter based on optim parameter values
}
