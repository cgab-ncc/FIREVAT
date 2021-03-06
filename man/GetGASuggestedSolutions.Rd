% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/firevat_brute_force.R
\name{GetGASuggestedSolutions}
\alias{GetGASuggestedSolutions}
\title{GetGASuggestedSolutions}
\usage{
GetGASuggestedSolutions(
  vcf.obj,
  bsg,
  config.obj,
  lower.upper.list,
  df.mut.pat.ref.sigs,
  target.mut.sigs,
  sequencing.artifact.mut.sigs,
  objective.fn,
  original.muts.seq.art.weights.sum,
  ga.preemptive.killing,
  verbose = TRUE
)
}
\arguments{
\item{vcf.obj}{A list from ReadVCF}

\item{bsg}{BSgenome.Hsapiens.UCSC object}

\item{config.obj}{A list from ParseConfigFile}

\item{lower.upper.list}{A list from GetParameterLowerUpperVector}

\item{df.mut.pat.ref.sigs}{A data.frame from MutPatParseRefMutSigs}

\item{target.mut.sigs}{A character vector of the target mutational signatures from reference mutational signatures.}

\item{sequencing.artifact.mut.sigs}{A character vector of the sequencing artifact mutational signatures from reference mutational signatures.}

\item{objective.fn}{Objective value derivation function.}

\item{original.muts.seq.art.weights.sum}{A numeric value. 'seq.art.sigs.weights.sum' from CheckIfVariantRefinementIsNecessary}

\item{ga.preemptive.killing}{If TRUE, then preemptively kills populations that yield greater sequencing artifact weights sum
compared to the original mutatational signatures analysis}

\item{verbose}{If TRUE, provides process detail. Default value is TRUE.}
}
\value{
A list with the following elements
\itemize{
 \item{judgment}{A boolean value}
 \item{seq.art.sigs.weights}{A numeric value. Sum of sequencing artifact weights.}
}
}
\description{
Computes suggested solutions
}
