# FIREVAT Mutational Patterns Functions
# Mutational Patterns (Blokzijl et al., Genome Medicine 2018; PMID 29695279)
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title MutPatParseVCFObj
#' @description Parses a vcf.obj and prepares it to run Mutational Patterns.
#'
#' @param vcf.obj A list from ReadVCF
#' @param bsg A BSgenome object
#' @param sample.id A string value
#'
#' @return A data.frame with the column sample.id and
#' row names corresponding to 96 substitution types
#'
#' @export
#' @importFrom deconstructSigs mut.to.sigs.input
MutPatParseVCFObj <- function(vcf.obj, bsg, sample.id = "sample") {
    vcf.obj$data$Sample <- rep(sample.id, nrow(vcf.obj$data))
    df.deconstructsigs.sigs.input <- mut.to.sigs.input(mut.ref = vcf.obj$data,
                                                       sample.id = "Sample",
                                                       chr = "CHROM",
                                                       pos = "POS",
                                                       ref = "REF",
                                                       alt = "ALT",
                                                       bsg = bsg)

    df.mutational.patterns.input <- t(df.deconstructsigs.sigs.input)
    df.mutational.patterns.input <- df.mutational.patterns.input + 0.0001
    return(df.mutational.patterns.input)
}


#' @title MutPatParseRefMutSigs
#' @description Parses a df.ref.mut.sigs and prepares it to run Mutational Patterns.
#'
#' @param df.ref.mut.sigs A data.frame of reference mutational signatures
#' @param target.mut.sigs A character vector of target mutational signatures names
#' @param signature.start.column.index = An integer value (e.g. column index corresponding to 'SBS1')
#' @param mutation.type.header = A string value
#'                               (name of header corresponding to column containing 'A[C>A]A' data))
#'
#' @return A data.frame of the format deconstructSigs::signatures.cosmic
#'
#' @export
MutPatParseRefMutSigs <- function(df.ref.mut.sigs,
                                  target.mut.sigs,
                                  signature.start.column.index = 4,
                                  mutation.type.header = "SomaticMutationType") {
    # Only consider target.mut.sigs
    if (length(target.mut.sigs) > 0) {
        df.ref.mut.sigs <- df.ref.mut.sigs[, c(mutation.type.header, target.mut.sigs)]
    }

    # Parse df.ref.mut.sigs
    row.names(df.ref.mut.sigs) <- df.ref.mut.sigs[, mutation.type.header]
    signature.names <- colnames(df.ref.mut.sigs)[2:ncol(df.ref.mut.sigs)]
    df.ref.mut.sigs <- df.ref.mut.sigs[,signature.names]
    df.ref.mut.sigs <- as.matrix(df.ref.mut.sigs)
    return(df.ref.mut.sigs)
}


#' @title RunMutPat
#' @description Identifies mutational signatures using Mutational Patterns
#'
#' @param mut.pat.input A list from \code{\link{MutPatParseVCFObj}}
#' @param df.mut.pat.ref.sigs A data.frame returned by \code{\link{MutPatParseRefMutSigs}}
#' @param target.mut.sigs A character vector of target mutational signatures names
#' @param verbose If true, provides process details
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{tumor.mutation.types.spectrum}{A numeric vector of length 96 - 'observed' spectrum}
#'  \item{identified.mutation.types.spectrum}{A numeric vector of length 96 - 'identified' spectrum}
#'  \item{residuals}{A numeric vector of length 96 - residuals}
#'  \item{mutation.types}{A character vector of length 96}
#'  \item{identified.mut.sigs}{A character vector where each element is a mutational signature identified}
#'  \item{identified.mut.sigs.contribution.weights}{A numeric vector where each element is the weight of mutational signature identified. The ordering follows identified.mut.sigs}
#'  \item{cosine.similarity.score}{A numeric value}
#' }
#'
#' @export
#' @importFrom MutationalPatterns fit_to_signatures
#' @examples
#' \dontrun{
#' vcf.obj <- ReadVCF(vcf.file = "../data/sample/HNT-082-BT.final.call.vcf", genome = "hg19")
#' df.ref.mut.sigs <- GetPCAWGMutSigs()
#' target.mut.sigs <- GetPCAWGMutSigsNames()
#' RunMutPat(vcf.obj = vcf.obj,
#' df.ref.mut.sigs = df.ref.mut.sigs,
#' target.mut.sigs = target.mut.sigs)
#' }
RunMutPat <- function(mut.pat.input,
                      df.mut.pat.ref.sigs,
                      target.mut.sigs,
                      verbose = TRUE) {
    if (verbose == TRUE) {
        PrintLog("* Started running Mutational Patterns")
    }

    # Take only unique mutational signatures
    target.mut.sigs <- unique(target.mut.sigs)

    # Mutational Patterns - identify mutational signatures
    mut.pat.results <- fit_to_signatures(mut.pat.input, df.mut.pat.ref.sigs)

    # Wrap Mutational Patterns results
    df.mut.pat.results <- data.frame(list(sig = rownames(mut.pat.results$contribution),
                                          weight = unname(mut.pat.results$contribution)),
                                     stringsAsFactors = F)
    df.mut.pat.results <- df.mut.pat.results[(df.mut.pat.results$weight > 0), ]

    # Compute cosine similarity score
    cosine.similarity.score <- lsa::cosine(as.numeric(mut.pat.results$reconstructed),
                                           as.numeric(mut.pat.input))

    # Prepare return data
    r <- list(tumor.mutation.types.spectrum = as.numeric(mut.pat.input),
              identified.mutation.types.spectrum = as.numeric(mut.pat.results$reconstructed),
              residuals = as.numeric(mut.pat.input) - as.numeric(mut.pat.results$reconstructed),
              mutation.types = row.names(mut.pat.results$reconstructed),
              identified.mut.sigs = df.mut.pat.results$sig,
              identified.mut.sigs.contribution.weights = df.mut.pat.results$weight,
              cosine.similarity.score = cosine.similarity.score)

    if (verbose == TRUE) {
        PrintLog("* Finished running Mutational Patterns")
    }

    return(r)
}
