# FIREVAT Mutational Patterns Functions
# Mutational Patterns (Blokzijl et al., Genome Medicine 2018; PMID 29695279)
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.k)rg
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title MutPatParseVCFObj
#' @description Parses a vcf.obj and prepares it to run Mutational Patterns.
#'
#' @param vcf.obj A list from ReadVCF
#' @param sample.id A string value
#'
#' @return A data.frame with the column sample.id and
#' row names corresponding to 96 substitution types
#'
#' @keywords internal
#' @importFrom deconstructSigs mut.to.sigs.input
MutPatParseVCFObj <- function(vcf.obj, sample.id = "sample") {
    if (vcf.obj$genome == "hg19") {
        bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    }
    if (vcf.obj$genome == "hg38") {
        bsg <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    }

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
#' @keywords internal
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
#' @param vcf.obj A list from ReadVCF
#' @param df.ref.mut.sigs A data.frame of reference mutational signatures
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
#' @keywords internal
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
RunMutPat <- function(vcf.obj,
                      df.ref.mut.sigs,
                      target.mut.sigs,
                      verbose = TRUE) {
    if (verbose == TRUE) {
        print("Started running Mutational Patterns")
    }

    target.mut.sigs <- unique(target.mut.sigs)
    mut.pat.input <- MutPatParseVCFObj(vcf.obj = vcf.obj)
    df.mut.pat.ref.sigs <- MutPatParseRefMutSigs(df.ref.mut.sigs = df.ref.mut.sigs,
                                                 target.mut.sigs = target.mut.sigs)

    # Mutational Patterns - identify mutational signatures
    mut.pat.results <- fit_to_signatures(mut.pat.input, df.mut.pat.ref.sigs)

    # Wrap Mutational Patterns results
    mut.sigs <- c()
    mut.sigs.contribution.weights <- c()
    for(i in 1:nrow(mut.pat.results$contribution)) {
        curr.signature <- rownames(mut.pat.results$contribution)[i]
        curr.signature.weight <- mut.pat.results$contribution[i,1]
        if (curr.signature.weight > 0) {
            mut.sigs <- c(mut.sigs, curr.signature)
            mut.sigs.contribution.weights <- c(mut.sigs.contribution.weights,
                                               curr.signature.weight)
        }
    }

    # Compute cosine similarity score
    cosine.similarity.score <- lsa::cosine(as.numeric(mut.pat.results$reconstructed),
                                           as.numeric(mut.pat.input))

    r <- list(tumor.mutation.types.spectrum = as.numeric(mut.pat.input),
              identified.mutation.types.spectrum = as.numeric(mut.pat.results$reconstructed),
              residuals = as.numeric(mut.pat.input) - as.numeric(mut.pat.results$reconstructed),
              mutation.types = row.names(mut.pat.results$reconstructed),
              identified.mut.sigs = mut.sigs,
              identified.mut.sigs.contribution.weights = mut.sigs.contribution.weights,
              cosine.similarity.score = cosine.similarity.score)

    if (verbose == TRUE) {
        print("Finished running Mutational Patterns")
    }

    return(r)
}
