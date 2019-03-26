# FIREVAT Strand Bias Analysis Functions
#
# Last revised date:
#   February 22, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title PerformStrandBiasAnalysis
#' @description Performs strand bias analysis
#'
#' @param vcf.obj \code{\link{ReadVCF}}
#' @param ref.forward.strand.var A string value.
#' @param ref.reverse.strand.var A string value.
#' @param alt.forward.strand.var A string value.
#' @param alt.reverse.strand.var A string value.
#' @param perform.fdr.correction If TRUE, then performs false discovery rate correction
#' @param fdr.correction.method A string value. FDR correction method (Refer to p.adjust() function)
#'
#' @return An updated vcf.obj
#'
#' @export
PerformStrandBiasAnalysis <- function(vcf.obj,
                                      ref.forward.strand.var,
                                      ref.reverse.strand.var,
                                      alt.forward.strand.var,
                                      alt.reverse.strand.var,
                                      perform.fdr.correction = TRUE,
                                      fdr.correction.method = "BH") {
    # Prepare data
    df <- data.frame(list("ref.forward" = vcf.obj$data[[ref.forward.strand.var]],
                          "ref.reverse" = vcf.obj$data[[ref.reverse.strand.var]],
                          "alt.forward" = vcf.obj$data[[alt.forward.strand.var]],
                          "alt.reverse" = vcf.obj$data[[alt.reverse.strand.var]]),
                     stringsAsFactors = F)

    # Perform Fisher's Exact Test
    Perform.Fisher.Exact.Test <- function(row) {
        test.mat <- rbind(c(row['ref.forward'], row['ref.reverse']),
                          c(row['alt.forward'], row['alt.reverse']))
        test.results <- fisher.test(test.mat)
        return(test.results$p.value)
    }
    fisher.exact.test.p.values <- apply(X = df[, colnames(df)], MARGIN = 1, FUN = Perform.Fisher.Exact.Test)
    vcf.obj$data$StrandBiasPValue <- fisher.exact.test.p.values

    if (perform.fdr.correction == TRUE) {
        fisher.exact.test.q.values <- p.adjust(p = fisher.exact.test.p.values,
                                               method = fdr.correction.method)
        vcf.obj$data$StrandBiasQValue <- fisher.exact.test.q.values
    }

    return(vcf.obj)
}
