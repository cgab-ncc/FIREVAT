# FIREVAT Pareto Optimality Functions
#
# Last revised date:
#   March 22, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.k)rg
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title GetGASuggestedSolutions
#' @description Computes suggested solutions
#'
#' @param vcf.obj A list from ReadVCF
#' @param bsg BSgenome.Hsapiens.UCSC object
#' @param config.obj A list from ParseConfigFile
#' @param lower.upper.list A list from GetParameterLowerUpperVector
#' @param df.mut.pat.ref.sigs A data.frame from MutPatParseRefMutSigs
#' @param target.mut.sigs A character vector of the target mutational signatures from reference mutational signatures.
#' @param sequencing.artifact.mut.sigs A character vector of the sequencing artifact mutational signatures from reference mutational signatures.
#' @param objective.fn Objective value derivation function.
#' @param original.muts.seq.art.weights.sum A numeric value. 'seq.art.sigs.weights.sum' from CheckIfVariantRefinementIsNecessary
#' @param verbose If TRUE, provides process detail. Default value is TRUE.
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{judgment}{A boolean value}
#'  \item{seq.art.sigs.weights}{A numeric value. Sum of sequencing artifact weights.}
#' }
#'
#' @export
GetGASuggestedSolutions <- function(vcf.obj,
                                    bsg,
                                    config.obj,
                                    lower.upper.list,
                                    df.mut.pat.ref.sigs,
                                    target.mut.sigs,
                                    sequencing.artifact.mut.sigs,
                                    objective.fn,
                                    original.muts.seq.art.weights.sum,
                                    verbose = TRUE) {
    vcf.filter <- MakeFilter(config.obj = config.obj)
    vcf.filter.params <- names(vcf.filter)

    lower.vector <- lower.upper.list$lower.vector
    upper.vector <- lower.upper.list$upper.vector

    # Turn off all filters
    param.values <- c()
    for (curr.filter.param.temp in vcf.filter.params) {
        curr.filter.param.temp.direction <- config.obj[[curr.filter.param.temp]]$direction
        if (curr.filter.param.temp.direction == "POS") {
            param.values <- c(param.values, lower.vector[curr.filter.param.temp])
            # param.values <- c(param.values, -Inf)
        } else { # NEG
            param.values <- c(param.values, upper.vector[curr.filter.param.temp])
            # param.values <- c(param.values, Inf)
        }
    }
    names(param.values) <- vcf.filter.params

    # Prepare data
    data.temp <- list(ga.type = "real-valued",
                      vcf.obj = vcf.obj,
                      bsg = bsg,
                      vcf.filter = vcf.filter,
                      config.obj = config.obj,
                      df.mut.pat.ref.sigs = df.mut.pat.ref.sigs,
                      target.mut.sigs = target.mut.sigs,
                      sequencing.artifact.mut.sigs = sequencing.artifact.mut.sigs,
                      original.muts.seq.art.weights.sum = original.muts.seq.art.weights.sum,
                      objective.fn = objective.fn,
                      verbose = verbose)

    df.suggested.solutions <- data.frame()
    for (curr.filter.param in vcf.filter.params) {
        curr.filter.param.min <- as.numeric(lower.vector[names(lower.vector) == curr.filter.param])
        curr.filter.param.max <- as.numeric(upper.vector[names(upper.vector) == curr.filter.param])
        print(paste0("Iterating through the filter parameter: ", curr.filter.param, " from ", curr.filter.param.min, " to ", curr.filter.param.max))

        df.suggested.solutions.temp <- data.frame()
        for (curr.filter.param.val in seq(curr.filter.param.min, curr.filter.param.max, by = 1)) {
            param.values.temp <- param.values
            param.values.temp[curr.filter.param] <- curr.filter.param.val
            results.temp <- GAOptimizationObjFnHelper(params.x = param.values.temp,
                                                      data = data.temp)
            param.values.temp$objective.value <- results.temp$objective.val
            df.results.temp <- data.frame(param.values.temp,
                                          stringsAsFactors = F)
            if (nrow(df.suggested.solutions.temp) == 0) {
                df.suggested.solutions.temp <- df.results.temp
            } else {
                df.suggested.solutions.temp <- rbind(df.suggested.solutions.temp, df.results.temp)
            }
        }
        # Take the row with the highest objective value
        df.suggested.solutions.temp <- df.suggested.solutions.temp[df.suggested.solutions.temp$objective.value > 0,]
        if (nrow(df.suggested.solutions.temp) == 0) {
            next
        }
        df.suggested.solutions.temp <- df.suggested.solutions.temp[which.max(df.suggested.solutions.temp$objective.value),]
        if (nrow(df.suggested.solutions) == 0) {
            df.suggested.solutions <- df.suggested.solutions.temp
        } else {
            df.suggested.solutions <- rbind(df.suggested.solutions, df.suggested.solutions.temp)
        }
    }

    # Get default values
    for (curr.filter.param in vcf.filter.params) {
        if ("default" %in% names(config.obj[[curr.filter.param]])) {
            param.values.temp <- param.values
            param.values.temp[curr.filter.param] <- config.obj[[curr.filter.param]]$default
            param.values.temp$objective.value <- NA
            df.results.temp <- data.frame(param.values.temp,
                                          stringsAsFactors = F)
            if (nrow(df.suggested.solutions) == 0) {
                df.suggested.solutions <- df.results.temp
            } else {
                df.suggested.solutions <- rbind(df.suggested.solutions, df.results.temp)
            }
        }
    }
    suggested.solutions.matrix <- as.matrix.data.frame(df.suggested.solutions[, !(names(df.suggested.solutions) %in% c("objective.value"))])
    return(list(df.suggested.solutions = df.suggested.solutions,
                suggested.solutions.matrix = suggested.solutions.matrix))
}
