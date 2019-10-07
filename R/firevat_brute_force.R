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
#' @param ga.preemptive.killing If TRUE, then preemptively kills populations that yield greater sequencing artifact weights sum
#' compared to the original mutatational signatures analysis
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
                                    ga.preemptive.killing,
                                    verbose = TRUE) {
    vcf.filter <- MakeFilter(config.obj = config.obj)
    vcf.filter.params <- names(vcf.filter)

    lower.vector <- lower.upper.list$lower.vector
    upper.vector <- lower.upper.list$upper.vector

    # Turn off all filters
    param.values.list <- list()
    for (curr.filter.param.name in vcf.filter.params) {
        curr.filter.param.direction <- config.obj[[curr.filter.param.name]]$direction
        if (curr.filter.param.direction == "POS") {
            param.values.list[[curr.filter.param.name]] <- lower.vector[[curr.filter.param.name]]
        } else { # NEG
            param.values.list[[curr.filter.param.name]] <- upper.vector[[curr.filter.param.name]]
        }
    }

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
                      ga.preemptive.killing = ga.preemptive.killing,
                      verbose = verbose)

    df.suggested.solutions <- data.frame()
    for (curr.filter.param in vcf.filter.params) {
        curr.filter.param.min <- as.numeric(lower.vector[names(lower.vector) == curr.filter.param])
        curr.filter.param.max <- as.numeric(upper.vector[names(upper.vector) == curr.filter.param])

        # Skip if min or max is either -Inf or Inf
        if (curr.filter.param.min == -Inf || curr.filter.param.min == Inf) {
            next
        }

        if (curr.filter.param.max == -Inf || curr.filter.param.max == Inf) {
            next
        }

        PrintLog(paste0("* Iterating through the filter parameter: ", curr.filter.param, " from ", curr.filter.param.min, " to ", curr.filter.param.max))

        df.suggested.solutions.temp <- data.frame()
        for (curr.filter.param.val in seq(curr.filter.param.min, curr.filter.param.max, by = 1)) {
            param.values.list.temp <- param.values.list
            param.values.list.temp[[curr.filter.param]] <- curr.filter.param.val
            results.temp <- GAOptimizationObjFnHelper(params.x = as.numeric(unlist(param.values.list.temp)),
                                                      data = data.temp)

            if (results.temp$valid == FALSE) {
                next
            }

            # Prepare log data
            log.list <- list()
            log.list[['iteration']] <- NA
            log.list[['iteration.ran.successfully']] <- TRUE
            log.list[['datetime']] <- Sys.time()
            log.list[['original.mutations.count']] <- as.numeric(results.temp$original.muts.count)
            log.list[['refined.mutations.count']] <- as.numeric(results.temp$refined.muts.count)
            log.list[['artifact.mutations.count']] <- as.numeric(results.temp$artifactual.muts.count)
            log.list[['refined.muts.cosine.similarity.score']] <- as.numeric(results.temp$refined.muts.cos.sim)
            log.list[['refined.muts.sequencing.artifact.signatures']] <- paste(as.character(results.temp$refined.muts.seq.art.sigs), collapse = ",")
            log.list[['refined.muts.sequencing.artifact.signatures.weights']] <- paste(as.character(results.temp$refined.muts.seq.art.sigs.weights), collapse = ",")
            log.list[['refined.muts.sequencing.artifact.signatures.weights.sum']] <- as.numeric(results.temp$refined.muts.seq.art.weights.sum)
            log.list[['refined.muts.proportion']] <- as.numeric(results.temp$refined.muts.proportion)
            log.list[['refined.muts.target.signatures']] <- paste(as.character(results.temp$refined.muts.target.sigs), collapse = ",")
            log.list[['refined.muts.target.signatures.weights']] <- paste(as.character(results.temp$refined.muts.target.sigs.weights), collapse = ",")
            log.list[['artifactual.muts.cosine.similarity.score']] <- as.numeric(results.temp$artifactual.muts.cos.sim)
            log.list[['artifactual.muts.sequencing.artifact.signatures']] <- paste(as.character(results.temp$artifactual.muts.seq.art.sigs), collapse = ",")
            log.list[['artifactual.muts.sequencing.artifact.signatures.weights']] <- paste(as.character(results.temp$artifactual.muts.seq.art.sigs.weights), collapse = ",")
            log.list[['artifactual.muts.sequencing.artifact.signatures.weights.sum']] <- as.numeric(results.temp$artifactual.muts.seq.art.weights.sum)
            log.list[['artifactual.muts.proportion']] <- as.numeric(results.temp$artifactual.muts.proportion)
            log.list[['artifactual.muts.target.signatures']] <- paste(as.character(results.temp$artifactual.muts.target.sigs), collapse = ",")
            log.list[['artifactual.muts.target.signatures.weights']] <- paste(as.character(results.temp$artifactual.muts.target.sigs.weights), collapse = ",")
            log.list[['objective.value']] <- as.numeric(results.temp$objective.val)

            # Get filter values
            log.list <- append(log.list, param.values.list.temp)

            # Append data to df.suggested.solutions.temp
            df.results.temp <- data.frame(log.list,
                                          stringsAsFactors = F)
            df.suggested.solutions.temp <- rbind(df.suggested.solutions.temp, df.results.temp)
        }

        # Take the row with the highest objective value
        df.suggested.solutions.temp <- df.suggested.solutions.temp[df.suggested.solutions.temp$objective.value > 0,]
        if (nrow(df.suggested.solutions.temp) == 0) {
            # If all zero then skip
            next
        } else {
            df.suggested.solutions.temp <- df.suggested.solutions.temp[which.max(df.suggested.solutions.temp$objective.value),]
            df.suggested.solutions <- rbind(df.suggested.solutions, df.suggested.solutions.temp)
        }
    }

    # Assign iteration based on increasing objective value
    df.suggested.solutions <- df.suggested.solutions[order(+df.suggested.solutions$objective.value),]
    df.suggested.solutions$iteration <- (((-1) * nrow(df.suggested.solutions)) + 1):0

    # Get a suggested solutions matrix with only the filter parameters
    suggested.solutions.matrix <- as.matrix.data.frame(df.suggested.solutions[, vcf.filter.params])

    return(list(df.suggested.solutions = df.suggested.solutions,
                suggested.solutions.matrix = suggested.solutions.matrix))
}
