# FIREVAT Optimization Functions
#
# Last revised date:
#   March 22, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.k)rg
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title CheckIfVariantRefinementIsNecessary
#' @description
#' Checks if variant refinement is necessary by identifying mutational signatures
#' related to sequencing artifact in the vcf.obj (set of original unrefined point mutations).
#'
#' @param vcf.obj A list from ReadVCF
#' @param bsg BSgenome.Hsapiens.UCSC object
#' @param df.mut.pat.ref.sigs A data.frame from MutPatParseRefMutSigs
#' @param target.mut.sigs A character vector of target mutational signatures from reference mutational signatures.
#' @param sequencing.artifact.mut.sigs A character vector of sequencing artifact mutational signatures from reference mutational signatures.
#' @param init.artifact.stop
#' Numeric value less than 1. If the sum of sequencing artifact weights in vcf.obj is less than or equal to this value then
#' this function returns judgment = FALSE, otherwise returns judgment = TRUE.
#' @param verbose If TRUE, provides process detail. Default value is TRUE.
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{judgment}{A boolean value}
#'  \item{seq.art.sigs.weights.sum}{A numeric value. Sum of sequencing artifact weights.}
#' }
#'
#' @export
CheckIfVariantRefinementIsNecessary <- function(vcf.obj,
                                                bsg,
                                                df.mut.pat.ref.sigs,
                                                target.mut.sigs,
                                                sequencing.artifact.mut.sigs,
                                                init.artifact.stop = 0.05,
                                                verbose = TRUE) {
    mut.pat.input <- MutPatParseVCFObj(vcf.obj, bsg = bsg)
    mut.pat.results <- RunMutPat(mut.pat.input = mut.pat.input,
                                 df.mut.pat.ref.sigs = df.mut.pat.ref.sigs,
                                 target.mut.sigs = unique(c(target.mut.sigs, sequencing.artifact.mut.sigs)),
                                 verbose = verbose)
    df.sigs <- data.frame(sig = mut.pat.results$identified.mut.sigs,
                          weight = mut.pat.results$identified.mut.sigs.contribution.weights,
                          stringsAsFactors = F)
    df.sigs.seq.art <- df.sigs[(df.sigs$sig %in% sequencing.artifact.mut.sigs), ]
    seq.art.sigs.weights <- df.sigs.seq.art$weight
    seq.art.sigs.weights <- seq.art.sigs.weights / sum(mut.pat.results$identified.mut.sigs.contribution.weights)
    seq.art.sigs.weights.sum <- sum(seq.art.sigs.weights)
    if (seq.art.sigs.weights.sum <= init.artifact.stop) {
        return(list(judgment = FALSE,
                    seq.art.sigs.weights.sum = seq.art.sigs.weights.sum))
    } else {
        return(list(judgment = TRUE,
                    seq.art.sigs.weights.sum = seq.art.sigs.weights.sum))
    }
}


#' @title ParameterToBits
#' @description Calculate the number of bits needed to conduct FIREVAT GA binary optimization.
#'
#' @param vcf.obj A list from ReadVCF
#' @param vcf.filter A list from MakeMuTect2Filter
#' @param config.obj A list from ParseConfigFile
#' @param multiplier A multiplier for convert fraction to integer (default = 100)
#'
#' @details
#' vcf.obj$data: if max(vcf.obj$data[[param]]) < 1, then multiply multiplier to the vector
#'
#' @return A list with the elements
#' \itemize{
#'  \item{params.bit.len}{A numeric vector. Each element is the bit length of each parameter value}
#'  \item{vcf.obj}{A vcf.obj (\code{\link{ReadVCF}}) with updated data}
#' }
#'
#' @export
#' @importFrom GA decimal2binary
ParameterToBits <- function(vcf.obj, config.obj, vcf.filter, multiplier = 100) {
    params.bit.len <- rep(0, length(vcf.filter))
    for (i in 1:length(vcf.filter)) {
        param <- names(vcf.filter)[i]

        # If custom range is given in config.obj,
        # calculate bit length with the given range.
        if ("range"%in% names(config.obj[[param]])) {
            max.range <- unlist(config.obj[[param]]["range"])
            params.bit.len[i] <- length(decimal2binary(max.range[[2]]))
        } else {
            if (max(vcf.obj$data[[param]]) > 1) {
                params.bit.len[i] <- length(decimal2binary(max(vcf.obj$data[[param]])))
            } else {
                # If max(vcf.obj$data[[param]]) <= 1,
                # then multiply multiplier to value
                params.bit.len[i] <- length(decimal2binary(max(multiplier * vcf.obj$data[[param]])))
                vcf.obj$data[[param]] <- multiplier * vcf.obj$data[[param]]
            }
        }
    }
    names(params.bit.len) <- names(vcf.filter)
    return(list(params.bit.len = params.bit.len, vcf.obj = vcf.obj))
}

#' @title GetParameterLowerUpperVector
#' @description Return a lower/upper vector needed to conduct FIREVAT GA real-valued optimization.
#'
#' @param vcf.obj A list from ReadVCF
#' @param vcf.filter A list from MakeMuTect2Filter
#' @param config.obj A list from ParseConfigFile
#' @param multiplier A multiplier for convert fraction to integer (default = 100)
#'
#' @details
#' vcf.obj$data: if max(vcf.obj$data[[param]]) < 1, then multiply multiplier to the vector
#'
#' @return A list with the elements
#' \itemize{
#'  \item{lower.vector}{ A numeric vector. Each element is the minimum value of each parameter}
#'  \item{upper.vector}{ A numeric vector. Each element is the maximum value of each parameter}
#'  \item{vcf.obj}{ vcf.obj with updated data}
#' }
#'
#' @importFrom IRanges IRanges intersect
#' @export
GetParameterLowerUpperVector <- function(vcf.obj, config.obj, vcf.filter, multiplier = 100) {
    # Make empty vectors
    lower.vector <- rep(0, length(vcf.filter))
    upper.vector <- rep(0, length(vcf.filter))

    for (i in 1:length(vcf.filter)) {
        param <- names(vcf.filter)[i]
        # If custom range is given in config.obj,
        # return lower/upper vectors with the given range.
        if ("range" %in% names(config.obj[[param]])) {
            param.range <- unlist(config.obj[[param]]["range"])
            if (max(vcf.obj$data[[param]]) <= 1) {
                vcf.obj$data[[param]] <- multiplier * vcf.obj$data[[param]]
            }
            user.specified.range <- IRanges(start = c(param.range[[1]]), end = c(param.range[[2]]))

            data.min <- min(vcf.obj$data[[param]])
            data.max <- max(vcf.obj$data[[param]])
            if (data.min == -Inf) {
                data.min <- param.range[[1]]
            }
            if (data.max == Inf) {
                data.max <- param.range[[2]]
            }
            data.range <- IRanges(start = c(data.min), end = c(data.max))

            range.intersection <- intersect(user.specified.range, data.range)
            if (length(range.intersection) > 0) {
                lower.vector[i] <- range.intersection@start
                upper.vector[i] <- (range.intersection@start + range.intersection@width) - 1
            } else {
                print(paste0("The range for the filter parameter ", param, " is invalid.\n",
                             "The actual data do not reflect this range.\n",
                             "config file min = ", param.range[[1]], "\n",
                             "vcf file min = ", min(vcf.obj$data[[param]]), "\n",
                             "config file max = ", param.range[[2]], "\n",
                             "vcf file max = ", max(vcf.obj$data[[param]]), "\n",
                             "Returning prematurely."))
                return(list(valid = FALSE))
            }
        } else {
            if (max(vcf.obj$data[[param]]) <= 1) {
                vcf.obj$data[[param]] <- multiplier * vcf.obj$data[[param]]

            }
            lower.vector[i] <- min(vcf.obj$data[[param]])
            upper.vector[i] <- max(vcf.obj$data[[param]])
        }
    }
    names(lower.vector) <- names(vcf.filter)
    names(upper.vector) <- names(vcf.filter)
    return(list(valid = TRUE,
                lower.vector = lower.vector,
                upper.vector = upper.vector, vcf.obj = vcf.obj))
}


#' @title GADecodeBinaryString
#' @description This function returns decimal vector decoded from binary string.
#'
#' @param string A binary expression of parameters
#' @param data A list which contains FIREVAT data
#'
#' @return A numeric vector corresponding to decoded string representation
#' of parameters
#'
#' @keywords internal
#' @importFrom GA binary2decimal gray2binary
GADecodeBinaryString <- function(string, data) {
    params.bit.len <- data$params.bit.len
    string <- gray2binary(string)
    decimal.vector <- rep(0, length(params.bit.len))
    names(decimal.vector) <- names(params.bit.len)

    start <- 1
    end <- params.bit.len[1]
    for (i in 1:length(params.bit.len)) {
        decimal.vector[i] <- binary2decimal(string[start:end])
        start <- start + params.bit.len[i]
        if (i != length(params.bit.len)) {
            end <- start + params.bit.len[i+1]-1
        }
    }
    return(decimal.vector)
}


#' @title GAOptimizationObjFnHelper
#' @description Objective helper function for FIREVAT optimization
#'
#' @param params.x = representation of parameter values
#' @param data A list with the following data
#' \itemize{
#'  \item{"vcf.obj"}{\code{\link{ReadVCF}}}
#'  \item{"vcf.filter"}{a list}
#'  \item{"vcf.calling.method"}{string value}
#'  \item{"df.ref.mut.sig"}{a data.frame}
#'  \item{"target.mut.sigs"}{a character vector}
#'  \item{"sequencing.artifact.mut.sigs"}{a charcter vector}
#'  \item{"config.obj"}{\code{\link{ParseConfigFile}}}
#' }
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{C.filtered}{numeric value}
#'  \item{A.filtered}{numeric value}
#'  \item{P.filtered}{numeric value}
#'  \item{C.artifact}{numeric value}
#'  \item{A.artifact}{numeric value}
#'  \item{P.artifact}{numeric value}
#'  \item{obj.val}{obj.val}
#'  \item{x}{x (numeric vector of GA parameter values)}
#'  \item{original mutations}{numeric value}
#'  \item{filtered.mutations}{numeric value}
#'  \item{artifact.mutations}{numeric value}
#'  \item{sequencing.artifact.signatures}{character vector (of signature names)}
#'  \item{sequencing.artifact.signatures.weights}{numeric vector}
#' }
#'
#' @keywords internal
GAOptimizationObjFnHelper <- function(params.x, data) {
    # Check ga.type and convert params.x
    if (data$ga.type == "binary") {
        # Convert binary string representation of filter parameters to integers
        x <- GADecodeBinaryString(string = params.x, data = data)
    } else if (data$ga.type == "real-valued") {
        # Convert numeric to integer
        x <- floor(params.x)
    }

    # Update filter
    vcf.filter <- UpdateFilter(vcf.filter = data$vcf.filter, param.values = x)

    # Filter vcf.data based on the updated vcf.filter
    vcf.objs <- FilterVCF(vcf.obj = data$vcf.obj, vcf.filter = vcf.filter,
                          config.obj = data$config.obj, verbose = FALSE)

    # Return if 50 or fewer point mutations left
    if ((nrow(vcf.objs$vcf.obj.filtered$data) <= 50) || (nrow(vcf.objs$vcf.obj.artifact$data) <= 50)) {
        return(list(valid = FALSE,
                    x = x,
                    original.muts.count = nrow(data$vcf.obj$data),
                    refined.mutations.count = nrow(vcf.objs$vcf.obj.filtered$data),
                    artifactual.muts.count = nrow(vcf.objs$vcf.obj.artifact$data),
                    objective.val = 0))
    }

    # Identify mutational signatures
    # C.value.refined
    # A.value.refined
    # P.value.refined
    refined.mut.pat.input <- MutPatParseVCFObj(vcf.obj = vcf.objs$vcf.obj.filtered,
                                               bsg = data$bsg)
    refined.muts.mut.pat.results <- RunMutPat(mut.pat.input = refined.mut.pat.input,
                                              df.mut.pat.ref.sigs = data$df.mut.pat.ref.sigs,
                                              target.mut.sigs = unique(c(data$target.mut.sigs, data$sequencing.artifact.mut.sigs)),
                                              verbose = FALSE)

    # Extract Mutational Patterns data
    df.refined.sigs <- data.frame(sig = refined.muts.mut.pat.results$identified.mut.sigs,
                                  weight = refined.muts.mut.pat.results$identified.mut.sigs.contribution.weights,
                                  stringsAsFactors = F)
    df.refined.sigs.seq.art <- df.refined.sigs[(df.refined.sigs$sig %in% data$sequencing.artifact.mut.sigs), ]
    refined.seq.art.sigs <- df.refined.sigs.seq.art$sig
    refined.seq.art.sigs.weights <- df.refined.sigs.seq.art$weight
    df.refined.sigs.target <- df.refined.sigs[!(df.refined.sigs$sig %in% data$sequencing.artifact.mut.sigs), ]
    refined.target.sigs <- df.refined.sigs.target$sig
    refined.target.sigs.weights <- df.refined.sigs.target$weight

    # Compute desired variables
    refined.seq.art.sigs.weights <- refined.seq.art.sigs.weights / sum(refined.muts.mut.pat.results$identified.mut.sigs.contribution.weights)
    refined.target.sigs.weights <- refined.target.sigs.weights / sum(refined.muts.mut.pat.results$identified.mut.sigs.contribution.weights)
    C.value.refined <- refined.muts.mut.pat.results$cosine.similarity.score
    A.value.refined <- sum(refined.seq.art.sigs.weights)
    P.value.refined <- nrow(vcf.objs$vcf.obj.filtered$data) / nrow(data$vcf.obj$data)

    # Identify mutational signatures
    # C.value.artifactual
    artifact.mut.pat.input <- MutPatParseVCFObj(vcf.obj = vcf.objs$vcf.obj.artifact,
                                                bsg = data$bsg)
    artifactual.muts.mut.pat.results <- RunMutPat(mut.pat.input = artifact.mut.pat.input,
                                                  df.mut.pat.ref.sigs = data$df.mut.pat.ref.sigs,
                                                  target.mut.sigs = data$sequencing.artifact.mut.sigs,
                                                  verbose = FALSE)
    C.value.artifactual <- artifactual.muts.mut.pat.results$cosine.similarity.score

    # Identify mutational signatures
    # A.value.artifactual
    # P.value.artifactual
    artifactual.muts.mut.pat.results <- RunMutPat(mut.pat.input = artifact.mut.pat.input,
                                                  df.mut.pat.ref.sigs = data$df.mut.pat.ref.sigs,
                                                  target.mut.sigs = unique(c(data$target.mut.sigs, data$sequencing.artifact.mut.sigs)),
                                                  verbose = FALSE)

    # Extract Mutational Patterns data
    df.artifact.sigs <- data.frame(list(sig = as.character(artifactual.muts.mut.pat.results$identified.mut.sigs),
                                        weight = as.numeric(artifactual.muts.mut.pat.results$identified.mut.sigs.contribution.weights)),
                                   stringsAsFactors = F)
    df.artifact.sigs.seq.art <- df.artifact.sigs[(df.artifact.sigs$sig %in% data$sequencing.artifact.mut.sigs), ]
    artifactual.muts.seq.art.sigs <- df.artifact.sigs.seq.art$sig
    artifactual.muts.seq.art.sigs.weights <- df.artifact.sigs.seq.art$weight
    df.artifact.sigs.target <- df.artifact.sigs[!(df.artifact.sigs$sig %in% data$sequencing.artifact.mut.sigs), ]
    artifactual.target.sigs <- df.artifact.sigs.target$sig
    artifactual.target.sigs.weights <- df.artifact.sigs.target$weight

    # Compute desired variables
    artifactual.muts.seq.art.sigs.weights <- artifactual.muts.seq.art.sigs.weights / sum(artifactual.muts.mut.pat.results$identified.mut.sigs.contribution.weights)
    artifactual.target.sigs.weights <- artifactual.target.sigs.weights / sum(artifactual.muts.mut.pat.results$identified.mut.sigs.contribution.weights)
    A.value.artifactual <- sum(artifactual.muts.seq.art.sigs.weights)
    P.value.artifactual <- nrow(vcf.objs$vcf.obj.artifact$data) / nrow(data$vcf.obj$data)

    # Objective value
    obj.val <- data$objective.fn(C.refined = C.value.refined,
                                 A.refined = A.value.refined,
                                 C.artifactual = C.value.artifactual,
                                 A.artifactual = A.value.artifactual)

    # A.value.refined must be lower than the sum of sequencing artifact weights in the original vcf file
    if (data$ga.preemptive.killing == TRUE) {
        if (A.value.refined > data$original.muts.seq.art.weights.sum) {
            obj.val <- 0
        }
    }

    return(list(valid = TRUE,
                x = x,
                objective.val = obj.val,
                refined.muts.cos.sim = C.value.refined,
                refined.muts.seq.art.weights.sum = A.value.refined,
                refined.muts.proportion = P.value.refined,
                artifactual.muts.cos.sim = C.value.artifactual,
                artifactual.muts.seq.art.weights.sum = A.value.artifactual,
                artifactual.muts.proportion = P.value.artifactual,
                original.muts.count = nrow(data$vcf.obj$data),
                refined.muts.count = nrow(vcf.objs$vcf.obj.filtered$data),
                artifactual.muts.count = nrow(vcf.objs$vcf.obj.artifact$data),
                refined.muts.seq.art.sigs = refined.seq.art.sigs,
                refined.muts.seq.art.sigs.weights = refined.seq.art.sigs.weights,
                refined.muts.target.sigs = refined.target.sigs,
                refined.muts.target.sigs.weights = refined.target.sigs.weights,
                artifactual.muts.seq.art.sigs = artifactual.muts.seq.art.sigs,
                artifactual.muts.seq.art.sigs.weights = artifactual.muts.seq.art.sigs.weights,
                artifactual.muts.target.sigs = artifactual.target.sigs,
                artifactual.muts.target.sigs.weights = artifactual.target.sigs.weights))
}


#' @title GAOptimizationObjFn
#' @description Objective function for FIREVAT optimization
#'
#' @param string = binary representation of parameter values
#' @param data = A list with the following data
#' \itemize{
#'  \item{"vcf.obj"}{\code{\link{ReadVCF}}}
#'  \item{"vcf.filter"}{a list}
#'  \item{"vcf.calling.method"}{string value}
#'  \item{"df.ref.mut.sig"}{a data.frame}
#'  \item{"target.mut.sigs"}{a character vector}
#'  \item{"sequencing.artifact.mut.sigs"}{a charcter vector}
#'  \item{"config.obj"}{\code{\link{ParseConfigFile}}}
#' }
#'
#' @return A numeric value corresponding to the optimization objective value
#'
#' @keywords internal
GAOptimizationObjFn <- function(params.x, data) {
    results <- GAOptimizationObjFnHelper(params.x = params.x, data = data)
    return(results$objective.val)
}


#' @title GAMonitorFn
#' @description Prints GA optimization iteration to console
#'
#' @param obj = A factor returned by GA package
#' @param data A list with the following data
#' \itemize{
#'  \item{"vcf.obj"}{\code{\link{ReadVCF}}}
#'  \item{"vcf.filter"}{a list}
#'  \item{"vcf.calling.method"}{string value}
#'  \item{"df.ref.mut.sig"}{a data.frame}
#'  \item{"target.mut.sigs"}{a character vector}
#'  \item{"sequencing.artifact.mut.sigs"}{a charcter vector}
#'  \item{"config.obj"}{\code{\link{ParseConfigFile}}}
#' }
#'
#' @keywords internal
#' @importFrom utils write.table
GAMonitorFn <- function(obj, data) {
    max.fitness <- max(obj@fitness)
    max.fitness.index <- which(obj@fitness == max.fitness)
    max.fitness.pop <- unique(as.data.frame(obj@population)[max.fitness.index,])
    results <- GAOptimizationObjFnHelper(params.x = as.numeric(max.fitness.pop[1,]), data = data)

    if (data$verbose == TRUE) {
        cat("\n")
        print(sprintf("%-50s%s", "Datetime", Sys.time()))
        print(sprintf("%-50s%i", "Iteration", obj@iter))
        print(sprintf("%-50s%f", "Max fitness", max.fitness))
        print(sprintf("%-50s%i", "Max fitness population counts", nrow(max.fitness.pop)))
    }

    for (i in 1:length(data$vcf.filter)) {
        param <- names(data$vcf.filter)[i]
        direction <- data$config.obj[[param]]["direction"]
        if (direction == "POS") {
            inequal.sign <- "(>=)"
        } else if (direction == "NEG") {
            inequal.sign <- "(<=)"
        }

        if (data$verbose == TRUE) {
            print(sprintf("%-50s%s", paste0(names(data$vcf.filter)[i], " ", inequal.sign), results$x[i]))
        }
    }

    if (results$valid) {
        if (data$verbose == TRUE) {
            print(sprintf("%-50s%f", "Refined mutations cosine sim score", results$refined.muts.cos.sim))
            print(sprintf("%-50s%f", "Refined mutations seq art sigs weights sum", results$refined.muts.seq.art.weights.sum))
            print(sprintf("%-50s%f", "Refined mutations proportion", results$refined.muts.proportion))
            print(sprintf("%-50s%f", "Artifactual mutations cosine sim score ", results$artifactual.muts.cos.sim))
            print(sprintf("%-50s%f", "Artifactual mutations seq art sigs weights sum", results$artifactual.muts.seq.art.weights.sum))
            print(sprintf("%-50s%f", "Artifactual mutations proportion", results$artifactual.muts.proportion))
            print(sprintf("%-50s%f", "Objective value", results$objective.val))
        }

        # Prepare data to write to log
        df.log <- data.frame("iteration" = c(obj@iter),
                             "iteration.ran.successfully" = c(results$valid),
                             "datetime" = c(Sys.time()),
                             "original.mutations.count" = c(results$original.muts.count),
                             "refined.mutations.count" = c(results$refined.muts.count),
                             "artifact.mutations.count" = c(results$artifactual.muts.count),
                             "refined.muts.cosine.similarity.score" = c(results$refined.muts.cos.sim),
                             "refined.muts.sequencing.artifact.signatures" = c(paste(as.character(results$refined.muts.seq.art.sigs), collapse = ",")),
                             "refined.muts.sequencing.artifact.signatures.weights" = c(paste(as.character(results$refined.muts.seq.art.sigs.weights), collapse = ",")),
                             "refined.muts.sequencing.artifact.signatures.weights.sum" = c(results$refined.muts.seq.art.weights.sum),
                             "refined.muts.proportion" = c(results$refined.muts.proportion),
                             "refined.muts.target.signatures" = c(paste(as.character(results$refined.muts.target.sigs), collapse = ",")),
                             "refined.muts.target.signatures.weights" = c(paste(as.character(results$refined.muts.target.sigs.weights), collapse = ",")),
                             "artifactual.muts.cosine.similarity.score" = c(results$artifactual.muts.cos.sim),
                             "artifactual.muts.sequencing.artifact.signatures" = c(paste(as.character(results$artifactual.muts.seq.art.sigs), collapse = ",")),
                             "artifactual.muts.sequencing.artifact.signatures.weights" = c(paste(as.character(results$artifactual.muts.seq.art.sigs.weights), collapse = ",")),
                             "artifactual.muts.sequencing.artifact.signatures.weights.sum" = c(results$artifactual.muts.seq.art.weights.sum),
                             "artifactual.muts.proportion" = c(results$artifactual.muts.proportion),
                             "artifactual.muts.target.signatures" = c(paste(as.character(results$artifactual.muts.target.sigs), collapse = ",")),
                             "artifactual.muts.target.signatures.weights" = c(paste(as.character(results$artifactual.muts.target.sigs.weights), collapse = ",")),
                             "objective.value" = c(results$objective.val))

        for (i in 1:length(data$vcf.filter)) {
            df.log[names(data$vcf.filter)[i]] <- c(results$x[i])
        }
    } else {
        if (data$verbose == TRUE) {
            print(sprintf("%-50s", "50 or fewer mutations. Skipping iteration"))
        }

        # Prepare data to write to log
        df.log <- data.frame("iteration" = c(obj@iter),
                             "iteration.ran.successfully" = c(results$valid),
                             "datetime" = c(Sys.time()),
                             "original.mutations.count" = c(results$original.muts.count),
                             "refined.mutations.count" = c(results$refined.mutations.count),
                             "artifact.mutations.count" = c(results$artifactual.muts.count),
                             "refined.muts.cosine.similarity.score" = c(""),
                             "refined.muts.sequencing.artifact.signatures" = c(""),
                             "refined.muts.sequencing.artifact.signatures.weights" = c(""),
                             "refined.muts.sequencing.artifact.signatures.weights.sum" = c(""),
                             "refined.muts.proportion" = c(""),
                             "refined.muts.target.signatures" = c(""),
                             "refined.muts.target.signatures.weights" = c(""),
                             "artifactual.muts.cosine.similarity.score" = c(""),
                             "artifactual.muts.sequencing.artifact.signatures" = c(""),
                             "artifactual.muts.sequencing.artifact.signatures.weights" = c(""),
                             "artifactual.muts.sequencing.artifact.signatures.weights.sum" = c(""),
                             "artifactual.muts.proportion" = c(""),
                             "artifactual.muts.target.signatures" = c(""),
                             "artifactual.muts.target.signatures.weights" = c(""),
                             "objective.value" = c(0))

        for (i in 1:length(data$vcf.filter)) {
            df.log[names(data$vcf.filter)[i]] <- c(results$x[i])
        }
    }

    # Write iteration to log file
    log.file <- paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_Optimization_Logs.tsv")
    if (!file.exists(log.file)) {
        write.table(df.log, log.file, row.names = F, sep = "\t")
    } else {
        write.table(df.log, log.file, row.names = F, sep = "\t", append = T, col.names = F)
    }
}
