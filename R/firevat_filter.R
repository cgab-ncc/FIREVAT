# FIREVAT Filtering Functions
#
# Last revised date:
#   March 22, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title MakeFilter
#' @description Creates a vcf filter from config.obj
#'
#' @param config.obj A list from ParseConfigFile
#' (any filter with "use_in_filter" value declared as FALSE is not considered)
#'
#' @return A list with the filter parameters
#'
#' @export
MakeFilter <- function(config.obj) {
    vcf.config.filter <- list()

    # If default tag is given in config file, save default value into filter list
    for (config.entry in config.obj) {
        # Any filter with FALSE "use_in_filter" value is not considered
        if (config.entry$use_in_filter == TRUE) {
            if ("default" %in% names(config.entry)) {
                vcf.config.filter[[config.entry$name]] <- config.entry$default
            } else {
                # Set default value according to the optimization direction
                # of the parameter: POS=0 / NEG=Inf
                if (config.entry$direction == "POS") {
                    vcf.config.filter[[config.entry$name]] <- 0
                } else if (config.entry$direction == "NEG") {
                    vcf.config.filter[[config.entry$name]] <- Inf
                } else {
                    # If direction is not in c("POS", "NEG"), invoke error
                    stop(paste("Unexpected direction of ", config.entry$name))
                }
            }
        }
    }
    return(vcf.config.filter)
}


#' @title UpdateFilter
#' @description Update filter based on optim parameter values
#'
#' @param vcf.filter A list from MakeFilterFromConfig
#' @param param.values A numeric vector contains filtering value
#' (same length with length(vcf.config.filter))
#'
#' @return Updated vcf.filter (list)
#'
#' @export
UpdateFilter <- function(vcf.filter, param.values) {
    vcf.filter.updated <- vcf.filter

    # Set the values in vcf.filter.updated with param.values
    for (i in 1:length(vcf.filter.updated)) {
        vcf.filter.updated[[i]] <- param.values[i]
    }

    return(vcf.filter.updated)
}


#' @title FilterVCF
#' @description
#' Filter vcf based on the filter
#' Filtering parameters are saved in config.obj
#' Split vcf.obj into vcf.obj.filtered & vcf.obj.artifact based on vcf.filter
#'
#' @param vcf.obj A list from ReadVCF
#' @param vcf.filter A list from MakeMuTect2Filter
#' @param config.obj A list from ParseConfigFile
#' @param include.array A boolean vector
#' @param force.include A boolean value. If TRUE, then uses 'include.array'
#' @param verbose If true, provides process detail
#'
#' @return A list with the following elements
#' \itemize{
#'   \item{1) Mutations which passed filtering}{vcf.obj.filtered = vcf.obj (list with data, header, genome)}
#'   \item{2) Mutations which did not pass filtering}{vcf.obj.artifact = vcf.obj (list with data, header, genome)}
#' }
#'
#' @export
FilterVCF <- function(vcf.obj,
                      vcf.filter,
                      config.obj,
                      include.array = NULL,
                      force.include = FALSE,
                      verbose = TRUE) {
    if (force.include == FALSE) {
        # Set up variables
        include <- rep(TRUE, nrow(vcf.obj$data))
        atcg.chars <- c("A", "T", "C", "G")

        condition.list <- list("include"=include)

        if (verbose == TRUE) {
            PrintLog(paste0("* Before applying filter: ", nrow(vcf.obj$data), " rows in VCF object"))
        }

        for (param in names(vcf.filter)) {
            # Get filtering direction from config.obj
            direction <- config.obj[[param]]["direction"]

            # "POS": pass values bigger than cutoff
            # "NEG": pass values smaller than cutoff
            if (direction == "POS") {
                condition.list[[param]] <- vcf.obj$data[[param]] >= vcf.filter[[param]]
            } else if (direction == "NEG") {
                condition.list[[param]] <- vcf.obj$data[[param]] <= vcf.filter[[param]]
            }
        }

        include <- AND.multiple(condition.list)

        # Split vcf.data with "INCLUDE"
        vcf.data.refined <- vcf.obj$data[include, ]
        vcf.data.artifact <- vcf.obj$data[!include, ]

        if (verbose == TRUE) {
            PrintLog("* After applying filter: ")
            PrintLog(paste0("** ", nrow(vcf.data.refined), " rows in vcf.data.filtered VCF object"))
            PrintLog(paste0("** ", nrow(vcf.data.artifact), " rows in vcf.data.artifact VCF object"))
        }

        # Return two vcf.obj
        return(list(vcf.obj.filtered = list(data = vcf.data.refined,
                                            header = vcf.obj$header,
                                            genome = vcf.obj$genome),
                    vcf.obj.artifact = list(data = vcf.data.artifact,
                                            header = vcf.obj$header,
                                            genome = vcf.obj$genome)))
    }
    else {
        if (is.null(include.array)) {
            stop("The parameter 'include.array' must not be NULL as force.include is TRUE.")
        }
        if (nrow(vcf.obj$data) != length(include.array)) {
            stop("The parameter 'include.array' length must be the same as the number of rows in vcf.obj$data")
        }
        vcf.data.include <- vcf.obj$data[include.array, ]
        vcf.data.exclude <- vcf.obj$data[!include.array, ]

        # Return two vcf.obj
        return(list(vcf.obj.filtered = list(data = vcf.data.include,
                                            header = vcf.obj$header,
                                            genome = vcf.obj$genome),
                    vcf.obj.artifact = list(data = vcf.data.exclude,
                                            header = vcf.obj$header,
                                            genome = vcf.obj$genome)))
    }
}


#' @title Transform default filtering parameters to a binary vector
#' @description
#' This function transforms default filtering parameter to binary vector
#' which can be used as a suggested solution in GA algorithm.
#'
#' @param vcf.filter A list generated in \code{\link{MakeFilter}}
#' @param params.bit.len A list with bit lengths of filtering parameters which is generated from \code{\link{ParameterToBits}}
#'
#' @return A binary vector
#'
#' @export
DefaultFilterToBinary <- function(vcf.filter, params.bit.len) {
    default.filter.decimal <- rep(0, length(vcf.filter))
    default.filter.binary <- c()

    for (i in 1:length(vcf.filter)) {
        param <- names(vcf.filter)[i]

        # max.binary = maximum value of parameter calculated from bit length.
        # max.binary is formatted in decimal format.
        max.binary <- 2^params.bit.len[[param]] - 1

        # Compare default filtering value with max.binary
        # If default filtering value > max.binary, return max.binary
        # Else, return default value in vcf.filter
        if (vcf.filter[[param]] > max.binary) {
            default.filter.decimal[i] <- max.binary
        } else {
            default.filter.decimal[i] <- vcf.filter[[param]]
        }
    }

    # Generate default.filter.binary with decimal values & bit length
    for (j in 1:length(default.filter.decimal)) {
        param.binary <- decimal2binary(default.filter.decimal[j], params.bit.len[j])
        default.filter.binary <- c(default.filter.binary, param.binary)
    }

    return(default.filter.binary)
}


#' @title AND.multiple
#' @description
#' An internal function for handling multiple filtering conditions.
#' Reduce conditions with "&" operations.
#'
#' @param ... Multiple filtering conditions
#'
#' @return A single condition reduced
#'
#' @keywords internal
AND.multiple <- function (...) {
    Reduce("&", ...)
}


#' @title OR.multiple
#' @description
#' An internal function for handling multiple filtering conditions.
#' Reduce conditions with "|" operations.
#'
#' @param ... Multiple filtering conditions
#'
#' @return A single condition reduced
#'
#' @keywords internal
OR.multiple <- function (...) {
    Reduce("|", ...)
}

