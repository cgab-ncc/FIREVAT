# FIREVAT Filtering Functions
#
# Last revised date:
#   February 18, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title MakeFilterFromConfig
#' @description Creates a vcf filter from config.obj
#'
#' @param config.obj A list from ParseConfigFile
#' (any filter with FALSE "use_in_filter" value is not considered)
#'
#' @return A list with the filter parameters
#' @export
MakeFilterFromConfig <- function(config.obj)  {

    vcf.config.filter <- list()
    
    # If default tag is given in config file, save default value into filter list
    for (config.entry in config.obj) {
        # Any filter with FALSE "use_in_filter" value is not considered
        if (config.entry$use_in_filter == TRUE){
            if ("default" %in% names(config.entry)) {
                vcf.config.filter[[config.entry$name]] <- config.entry$default
            } else {
                # Set default value according to the optimization direction
                # of the parameter: POS=0 / NEG=Inf
                if (config.entry$direction == "POS"){
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

#' @title UpdateFilterFromConfig
#' @description Update filter based on optim parameter values
#'
#' @param vcf.filter A list from MakeFilterFromConfig
#' @param param.values A numeric vector contains filtering value 
#'  (same length with length(vcf.config.filter))
#'
#' @return Updated vcf.filter (list)
#' @export
UpdateFilterFromConfig <- function(vcf.filter, param.values){

    vcf.filter.updated <- vcf.filter

    # Set the values in vcf.filter.updated with param.values
    for (i in 1:length(vcf.filter.updated)) {
        vcf.filter.updated[[i]] <- param.values[i]
    }

    return(vcf.filter.updated)
}

#' @title FilterVCFFromConfig
#' @description 
#' Filter vcf based on the filter
#' Filtering parameters are saved in config.obj
#' Split vcf.obj into vcf.obj.filtered & vcf.obj.artifact based on vcf.filter
#'
#' @param vcf.obj A list from ReadVCF
#' @param vcf.filter A list from MakeMuTect2Filter
#' @param config.obj A list from ParseConfigFile
#' @param verbose If true, provides process detail
#' 
#' @return A list with the following elements 
#' \itemize{
#'   \item{1) Mutations which passed filtering}{vcf.obj.filtered = vcf.obj (list with data, header, genome)}
#'   \item{2) Mutations which did not pass filtering}{vcf.obj.artifact = vcf.obj (list with data, header, genome)}
#' }
FilterVCFFromConfig <- function(vcf.obj, vcf.filter, config.obj, verbose = FALSE){

    # Set up variables
    vcf.data.temp <- vcf.obj$data
    include <- rep(TRUE, nrow(vcf.data.temp))
    atcg.chars <- c("A", "T", "C", "G")
    
    condition.list <- list("include"=include)

    if (verbose) { 
        print(paste0("Before applying filter: ", nrow(vcf.data.temp), " rows")) 
    }

    for (param in names(vcf.filter)) {

        # Get filtering direction from config.obj
        direction <- config.obj[[param]]["direction"]

        # "POS": pass values bigger than cutoff 
        # "NEG": pass values smaller than cutoff 
        if (direction == "POS") {
            condition.list[[param]] <- vcf.obj$data[[param]] > vcf.filter[[param]]
        } else if (direction == "NEG") {
            condition.list[[param]] <- vcf.obj$data[[param]] < vcf.filter[[param]]
        }
    }

    include <- AND.multiple(condition.list)

    # Split vcf.data with "INCLUDE"
    vcf.data.temp <- vcf.data.temp[include, ]
    vcf.data.artifact <- vcf.obj$data[!include, ]
    
    if (verbose)  {
        print(paste0("After applying filter: ", nrow(vcf.data.temp), 
                     " rows in vcf.data.filtered"))
        print(paste0("After applying filter: ", nrow(vcf.data.artifact), 
                     " rows in vcf.data.artifact"))
    }
    
    # Return two vcf.obj
    return(list(vcf.obj.filtered = list(data = vcf.data.temp, 
                                        header = vcf.obj$header, 
                                        genome = vcf.obj$genome),
                vcf.obj.artifact = list(data = vcf.data.artifact,
                                        header = vcf.obj$header, 
                                        genome = vcf.obj$genome)))
}