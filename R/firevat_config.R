# FIREVAT Configure Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title ParseConfigFile
#' @description 
#' This function returns config.obj from JSON or YAML config file.
#' - Check if the config file is in JSON format or YAML format
#' - Return config.obj 
#'
#' @param config.path A string for config file path
#' @param verbose If true, provides process detail
#'
#' @return config.obj: list of parameters
#'
#' @examples
#' \dontrun{
#' ParseConfigFile("example.variant.caller.json")
#' ParseConfigFile("example.variant.caller.json", verbose=False)
#' }
#' @export
#' @importFrom jsonlite validate read_json
#' @importFrom yaml read_yaml
ParseConfigFile <- function(config.path, verbose=TRUE) {
    if(verbose) {
        print(paste("Parsing",config.path))
    }

    suppressWarnings(
        config.string <- readChar(config.path, file.info(config.path)$size)
    )

    # Check if config file is in json format
    if (jsonlite::validate(config.string)==TRUE) {
        config.obj <- read_json(config.path)
    } else {
        tryCatch ({
            config.obj <- read_yaml(config.path)
        }, error = function(e) {
            # Error message if config file is in neither json nor yaml
            read.error.message <- paste("ERROR: Cannot read config file:",
                                        config.path)
            stop(read.error.message)
        })
    }

    # Set names for each parameter. Names should be given in config file.
    names(config.obj)<-lapply(config.obj, function(x) return(x$name))
    
    attr(config.obj, "FORMAT_value_index") <- rep(FALSE, length(config.obj))
    attr(config.obj, "INFO_value_index") <- rep(FALSE, length(config.obj))
    attr(config.obj, "OP_value_index") <- rep(FALSE, length(config.obj))
    
    # Evaluate custom functions
    for(i in 1:length(config.obj)) {
        config.entry <- config.obj[[names(config.obj)[i]]]
        
        if ("op" %in% names(config.entry)) {
            tryCatch ({
                eval_function <- eval(parse(text = config.entry$op$function_string))
                config.obj[[names(config.obj)[i]]][["eval_function"]] <- eval_function
                attr(config.obj,"OP_value_index")[i] <- TRUE
            }, error = function(f) {
                # Error message if eval does not work
                eval.error.message <- paste(
                    "ERROR: Cannot evaluate given function:",
                    config.entry$op$function_string)
                stop(eval.error.message)
            })
        } else if ("field" %in% names(config.entry)) {
            if (config.entry[["field"]][["field_type"]] == "FORMAT") {
                attr(config.obj,"FORMAT_value_index")[i] <- TRUE
            } else if (config.entry[["field"]][["field_type"]] == "INFO") {
                attr(config.obj,"INFO_value_index")[i] <- TRUE
            }
        }   
    }

    if (verbose) {
        print(paste("Completed parsing",config.path))
    }

    return(config.obj)
}


#' @title TransformValueType
#' @description This function checks the type of value and returns type-transformed value  
#' 
#' @param curr.value A value which need type conversion
#' @param value_type A character value
#'
#' @return curr.value. Type-converted value (integer or numeric)
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' TransformValueType("0.98", "double")
#' TransformValueType("3.0", "integer")
#' }
TransformValueType <- function(curr.value, value_type) {
    # Alias for each value type
    # Can be modified later
    integer.alias <- c("int","integer","i")
    numeric.alias <- c("float","double","numeric","decimal","f","d","n")
    boolean.alias <- c("bool","bin","boolean","binary","logical","b","l")
    character.alias <- c("char","str","character","string","c","s")

    value_type <- tolower(value_type)

    # Change types
    if (value_type %in% integer.alias) {
        curr.value <- as.integer(curr.value)
    } else if (value_type %in% numeric.alias) {
        curr.value <- as.numeric(curr.value)
    } else if (value_type %in% boolean.alias) {
        curr.value <- as.logical(curr.value)
    } else if (value_type %in% character.alias) {
        curr.value <- as.character(curr.value)
    }

    return(curr.value)
}


#' @title An internal function for parsing INFO column.
#' @description
#' Check whether an INFO entry is given with single key or key-value pair.
#' Return values accordingly.
#' 
#' @param INFO.list A list of INFO values
#'
#' @return A vector of matching INFO values
#' 
#' @keywords internal
MatchINFOKeyValue <- function(INFO.list) {
    out.char <- sapply(INFO.list, function(INFO.entry) {
        if (length(INFO.entry) == 1) {
            names(INFO.entry) <- INFO.entry
            return(INFO.entry)
        } else {
            out.INFO.entry <- INFO.entry[2]
            names(out.INFO.entry) <- INFO.entry[1]
            return(out.INFO.entry)
        }
    })

    return(out.char)
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