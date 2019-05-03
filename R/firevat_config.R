# FIREVAT Configure Functions
#
# Last revised date:
#   February 22, 2019
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
ParseConfigFile <- function(config.path, verbose = TRUE) {
    if (verbose == TRUE) {
        PrintLog("* Started parsing config file: ")
        PrintLog(paste0("** ", config.path))
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
            stop(paste0("Cannot read config file: ", config.path))
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

    if (verbose == TRUE) {
        PrintLog("* Finished parsing config file:")
        PrintLog(paste0("** ", config.path))
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


#' @title Generate config.obj by checking vcf header
#' @description
#' This function generate config.obj by checking vcf header.
#' Users should fill in the information needed in console.
#' In current version, only Integers & Float values can be used in
#' config.obj for running FIREVAT.
#'
#' @param vcf.obj A list from \code{\link{ReadVCF}}
#' @param save.config If true, save config.obj to config.path
#' @param config.path File path to write config.obj (json or yaml)
#'
#' @return config.obj
#'
#' @export
#' @importFrom jsonlite write_json
#' @importFrom yaml write_yaml
#' @importFrom tools file_ext
GenerateConfigObj <- function(vcf.obj, save.config = TRUE,
                              config.path = "../temp/FIREVAT_configure.json") {
    # FORMAT
    ## mandatory fields: name/direction/field_type/column_header/key/index/type
    ## optional fields: default/range
    # INFO
    ## mandatory fields: name/direction/field_type/key/index/type
    ## optional fields: default/range

    # Define config.obj
    print(paste("VCF Spec:", vcf.obj$header$fileformat))
    print("Checklist protocol, initiated.")

    config.obj <- list()

    column.names <- c("FORMAT", "INFO")

    # Check FORMAT & INFO columns.
    for (col.name in column.names) {
        col.name.info <- vcf.obj$header[[col.name]]

        # Show fields
        print(paste(length(col.name.info[,"ID"]),col.name,"fields identified."))
        print(paste0(col.name.info[,"ID"], ": ",
                     col.name.info[,"Description"],
                     ", Type: ", col.name.info[,"Type"]))

        i <- 1

        repeat{
            if (i >  length(col.name.info[,"ID"])){
                break
            }

            # Get values
            key <- col.name.info[i,"ID"]
            print(paste("Checking", key))
            type <- col.name.info[i,"Type"]
            col.samples <- vcf.obj$header$SAMPLE[,"ID"]

            key.count <- as.integer(col.name.info[i,"Number"])

            if (is.na(key.count)){
                i <- i + 1
                next
            }

            # Skip if type == String or Flag
            if (type %in% c("String", "Flag")) {
                i <- i + 1
                next
            }
            # Check key.count: 1 or more
            if (key.count == 1) {
                print(paste("There is a value in",key,"field"))
                use.in.filter <- as.logical(readline(
                    prompt="Include this in config.obj (T/F): "))
                if (use.in.filter==TRUE) {
                    field.index <- 1
                } else {
                    i <- i + 1
                    next
                }
            } else {
                print(paste("There are",key.count,"values in",key,"field"))
                use.in.filter <- as.logical(readline(
                    prompt="Include this in config.obj (T/F): "))
                if (use.in.filter==TRUE) {
                    field.index <- NA
                } else {
                    i <- i + 1
                    next
                }
            }

            if (use.in.filter==TRUE) {
                param.name <- readline(prompt="Parameter name: ")
                field.type <- col.name
                field.key <- key
                # if key.count > 1, get field_index here
                if (key.count > 1) {
                    field.index <- as.integer(readline(prompt=
                        paste("Index of the value to use (1-based): ")
                    ))
                }
                # if col.name == "FORMAT", get column names to extract value.
                if (col.name == "FORMAT") {

                    field.column.header <- readline(prompt=
                        paste0("Get values from column (",
                               paste(col.samples, collapse=", "),
                               "): ")
                    )
                    # Check whether given column header is present
                    if (!(field.column.header %in% col.samples)) {
                        stop(paste("Fail to find column",field.column.header))
                    }
                } else if (col.name == "INFO") {
                    field.column.header <- ""
                }

                direction <- readline(prompt=
                    paste0("Filtering direction (",
                           "POS: Parameter of Refined variants > CUTOFF, ",
                           "NEG: Parameter of Refined variants < CUTOFF): ")
                )
                # Check whether given direction is valid
                if (!(direction %in% c("POS","NEG"))) {
                    stop(paste("Fail to recognize direction:", direction))
                }
                # (Optional field) Get default cutoff & filtering range
                default.cutoff <- readline(prompt=
                    "Default cutoff (Optional: Press ENTER to ignore): "
                )
                filtering.range <- readline(prompt=
                    "Filtering range (min,max; Optional: Press ENTER to ignore): "
                )

                # Check whether given information is right.
                print(paste0("Check this information.",
                             "name: ", param.name, ";",
                             "direction: ", direction, ";",
                             "field_type: ", field.type, ";",
                             "column_header: ", field.column.header, ";",
                             "key: ", field.key, ";",
                             "index: ", field.index, ";",
                             "type: ", type, ";",
                             "default: ", default.cutoff, ";",
                             "range: ", filtering.range, ";",
                             "use_in_filter: ", use.in.filter))
                check.before.save <- as.logical(readline(prompt=
                                     "Save this information to config.obj (T/F): "))

                suppressWarnings({
                    if (check.before.save) {
                        config.obj[[param.name]] <- list(
                            "name" = param.name,
                            "direction" = direction,
                            "field" = list(
                                "field_type" = field.type,
                                "key" = field.key,
                                "index" = field.index
                            ),
                            "type" = type,
                            "use_in_filter" = use.in.filter
                        )
                        # If col.name == "FORMAT", save column header
                        if (col.name == "FORMAT"){
                            config.obj[[param.name]]["field"]["column_header"] <- field.column.header
                        }
                        # default: Optional
                        if (default.cutoff != "") {
                            default.cutoff <- as.integer(default.cutoff)
                            config.obj[[param.name]]["default"] <- default.cutoff
                        }
                        # range: Optional
                        if (filtering.range != "") {
                            filtering.range <- as.integer(unlist(strsplit(
                                gsub(" ","",filtering.range), ",")))
                            config.obj[[param.name]]["range"] <- filtering.range
                        }
                        print(paste("Saved",param.name,"to config.obj"))
                    }
                })
            }

            # Check whether to repeat generatation parameters from this key.
            print(paste("Current field:",col.name))
            print(paste("Current key:", key))
            if(col.name=="FORMAT") {
                print(paste("Sample columns:", paste0(col.samples,collapse = ", ")))
                print(paste("Current sample column:", field.column.header))
            }
            print(paste("Current index of value:", field.index))
            print(paste("Number of values in this key:", key.count))
            # Check whether to finish or not
            if (i == length(col.name.info[,"ID"])) {
                end.flag <- as.logical(readline(prompt=
                    paste0("You reached last values of ", col.name, ".",
                           "Proceed to next step (T/F): ")
                ))
                if (end.flag==TRUE) {
                    break
                } else {
                    next
                }
            }
            # Jump to next parameter
            key.next <- as.logical(readline(prompt="Jump to next param (T/F): "))

            if (key.next) {
                i <- i + 1
                next
            }
        }
    }

    # Save config.obj if save.config == TRUE
    if (save.config) {
        config.extension <- file_ext(config.path)
        if (config.extension == "json") {
            write_json(config.obj, config.path,
                       auto_unbox = TRUE)
        } else if (config.extension == "yaml") {
            write_yaml(config.obj, config.path)
        } else {
            stop("FIREVAT config file should be in json or yaml format.")
        }
    }

    # Jobs Done!
    print("Completed.")
    return(config.obj)
}
