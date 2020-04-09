# FIREVAT VCF Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title ReadVCF
#' @description Reads a .vcf file
#'
#' @param vcf.file (full path of a .vcf file)
#' @param genome A genome name for BSgenome (default: hg19)
#' @param split.info A boolean value. If TRUE, then makes the INFO column in the vcf
#' as a separate column. Default value is FALSE.
#' @param check.chromosome.name A boolean value. If TRUE, then check whether converts
#' @param match.bsg A boolean value. If TRUE, check whether chromosome matches with those in BSgenome
#' 'MT' to 'M' and adds 'chr' to the CHROM column. Default value is TRUE.
#'
#' @return A list with elements 'data', 'header', 'genome'
#'
#' @export
#' @importFrom bedr read.vcf
ReadVCF <- function(vcf.file,
                    genome = "hg19",
                    split.info = FALSE,
                    check.chromosome.name = TRUE,
                    match.bsg = FALSE) {
    vcf.temp <- read.vcf(x = vcf.file,
                         split.info = split.info,
                         verbose = FALSE)

    raw.data <- readLines(vcf.file)
    raw.header <- raw.data[(1:grep("#CHROM*", raw.data) - 1)]

    vcf.obj <- list(data = vcf.temp$vcf,
                    header = vcf.temp$header,
                    header.raw = raw.header,
                    genome = genome)

    if (check.chromosome.name == TRUE) {
        # Check "CHROM" column.
        # If name of mitochondrial chromosome is given as "MT",
        # change it into "M"
        chrom.MT.check <- grepl("MT", vcf.obj$data$CHROM)
        vcf.obj$data$CHROM[chrom.MT.check] <- "M"

        # Check if chromosome names start with "chr"
        # If not, paste "chr"
        chrom.chr.check <- grepl("chr", vcf.obj$data$CHROM)

        vcf.obj$data$CHROM[chrom.chr.check == FALSE] <- paste0(
            "chr", vcf.obj$data$CHROM[chrom.chr.check == FALSE]
        )
    }

    if (match.bsg == TRUE) {
        # If match.bsg is TRUE,
        # compare chromosome names of vcf with names in BSgenome
        bsg <- BSgenome::getBSgenome(genome)
        chrom.name.match <- vcf.obj$data$CHROM %in% bsg@seqinfo@seqnames
        vcf.obj$data <- vcf.obj$data[chrom.name.match,]
    }

    return(vcf.obj)
}


#' @title WriteVCF
#' @description Writes a vcf.obj to a .vcf file
#'
#' @param vcf.obj (from the function ReadVCF)
#' @param save.file (full path including filename)
#'
#' @export
#' @importFrom bedr write.vcf
WriteVCF <- function(vcf.obj, save.file) {
    vcf.obj.temp <- list(header = vcf.obj$header,
                         vcf = vcf.obj$data)
    attr(vcf.obj.temp, "vcf") <- TRUE

    result = tryCatch({
        write.vcf(vcf.obj.temp, filename = save.file, verbose = FALSE)
    }, error = function(error_condition) {
        print("There appears to be a format issue with the VCF file. Writing the file header and data separately.")
        # Write header and data separately
        write(vcf.obj$header.raw, file = save.file)
        write.table(vcf.obj$data,
                    sep = "\t",
                    col.names = F,
                    row.names = F,
                    quote = F,
                    append = T,
                    file = save.file)
    })
}


#' @title InitializeVCF
#' @description
#' Initialize VCF with FIREVAT config file
#' This functions selects point mutations and
#' appends filter values to vcf.obj$data
#'
#' @param vcf.obj A list from ReadVCF
#' @param config.obj A list from ParseConfigFile
#' @param verbose If true, provides process detail
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{vcf.obj.filtered}{vcf.obj (high-quality vcf)}
#'  \item{vcf.obj.artifact}{vcf.obj (low-quality vcf)}
#' }
#'
#' @export
#' @import data.table
InitializeVCF <- function(vcf.obj, config.obj, verbose = TRUE) {
    # Boolean vector that decides whether to keep each each row
    # If point mutation then insert the row into vcf.obj.filtered
    # else insert the row into vcf.obj.artifact
    atcg.chars <- c("A", "T", "C", "G")
    nrows <- nrow(vcf.obj$data)

    if (verbose == TRUE) {
        PrintLog(paste0("* Before initialization: ", nrow(vcf.obj$data), " rows in the original VCF file."))
        PrintLog("* Reading VCF file FORMAT column.")
    }

    # Parse FORMAT keys
    FORMAT.keys <- strsplit(vcf.obj$data$FORMAT, ":")

    if (verbose == TRUE) {
        PrintLog("* Reading VCF file INFO column.")
    }
    # Parse INFO column
    INFO.vals <- sapply(strsplit(vcf.obj$data$INFO,";"), function(x){
        return(strsplit(x, "="))}, simplify = FALSE)
    INFO.vals <- sapply(INFO.vals, MatchINFOKeyValue, simplify = FALSE)

    if (verbose == TRUE) {
        PrintLog(paste0("** Parsing VCF column: ", colnames(vcf.obj$data)[10:ncol(vcf.obj$data)]))
    }

    # Parse FORMAT values
    # Return values as matrix
    FORMAT.vals <- mapply(function(x) strsplit(x,":"),
                          vcf.obj$data[,10:ncol(vcf.obj$data)], SIMPLIFY = T)

    for (col in colnames(FORMAT.vals)) {
        FORMAT.vals[,col] <- mapply(function(x, y){names(x) <- y; return(x)},
                                    FORMAT.vals[,col], FORMAT.keys, SIMPLIFY = F)
    }

    # Save parameter value vectors to vcf.obj$data
    for (j in 1:length(config.obj)) {
        # Get FORMAT values
        if (attr(config.obj, "FORMAT_value_index")[j]) {
            name <- as.character(config.obj[[j]]$name)

            if (verbose == TRUE) {
                PrintLog(paste0("*** Getting VCF attribute: ", name))
            }

            column_header <- as.character(config.obj[[j]]$field$column_header)

            # If column_header was given with "$",
            # use the integer behind "$" as the column index
            if (startsWith(column_header,"$")) {
                column_index <- as.integer(substring(column_header,2))
                column_header <- colnames(vcf.obj$data)[column_index]
            }

            key <- as.character(config.obj[[j]]$field$key)
            index <- as.integer(config.obj[[j]]$field$index)
            type <- as.character(config.obj[[j]]$type)
            vcf.obj$data[[name]] <-sapply(sapply(FORMAT.vals[,column_header],
                function(x) return(strsplit(x[key],","))),
                function(y) return(y[index]))

            # Transform type of acquired values
            vcf.obj$data[[name]] <- TransformValueType(vcf.obj$data[[name]], type)

        # Get INFO values
        } else if (attr(config.obj, "INFO_value_index")[j]) {
            name <- as.character(config.obj[[j]]$name)

            if (verbose == TRUE) {
                PrintLog(paste0("*** Getting VCF attribute: ", name))
            }

            key <- as.character(config.obj[[j]]$field$key)
            index <- as.integer(config.obj[[j]]$field$index)
            type <- as.character(config.obj[[j]]$type)


            vcf.obj$data[[name]] <- sapply(sapply(INFO.vals,
            function(x) return(strsplit(x[key],","))),
            function(y) return(y[index]))

            # Transform type of acquired values
            vcf.obj$data[[name]] <- TransformValueType(vcf.obj$data[[name]], type)
        }
    }

    # Generate values with custom operations
    for (op.entry in config.obj[attr(config.obj,"OP_value_index")]) {
        name <- as.character(op.entry$name)

        if (verbose) {
            PrintLog(paste0("*** Calculating VCF attribute: ", name))
        }

        function_string <- as.character(op.entry$op$function_string)
        eval_function <- op.entry$eval_function
        vcf.args <- as.character(op.entry$op$args)
        type <- as.character(op.entry$type)

        args.names <- as.character(formalArgs(eval_function))
        args.values <- vcf.obj$data[, vcf.args, with=FALSE]
        args.for.function <- as.list(args.values)

        names(args.for.function) <- args.names

        vcf.obj$data[[name]] <- do.call(eval_function, args.for.function)
    }

    # Filter data
    # Only consider point mutations for downstream analysis
    if (verbose == TRUE) {
        PrintLog("* Checking if given mutations are point mutations.")
    }

    condition.list <- list(
        "REF.atcg" = vcf.obj$data$REF %in% atcg.chars,
        "ALT.atcg" = vcf.obj$data$ALT %in% atcg.chars
    )

    if (verbose == TRUE) {
        PrintLog("* Removing VCF rows with NA parameters.")
    }

    for (name in names(config.obj)) {
        condition.list[[name]] <- !is.na(vcf.obj$data[[name]])
    }

    # Filter vcf with generated conditions
    include <- AND.multiple(condition.list)

    vcf.data.filtered <- vcf.obj$data[include, ]
    vcf.data.artifact <- vcf.obj$data[!include, ]

    if (verbose == TRUE) {
        PrintLog("* After initialization:")
        PrintLog(paste0("** ", nrow(vcf.data.filtered), " rows in vcf.data.filtered VCF object."))
        PrintLog(paste0("** ", nrow(vcf.data.artifact), " rows in vcf.data.artifact VCF object."))
    }

    return(list(vcf.obj.filtered = list(data = vcf.data.filtered,
                                        header = vcf.obj$header,
                                        genome = vcf.obj$genome),
                vcf.obj.artifact = list(data = vcf.data.artifact,
                                        header = vcf.obj$header,
                                        genome = vcf.obj$genome)))
}


#' @title GetVCFValues
#' @description Get values of filtering parameters from vcf.obj
#'
#' @param vcf.obj A list of vcf data
#' @param vcf.filter  A list of vcf filtering information
#'
#' @return A list with filtering parameter values
#'
#' @keywords internal
#' @import data.table
GetVCFValues <- function(vcf.obj, vcf.filter) {
    vcf.params <- names(vcf.filter)
    output.list <- vcf.obj$data[, vcf.params, with=FALSE]
    return(output.list)
}


#' @title FilterByStrandBiasAnalysis
#' @description Filters refined.vcf.obj by strand bias analysis and
#' moves these filtered variants to artifactual.vcf.obj
#'
#' @param refined.vcf.obj A list of vcf data
#' @param artifactual.vcf.obj A list of vcf data
#' @param perform.fdr.correction A boolean value.
#' @param filter.by.strand.bias.analysis.cutoff A numeric value.
#'
#' @return A list with filtering parameter values
#' \itemize{
#'  \item{refined.vcf.obj}{ updated refined.vcf.obj}
#'  \item{artifactual.vcf.obj}{ updated artifactual.vcf.obj}
#' }
#'
#' @importFrom dplyr select
#' @export
FilterByStrandBiasAnalysis <- function(refined.vcf.obj,
                                       artifactual.vcf.obj,
                                       perform.fdr.correction,
                                       filter.by.strand.bias.analysis.cutoff) {
    if (perform.fdr.correction == TRUE) { # filter by q value
        refined.vcf.obj.temp <- refined.vcf.obj
        refined.vcf.obj.temp$data <- refined.vcf.obj$data[refined.vcf.obj$data$StrandBiasQValue < filter.by.strand.bias.analysis.cutoff,]
        refined.vcf.obj.temp$data <- select (refined.vcf.obj.temp$data, -c('StrandBiasQValue', 'StrandBiasPValue'))
        # keep.columns <- colnames(refined.vcf.obj.temp$data)[!(colnames(refined.vcf.obj.temp$data) %in% c('StrandBiasQValue', 'StrandBiasPValue'))]
        # refined.vcf.obj.temp$data <- refined.vcf.obj.temp$data[, ..keep.columns]
        artifactual.vcf.obj$data <- rbind(artifactual.vcf.obj$data, refined.vcf.obj.temp$data)
        refined.vcf.obj$data <- refined.vcf.obj$data[refined.vcf.obj$data$StrandBiasQValue >= filter.by.strand.bias.analysis.cutoff,]
    } else { # filter by p value
        refined.vcf.obj.temp <- refined.vcf.obj
        refined.vcf.obj.temp$data <- refined.vcf.obj$data[refined.vcf.obj$data$StrandBiasPValue < filter.by.strand.bias.analysis.cutoff,]
        refined.vcf.obj.temp$data <- select (refined.vcf.obj.temp$data, -c('StrandBiasPValue'))
        # keep.columns <- colnames(refined.vcf.obj$data)[!(colnames(refined.vcf.obj$data) %in% c('StrandBiasPValue'))]
        # refined.vcf.obj.temp$data <- refined.vcf.obj.temp$data[, ..keep.columns]
        artifactual.vcf.obj$data <- rbind(artifactual.vcf.obj$data, refined.vcf.obj.temp$data)
        refined.vcf.obj$data <- refined.vcf.obj$data[refined.vcf.obj$data$StrandBiasPValue >= filter.by.strand.bias.analysis.cutoff,]
    }
    return(list(refined.vcf.obj = refined.vcf.obj,
                artifactual.vcf.obj = artifactual.vcf.obj))
}


