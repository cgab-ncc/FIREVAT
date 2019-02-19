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
#' @param genome ('hg19' or 'hg38')
#'
#' @return A list with elements 'data', 'header', 'genome'
#'
#' @export
#' @importFrom bedr read.vcf
ReadVCF <- function(vcf.file, genome = "hg19") {
    vcf.temp <- read.vcf(x = vcf.file)
    return(list(data = vcf.temp$vcf,
                header = vcf.temp$header,
                genome = genome))
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
    write.vcf(vcf.obj.temp, filename = save.file)
}


#' @title InitializeVCFFromConfig
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
InitializeVCFFromConfig <- function(vcf.obj, config.obj, verbose=TRUE) {
    if (verbose) {
        print(Sys.time())
    }

    # Boolean vector that decides whether to keep each each row
    # If point mutation then insert the row into vcf.obj.filtered
    # else insert the row into vcf.obj.artifact
    atcg.chars <- c("A", "T", "C", "G")
    nrows <- nrow(vcf.obj$data)

    if (verbose) {
        print(paste0("Before initialization, there are ",
                     nrow(vcf.obj$data), " rows in the original vcf"))
        print(paste("Reading FORMAT column"))
    }

    # Parse FORMAT keys
    FORMAT.keys <- strsplit(vcf.obj$data$FORMAT, ":")

    if (verbose) {
        print(paste("Reading INFO column"))
    }
    # Parse INFO column
    INFO.vals <- sapply(strsplit(vcf.obj$data$INFO,";"), function(x){
        return(strsplit(x, "="))})
    INFO.vals <- sapply(INFO.vals, MatchINFOKeyValue)

    if (verbose) {
        print(paste("Parsing column",
        colnames(vcf.obj$data)[10:ncol(vcf.obj$data)]))
    }

    # Parse FORMAT values
    # Return values as matrix
    FORMAT.vals <- mapply(function(x) strsplit(x,":"),
                          vcf.obj$data[,10:ncol(vcf.obj$data)],SIMPLIFY = T)

    for (col in colnames(FORMAT.vals)) {
        FORMAT.vals[,col] <- mapply(function(x, y){names(x) <- y; return(x)},
                                    FORMAT.vals[,col], FORMAT.keys, SIMPLIFY = F)
    }

    # Save parameter value vectors to vcf.obj$data
    for (j in 1:length(config.obj)) {
        # Get FORMAT values
        if (attr(config.obj, "FORMAT_value_index")[j]) {
            name <- as.character(config.obj[[j]]$name)

            if (verbose) {
                print(paste("Getting", name))
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

            if (verbose) {
                print(paste("Getting", name))
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
            print(paste("Calculating", name))
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
    if (verbose) {
        print(paste("Checking if given mutations are point mutations"))
    }

    condition.list <- list(
        "REF.atcg" = vcf.obj$data$REF %in% atcg.chars,
        "ALT.atcg" = vcf.obj$data$ALT %in% atcg.chars
    )

    # Check "CHROM" column
    ## If name of mitochondrial chromosome is given as "MT",
    ## change it into "M"
    chrom.MT.check <- grepl("MT",vcf.obj$data$CHROM)
    vcf.obj$data$CHROM[chrom.MT.check] <- "M"

    ## Check if chromosome names start with "chr"
    ## If not, paste "chr"
    chrom.chr.check <- grepl("chr",vcf.obj$data$CHROM)

    vcf.obj$data$CHROM[chrom.chr.check==F] <- paste0(
        "chr", vcf.obj$data$CHROM[chrom.chr.check==F]
    )

    if (verbose) {
        print(paste("Removing lines with NA parameters"))
    }

    for (name in names(config.obj)) {
        condition.list[[name]] <- !is.na(vcf.obj$data[[name]])
    }

    # Filter vcf with generated conditions
    include <- AND.multiple(condition.list)

    vcf.data.filtered <- vcf.obj$data[include, ]
    vcf.data.artifact <- vcf.obj$data[!include, ]

    if (verbose) {
        print("After initialization,")
        print(paste0("There are ", nrow(vcf.data.filtered), " rows in vcf.data.filtered vcf"))
        print(paste0("There has ", nrow(vcf.data.artifact), " rows in vcf.data.artifact vcf"))
        print(Sys.time())
    }

    return(list(vcf.obj.filtered = list(data = vcf.data.filtered,
                                        header = vcf.obj$header,
                                        genome = vcf.obj$genome),
                vcf.obj.artifact = list(data = vcf.data.artifact,
                                        header = vcf.obj$header,
                                        genome = vcf.obj$genome)))
}