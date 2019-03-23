# FIREVAT Annotation Functions
#
# Last revised date:
#   February 21, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title PrepareAnnotationDB
#' @description
#' Prepares df.genes.of.interest from a vcf.obj (\code{\link{ReadVCF}})
#' of COSMIC or ClinVar vcf for \code{\link{AnnotateVCFObj}}
#'
#' @param annotation.vcf.obj vcf.obj of COSMIC or ClinVar vcf file
#'
#' @return A data.frame of annotation.vcf.obj
#'
#' @export
PrepareAnnotationDB <- function(annotation.vcf.obj) {
    df.annotation.db <- data.frame(annotation.vcf.obj$data)
    return(df.annotation.db)
}


#' @title AnnotateVCFObj
#' @description
#' Annotates a vcf.obj using df.variants.of.interest (from \code{\link{PrepareAnnotationDB}})
#'
#' @param vcf.obj \code{\link{ReadVCF}}
#' @param df.annotation.db A data.frame from \code{\link{PrepareAnnotationDB}}.
#' This data.frame must have the columns 'CHROM', 'POS', 'REF', 'ALT'
#' @param columns.to.include A character vector of columns to include.
#' Note that existing columns in vcf.obj will not be affected.
#' @param include.all.columns A boolean value. If TRUE, then annotates vcf.obj with
#' all columns present in df.variants.of.interest. If FALSE, columns.to.include must be
#' supplied.
#'
#' @return An annotated vcf.obj
#'
#' @importFrom dplyr left_join
#' @export
AnnotateVCFObj <- function(vcf.obj,
                           df.annotation.db,
                           columns.to.include,
                           include.all.columns = FALSE) {
    if (include.all.columns) {
        columns.to.include <- colnames(df.annotation.db)
    }

    query.cols <- c("CHROM", "POS", "REF", "ALT")
    vcf.obj.data.cols <- names(vcf.obj$data)
    additional.cols <- setdiff(columns.to.include, vcf.obj.data.cols)
    df.annotated.vcf.obj <- left_join(vcf.obj$data, df.annotation.db, by = query.cols)

    # Annotate vcf.obj
    for (i in additional.cols) {
        vcf.obj$data[[`i`]] <- df.annotated.vcf.obj[,i]
    }
    return(vcf.obj)
}


#' @title FilterAnnotatedVCF
#' @description
#' Annotates a vcf.obj using df.variants.of.interest (from (\code{\link{PrepareAnnotationDB}})
#'
#' @param vcf.obj.annotated \code{\link{AnnotateVCFObj}}
#' @param filter.key.value.pairs A list with the key as the column name and value as the
#' filtering values. E.g. list("CLNSIG" = c("Pathogenic", "Pathogenic/Likely_pathogenic"))
#' @param filter.condition 'AND' or 'OR'.
#'
#' @return A vcf.obj
#'
#' @export
QueryAnnotatedVCF <- function(vcf.obj.annotated,
                              filter.key.value.pairs,
                              filter.condition = 'AND') {
    # Query vcf.obj
    condition.list <- list()
    for (i in names(filter.key.value.pairs)) {
        curr.key <- i
        curr.values <- filter.key.value.pairs[[curr.key]]
        condition.list[[curr.key]] <- vcf.obj.annotated$data[[`curr.key`]] %in% curr.values
    }

    # Apply filter condition
    if (filter.condition == "AND") {
        include <- AND.multiple(condition.list)
    }
    else if (filter.condition == "OR") {
        include <- OR.multiple(condition.list)
    }
    else {
        stop("The parameter 'filter.condition' must be either 'AND' or 'OR'")
    }

    # Prepare query data
    vcf.data.survived <- vcf.obj.annotated$data[include, ]
    queried.vcf.obj <- list(data = vcf.data.survived,
                            header = vcf.obj.annotated$header,
                            genome = vcf.obj.annotated$genome)
    return(queried.vcf.obj)
}

