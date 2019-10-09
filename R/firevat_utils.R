# FIREVAT Utility Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title GetPCAWGMutSigs
#' @description Returns the PCAWG mutational signatures data
#'
#' @param sequencing.type A string value.
#' It can be either 'wes' for whole-exome sequencing or 'wgs' for whole-genome sequencing
#'
#' @return a data.frame of the PCAWG mutatioanl signatures
#'
#' @export
#' @importFrom utils read.csv
GetPCAWGMutSigs <- function(sequencing.type = "wes") {
    if (sequencing.type == "wes") {
        f <- system.file("extdata", "PCAWG_sigProfiler_exome_SBS_signatures_20180321.csv", package = "FIREVAT")
    } else if (sequencing.type == "wgs") {
        f <- system.file("extdata", "PCAWG_sigProfiler_genome_SBS_Signatures_20180328.csv", package = "FIREVAT")
    } else {
        print("Invalid sequencing.type parameter. It must be either 'wes' or 'wgs'")
        return(NULL)
    }
    df.pcawg.mut.sigs <- read.csv(f,
                                  header = T,
                                  check.names = F,
                                  stringsAsFactors = F,
                                  dec = ".")
    colnames(df.pcawg.mut.sigs)[1] <- "Type"
    colnames(df.pcawg.mut.sigs)[2] <- "SubType"
    colnames(df.pcawg.mut.sigs)[3] <- "SomaticMutationType"
    return(df.pcawg.mut.sigs)
}


#' @title GetPCAWGMutSigsNames
#' @description Returns the PCAWG mutational signatures names
#'
#' @return a character vector of the PCAWG mutational signatures names
#'
#' @export
GetPCAWGMutSigsNames <- function() {
    df.pcawg.mut.sigs <- GetPCAWGMutSigs()
    pcawg.mut.sig.names <- as.character(colnames(df.pcawg.mut.sigs)[4:ncol(df.pcawg.mut.sigs)])
    return(pcawg.mut.sig.names)
}


#' @title GetPCAWGMutSigsEtiologiesColors
#' @description Returns the PCAWG mutational signatures etiologies and colors
#'
#' @return a data.frame with the columns 'signature', 'group', 'color'
#'
#' @export
#' @importFrom utils read.csv
GetPCAWGMutSigsEtiologiesColors <- function() {
    f <- system.file("extdata", "PCAWG_Signatures_Groups_Colors.csv", package = "FIREVAT")
    df.pcawg.mut.sigs <- read.csv(f,
                                  header = T,
                                  check.names = F,
                                  stringsAsFactors = F)
    return(df.pcawg.mut.sigs)
}


#' @title GetPCAWGPlatinumMutSigs
#' @description Returns the PCAWG platinum mutational signatures data
#'
#' @return a data.frame of the PCAWG platinum mutatioanl signatures
#'
#' @export
#' @importFrom utils read.csv
GetPCAWGPlatinumMutSigs <- function() {
    f <- system.file("extdata", "PCAWG_Platinum_Signatures.csv", package = "FIREVAT")
    df.pcawg.mut.sigs <- read.csv(f,
                                  header = T,
                                  check.names = F,
                                  stringsAsFactors = F,
                                  dec = ".")
    colnames(df.pcawg.mut.sigs)[1] <- "Type"
    colnames(df.pcawg.mut.sigs)[2] <- "SubType"
    colnames(df.pcawg.mut.sigs)[3] <- "SomaticMutationType"
    return(df.pcawg.mut.sigs)
}



#' @title GetPCAWGPlatinumMutSigsNames
#' @description Returns the PCAWG platinum mutational signatures names
#'
#' @return a character vector of the PCAWG platinum mutational signatures names
#'
#' @export
GetPCAWGPlatinumMutSigsNames <- function() {
    df.pcawg.mut.sigs <- GetPCAWGPlatinumMutSigs()
    pcawg.mut.sig.names <- as.character(colnames(df.pcawg.mut.sigs)[4:ncol(df.pcawg.mut.sigs)])
    return(pcawg.mut.sig.names)
}


#' @title GetPCAWGPlatinumMutSigsEtiologiesColors
#' @description Returns the PCAWG platinum mutational signatures etiologies and colors
#'
#' @return a data.frame with the columns 'signature', 'group', 'color'
#'
#' @export
#' @importFrom utils read.csv
GetPCAWGPlatinumMutSigsEtiologiesColors <- function() {
    f <- system.file("extdata", "PCAWG_Platinum_Signatures_Groups_Colors.csv", package = "FIREVAT")
    df.pcawg.mut.sigs <- read.csv(f,
                                  header = T,
                                  check.names = F,
                                  stringsAsFactors = F)
    return(df.pcawg.mut.sigs)
}


#' @title GetCOSMICMutSigs
#' @description Returns a data.frame of the COSMIC mutational signature reference file
#' from the data directory
#'
#' @return a data.frame of the COSMIC reference mutational signatures
#'
#' @export
#' @importFrom utils read.csv
GetCOSMICMutSigs <- function() {
    f <- system.file("extdata", "COSMIC_Signatures.csv", package = "FIREVAT")
    df.cosmic.mut.sigs <- read.csv(f,
                                   header = T,
                                   check.names = F,
                                   stringsAsFactors = F,
                                   sep = "\t",
                                   dec = ".")

    colnames(df.cosmic.mut.sigs)[1] <- "Type"
    colnames(df.cosmic.mut.sigs)[2] <- "SubType"
    colnames(df.cosmic.mut.sigs)[3] <- "SomaticMutationType"

    return(df.cosmic.mut.sigs)
}


#' @title GetCOSMICMutSigsNames
#' @description Returns all COSMIC mutational signature names
#'
#' @return a character vector
#'
#' @export
GetCOSMICMutSigsNames <- function() {
    df.cosmic.mut.sigs <- GetCOSMICMutSigs()
    cosmic.mut.sig.names <- as.character(colnames(df.cosmic.mut.sigs)[4:ncol(df.cosmic.mut.sigs)])
    return(cosmic.mut.sig.names)
}


#' @title GetCOSMICMutSigsNames
#' @description Returns all COSMIC mutational signature etiologies and colors
#'
#' @return data.frame with following columns: signature, group and color.
#'
#' @export
#' @importFrom utils read.csv
GetCOSMICMutSigsEtiologiesColors <- function() {
    f <- system.file("extdata", "COSMIC_Signatures_Groups_Colors.csv", package = "FIREVAT")
    df.cosmic.mut.sigs <- read.csv(f,
                                   header = T,
                                   check.names = F,
                                   stringsAsFactors = F)
    return(df.cosmic.mut.sigs)
}


#' @title EnumerateTriNucCounts
#' @description Returns C>A, C>G, C>T, T>A, T>C, T>G counts
#'
#' @param spectrum a numeric vector with 96 numeric values
#'
#' @details
#' Please note that this function assumes that 'spectrum' is sorted
#' (i.e. 1:16  --> C>A;
#'       17:32 --> C>G;
#'       33:48 --> C>T;
#'       49:64 --> T>A;
#'       65:80 --> T>C;
#'       81:96 --> T>G)
#'
#' @return a numeric vector of length 6 corresponding to the counts of each trinucleotide change (C>A, C>G, C>T, T>A, T>C, T>G)
#'
#' @export
EnumerateTriNucCounts <- function(spectrum) {
    C.to.A <- sum(spectrum[1:16])
    C.to.G <- sum(spectrum[17:32])
    C.to.T <- sum(spectrum[33:48])
    T.to.A <- sum(spectrum[49:64])
    T.to.C <- sum(spectrum[65:80])
    T.to.G <- sum(spectrum[81:96])
    return(c(C.to.A,
             C.to.G,
             C.to.T,
             T.to.A,
             T.to.C,
             T.to.G))
}


#' @title WriteFIREVATResultsToTSV
#' @description Writes FIREVAT results to a csv file
#'
#' @param firevat.results List returned from RunFIREVAT
#'
#' @export
WriteFIREVATResultsToTSV <- function(firevat.results) {
    save.list <- list(vcf.file = NA,
                      vcf.file.basename = NA,
                      config.file = NA,
                      start.datetime = NA,
                      end.datetime = NA,
                      variant.refinement.performed = NA,
                      variant.refinement.terminiation.log = NA,
                      original.muts.seq.art.weights.sum = NA,
                      target.mut.sigs = NA,
                      sequencing.artifact.mut.sigs = NA,
                      output.dir = NA,
                      ga.pop.size = NA,
                      ga.max.iter = NA,
                      ga.pmutation = NA,
                      ga.run = NA,
                      mutalisk.method = NA,
                      mutalisk.random.sampling.count = NA,
                      mutalisk.random.sampling.max.iter = NA,
                      mode = NA,
                      perform.strand.bias.analysis = NA,
                      strand.bias.fdr.correction.method = NA,
                      strand.bias.perform.fdr.correction = NA,
                      raw.mutational.signatures = NA,
                      raw.mutational.signatures.probs = NA,
                      raw.mutational.signatures.cos.sim = NA,
                      raw.mutational.signatures.rss = NA,
                      refined.mutational.signatures = NA,
                      refined.mutational.signatures.probs = NA,
                      refined.mutational.signatures.cos.sim = NA,
                      refined.mutational.signatures.rss = NA,
                      artifactual.mutational.signatures = NA,
                      artifactual.mutational.signatures.probs = NA,
                      artifactual.mutational.signatures.cos.sim = NA,
                      artifactual.mutational.signatures.rss = NA,
                      refinement.filters = NA,
                      refinement.filter.cutoffs = NA,
                      refinement.filter.directions = NA)

    if (is.null(firevat.results$vcf.file) == FALSE) {
        save.list$vcf.file <- firevat.results$vcf.file
    }
    if (is.null(firevat.results$vcf.file.basename) == FALSE) {
        save.list$vcf.file.basename <- firevat.results$vcf.file.basename
    }
    if (is.null(firevat.results$config.file) == FALSE) {
        save.list$config.file <- firevat.results$config.file
    }
    if (is.null(firevat.results$start.datetime) == FALSE) {
        save.list$start.datetime <- firevat.results$start.datetime
    }
    if (is.null(firevat.results$end.datetime) == FALSE) {
        save.list$end.datetime <- firevat.results$end.datetime
    }
    if (is.null(firevat.results$variant.refinement.performed) == FALSE) {
        save.list$variant.refinement.performed <- firevat.results$variant.refinement.performed
    }
    if (is.null(firevat.results$variant.refinement.terminiation.log) == FALSE) {
        save.list$variant.refinement.terminiation.log <- firevat.results$variant.refinement.terminiation.log
    }
    if (is.null(firevat.results$original.muts.seq.art.weights.sum) == FALSE) {
        save.list$original.muts.seq.art.weights.sum <- firevat.results$original.muts.seq.art.weights.sum
    }
    if (is.null(firevat.results$target.mut.sigs) == FALSE) {
        save.list$target.mut.sigs <- paste0(firevat.results$target.mut.sigs, collapse = ",")
    }
    if (is.null(firevat.results$sequencing.artifact.mut.sigs) == FALSE) {
        save.list$sequencing.artifact.mut.sigs <- paste0(firevat.results$sequencing.artifact.mut.sigs, collapse = ",")
    }
    if (is.null(firevat.results$output.dir) == FALSE) {
        save.list$output.dir <- firevat.results$output.dir
    }
    if (is.null(firevat.results$ga.pop.size) == FALSE) {
        save.list$ga.pop.size <- firevat.results$ga.pop.size
    }
    if (is.null(firevat.results$ga.max.iter) == FALSE) {
        save.list$ga.max.iter <- firevat.results$ga.max.iter
    }
    if (is.null(firevat.results$ga.pmutation) == FALSE) {
        save.list$ga.pmutation <- firevat.results$ga.pmutation
    }
    if (is.null(firevat.results$ga.run) == FALSE) {
        save.list$ga.run <- firevat.results$ga.run
    }
    if (is.null(firevat.results$mutalisk.method) == FALSE) {
        save.list$mutalisk.method <- firevat.results$mutalisk.method
    }
    if (is.null(firevat.results$mutalisk.random.sampling.count) == FALSE) {
        save.list$mutalisk.random.sampling.count <- firevat.results$mutalisk.random.sampling.count
    }
    if (is.null(firevat.results$mutalisk.random.sampling.max.iter) == FALSE) {
        save.list$mutalisk.random.sampling.max.iter <- firevat.results$mutalisk.random.sampling.max.iter
    }
    if (is.null(firevat.results$mode) == FALSE) {
        save.list$mode <- firevat.results$mode
    }
    if (is.null(firevat.results$perform.strand.bias.analysis) == FALSE) {
        save.list$perform.strand.bias.analysis <- firevat.results$perform.strand.bias.analysis
    }
    if (is.null(firevat.results$strand.bias.fdr.correction.method) == FALSE) {
        save.list$strand.bias.fdr.correction.method <- firevat.results$strand.bias.fdr.correction.method
    }
    if (is.null(firevat.results$strand.bias.perform.fdr.correction) == FALSE) {
        save.list$strand.bias.perform.fdr.correction <- firevat.results$strand.bias.perform.fdr.correction
    }

    # Raw/original mutations mutational signatures
    if (is.null(firevat.results$raw.muts.mutalisk.results) == FALSE) {
        if (is.null(firevat.results$raw.muts.mutalisk.results$identified.mut.sigs) == FALSE) {
            save.list$raw.mutational.signatures <- paste0(firevat.results$raw.muts.mutalisk.results$identified.mut.sigs, collapse = ",")
        }
        if (is.null(firevat.results$raw.muts.mutalisk.results$identified.mut.sigs.probs) == FALSE) {
            save.list$raw.mutational.signatures.probs <- paste0(firevat.results$raw.muts.mutalisk.results$identified.mut.sigs.probs, collapse = ",")
        }
        if (is.null(firevat.results$raw.muts.mutalisk.results$cos.sim.score) == FALSE) {
            save.list$raw.mutational.signatures.cos.sim <- firevat.results$raw.muts.mutalisk.results$cos.sim.score
        }
        if (is.null(firevat.results$raw.muts.mutalisk.results$rss) == FALSE) {
            save.list$raw.mutational.signatures.rss <- firevat.results$raw.muts.mutalisk.results$rss
        }
    }

    # Refined mutations mutational signatures
    if (is.null(firevat.results$refined.muts.mutalisk.results) == FALSE) {
        if (is.null(firevat.results$refined.muts.mutalisk.results$identified.mut.sigs) == FALSE) {
            save.list$refined.mutational.signatures <- paste0(firevat.results$refined.muts.mutalisk.results$identified.mut.sigs, collapse = ",")
        }
        if (is.null(firevat.results$refined.muts.mutalisk.results$identified.mut.sigs.probs) == FALSE) {
            save.list$refined.mutational.signatures.probs <- paste0(firevat.results$refined.muts.mutalisk.results$identified.mut.sigs.probs, collapse = ",")
        }
        if (is.null(firevat.results$refined.muts.mutalisk.results$cos.sim.score) == FALSE) {
            save.list$refined.mutational.signatures.cos.sim <- firevat.results$refined.muts.mutalisk.results$cos.sim.score
        }
        if (is.null(firevat.results$refined.muts.mutalisk.results$rss) == FALSE) {
            save.list$refined.mutational.signatures.rss <- firevat.results$refined.muts.mutalisk.results$rss
        }
    }

    # Artifactual mutations mutational signatures
    if (is.null(firevat.results$artifactual.muts.mutalisk.results) == FALSE) {
        if (is.null(firevat.results$artifactual.muts.mutalisk.results$identified.mut.sigs) == FALSE) {
            save.list$artifactual.mutational.signatures <- paste0(firevat.results$artifactual.muts.mutalisk.results$identified.mut.sigs, collapse = ",")
        }
        if (is.null(firevat.results$artifactual.muts.mutalisk.results$identified.mut.sigs.probs) == FALSE) {
            save.list$artifactual.mutational.signatures.probs <- paste0(firevat.results$artifactual.muts.mutalisk.results$identified.mut.sigs.probs, collapse = ",")
        }
        if (is.null(firevat.results$artifactual.muts.mutalisk.results$cos.sim.score) == FALSE) {
            save.list$artifactual.mutational.signatures.cos.sim <- firevat.results$artifactual.muts.mutalisk.results$cos.sim.score
        }
        if (is.null(firevat.results$artifactual.muts.mutalisk.results$rss) == FALSE) {
            save.list$artifactual.mutational.signatures.rss <- firevat.results$artifactual.muts.mutalisk.results$rss
        }
    }

    refinement.filters <- c()
    refinement.filter.cutoffs <- c()
    refinement.filter.directions <- c()
    for (curr.filter in names(firevat.results$vcf.filter)) {
        curr.filter.cutoff <- as.numeric(firevat.results$x.solution.decimal[[curr.filter]])
        curr.filter.direction <- firevat.results$config.obj[[curr.filter]]$direction
        refinement.filters <- c(refinement.filters, curr.filter)
        refinement.filter.cutoffs <- c(refinement.filter.cutoffs, curr.filter.cutoff)
        if (curr.filter.direction == "POS") {
            refinement.filter.directions <- c(refinement.filter.directions, ">=")
        } else {
            refinement.filter.directions <- c(refinement.filter.directions, "<=")
        }
    }
    save.list$refinement.filters = paste0(refinement.filters, collapse = ",")
    save.list$refinement.filter.cutoffs = paste0(refinement.filter.cutoffs, collapse = ",")
    save.list$refinement.filter.directions = paste0(refinement.filter.directions, collapse = ",")

    df.save <- data.frame(save.list,
                          stringsAsFactors = F)
    write.table(df.save,
                paste0(firevat.results$output.dir, firevat.results$vcf.file.basename, "_FIREVAT_data.tsv"),
                sep = "\t", row.names = F)
}


#' @title GetSeedForVCF
#' @description Returns a seed integer based on VCF file size
#'
#' @param vcf.file (full path of a .vcf file)
#'
#' @details
#' Returns the same seed integer for the same VCF file (based on file size)
#'
#' @return an integer value
#'
#' @export
GetSeedForVCF <- function(vcf.file) {
    file.size <- file.info(vcf.file)$size
    file.size <- round(file.size / 1000)
    return(file.size)
}
