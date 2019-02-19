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
#' @return a data.frame of the PCAWG mutatioanl signatures
#'
#' @export
#' @importFrom utils read.csv
GetPCAWGMutSigs <- function() {
    f <- system.file("extdata", "PCAWG_sigProfiler_SBS_Signatures_2018_03_28.csv", package = "FIREVAT")
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