# FIREVAT Constants
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title Constant
#' @description
#' PCAWG target mutational signatures reported to be unrelated to sequencing artifacts
#' @export
PCAWG.Target.Mutational.Signatures <- c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6",
                                        "SBS7a","SBS7b","SBS7c","SBS7d","SBS8","SBS9",
                                        "SBS10a","SBS10b","SBS11","SBS12","SBS13","SBS14",
                                        "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19",
                                        "SBS20","SBS21","SBS22","SBS23","SBS24","SBS25","SBS26",
                                        "SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
                                        "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40",
                                        "SBS41","SBS42","SBS44")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.Known.Sequencing.Artifact.Signatures <- c("SBS60")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.Possible.Sequencing.Artifact.Signatures <- c("SBS27","SBS43","SBS45","SBS46",
                                                   "SBS47","SBS48","SBS49","SBS50",
                                                   "SBS51","SBS52","SBS53","SBS54",
                                                   "SBS55","SBS56","SBS57","SBS58",
                                                   "SBS59")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.All.Sequencing.Artifact.Signatures <- c(PCAWG.Known.Sequencing.Artifact.Signatures,
                                              PCAWG.Possible.Sequencing.Artifact.Signatures)


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.Platinum.All.Technology.Related.Artifact.Signatures <- c("SBSR1",
                                                               "SBSR2",
                                                               "SBSR3",
                                                               "SBSR4",
                                                               "SBSR5",
                                                               "SBSR6",
                                                               "SBSR7",
                                                               "SBSR8",
                                                               "SBSR9")


#' @title Constant
#' @description
#' Hex codes for the mutation types (for plotting purposes)
#' @export
TriNuc.Mutation.Type.Hex.Colors <- c("#029ACE", # C>A
                                     "#231F20", # C>G
                                     "#CD5D5D", # C>T
                                     "#B4B3B4", # T>A
                                     "#9ACA3C", # T>C
                                     "#D7BFD7") # T>G


#' @title Constant
#' @description
#' Chromosome names for FIREVAT.
#' Chromosome names should be given in the format of "chr" + chromosome number.
#' @export
Chromosome.Names <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                      "chr8","chr9","chr10","chr11","chr12","chr13",
                      "chr14","chr15","chr16","chr17","chr18","chr19",
                      "chr20","chr21","chr22","chrX","chrY","chrM")
