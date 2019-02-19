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
                                        "SBS27","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
                                        "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40",
                                        "SBS41","SBS42","SBS44","SBS45")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.Known.Sequencing.Artifact.Signatures <- c("SBS60")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.Likely.Sequencing.Artifact.Signatures <- c("SBS46","SBS47", "SBS48","SBS50","SBS53")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.Possible.Sequencing.Artifact.Signatures <- c("SBS27","SBS43","SBS49","SBS51","SBS52",
                                                   "SBS54","SBS55","SBS56","SBS57","SBS58",
                                                   "SBS59")


#' @title Constant
#' @description
#' PCAWG mutational signatures reported to be associated with sequencing artifacts
#' @export
PCAWG.All.Sequencing.Artifact.Signatures <- c(PCAWG.Known.Sequencing.Artifact.Signatures,
                                              PCAWG.Likely.Sequencing.Artifact.Signatures,
                                              PCAWG.Possible.Sequencing.Artifact.Signatures)


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