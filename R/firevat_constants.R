# FIREVAT Constants
#
# Last revised date:
#   February 18, 2019
#
# Authors: 
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


# PCAWG Target Mutational Signatures Reported to be Unrelated to Sequencing Artifacts)
PCAWG.Target.Mutational.Signatures <- c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6",
                                        "SBS7a","SBS7b","SBS7c","SBS7d","SBS8","SBS9",
                                        "SBS10a","SBS10b","SBS11","SBS12","SBS13","SBS14",
                                        "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19",
                                        "SBS20","SBS21","SBS22","SBS23","SBS24","SBS25","SBS26",
                                        "SBS27","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
                                        "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40",
                                        "SBS41","SBS42","SBS44","SBS45")

# PCAWG Mutational Signatures Reported to be Associated with Sequencing Artifacts
PCAWG.Known.Sequencing.Artifact.Signatures <- c("SBS60")
PCAWG.Likely.Sequencing.Artifact.Signatures <- c("SBS46","SBS47", "SBS48","SBS50","SBS53")
PCAWG.Possible.Sequencing.Artifact.Signatures <- c("SBS27","SBS43","SBS49","SBS51","SBS52",
                                                   "SBS54","SBS55","SBS56","SBS57","SBS58",
                                                   "SBS59")
PCAWG.All.Sequencing.Artifact.Signatures <- c(PCAWG.Known.Sequencing.Artifact.Signatures,
                                              PCAWG.Likely.Sequencing.Artifact.Signatures,
                                              PCAWG.Possible.Sequencing.Artifact.Signatures)

# Mutation Calling Methods
Mutation.Calling.Methods <- c("MuTect2")

# For Plotting
TriNuc.Mutation.Type.Hex.Colors <- c("#029ACE", # C>A
                                     "#231F20", # C>G
                                     "#CD5D5D", # C>T
                                     "#B4B3B4", # T>A
                                     "#9ACA3C", # T>C
                                     "#D7BFD7") # T>G

# Constants for TriNucSubCounts
subtype.levels=c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C",
                 "C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T",
                 "T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C",
                 "A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
                 "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C",
                 "T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T",
                 "C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C",
                 "G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
                 "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C",
                 "C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T",
                 "T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C",
                 "A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T",
                 "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C",
                 "T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
                 "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C",
                 "G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

sub.motif.levels = c("C[G>A]A", "T[C>G]C", "A[G>A]C", "T[C>T]A", "T[G>C]A", "C[G>A]C",
                     "A[A>G]C", "T[G>A]G", "G[G>A]A", "G[C>T]G", "T[G>A]A", "T[C>A]T",
                     "C[C>T]A", "G[G>C]A", "T[C>A]C", "T[C>T]G", "A[G>C]A", "G[C>T]T",
                     "G[A>G]G", "T[C>T]T", "C[C>T]C", "C[G>A]G", "T[T>A]C", "C[G>A]T",
                     "C[C>T]G", "A[G>A]A", "T[G>T]G", "A[G>T]A", "C[C>G]T", "A[T>C]A",
                     "G[G>T]A", "T[C>T]C", "T[C>A]A", "T[C>G]A", "A[T>G]C", "G[G>A]C",
                     "A[G>T]T", "A[G>A]G", "T[C>G]T", "C[C>A]A", "C[T>G]A", "C[A>G]C",
                     "C[C>A]C", "T[C>G]G", "A[C>G]G", "C[C>G]C", "T[G>T]T", "T[G>T]C",
                     "G[G>T]C", "A[T>C]G", "C[T>G]T", "A[C>T]G", "G[T>C]A", "G[A>C]G",
                     "A[G>C]T", "T[G>C]G", "A[C>G]A", "C[G>C]A", "T[G>C]T", "T[G>A]C",
                     "A[T>C]C", "T[G>A]T", "A[T>G]T", "C[T>A]T", "A[A>G]T", "A[T>A]G",
                     "G[G>T]G", "G[G>C]C", "G[C>A]T", "T[G>T]A", "A[C>A]A", "C[C>A]T",
                     "C[T>C]C", "G[G>A]G", "G[A>C]A", "A[T>C]T", "T[A>T]T", "C[A>T]G",
                     "A[C>T]C", "T[T>C]C", "G[C>T]C", "A[A>G]A", "A[C>A]T", "G[T>C]C",
                     "A[G>T]G", "A[T>A]T", "G[G>C]T", "C[G>T]G", "A[C>G]C", "C[G>T]A",
                     "A[C>T]A", "C[A>T]T", "C[G>T]T", "A[G>T]C", "G[C>T]A", "C[C>A]G",
                     "G[C>A]A", "G[C>A]G", "G[T>G]C", "G[G>T]T", "G[A>G]C", "T[C>A]G",
                     "C[C>T]T", "A[G>A]T", "A[A>T]C", "C[T>C]A", "C[C>G]A", "T[T>A]T",
                     "C[T>C]G", "A[G>C]C", "G[C>G]C", "A[C>T]T", "T[G>C]C", "G[C>G]T",
                     "G[C>A]C", "C[T>A]C", "C[A>G]T", "A[A>C]T", "G[G>C]G", "A[G>C]G",
                     "T[T>C]T", "A[A>T]G", "C[G>C]C", "C[T>A]G", "T[T>A]A", "T[A>C]T",
                     "A[C>A]C", "G[T>C]T", "A[T>A]C", "C[A>C]T", "T[A>C]C", "G[A>T]G",
                     "T[A>C]A", "G[A>G]T", "G[A>T]T", "C[T>C]T", "G[T>A]T", "G[A>G]A",
                     "T[T>G]A", "T[T>G]T", "T[A>G]T", "G[A>T]C", "A[A>C]A", "G[G>A]T",
                     "C[A>C]G", "C[A>G]A", "G[T>G]A", "G[A>T]A", "T[A>C]G", "A[C>G]T",
                     "G[T>A]A", "T[A>G]C", "C[G>C]G", "T[T>C]G", "T[A>G]G", "A[A>T]T",
                     "T[T>G]G", "C[A>T]A", "C[A>C]C", "C[G>T]C", "G[T>A]C", "C[A>G]G",
                     "C[A>C]A", "T[A>T]A", "G[T>C]G", "A[A>C]C", "C[A>T]C", "C[T>G]G",
                     "A[A>G]G", "G[C>G]A", "A[A>C]G", "T[A>T]G", "C[G>C]T", "T[T>G]C",
                     "A[T>G]G", "G[T>G]T", "G[T>A]G", "C[C>G]G", "A[C>A]G", "A[A>T]A",
                     "A[T>G]A", "G[A>C]T", "G[T>G]G", "G[A>C]C", "A[T>A]A", "C[T>G]C",
                     "T[A>G]A", "G[C>G]G", "T[T>C]A", "C[T>A]A", "T[T>A]G", "T[A>T]C")

