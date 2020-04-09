#!/usr/bin/env Rscript --vanilla

# FIREVAT Run
#
# Last revised date:
#   February 21, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea

# Scripts for FIREVAT command line usage
# Rscript firevat_run.R --configure-file example.json --vcf query.vcf --output output.vcf --report

# Load FIREVAT package
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(FIREVAT))


# Specify necessary FIREVAT options in firevat.run.option.list
firevat.run.option.list <- list(
    # vcf file
    make_option(c("-v", "--vcf"),
                type="character",
                help="Input vcf file for FIREVAT pipeline"),
    # config file
    make_option(c("-c", "--config"),
                type="character",
                help="Configuration JSON/YAML for FIREVAT pipeline"),
    # output directory
    make_option(c("-o", "--output"),
                type="character",
                default=paste0("../experiment/",format(Sys.time(), "%Y%m%d_%H%M%S/")),
                help="Output directory [default \"%default\"]"),
    # reference genome version
    make_option(c("-g","--genome"),
                type="character",
                help="Genome version of input vcf"),
    # running mode of firevat: v0.1.2 - manual / ga
    make_option(c("-m","--mode"),
                type="character",
                default="ga",
                help="Running mode of FIREVAT (manual or ga)"),
    # annotate
    make_option(c("--annotate"),
                default= TRUE,
                help="Annotate [default]"),
    # perform strand bias analysis
    make_option(c("--strandbias"),
                default= FALSE,
                help="Perform strand bias analysis [default: FALSE]"),
    # number of cores
    make_option(c("-n","--numcores"),
                type="integer",
                default=1,
                dest="num.cores",
                help="Number of cores using in FIREVAT pipeline [default\"%default\"]"),
    # GA: population size
    make_option(c("-p","--popsize"),
                type="integer",
                default= 200,
                dest="pop.size",
                help="Population size for GA algorithm [default\"%default\"]"),
    # GA: maximum iteration
    make_option(c("-i","--maxiter"),
                type="integer",
                default= 50,
                dest="max.iter",
                help="Maximum iteration for GA algorithm [default\"%default\"]"),
    make_option(c("-r","--run"),
                type="integer",
                default= 25,
                help="Limit for consecutive generations without changes [default\"%default\"]"),
    # GA: mutation probability
    make_option("--pmutation",
                type="numeric",
                default= 0.25,
                dest="p.mutation",
                help="Mutation probability for GA algorithm [default\"%default\"]"),
    # verbose / quiet option
    make_option("--verbose",
                action="store_true",
                default= TRUE,
                help="Print logs [default]"),
    make_option("--quiet",
                action="store_false",
                dest="verbose",
                help="Print less logs [default]")
)

parser <- OptionParser(usage="%prog [options]",
                       option_list = firevat.run.option.list)

run.args <- parse_args(parser, args = commandArgs(trailingOnly=TRUE))


if (!interactive()){

    # Execution
    verbose <- run.args$verbose
    if (verbose) print(paste0("Started datetime ", Sys.time()))

    if (run.args$annotate == TRUE) {
        # Annotation DB
        clinvar.vcf.file <- system.file("extdata", "clinvar_hg19_20190212.vcf", package = "FIREVAT")
        clinvar.vcf.obj <- ReadVCF(vcf.file = clinvar.vcf.file, genome = "hg19", split.info = TRUE)
        df.annotation.db <- PrepareAnnotationDB(annotation.vcf.obj = clinvar.vcf.obj)

        # Annotation parameters
        cols.to.display = c("GENEINFO", "CLNSIG")
        filter.key.value.pairs <- list("CLNSIG" = c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"))
    } else {
        df.annotation.db <- NULL
        cols.to.display <- c()
        filter.key.value.pairs <- list()
    }
    # Run FIREVAT
    results <- RunFIREVAT(vcf.file = run.args$vcf,
                          vcf.file.genome = run.args$genome,
                          config.file = run.args$config,
                          df.ref.mut.sigs = GetPCAWGMutSigs(),
                          target.mut.sigs = GetPCAWGMutSigsNames(),
                          sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
                          output.dir = run.args$output,
                          num.cores = run.args$num.cores,
                          ga.pop.size = run.args$pop.size,
                          ga.max.iter = run.args$max.iter,
                          ga.run = run.args$run,
                          ga.pmutation = run.args$p.mutation,
                          mutalisk.method = "all",
                          perform.strand.bias.analysis = run.args$strand_bias,
                          annotate = run.args$annotate,
                          df.annotation.db = df.annotation.db,
                          annotated.columns.to.display = cols.to.display,
                          annotation.filter.key.value.pairs = filter.key.value.pairs,
                          annotation.filter.condition = "AND")

    if (verbose) print(paste0("Finished datetime ", Sys.time()))

    x.solution.decimal <- results$x.solution.decimal
    data <- results$data

    print("Finished.")
    print(x.solution.decimal)
}
