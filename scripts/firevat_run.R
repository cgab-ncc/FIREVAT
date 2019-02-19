#!/usr/bin/env Rscript --vanilla

# FIREVAT Run
#
# Last revised date:
#   February 19, 2019
#
# Authors: 
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea

## Scripts for FIREVAT command line usage
## $FIREVAT_path/firevat_run.R --configure-file example.json --vcf query.vcf --output output.vcf --report

# Load FIREVAT package
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(FIREVAT))

# Specify necessary FIREVAT options in firevat.run.option.list

firevat.run.option.list <- list( 
    # vcf file
    make_option(c("-v", "--vcf"), type="character",
        help="Input vcf file for FIREVAT pipeline"),
    # config file
    make_option(c("-c", "--config"), type="character", 
        help="Configuration JSON/YAML for FIREVAT pipeline"),
    # output directory
    make_option(c("-o", "--output"), type="character",
        default=paste0("../experiment/",format(Sys.time(), "%Y%m%d_%H%M%S/")),
        help="Output directory [default \"%default\"]"),
    # reference genome version
    make_option(c("-g","--genome"), type="character",
        help="Genome version of input vcf (hg19/hg38)"),
    # number of cores
    make_option(c("-n","--numcores"), type="integer", default=1, dest="num.cores",
        help="Number of cores using in FIREVAT pipeline [default\"%default\"]"),
    # GA: population size
    make_option(c("-p","--popsize"), type="integer", 
        default= 200,dest="pop.size",
        help="Population size for GA algorithm [default\"%default\"]"),
    # GA: maximum iteration
    make_option(c("-i","--maxiter"), type="integer", default= 50, dest="max.iter",
        help="Maximum iteration for GA algorithm [default\"%default\"]"),
    make_option(c("-r","--run"), type="integer", default= 25,
        help="Limit for consecutive generations without changes [default\"%default\"]"),
    # GA: mutation probability
    make_option(c("-m","--pmutation"), type="numeric",
        default= 0.25, dest="p.mutation",
        help="Mutation probability for GA algorithm [default\"%default\"]"),
    # verbose / quiet option
    make_option("--verbose", action="store_true", default= TRUE,
        help="Print logs [default]"),
    make_option("--quiet", action="store_false", dest="verbose",
        help="Print less logs [default]")
)

parser <- OptionParser(usage="%prog [options]",
                       option_list = firevat.run.option.list)

firevat.run.args <- parse_args(parser, args = commandArgs(trailingOnly=TRUE))


if (!interactive()){

    # Execution
    verbose <- firevat.run.args$verbose
    if (verbose) print(paste0("Started datetime ", Sys.time()))
    df.ref.mut.sigs <- GetPCAWGMutSigs()
    target.mut.sigs <- GetPCAWGMutSigsNames()
    results <- RunFIREVAT(
        vcf.file = firevat.run.args$vcf,
        vcf.file.genome = firevat.run.args$genome,
        config.file = firevat.run.args$config,
        df.ref.mut.sigs = df.ref.mut.sigs,
        target.mut.sigs = target.mut.sigs,
        sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
        num.cores = firevat.run.args$num.cores,
        output.dir = firevat.run.args$output,
        pop.size = firevat.run.args$pop.size,
        max.iter = firevat.run.args$max.iter,
        run = firevat.run.args$run,
        pmutation = firevat.run.args$p.mutation,
        verbose = verbose
        )
    if (verbose) print(paste0("Finished datetime ", Sys.time()))

    x.solution.decimal <- results$x.solution.decimal
    data <- results$data

    print("Finished.")
    print(x.solution.decimal)
}