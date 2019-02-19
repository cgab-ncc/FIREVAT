# FIREVAT Main Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title RunFIREVAT
#' @description
#' Runs FIREVAT using configuration data.Filters point mutations in the specified vcf.file based on mutational signature
#' decomposition and outputs the refined and artifact vcf as well as metadata related to
#' the refinement process. 
#'
#' @param vcf.file String value corresponding to input .vcf file (full path)
#' @param vcf.file.genome Genome assembly of the input .vcf file
#' @param config.file String value corresponding to input configuration file (refer to ...)
#' @param df.ref.mut.sigs A data.frame of the reference mutational signatures
#' @param target.mut.sigs A character vector of the target mutational signatures (from reference mutational signatures)
#' @param sequencing.artifact.mut.sigs A character vector of the sequencing artifact mutational signatures (from reference mutational signatures)
#' @param num.cores Number of cores to allocate
#' @param output.dir String value of the desired output directory
#' @param pop.size Integer value of the Genetic Algorithm 'population size' parameter
#' @param max.iter Integer value of the Genetic Algorithm 'maximum iterations' parameter
#' @param run Integer value of the Genetic Algorithm 'run' parameter
#' @param pmutation Float value of the Genetic Algorithm 'mutation probability' parameter
#' @param verbose If true, provides process detail
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{x.solution.decimal}{A numeric vector of optimized parameter values}
#'  \item{data}{}
#' }
#'
#' @importFrom GA ga
#' @export
RunFIREVAT <- function(vcf.file,
                       vcf.file.genome,
                       config.file,
                       df.ref.mut.sigs,
                       target.mut.sigs,
                       sequencing.artifact.mut.sigs,
                       num.cores,
                       output.dir,
                       pop.size = 200,
                       max.iter = 200,
                       run = 50,
                       pmutation = 0.25,
                       verbose = TRUE) {
    # Check input parameters
    if (is.character(vcf.file) == F) {
        stop("The parameter 'vcf.file' must be a string")
    }
    if (vcf.file.genome != 'hg19' && vcf.file.genome != 'hg38') {
        stop("The parameter 'vcf.file.genome' must be either 'hg19' or 'hg38'")
    }
    if (is.character(config.file) == F) {
        stop("The parameter 'config.file' must be a string")
    }
    if (is.data.frame(df.ref.mut.sigs) == F) {
        stop("The parameter 'df.ref.mut.sigs' must be a data.frame")
    }

    start.datetime <- Sys.time()
    if (verbose) {
        print(start.datetime)
        print("Initializing FIREVAT variant filtering pipeline")
    }

    # Create the output directory
    if (!dir.exists(output.dir)) {
        dir.create(output.dir, recursive = T)
    }
    else {
        print("output.dir already exists")
    }

    # Read the vcf file
    vcf.obj <- ReadVCF(vcf.file = vcf.file, genome = vcf.file.genome)

    # Read the config file
    config.obj <- ParseConfigFile(config.file, verbose = verbose)

    # Initially select point mutations that satisfy the config file and
    # consider these the initial set of mutations
    vcf.objs <- InitializeVCFFromConfig(vcf.obj = vcf.obj,
                                        config.obj = config.obj,
                                        verbose = verbose)
    vcf.obj <- vcf.objs$vcf.obj.filtered

    # Make filter from config file
    vcf.filter <- MakeFilterFromConfig(config.obj)

    # FIREVAT can only be run if there are more than 50 point mutations in the initial vcf file
    if (nrow(vcf.obj$data) <= 50) {
        print("FIREVAT must have at least 50 mutations to run.")
        return()
    }

    print(paste0("Starting with ", nrow(vcf.obj$data), " point mutations"))

    bits.list <- ParameterToBits(vcf.obj, config.obj, vcf.filter)
    params.bit.len <- bits.list$params.bit.len
    vcf.obj <- bits.list$vcf.obj

    n.bits <- sum(params.bit.len)

    # Prepare data for optimization
    data = list(n.bits = n.bits,
                params.bit.len = params.bit.len,
                vcf.file = vcf.file,
                vcf.file.basename = gsub("\\.vcf", "", basename(vcf.file)),
                vcf.obj = vcf.obj,
                config.file = config.file,
                config.obj = config.obj,
                vcf.filter = vcf.filter,
                df.ref.mut.sigs = df.ref.mut.sigs,
                df.ref.mut.sigs.cosmic = GetCOSMICMutSigs(),
                target.mut.sigs = target.mut.sigs,
                sequencing.artifact.mut.sigs = sequencing.artifact.mut.sigs,
                output.dir = output.dir,
                start.datetime = start.datetime,
                pmutation = pmutation,
                pop.size = pop.size,
                max.iter = max.iter,
                run = run)

    # Optimize filter parameters
    print("Started running optimization part")
    ga.results <- ga(type = "binary",
                     fitness =  function(string) GAOptimizationObjFn(string, data),
                     nBits = data$n.bits,
                     popSize = pop.size,
                     maxiter = max.iter,
                     parallel = num.cores,
                     run = run,
                     pmutation = pmutation,
                     monitor = function(obj) GAMonitorFn(obj, data),
                     keepBest = T)
    print("Finished running optimization part")

    # Parse optimized results
    x.solution.binary <- as.numeric(ga.results@solution[1,])
    x.solution.decimal <- GADecodeBinaryString(string = x.solution.binary,
                                               data = data)

    print("FIREVAT results")
    print(summary(ga.results))

    print("FIREVAT optimized binary values of filter parameters")
    print(x.solution.binary)

    print("FIREVAT Optimized decimal values of filter parameters")
    print(x.solution.decimal)

    data$end.datetime <- Sys.time()

    # Report results
    ReportFIREVATResults(x.solution.decimal = x.solution.decimal, data = data)

    return(list(x.solution.decimal = x.solution.decimal, data = data))
}