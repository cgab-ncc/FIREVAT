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
#' @param vcf.file String value corresponding to input .vcf file. Please provide the full path.
#' @param vcf.file.genome Genome assembly of the input .vcf file. The value should be eitehr 'hg19' or 'hg38'.
#' @param config.file String value corresponding to input configuration file. For more details please refer to ...
#' @param df.ref.mut.sigs A data.frame of the reference mutational signatures
#' @param target.mut.sigs A character vector of the target mutational signatures from reference mutational signatures.
#' @param sequencing.artifact.mut.sigs A character vector of the sequencing artifact mutational signatures from reference mutational signatures.
#' @param num.cores Number of cores to allocate
#' @param output.dir String value of the desired output directory
#' @param mode String value. The value should be either 'ga' or 'manual'.
#' @param objective.fn Objective value derivation function. Default: Default.Obj.Fn.
#' @param use.suggested.soln Boolean value. If TRUE, then FIREVAT passes the default values
#' of filter variables declared as 'use_in_filter' in the config file to the 'suggestions' parameter of
#' the Genetic Algorithm package. If FALSE, then FIREVAT supplies NULL to the GA package 'suggestions' parameter.
#' @param ga.pop.size Integer value of the Genetic Algorithm 'population size' parameter. Default: 200.
#' This value should be set based on the number of filter parameters. Recommendation: 40 per filter parameter.
#' @param ga.max.iter Integer value of the Genetic Algorithm 'maximum iterations' parameter. Ddefault: 200.
#' This value should be set based on the number of filter parameters. Recommendation: same as 'ga.pop.size'.
#' @param ga.run Integer value of the Genetic Algorithm 'run' parameter. Default: 50.
#' This value should be set based on the 'ga.max.iter' parameter. Recommendation: 25 percent of 'ga.max.iter'.
#' @param ga.pmutation Float value of the Genetic Algorithm 'mutation probability' parameter. Default: 0.25.
#' @param mutalisk.method Mutalisk signature identification method. Default: 'random.sampling'.
#' The value can be either 'all' or 'random.sampling'.
#' 'all' uses all target.mut.sigs to identify mutational signatures.
#' 'random.sampling' randomly samples from target.mut.sigs to identify mutational signatures.
#' @param mutalisk.random.sampling.count Mutalisk random sampling count. Default: 20.
#' The number of signatures to sample from target.mut.sigs
#' @param mutalisk.random.sampling.max.iter Mutalisk random sampling maximum iteration. Default: 10.
#' The number of times Mutalisk randomly samples from target.mut.sigs before determining the candidate signatures.
#' @param report.format The format of FIREVAT report. We currently only support 'html'.
#' @param perform.strand.bias.analysis If TRUE, then performs strand bias analysis.
#' @param strand.bias.perform.fdr.correction If TRUE, then performs false discovery rate
#' correction for strand bias analysis.
#' @param strand.bias.fdr.correction.method A string value. Default value is 'BH'.
#' Refer to 'p.adjust()' function method.
#' @param ref.forward.strand.var A string value.
#' @param ref.reverse.strand.var A string value,
#' @param alt.forward.strand.var A string value,
#' @param alt.reverse.strand.var A string value,
#' @param verbose If TRUE, provides process detail. Default: TRUE.
#' @param annotate A boolean value. Default value is TRUE.
#' @param df.annotation.db A data.frame. Please refer to \code{\link{PrepareAnnotationDB}}
#' @param annotated.columns.to.display A character vector.
#' @param annotation.filter.key.value.pairs A list.
#' @param annotation.filter.condition 'AND' or 'OR'.
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{f}{ = A ggarrange object}
#'  \item{graphs}{ = A list of length 3; each element is a ggplot histogram}
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
                       mode = "ga",
                       objective.fn = Default.Obj.Fn,
                       use.suggested.soln = TRUE,
                       ga.pop.size = 200,
                       ga.max.iter = 200,
                       ga.run = 50,
                       ga.pmutation = 0.25,
                       mutalisk.method = "all",
                       mutalisk.random.sampling.count = 20,
                       mutalisk.random.sampling.max.iter = 10,
                       perform.strand.bias.analysis = TRUE,
                       strand.bias.perform.fdr.correction = TRUE,
                       strand.bias.fdr.correction.method = "BH",
                       ref.forward.strand.var = NULL,
                       ref.reverse.strand.var = NULL,
                       alt.forward.strand.var = NULL,
                       alt.reverse.strand.var = NULL,
                       report.format = "html",
                       annotate = TRUE,
                       df.annotation.db = NULL,
                       annotated.columns.to.display = NULL,
                       annotation.filter.key.value.pairs = NULL,
                       annotation.filter.condition = "AND",
                       verbose = TRUE) {
    # Check input parameters
    if (is.character(vcf.file) == FALSE) {
        stop("The parameter 'vcf.file' must be a string")
    }
    if (vcf.file.genome != 'hg19' && vcf.file.genome != 'hg38') {
        stop("The parameter 'vcf.file.genome' must be either 'hg19' or 'hg38'")
    }
    if (is.character(config.file) == FALSE) {
        stop("The parameter 'config.file' must be a string")
    }
    if (is.data.frame(df.ref.mut.sigs) == FALSE) {
        stop("The parameter 'df.ref.mut.sigs' must be a data.frame")
    }
    if (all(target.mut.sigs %in% colnames(df.ref.mut.sigs)) == FALSE)  {
        stop("The parameter 'target.mut.sigs' must be present in df.ref.mut.sigs columns")
    }
    if (all(sequencing.artifact.mut.sigs %in% colnames(df.ref.mut.sigs)) == FALSE)  {
        stop("The parameter 'sequencing.artifact.mut.sigs' must be present in df.ref.mut.sigs columns")
    }
    if (is.numeric(num.cores) == FALSE || num.cores <= 0)  {
        stop("The parameter 'num.cores' must be an integer and be greater than 0")
    }
    if (is.character(output.dir) == FALSE)  {
        stop("The parameter 'output.dir' must be a string")
    }
    if (mode != "ga" && mode != "manual")  {
        stop("The parameter 'mode' must be either 'ga' or 'manual'")
    }
    if (is.numeric(ga.pop.size) == FALSE || ga.pop.size <= 0)  {
        stop("The parameter 'ga.pop.size' must be an integer and be greater than 0")
    }
    if (is.numeric(ga.max.iter) == FALSE || ga.max.iter <= 0)  {
        stop("The parameter 'ga.max.iter' must be an integer and be greater than 0")
    }
    if (is.numeric(ga.run) == FALSE || ga.run <= 0)  {
        stop("The parameter 'ga.run' must be an integer and be greater than 0")
    }
    if (is.numeric(ga.pmutation) == FALSE || ga.pmutation < 0 || ga.pmutation > 1)  {
        stop("The parameter 'ga.pmutation' must be a value between 0 and 1")
    }
    if (verbose != TRUE && verbose != FALSE)  {
        stop("The parameter 'verbose' must be a boolean")
    }

    start.datetime <- Sys.time()
    if (verbose == TRUE) {
        print(start.datetime)
        print("Initializing FIREVAT variant filtering pipeline")
    }

    # Create the output directory
    if (dir.exists(output.dir) == FALSE) {
        dir.create(output.dir, recursive = T)
    } else {
        if (verbose == TRUE) {
            print(paste0(output.dir, " already exists"))
        }
    }

    # 1. Read the vcf file
    vcf.obj <- ReadVCF(vcf.file = vcf.file, genome = vcf.file.genome)

    # 2. Read the config file
    config.obj <- ParseConfigFile(config.file, verbose = verbose)

    # 3.Initially select point mutations that satisfy the config file and
    #   consider these the initial set of mutations
    vcf.objs <- InitializeVCF(vcf.obj = vcf.obj,
                              config.obj = config.obj,
                              verbose = verbose)
    vcf.obj <- vcf.objs$vcf.obj.filtered
    # FIREVAT can only be run if there are more than 50 point mutations in the initial vcf file
    if (nrow(vcf.obj$data) <= 50) {
        warning("FIREVAT must have at least 50 mutations to run. Returning without running FIREVAT.")
        return()
    } else {
        if (verbose == TRUE) {
            print(paste0("Starting with ", nrow(vcf.obj$data), " point mutations"))
        }
    }

    # 4. Make filter from config file
    vcf.filter <- MakeFilter(config.obj)

    # 5. Convert filter parameters to bits
    bits.list <- ParameterToBits(vcf.obj, config.obj, vcf.filter)
    params.bit.len <- bits.list$params.bit.len
    vcf.obj <- bits.list$vcf.obj
    n.bits <- sum(params.bit.len)

    # 6. Prepare data for optimization
    data = list(start.datetime = start.datetime,
                objective.fn = objective.fn,
                n.bits = n.bits,
                params.bit.len = params.bit.len,
                vcf.file = vcf.file,
                vcf.file.basename = gsub("\\.vcf", "", basename(vcf.file)),
                vcf.obj = vcf.obj,
                config.file = config.file,
                config.obj = config.obj,
                vcf.filter = vcf.filter,
                df.ref.mut.sigs = df.ref.mut.sigs,
                target.mut.sigs = target.mut.sigs,
                sequencing.artifact.mut.sigs = sequencing.artifact.mut.sigs,
                output.dir = output.dir,
                ga.pmutation = ga.pmutation,
                ga.pop.size = ga.pop.size,
                ga.max.iter = ga.max.iter,
                ga.run = ga.run,
                mutalisk.method = mutalisk.method,
                mutalisk.random.sampling.count = mutalisk.random.sampling.count,
                mutalisk.random.sampling.max.iter = mutalisk.random.sampling.max.iter,
                report.format = report.format,
                mode = mode,
                use.suggested.soln = use.suggested.soln,
                num.cores = num.cores,
                verbose = verbose)

    # 7. Optimize filter parameters
    if (mode == "ga") {
        if (verbose == TRUE) {
            print("Started running FIREVAT Genetic Algorithm optimization")
        }
        # Suggestions
        if (use.suggested.soln == TRUE) { # use default parameters as suggestions
            suggestions =  NULL # TODO use the real GetDefaultValues() function
        } else {
            suggestions = NULL
        }

        # Prepare data for Mutational Patterns (memoization)
        data$df.mut.pat.ref.sigs <- MutPatParseRefMutSigs(df.ref.mut.sigs = df.ref.mut.sigs,
                                                          target.mut.sigs = target.mut.sigs)
        if (data$vcf.obj$genome == "hg19") {
            data$bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
        }
        if (data$vcf.obj$genome == "hg38") {
            data$bsg <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        }

        ga.results <- ga(type = "binary",
                         fitness =  function(string) GAOptimizationObjFn(string, data),
                         nBits = n.bits,
                         popSize = ga.pop.size,
                         maxiter = ga.max.iter,
                         parallel = num.cores,
                         run = ga.run,
                         pmutation = ga.pmutation,
                         suggestions = suggestions,
                         monitor = function(obj) GAMonitorFn(obj, data),
                         keepBest = TRUE)
        if (verbose == TRUE) {
            print("Finished running FIREVAT Genetic Algorithm optimization")
        }

        # Parse optimized results
        data$x.solution.binary <- as.numeric(ga.results@solution[1,])
        data$x.solution.decimal <- GADecodeBinaryString(string = data$x.solution.binary,
                                                        data = data)
        if (verbose == TRUE) {
            print("FIREVAT results")
            print(summary(ga.results))
            print("FIREVAT optimized binary values of filter parameters")
            print(data$x.solution.binary)
            print("FIREVAT Optimized decimal values of filter parameters")
            print(data$x.solution.decimal)
        }
    } else if (mode == "manual") {
        if (verbose == TRUE) {
            print("Running FIREVAT manual optimization")
        }
        data$x.solution.decimal <- unlist(data$vcf.filter) # default values
    }

    # 8. Filter vcf.data based on the updated vcf.filter
    data$vcf.filter <- UpdateFilter(vcf.filter = data$vcf.filter,
                                    param.values = data$x.solution.decimal)
    optimized.vcf.objs <- FilterVCF(vcf.obj = data$vcf.obj,
                                    config.obj = data$config.obj,
                                    vcf.filter = data$vcf.filter,
                                    verbose = data$verbose)

    data$refined.vcf.obj <- optimized.vcf.objs$vcf.obj.filtered
    data$artifactual.vcf.obj <- optimized.vcf.objs$vcf.obj.artifact

    # 9. Strand bias analysis
    data$perform.strand.bias.analysis <- perform.strand.bias.analysis
    data$ref.forward.strand.var <- ref.forward.strand.var
    data$ref.reverse.strand.var <- ref.reverse.strand.var
    data$alt.forward.strand.var <- alt.forward.strand.var
    data$alt.reverse.strand.var <- alt.reverse.strand.var
    data$strand.bias.fdr.correction.method <- strand.bias.fdr.correction.method
    data$strand.bias.perform.fdr.correction <- strand.bias.perform.fdr.correction

    if (data$perform.strand.bias.analysis == TRUE) {
        data$refined.vcf.obj <- PerformStrandBiasAnalysis(
            vcf.obj = data$refined.vcf.obj,
            ref.forward.strand.var = data$ref.forward.strand.var,
            ref.reverse.strand.var = data$ref.reverse.strand.var,
            alt.forward.strand.var = data$alt.forward.strand.var,
            alt.reverse.strand.var = data$alt.reverse.strand.var,
            perform.fdr.correction = data$strand.bias.perform.fdr.correction,
            fdr.correction.method = data$strand.bias.fdr.correction.method)
        data$artifactual.vcf.obj <- PerformStrandBiasAnalysis(
            vcf.obj = data$artifactual.vcf.obj,
            ref.forward.strand.var = data$ref.forward.strand.var,
            ref.reverse.strand.var = data$ref.reverse.strand.var,
            alt.forward.strand.var = data$alt.forward.strand.var,
            alt.reverse.strand.var = data$alt.reverse.strand.var,
            perform.fdr.correction = data$strand.bias.perform.fdr.correction,
            fdr.correction.method = data$strand.bias.fdr.correction.method)
    }

    # 10. Annotate
    data$annotate = annotate
    data$df.annotation.db = df.annotation.db
    data$annotated.columns.to.display = annotated.columns.to.display
    data$annotation.filter.key.value.pairs = annotation.filter.key.value.pairs
    data$annotation.filter.condition = annotation.filter.condition

    if (data$annotate == TRUE) {
        # Annotate VCFs
        data$vcf.obj.annotated <- AnnotateVCFObj(vcf.obj = data$vcf.obj,
                                                 df.annotation.db = data$df.annotation.db,
                                                 include.all.columns = TRUE)
        data$refined.vcf.obj.annotated <- AnnotateVCFObj(vcf.obj = data$refined.vcf.obj,
                                                         df.annotation.db = data$df.annotation.db,
                                                         include.all.columns = TRUE)
        data$artifactual.vcf.obj.annotated <- AnnotateVCFObj(vcf.obj = data$artifactual.vcf.obj,
                                                             df.annotation.db = data$df.annotation.db,
                                                             include.all.columns = TRUE)

        # Query annotated VCFs
        data$refined.vcf.obj.annotated.queried <- QueryAnnotatedVCF(vcf.obj.annotated = data$refined.vcf.obj.annotated,
                                                                    filter.key.value.pairs = data$annotation.filter.key.value.pairs,
                                                                    filter.condition = data$annotation.filter.condition)
        data$artifactual.vcf.obj.annotated.queried <- QueryAnnotatedVCF(vcf.obj.annotated = data$artifactual.vcf.obj.annotated,
                                                                        filter.key.value.pairs = data$annotation.filter.key.value.pairs,
                                                                        filter.condition = data$annotation.filter.condition)
    }

    # 11. Run Mutalisk
    # Identify target mutational signatures
    # Here we fetch all signatures ever identified by Mutational Patterns
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)
    Split.Sigs <- function(sigs, weights, cutoff = 0.05) {
        include <- rep(TRUE, length(sigs))
        include[which(sigs == "")] <- FALSE
        sigs <- as.character(sigs[include])
        weights <- as.character(weights[include])
        sigs <- lapply(sigs, function(x) strsplit(x, ',')[[1]])
        weights <- lapply(weights, function(x) strsplit(x, ',')[[1]])

        candidate.sigs <- c()
        for (i in 1:length(sigs)) {
            df <- data.frame(sig = sigs[[i]],
                             weight = weights[[i]],
                             stringsAsFactors = F)
            df <- df[df$weight >= cutoff, ]
            candidate.sigs <- c(candidate.sigs, df$sig)
        }
        return(candidate.sigs)
    }
    sigs1 <- Split.Sigs(sigs = df.optimization.logs$refined.muts.target.signatures,
                        weights = df.optimization.logs$refined.muts.target.signatures.weights)
    sigs2 <- Split.Sigs(sigs = df.optimization.logs$refined.muts.sequencing.artifact.signatures,
                        weights = df.optimization.logs$refined.muts.sequencing.artifact.signatures.weights)
    sigs3 <- Split.Sigs(sigs = df.optimization.logs$artifactual.muts.target.signatures,
                        weights = df.optimization.logs$artifactual.muts.target.signatures.weights)
    sigs4 <- Split.Sigs(sigs = df.optimization.logs$artifactual.muts.sequencing.artifact.signatures,
                        weights = df.optimization.logs$artifactual.muts.sequencing.artifact.signatures.weights)
    data$mut.pat.target.sigs <- unique(c(sigs1, sigs2, sigs3, sigs4, data$sequencing.artifact.mut.sigs))

    # Original vcf
    data$raw.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$vcf.obj,
                                                  df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                  target.mut.sigs = data$mut.pat.target.sigs,
                                                  method = data$mutalisk.method,
                                                  n.sample = data$mutalisk.random.sampling.count,
                                                  n.iter = data$mutalisk.random.sampling.max.iter,
                                                  verbose = data$verbose)
    # Refined vcf
    data$refined.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$refined.vcf.obj,
                                                      df.ref.mut.sigs = df.ref.mut.sigs,
                                                      target.mut.sigs = data$mut.pat.target.sigs,
                                                      method = data$mutalisk.method,
                                                      n.sample = data$mutalisk.random.sampling.count,
                                                      n.iter = data$mutalisk.random.sampling.max.iter,
                                                      verbose = data$verbose)
    # Artifact vcf
    data$artifactual.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$artifactual.vcf.obj,
                                                          df.ref.mut.sigs = df.ref.mut.sigs,
                                                          target.mut.sigs = data$mut.pat.target.sigs,
                                                          method = data$mutalisk.method,
                                                          n.sample = data$mutalisk.random.sampling.count,
                                                          n.iter = data$mutalisk.random.sampling.max.iter,
                                                          verbose = data$verbose)

    # 12. Write VCF Files
    WriteVCF(vcf.obj = data$vcf.obj,
             save.file = paste0(data$output.dir, data$vcf.file.basename, "_Original.vcf"))
    WriteVCF(vcf.obj = data$refined.vcf.obj,
             save.file = paste0(data$output.dir, data$vcf.file.basename, "_Refined.vcf"))
    WriteVCF(vcf.obj = data$artifactual.vcf.obj,
             save.file = paste0(data$output.dir, data$vcf.file.basename, "_Artifact.vcf"))
    if (data$annotate == TRUE) {
        WriteVCF(vcf.obj = data$data$vcf.obj.annotated,
                 save.file = paste0(data$output.dir, data$vcf.file.basename, "_Original_Annotated.vcf"))
        WriteVCF(vcf.obj = data$refined.vcf.obj.annotated,
                 save.file = paste0(data$output.dir, data$vcf.file.basename, "_Refined_Annotated.vcf"))
        WriteVCF(vcf.obj = data$artifactual.vcf.obj.annotated,
                 save.file = paste0(data$output.dir, data$vcf.file.basename, "_Artifact_Annotated.vcf"))
    }

    data$end.datetime <- Sys.time()

    # 13. Report results
    data <- ReportFIREVATResults(data = data)

    # 14. Save data
    save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))

    return(data)
}
