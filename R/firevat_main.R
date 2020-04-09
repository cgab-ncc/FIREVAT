# FIREVAT Main Functions
#
# Last revised date:
#   Oct 08, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title RunGAMode
#' @description Runs FIREVAT ga mode
#'
#' @param data A list from RunFIREVAT
#'
#' @return A list
#'
#' @importFrom BSgenome getBSgenome
#' @export
RunGAMode <- function(data) {
    if (data$verbose == TRUE) {
        PrintLog("Step 02. Run FIREVAT variant refinement optimization.")
    }

    # Prepare data for Mutational Patterns (memoization)
    data$df.mut.pat.ref.sigs <- MutPatParseRefMutSigs(df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                      target.mut.sigs = data$target.mut.sigs)

    data$bsg <- BSgenome::getBSgenome(data$vcf.obj$genome)

    # 02-1. Check if refinement is necessary based on the sum of sequencing artifact weights in the original VCF file
    if (data$verbose == TRUE) {
        PrintLog("Step 02-1. Check if FIREVAT variant refinement is necessary [firevat_optimization::CheckIfVariantRefinementIsNecessary]")
    }

    is.variant.refinement.necessary <- CheckIfVariantRefinementIsNecessary(vcf.obj = data$vcf.obj,
                                                                           bsg = data$bsg,
                                                                           df.mut.pat.ref.sigs = data$df.mut.pat.ref.sigs,
                                                                           target.mut.sigs = data$target.mut.sigs,
                                                                           sequencing.artifact.mut.sigs = data$sequencing.artifact.mut.sigs,
                                                                           init.artifact.stop = data$init.artifact.stop,
                                                                           verbose = data$verbose)
    data$original.muts.seq.art.weights.sum <- is.variant.refinement.necessary$seq.art.sigs.weights.sum

    if (is.variant.refinement.necessary$judgment == TRUE) {
        PrintLog(paste0("* Sum of sequencing artifact weights: ", is.variant.refinement.necessary$seq.art.sigs.weights.sum))
        PrintLog(paste0("** This value is greater than 'init.artifact.stop' (", data$init.artifact.stop, ")"))
        PrintLog("** FIREVAT will now begin performing variant refinement.")
        data$variant.refinement.performed <- TRUE
    } else {
        PrintLog(paste0("* Sum of sequencing artifact weights: ", is.variant.refinement.necessary$seq.art.sigs.weights.sum))
        PrintLog(paste0("** This value is equal to or smaller than 'init.artifact.stop' (", data$init.artifact.stop, ")"))
        PrintLog("** FIREVAT will return without performing variant refinement.", type = "WARNING")
        data$variant.refinement.performed <- FALSE
        data$end.datetime <- Sys.time()
        data$variant.refinement.terminiation.log <- "Initial sequencing artifact weights sum is less than or equal to init.artifact.stop"
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }

    # Get lower / upper vectors from config file
    lower.upper.list <- GetParameterLowerUpperVector(data$vcf.obj, data$config.obj, data$vcf.filter)
    if (lower.upper.list$valid == FALSE) {
        data$end.datetime <- Sys.time()
        data$variant.refinement.performed <- FALSE
        data$variant.refinement.terminiation.log <- "Config file parameter filter range is invalid."
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }
    data$vcf.obj <- lower.upper.list$vcf.obj
    lower.vector <- lower.upper.list$lower.vector
    upper.vector <- lower.upper.list$upper.vector
    data$lower.vector <- lower.vector
    data$upper.vector <- upper.vector

    # 02-2. Generate candidate suggestions
    if (data$use.suggested.soln == TRUE) { # use default parameters as suggestions
        PrintLog("Step 02-2. Compute suggested solutions [firevat_brute_force::GetGASuggestedSolutions]")
        suggested.solutions <- GetGASuggestedSolutions(vcf.obj = data$vcf.obj,
                                                       bsg = data$bsg,
                                                       config.obj = data$config.obj,
                                                       lower.upper.list = lower.upper.list,
                                                       df.mut.pat.ref.sigs = data$df.mut.pat.ref.sigs,
                                                       target.mut.sigs = data$target.mut.sigs,
                                                       sequencing.artifact.mut.sigs = data$sequencing.artifact.mut.sigs,
                                                       objective.fn = data$objective.fn,
                                                       original.muts.seq.art.weights.sum = data$original.muts.seq.art.weights.sum,
                                                       ga.preemptive.killing = data$ga.preemptive.killing,
                                                       verbose = data$verbose)
        data$df.suggested.solutions <- suggested.solutions$df.suggested.solutions
        data$suggested.solutions.matrix <- suggested.solutions$suggested.solutions.matrix
        suggestions <- data$suggested.solutions.matrix

        # Write to log file
        log.file <- paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_Optimization_Logs.tsv")
        write.table(data$df.suggested.solutions, log.file, row.names = F, sep = "\t")
    } else {
        suggestions <- NULL
    }

    # 02-3. Run variant refinement optimization
    PrintLog("Step 02-3. Run variant optimization guided by mutational signatures [firevat_optimization::GAOptimizationObjFn]")
    if (data$ga.type == "binary") {
        # Convert filter parameters to bits
        bits.list <- ParameterToBits(data$vcf.obj, data$config.obj, data$vcf.filter)
        params.bit.len <- bits.list$params.bit.len
        data$vcf.obj <- bits.list$vcf.obj
        n.bits <- sum(params.bit.len)
        data$n.bits <- n.bits
        data$params.bit.len <- params.bit.len
        # Run GA with binary type
        ga.results <- ga(type = "binary",
                         fitness =  function(string) GAOptimizationObjFn(string, data),
                         nBits = n.bits,
                         popSize = data$ga.pop.size,
                         maxiter = data$ga.max.iter,
                         parallel = data$num.cores,
                         run = data$ga.run,
                         pmutation = data$ga.pmutation,
                         suggestions = NULL,
                         monitor = function(obj) GAMonitorFn(obj, data),
                         keepBest = TRUE,
                         seed = data$ga.seed)

    } else if (data$ga.type == "real-valued") {
        # Run GA with real-valued type
        if (packageVersion("GA") >= '3.1') {
            ga.results <- ga(type = "real-valued",
                             fitness =  function(cutoffs) GAOptimizationObjFn(cutoffs, data),
                             lower = lower.vector, # <-- vector from config
                             upper = upper.vector, # <-- vector from config
                             popSize = data$ga.pop.size,
                             maxiter = data$ga.max.iter,
                             parallel = data$num.cores,
                             run = data$ga.run,
                             pmutation = data$ga.pmutation,
                             suggestions = suggestions,
                             monitor = function(obj) GAMonitorFn(obj, data),
                             keepBest = TRUE,
                             seed = data$ga.seed)
        } else {
            ga.results <- ga(type = "real-valued",
                             fitness =  function(cutoffs) GAOptimizationObjFn(cutoffs, data),
                             min = lower.vector, # <-- vector from config
                             max = upper.vector, # <-- vector from config
                             popSize = data$ga.pop.size,
                             maxiter = data$ga.max.iter,
                             parallel = data$num.cores,
                             run = data$ga.run,
                             pmutation = data$ga.pmutation,
                             suggestions = suggestions,
                             monitor = function(obj) GAMonitorFn(obj, data),
                             keepBest = TRUE,
                             seed = data$ga.seed)
        }
    }

    # Parse optimized results
    if (data$ga.type == "binary") {
        data$x.solution.binary <- as.numeric(ga.results@solution[1,])
        data$x.solution.decimal <- GADecodeBinaryString(string = data$x.solution.binary,
                                                        data = data)
    } else if (data$ga.type == "real-valued"){
        data$x.solution.decimal <- floor(as.numeric(ga.results@solution[1,]))
        names(data$x.solution.decimal) <- names(data$vcf.filter)
    }

    if (data$verbose == TRUE) {
        print("FIREVAT results")
        print(summary(ga.results))
        if (data$ga.type == "binary") {
            PrintLog("* FIREVAT optimized binary values of filter parameters:")
            print(data$x.solution.binary)
        }
        PrintLog("* FIREVAT optimized integer values of filter parameters:")
        print(data$x.solution.decimal)
    }

    # 02-4. Filter vcf.data based on the updated vcf.filter
    PrintLog("Step 02-4. Filter VCF based on optmized filter parameters.")
    data$vcf.filter <- UpdateFilter(vcf.filter = data$vcf.filter,
                                    param.values = data$x.solution.decimal)
    optimized.vcf.objs <- FilterVCF(vcf.obj = data$vcf.obj,
                                    config.obj = data$config.obj,
                                    vcf.filter = data$vcf.filter,
                                    verbose = data$verbose)
    data$refined.vcf.obj <- optimized.vcf.objs$vcf.obj.filtered
    data$artifactual.vcf.obj <- optimized.vcf.objs$vcf.obj.artifact

    # Check if point mutations remaining in either refined.vcf.obj or artifactual.vcf.obj.
    # Either refinement is too liberal or too stringent. Return with a message.
    if (nrow(data$refined.vcf.obj$data) == 0) {
        PrintLog("After performing variant refinement there are no mutations remaining in the refined set.", type = "WARNING")
        data$end.datetime <- Sys.time()
        data$variant.refinement.terminiation.log <- "Successful but after performing variant refinement there are no mutations remaining in the refined set"
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }
    if (nrow(data$artifactual.vcf.obj$data) == 0) {
        PrintLog("After performing variant refinement there are no mutations remaining in the artifactual set.", type = "WARNING")
        data$end.datetime <- Sys.time()
        data$variant.refinement.terminiation.log <- "Successful but after performing variant refinement there are no mutations remaining in the artifactual set"
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }

    # Check if variant refinement failed for all iterations
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)
    if (tail(df.optimization.logs$iteration.ran.successfully, 1) == FALSE) {
        PrintLog("Variant refinement is unsuccessful because there are not enough mutations.", type = "WARNING")
        data$end.datetime <- Sys.time()
        data$variant.refinement.terminiation.log <- "Unsuccessful because there are not enough mutations for variant refinement"
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }

    # 03. Additional analysis
    # 03-1. Strand bias analysis
    PrintLog("Step 03. Additional analysis.")
    if (data$perform.strand.bias.analysis == TRUE) {
        PrintLog("Step 03-1. Perform strand bias analysis [firevat_strand_bias::PerformStrandBiasAnalysis]")

        data$refined.vcf.obj <- PerformStrandBiasAnalysis(
            vcf.obj = data$refined.vcf.obj,
            ref.forward.strand.var = data$ref.forward.strand.var,
            ref.reverse.strand.var = data$ref.reverse.strand.var,
            alt.forward.strand.var = data$alt.forward.strand.var,
            alt.reverse.strand.var = data$alt.reverse.strand.var,
            perform.fdr.correction = data$strand.bias.perform.fdr.correction,
            fdr.correction.method = data$strand.bias.fdr.correction.method)

        if (data$filter.by.strand.bias.analysis == TRUE) {
            PrintLog("* Filter by strand bias analysis results [firevat_strand_bias::FilterByStrandBiasAnalysis]")
            filtered.vcf.objs <- FilterByStrandBiasAnalysis(
                refined.vcf.obj = data$refined.vcf.obj,
                artifactual.vcf.obj = data$artifactual.vcf.obj,
                perform.fdr.correction = data$strand.bias.perform.fdr.correction,
                filter.by.strand.bias.analysis.cutoff = data$filter.by.strand.bias.analysis.cutoff)
            data$refined.vcf.obj <- filtered.vcf.objs$refined.vcf.obj
            data$artifactual.vcf.obj <- filtered.vcf.objs$artifactual.vcf.obj
        }

        data$artifactual.vcf.obj <- PerformStrandBiasAnalysis(
            vcf.obj = data$artifactual.vcf.obj,
            ref.forward.strand.var = data$ref.forward.strand.var,
            ref.reverse.strand.var = data$ref.reverse.strand.var,
            alt.forward.strand.var = data$alt.forward.strand.var,
            alt.reverse.strand.var = data$alt.reverse.strand.var,
            perform.fdr.correction = data$strand.bias.perform.fdr.correction,
            fdr.correction.method = data$strand.bias.fdr.correction.method)
    }

    # 03-2. Annotate
    if (data$annotate == TRUE) {
        PrintLog("Step 03-2. Annotate variants [firevat_annotation::AnnotateVCFObj]")

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

    # 03-3. Run Mutalisk
    # Identify target mutational signatures
    # Here we fetch all signatures ever identified by Mutational Patterns
    if (data$mutalisk == TRUE) {
        PrintLog("Step 03-3. Perform Mutalisk mutational signature analysis [firevat_mutalisk::RunMutalisk]")
        PrintLog("* Preparing data")

        df.optimization.logs <- ReadOptimizationIterationReport(data = data)
        Split.Sigs <- function(sigs, weights, cutoff = 0.05) {
            include <- rep(TRUE, length(sigs))
            include[which(sigs == "")] <- FALSE
            include[is.na(sigs)] <- FALSE

            sigs <- as.character(sigs[include])
            weights <- as.character(weights[include])
            sigs <- lapply(sigs, function(x) strsplit(x, ',')[[1]])
            weights <- lapply(weights, function(x) strsplit(x, ',')[[1]])

            if (length(sigs) == 0) {
                return(c())
            }

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
        if (is.null(data$mutalisk.must.include.sigs) == FALSE) {
            data$mut.pat.target.sigs <- unique(c(sigs1, sigs2, sigs3, sigs4, data$sequencing.artifact.mut.sigs, data$mutalisk.must.include.sigs))
        } else {
            data$mut.pat.target.sigs <- unique(c(sigs1, sigs2, sigs3, sigs4, data$sequencing.artifact.mut.sigs))
        }

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
                                                          df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                          target.mut.sigs = data$mut.pat.target.sigs,
                                                          method = data$mutalisk.method,
                                                          n.sample = data$mutalisk.random.sampling.count,
                                                          n.iter = data$mutalisk.random.sampling.max.iter,
                                                          verbose = data$verbose)
        # Artifact vcf
        data$artifactual.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$artifactual.vcf.obj,
                                                              df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                              target.mut.sigs = data$mut.pat.target.sigs,
                                                              method = data$mutalisk.method,
                                                              n.sample = data$mutalisk.random.sampling.count,
                                                              n.iter = data$mutalisk.random.sampling.max.iter,
                                                              verbose = data$verbose)
    }

    # 04. Write VCF Files
    if (data$write.vcf == TRUE) {
        PrintLog("Step 04. Write refined and artifactual VCF files [firevat_vcf::WriteVCF]")

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
    }

    data$end.datetime <- Sys.time()
    data$variant.refinement.terminiation.log <- "Successful"

    # 05. Report results
    if (data$report == TRUE) {
        PrintLog("Step 05. Generate FIREVAT report")
        data <- ReportFIREVATResults(data = data)
    }

    # 06. Save data
    if (data$save.rdata == TRUE) {
        PrintLog("Step 06. Write .RData")
        save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
    }

    # 07. Save tsv
    if (data$save.tsv == TRUE) {
        PrintLog("Step 07. Write FIREVAT results to .tsv file")
        WriteFIREVATResultsToTSV(firevat.results = data)
    }

    return(data)
}

#' @title RunManualMode
#' @description Runs FIREVAT manual mode
#'
#' @param data A list from RunFIREVAT
#'
#' @return A list
#'
#' @export
RunManualMode <- function(data) {
    if (data$verbose == TRUE) {
        PrintLog("Step 02. Run FIREVAT manual variant refinement.")
    }
    data$x.solution.decimal <- unlist(data$vcf.filter) # default values

    # 02-1. Filter vcf.data based on the updated vcf.filter
    PrintLog("Step 02-1. Filter VCF based on optmized filter parameters.")
    data$vcf.filter <- UpdateFilter(vcf.filter = data$vcf.filter,
                                    param.values = data$x.solution.decimal)
    optimized.vcf.objs <- FilterVCF(vcf.obj = data$vcf.obj,
                                    config.obj = data$config.obj,
                                    vcf.filter = data$vcf.filter,
                                    verbose = data$verbose)
    data$refined.vcf.obj <- optimized.vcf.objs$vcf.obj.filtered
    data$artifactual.vcf.obj <- optimized.vcf.objs$vcf.obj.artifact

    # Check if point mutations remaining in either refined.vcf.obj or artifactual.vcf.obj.
    # Either refinement is too liberal or too stringent. Return with a message.
    if (nrow(data$refined.vcf.obj$data) == 0) {
        PrintLog("After performing variant refinement there are no mutations remaining in the refined set.", type = "WARNING")
        data$end.datetime <- Sys.time()
        data$variant.refinement.terminiation.log <- "Successful but after performing variant refinement there are no mutations remaining in the refined set"
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }
    if (nrow(data$artifactual.vcf.obj$data) == 0) {
        PrintLog("After performing variant refinement there are no mutations remaining in the artifactual set.", type = "WARNING")
        data$end.datetime <- Sys.time()
        data$variant.refinement.terminiation.log <- "Successful but after performing variant refinement there are no mutations remaining in the artifactual set"
        if (data$save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (data$save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    }

    # 03. Additional analysis
    # 03-1. Strand bias analysis
    PrintLog("Step 03. Additional analysis.")
    if (data$perform.strand.bias.analysis == TRUE) {
        PrintLog("Step 03-1. Perform strand bias analysis [firevat_strand_bias::PerformStrandBiasAnalysis]")

        data$refined.vcf.obj <- PerformStrandBiasAnalysis(
            vcf.obj = data$refined.vcf.obj,
            ref.forward.strand.var = data$ref.forward.strand.var,
            ref.reverse.strand.var = data$ref.reverse.strand.var,
            alt.forward.strand.var = data$alt.forward.strand.var,
            alt.reverse.strand.var = data$alt.reverse.strand.var,
            perform.fdr.correction = data$strand.bias.perform.fdr.correction,
            fdr.correction.method = data$strand.bias.fdr.correction.method)

        if (data$filter.by.strand.bias.analysis == TRUE) {
            PrintLog("* Filter by strand bias analysis results [firevat_strand_bias::FilterByStrandBiasAnalysis]")
            filtered.vcf.objs <- FilterByStrandBiasAnalysis(
                refined.vcf.obj = data$refined.vcf.obj,
                artifactual.vcf.obj = data$artifactual.vcf.obj,
                perform.fdr.correction = data$strand.bias.perform.fdr.correction,
                filter.by.strand.bias.analysis.cutoff = data$filter.by.strand.bias.analysis.cutoff)
            data$refined.vcf.obj <- filtered.vcf.objs$refined.vcf.obj
            data$artifactual.vcf.obj <- filtered.vcf.objs$artifactual.vcf.obj
        }

        data$artifactual.vcf.obj <- PerformStrandBiasAnalysis(
            vcf.obj = data$artifactual.vcf.obj,
            ref.forward.strand.var = data$ref.forward.strand.var,
            ref.reverse.strand.var = data$ref.reverse.strand.var,
            alt.forward.strand.var = data$alt.forward.strand.var,
            alt.reverse.strand.var = data$alt.reverse.strand.var,
            perform.fdr.correction = data$strand.bias.perform.fdr.correction,
            fdr.correction.method = data$strand.bias.fdr.correction.method)
    }

    # 03-2. Annotate
    if (data$annotate == TRUE) {
        PrintLog("Step 03-2. Annotate variants [firevat_annotation::AnnotateVCFObj]")

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

    # 03-3. Run Mutalisk
    if (data$mutalisk == TRUE) {
        PrintLog("Step 03-3. Perform Mutalisk mutational signature analysis [firevat_mutalisk::RunMutalisk]")

        # Original vcf
        data$raw.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$vcf.obj,
                                                      df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                      target.mut.sigs = data$target.mut.sigs,
                                                      method = data$mutalisk.method,
                                                      n.sample = data$mutalisk.random.sampling.count,
                                                      n.iter = data$mutalisk.random.sampling.max.iter,
                                                      verbose = data$verbose)
        # Refined vcf
        data$refined.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$refined.vcf.obj,
                                                          df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                          target.mut.sigs = data$target.mut.sigs,
                                                          method = data$mutalisk.method,
                                                          n.sample = data$mutalisk.random.sampling.count,
                                                          n.iter = data$mutalisk.random.sampling.max.iter,
                                                          verbose = data$verbose)
        # Artifact vcf
        data$artifactual.muts.mutalisk.results <- RunMutalisk(vcf.obj = data$artifactual.vcf.obj,
                                                              df.ref.mut.sigs = data$df.ref.mut.sigs,
                                                              target.mut.sigs = data$target.mut.sigs,
                                                              method = data$mutalisk.method,
                                                              n.sample = data$mutalisk.random.sampling.count,
                                                              n.iter = data$mutalisk.random.sampling.max.iter,
                                                              verbose = data$verbose)
    }

    # 04. Write VCF Files
    if (data$write.vcf == TRUE) {
        PrintLog("Step 04. Write refined and artifactual VCF files [firevat_vcf::WriteVCF]")

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
    }

    data$end.datetime <- Sys.time()

    # 05. Report results
    if (data$report == TRUE) {
        PrintLog("Step 05. Generate FIREVAT report")
        data <- ReportFIREVATResults(data = data)
    }

    # 06. Save data
    if (data$save.rdata == TRUE) {
        PrintLog("Step 06. Write .RData")
        save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
    }

    # 07. Save tsv
    if (data$save.tsv == TRUE) {
        PrintLog("Step 07. Write FIREVAT results to .tsv file")
        WriteFIREVATResultsToTSV(firevat.results = data)
    }

    return(data)
}


#' @title RunFIREVAT
#' @description
#' Runs FIREVAT using configuration data. Filters point mutations in the user-specified vcf file based on mutational signature
#' identification and outputs the refined and artifact vcf files as well as metadata related to the refinement process.
#'
#' @param vcf.file String value corresponding to input .vcf file. Please provide the full path.
#' @param vcf.file.genome Genome assembly of the input .vcf file. The genome should be supported by BSgenome.
#' @param check.chromosome.name Boolean value. If TRUE, FIREVAT checks chromosome names
#' @param config.file String value corresponding to input configuration file. For more details please refer to ...
#' @param df.ref.mut.sigs A data.frame of the reference mutational signatures
#' @param target.mut.sigs A character vector of the target mutational signatures from reference mutational signatures.
#' @param df.ref.mut.sigs.groups.colors A data.frame of the reference mutational signatures groups and colors
#' @param sequencing.artifact.mut.sigs A character vector of the sequencing artifact mutational signatures from reference mutational signatures.
#' @param num.cores Number of cores to allocate
#' @param output.dir String value of the desired output directory
#' @param mode String value. The value should be either 'ga' or 'manual'.
#' @param init.artifact.stop
#' Numeric value  less than 1. If the sum of sequencing artifact weights in the user-specified original VCF file (i.e. vcf.file)
#' is less than or equal to this value then FIREVAT does not perform variant refinement.
#' Default value is 0.05. Note that this option does not apply if 'mode' is 'manual'.
#' @param ga.type String value. The value should be either 'binray' or 'real-valued'.
#' @param objective.fn Objective value derivation function. Default: Default.Obj.Fn.
#' @param use.suggested.soln Boolean value. If TRUE, then FIREVAT passes the default values
#' of filter variables declared as 'use_in_filter' in the config file to the 'suggestions' parameter of
#' the Genetic Algorithm package. If FALSE, then FIREVAT supplies NULL to the GA package 'suggestions' parameter.
#' FIREVAT also computes baseline performance of each filter variable and uses fittest population from each variable
#' as a suggested solution.
#' @param ga.pop.size Integer value of the Genetic Algorithm 'population size' parameter. Default: 100.
#' This value should be set based on the number of filter parameters. Recommendation: 40 per filter parameter.
#' @param ga.max.iter Integer value of the Genetic Algorithm 'maximum iterations' parameter. Default: 100.
#' This value should be set based on the number of filter parameters. Recommendation: same as 'ga.pop.size'.
#' @param ga.run Integer value of the Genetic Algorithm 'run' parameter. Default: 50.
#' This value should be set based on the 'ga.max.iter' parameter. Recommendation: 25 percent of 'ga.max.iter'.
#' @param ga.pmutation Float value of the Genetic Algorithm 'mutation probability' parameter. Default: 0.1.
#' @param ga.preemptive.killing If TRUE, then preemptively kills populations that yield greater sequencing artifact weights sum
#' compared to the original mutatational signatures analysis
#' @param ga.seed Integer value of the Genetic Algorithm 'seed' parameter. Default: NULL.
#' @param mutalisk If TRUE, confirm mutational signature analysis with Mutalisk. Default: TRUE.
#' @param mutalisk.method Mutalisk signature identification method. Default: 'random.sampling'.
#' The value can be either 'all' or 'random.sampling'.
#' 'all' uses all target.mut.sigs to identify mutational signatures.
#' 'random.sampling' randomly samples from target.mut.sigs to identify mutational signatures.
#' @param mutalisk.must.include.sigs Signatures that must be included in the Mutalisk signature identification
#' A character vector corresponding to the signature names.
#' @param mutalisk.random.sampling.count Mutalisk random sampling count. Default: 20.
#' The number of signatures to sample from target.mut.sigs
#' @param mutalisk.random.sampling.max.iter Mutalisk random sampling maximum iteration. Default: 10.
#' The number of times Mutalisk randomly samples from target.mut.sigs before determining the candidate signatures.
#' @param perform.strand.bias.analysis If TRUE, then performs strand bias analysis.
#' @param filter.by.strand.bias.analysis If TRUE, then filters out variants in refined vcf based on strand bias analysis results
#' @param filter.by.strand.bias.analysis.cutoff The p.value or q value cutoff for filtering out variants.
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
#' @param write.vcf If TRUE, write original/refined/artifact vcfs. Default: TRUE.
#' @param report If TRUE, generate report. Default: TRUE.
#' @param report.format The format of FIREVAT report. We currently only support 'html'.
#' @param save.rdata If TRUE, save rdata. Default: TRUE.
#' @param save.tsv If TRUE, save tsv. Default: TRUE.
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{f}{ = A ggarrange object}
#'  \item{graphs}{ = A list of length 3; each element is a ggplot histogram}
#' }
#'
#' @importFrom GA ga
#' @importFrom BSgenome available.genomes
#' @export
RunFIREVAT <- function(vcf.file,
                       vcf.file.genome = 'hg19',
                       check.chromosome.name = TRUE,
                       config.file,
                       df.ref.mut.sigs = GetPCAWGMutSigs(),
                       target.mut.sigs = GetPCAWGMutSigsNames(),
                       df.ref.mut.sigs.groups.colors = GetPCAWGMutSigsEtiologiesColors(),
                       sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
                       num.cores = 2,
                       output.dir,
                       mode = "ga",
                       init.artifact.stop = 0.05,
                       # GA parameters
                       objective.fn = Default.Obj.Fn,
                       use.suggested.soln = TRUE,
                       ga.type = "real-valued",
                       ga.pop.size = 100,
                       ga.max.iter = 100,
                       ga.run = 50,
                       ga.pmutation = 0.1,
                       ga.preemptive.killing = FALSE,
                       ga.seed = NULL,
                       # Mutalisk parameters
                       mutalisk = TRUE,
                       mutalisk.method = "all",
                       mutalisk.must.include.sigs = NULL,
                       mutalisk.random.sampling.count = 20,
                       mutalisk.random.sampling.max.iter = 10,
                       # Strand bias analysis parameters
                       perform.strand.bias.analysis = FALSE,
                       filter.by.strand.bias.analysis = TRUE,
                       filter.by.strand.bias.analysis.cutoff = 0.25,
                       strand.bias.perform.fdr.correction = TRUE,
                       strand.bias.fdr.correction.method = "BH",
                       ref.forward.strand.var = NULL,
                       ref.reverse.strand.var = NULL,
                       alt.forward.strand.var = NULL,
                       alt.reverse.strand.var = NULL,
                       # Annotation parameters
                       annotate = FALSE,
                       df.annotation.db = NULL,
                       annotated.columns.to.display = NULL,
                       annotation.filter.key.value.pairs = NULL,
                       annotation.filter.condition = "AND",
                       write.vcf = TRUE,
                       report = TRUE,
                       save.rdata = TRUE,
                       save.tsv = TRUE,
                       report.format = "html",
                       verbose = TRUE) {
    # 00. Check input parameters
    if (is.character(vcf.file) == FALSE) {
        stop("The parameter 'vcf.file' must be a string.")
    }
    bsg.available <- BSgenome::available.genomes(splitNameParts = TRUE)
    if (!(vcf.file.genome %in% c(bsg.available$pkgname, bsg.available$provider_version))) {
        stop("The parameter 'vcf.file.genome' should be supported by BSgenome.")
    }
    if (is.character(config.file) == FALSE) {
        stop("The parameter 'config.file' must be a string.")
    }
    if (is.data.frame(df.ref.mut.sigs) == FALSE) {
        stop("The parameter 'df.ref.mut.sigs' must be a data.frame.")
    }
    if (all(target.mut.sigs %in% colnames(df.ref.mut.sigs)) == FALSE)  {
        stop("All elements of the parameter 'target.mut.sigs' must be present in 'df.ref.mut.sigs' columns.")
    }
    if (all(sequencing.artifact.mut.sigs %in% colnames(df.ref.mut.sigs)) == FALSE)  {
        stop("All elements of the parameter 'sequencing.artifact.mut.sigs' must be present in 'df.ref.mut.sigs' columns.")
    }
    if (is.numeric(num.cores) == FALSE || num.cores <= 0)  {
        stop("The parameter 'num.cores' must be an integer and be greater than 0.")
    }
    if (is.character(output.dir) == FALSE)  {
        stop("The parameter 'output.dir' must be a string.")
    }
    if (mode != "ga" && mode != "manual")  {
        stop("The parameter 'mode' must be either 'ga' or 'manual'.")
    }
    if (is.numeric(init.artifact.stop) == FALSE || init.artifact.stop >= 1) {
        stop("The parameter 'init.artifact.stop' must be a numeric value less than 1.")
    }
    if (is.numeric(ga.pop.size) == FALSE || ga.pop.size <= 0)  {
        stop("The parameter 'ga.pop.size' must be an integer and be greater than 0.")
    }
    if (is.numeric(ga.max.iter) == FALSE || ga.max.iter <= 0)  {
        stop("The parameter 'ga.max.iter' must be an integer and be greater than 0.")
    }
    if (is.numeric(ga.run) == FALSE || ga.run <= 0)  {
        stop("The parameter 'ga.run' must be an integer and be greater than 0.")
    }
    if (is.numeric(ga.pmutation) == FALSE || ga.pmutation < 0 || ga.pmutation > 1)  {
        stop("The parameter 'ga.pmutation' must be a value between 0 and 1.")
    }
    if (verbose != TRUE && verbose != FALSE)  {
        stop("The parameter 'verbose' must be a boolean.")
    }

    start.datetime <- Sys.time()
    if (verbose == TRUE) {
        PrintLog("Started running FIREVAT.")
        PrintLog("Step 01. Initialize FIREVAT input parameters.")
    }

    # 01. Create the output directory
    if (dir.exists(output.dir) == FALSE) {
        dir.create(output.dir, recursive = T)
        PrintLog(paste0("Step 01-1. Created output directory at ", output.dir))
    } else {
        if (verbose == TRUE) {
            PrintLog(paste0("Step 01-1. Output directory ", output.dir, " already exists."))
        }
    }

    # Remove existing '*FIREVAT_Optimization_Logs.tsv'
    existing.firevat.optimization.log.tsv.file <- list.files(output.dir, pattern = "*FIREVAT_Optimization_Logs.tsv")
    if (identical(existing.firevat.optimization.log.tsv.file, character(0)) == FALSE) {
        # Remove
        PrintLog(paste0("Removing existing FIREVAT optimization logs (tsv) file: ",
                        existing.firevat.optimization.log.tsv.file))
        file.remove(paste0(output.dir, existing.firevat.optimization.log.tsv.file))
    }

    # 02. Read the vcf file
    PrintLog("Step 01-2. Read VCF file [firevat_vcf::ReadVCF].")
    vcf.obj <- ReadVCF(vcf.file = vcf.file, genome = vcf.file.genome,
                       check.chromosome.name = check.chromosome.name)

    # 03. Read the config file
    PrintLog("Step 01-3. Read config file [firevat_config::ParseConfigFile]")
    config.obj <- ParseConfigFile(config.file, verbose = verbose)

    # 04.Initially select point mutations that satisfy the config file and
    # consider these the initial set of mutations
    PrintLog("Step 01-4. Initialize VCF file [firevat_vcf::InitializeVCF]")
    vcf.objs <- InitializeVCF(vcf.obj = vcf.obj,
                              config.obj = config.obj,
                              verbose = verbose)
    vcf.obj <- vcf.objs$vcf.obj.filtered

    # Write the unidentifiable data
    WriteVCF(vcf.obj = vcf.objs$vcf.obj.artifact,
             save.file = paste0(output.dir, gsub("\\.vcf", "", basename(vcf.file)), "_Unknown.vcf"))

    # 05. Make filter from config file
    PrintLog("Step 01-5. Make VCF filter [firevat_filter::MakeFilter]")
    vcf.filter <- MakeFilter(config.obj)

    # 06. Prepare data for optimization
    PrintLog("Step 01-6. Prepare data for optimization")
    data = list(vcf.file = vcf.file,
                vcf.file.basename = gsub("\\.vcf", "", basename(vcf.file)),
                vcf.obj = vcf.obj,
                vcf.filter = vcf.filter,
                config.file = config.file,
                config.obj = config.obj,
                df.ref.mut.sigs = df.ref.mut.sigs,
                target.mut.sigs = target.mut.sigs,
                df.ref.mut.sigs.groups.colors = df.ref.mut.sigs.groups.colors,
                sequencing.artifact.mut.sigs = sequencing.artifact.mut.sigs,
                output.dir = output.dir,
                mode = mode,
                num.cores = num.cores,
                init.artifact.stop = init.artifact.stop,
                # GA parameters
                objective.fn = objective.fn,
                ga.type = ga.type,
                ga.pmutation = ga.pmutation,
                ga.pop.size = ga.pop.size,
                ga.max.iter = ga.max.iter,
                ga.run = ga.run,
                ga.preemptive.killing = ga.preemptive.killing,
                ga.seed = ga.seed,
                use.suggested.soln = use.suggested.soln,
                # Mutalisk parameters
                mutalisk = mutalisk,
                mutalisk.method = mutalisk.method,
                mutalisk.must.include.sigs = mutalisk.must.include.sigs,
                mutalisk.random.sampling.count = mutalisk.random.sampling.count,
                mutalisk.random.sampling.max.iter = mutalisk.random.sampling.max.iter,
                report.format = report.format,
                # Strand bias analysis parameters
                perform.strand.bias.analysis = perform.strand.bias.analysis,
                filter.by.strand.bias.analysis = filter.by.strand.bias.analysis,
                filter.by.strand.bias.analysis.cutoff = filter.by.strand.bias.analysis.cutoff,
                ref.forward.strand.var = ref.forward.strand.var,
                ref.reverse.strand.var = ref.reverse.strand.var,
                alt.forward.strand.var = alt.forward.strand.var,
                alt.reverse.strand.var = alt.reverse.strand.var,
                strand.bias.fdr.correction.method = strand.bias.fdr.correction.method,
                strand.bias.perform.fdr.correction = strand.bias.perform.fdr.correction,
                # Annotation parameters
                annotate = annotate,
                df.annotation.db = df.annotation.db,
                annotated.columns.to.display = annotated.columns.to.display,
                annotation.filter.key.value.pairs = annotation.filter.key.value.pairs,
                annotation.filter.condition = annotation.filter.condition,
                # Miscellaneous
                start.datetime = start.datetime,
                firevat.version = packageVersion("FIREVAT"),
                write.vcf = write.vcf,
                report = report,
                save.rdata = save.rdata,
                save.tsv = save.tsv,
                report.format = report.format,
                verbose = verbose)

    # 07. FIREVAT can only be run if there are more than 50 point mutations in the initial vcf file
    PrintLog("Step 01-7. Check mutations count")
    if (nrow(vcf.obj$data) <= 50) {
        PrintLog("FIREVAT must have at least 50 mutations to run. Returning without running FIREVAT.", type = "WARNING")
        data$end.datetime <- Sys.time()
        data$variant.refinement.performed <- FALSE
        data$variant.refinement.terminiation.log <- "Not enough point mutations (less than or equal 50)"
        if (save.rdata == TRUE) {
            save(data, file = paste0(data$output.dir, data$vcf.file.basename, "_FIREVAT_data.RData"))
        }
        if (save.tsv == TRUE) {
            WriteFIREVATResultsToTSV(firevat.results = data)
        }
        return(data)
    } else {
        if (verbose == TRUE) {
            PrintLog(paste0("* ", nrow(vcf.obj$data), " point mutations identified for potential FIREVAT variant refinement."))
        }
    }

    # 08. Optimize filter parameters
    if (mode == "ga") {
        data <- RunGAMode(data = data)
    } else if (mode == "manual") {
        data <- RunManualMode(data = data)
    }

    return(data)
}
