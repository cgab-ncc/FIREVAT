# FIREVAT Report Functions
#
# Last revised date:
#   February 18, 2019
#
# Authors: 
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title ReadOptimizationIterationReport
#' @description Read optimization iteration report
#'
#' @param data A list of main data from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame of FIREVAT optimization logs
#' @export
ReadOptimizationIterationReport <- function(data)  {
    log.file <- paste0(data$output.dir,
                       data$vcf.file.basename, "_FIREVAT_Optimization_Logs.tsv")
    df <- read.table(log.file, sep = "\t",
                     header = T, check.names = F, stringsAsFactors = F)
    df <- df[complete.cases(df),]
    return(df)
}


#' @title GetOptimizedSignatures
#' @description 
#' This function fetches the last row from the optimization iteration log
#' and returns the target and artifactual mutational signatures
#' for the type of mutations ('refined' or 'artifactual')
#'
#' @param data A list of main data from \code{\link{RunFIREVAT}}
#' @param mutations.type A string for type of mutations ('refined' or 'artifact')
#' @param signatures A string ('all', 'target', 'artifact')
#'
#' @return A data.frame with the columns 'signature' and 'weight'
#' @export
GetOptimizedSignatures <- function(data, mutations.type = "refined",
                                   signatures = "all")  {

    if (mutations.type == "refined")  {
        target.signatures.str <- "refined.muts.target.signatures"
        target.signatures.weights.str <- "refined.muts.target.signatures.weights"
        artifact.signatures.str <- "refined.muts.sequencing.artifact.signatures"
        artifact.signatures.weights.str <- "refined.muts.sequencing.artifact.signatures.weights"
    }
    else if (mutations.type == "artifact") {
        target.signatures.str <- "artifactual.muts.target.signatures"
        target.signatures.weights.str <- "artifactual.muts.target.signatures.weights"
        artifact.signatures.str <- "artifactual.muts.sequencing.artifact.signatures"
        artifact.signatures.weights.str <- "artifactual.muts.sequencing.artifact.signatures.weights"
    }
    else  {
        print(paste0("ERROR - unknown parameter passed.",
                     "'mutations.type' must be one of the following: ",
                     "'refined' or 'artifact'"))
        return(NULL)
    }
    
    # Get the signatures from the last optimization iteration from Mutational Patterns 
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)
    
    target.signatures <- tail(df.optimization.logs[,target.signatures.str], 1)
    target.signatures <- strsplit(target.signatures, ",")[[1]]
    target.signatures.weights <- tail(df.optimization.logs[,target.signatures.weights.str], 1)
    target.signatures.weights <- strsplit(target.signatures.weights, ",")[[1]]
    
    artifact.signatures <- tail(df.optimization.logs[,artifact.signatures.str], 1)
    artifact.signatures <- strsplit(artifact.signatures, ",")[[1]]
    artifact.signatures.weights <- tail(df.optimization.logs[,artifact.signatures.weights.str], 1)
    artifact.signatures.weights <- strsplit(artifact.signatures.weights, ",")[[1]]
    
    if (signatures == "all")  {
        df <- data.frame(list(signature = c(target.signatures,artifact.signatures),
                              weight = c(target.signatures.weights, artifact.signatures.weights)),
                         stringsAsFactors = F,
                         check.names = F)
    }
    else if (signatures == "target")  {
        df <- data.frame(list(signature = target.signatures,
                              weight = target.signatures.weights),
                         stringsAsFactors = F,
                         check.names = F)
    }
    else if (signatures == "artifact")  {
        df <- data.frame(list(signature = artifact.signatures,
                              weight = artifact.signatures.weights),
                         stringsAsFactors = F,
                         check.names = F)
    }
    else  {
        print(paste0("ERROR - unknown parameter passed.",
                     "'signatures' must be one of the following: ",
                     "'all', 'target' or 'artifact'"))
        return(NULL)
    }
    return(df)
}


#' @title PrepareReportDatafromConfig
#' @description Prepares report data (table data)
#'
#' @param x.solution.decimal A numeric vector
#' @param data A list of main data from \code{\link{RunFIREVAT}}
#'
#' @return A list with following items:
#' \itemize{
#'  \item{df.parameters}{A data.frame of FIREVAT Genetic Algorithm (GA) Parameters}
#'  \item{df.filter.cutoffs}{A data.frame of variant filtering cutoffs}
#'  \item{df.results}{A data.frame of objective values}
#'  \item{df.mutalisk.results}{A data.frame of mutational signature analysis results from Mutalisk}
#' }
#' @export
PrepareReportDatafromConfig <- function(x.solution.decimal, data)  {

    # 0. Write VCF Files
    # Update filter
    vcf.filter <- UpdateFilterFromConfig(vcf.filter = data$vcf.filter, 
                                         param.values = x.solution.decimal)
    # Filter vcf.data based on the updated vcf.filter
    vcf.objs <- FilterVCFFromConfig(vcf.obj = data$vcf.obj,
                                    config.obj = data$config.obj,
                                    vcf.filter = vcf.filter, 
                                    verbose = T)
    # Write refined and artifact vcf files
    WriteVCF(vcf.obj = data$vcf.obj,
             save.file = paste0(data$output.dir, data$vcf.file.basename, "_Original.vcf"))
    WriteVCF(vcf.obj = vcf.objs$vcf.obj.filtered, 
             save.file = paste0(data$output.dir, data$vcf.file.basename, "_Refined.vcf"))
    WriteVCF(vcf.obj = vcf.objs$vcf.obj.artifact, 
             save.file = paste0(data$output.dir, data$vcf.file.basename, "_Artifact.vcf"))

    # 1. FIREVAT Genetic Algorithm (GA) Parameters table
    #
    # FIREVAT Genetic Algorithm (GA) Parameters
    # GA Population Size                                value
    # GA Maximum Iteration                              value
    # GA Run                                            value
    # GA Mutation Probability                           value
    #
    df.parameters <- data.frame(
        list("FIREVAT Genetic Algorithm (GA) Parameters" = c(
                "GA Population Size", "GA Maximum Iteration", "GA Run",
                "GA Mutation Probability"),
            " " = c(format(x = data$pop.size, nsmall = 0, big.mark = ","),
                    format(x = data$max.iter, nsmall = 0, big.mark = ","),
                    format(x = data$run, nsmall = 0, big.mark = ","),
                    format(x = data$pmutation, nsmall = 3))),
        stringsAsFactors = F,
        check.names = F)
    
    # 2. Filter Cutoffs
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)

    filter.param.values <- sapply(
        names(vcf.filter), function(x) return(tail(df.optimization.logs[x], 1))
    )
    
    filter.direction <- sapply(
        names(vcf.filter), function(x) return(data$config.obj[[x]]["direction"])
    )
    filter.direction <- replace(filter.direction, filter.direction=="POS", ">")
    filter.direction <- replace(filter.direction, filter.direction=="NEG", "<")
    #
    # Filter Variable       Filter Direction    Optimized Cutoff
    # Variable.1                 > or <             value
    #
    df.filter.cutoffs <- data.frame(
        list(
            "Filter Variable" = names(vcf.filter),
            "Filter Direction" = as.vector(unlist(filter.direction)),
            "Optimized Cutoff" = as.vector(unlist(filter.param.values))
            ),
        stringsAsFactors = F,
        check.names = F)

    # 3. FIREVAT Results
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)
    refined.muts.cos.sim <- tail(
        df.optimization.logs$refined.muts.cosine.similarity.score, 1)
    refined.muts.count <- tail(
        df.optimization.logs$refined.mutations.count, 1)
    refined.muts.proportion <- tail(
        df.optimization.logs$refined.muts.proportion, 1)
    refined.muts.seq.art.sigs.weights.sum <- tail(
        df.optimization.logs$refined.muts.sequencing.artifact.signatures.weights.sum, 1)
    artifactual.muts.cos.sim <- tail(
        df.optimization.logs$artifactual.muts.cosine.similarity.score, 1)
    artifactual.muts.count <- tail(
        df.optimization.logs$artifact.mutations.count, 1)
    artifactual.muts.proportion <- tail(
        df.optimization.logs$artifactual.muts.proportion, 1)
    artifactual.muts.seq.art.sigs.weights.sum <- tail(
        df.optimization.logs$artifactual.muts.sequencing.artifact.signatures.weights.sum, 1)
    objective.value <- tail(df.optimization.logs$objective.value, 1)
    #
    # Objective Value   C.refined   W.refined   C.artifact  W.artifacdt
    #       value         value       value       value       value
    #
    df.results <- data.frame(
        list("Objective Value " = c(format(x = objective.value, digits = 3)),
            "C.refined" = c(format(x = refined.muts.cos.sim, digits = 3)),
            "W.refined" = c(format(x = refined.muts.seq.art.sigs.weights.sum,
                digits = 3)),
            "C.artifact" = c(format(x = artifactual.muts.cos.sim, digits = 3)),
            "W.artifact" = c(format(x = artifactual.muts.seq.art.sigs.weights.sum,
                digits = 3))),
        stringsAsFactors = F,
        check.names = F)

    # 4. Mutalisk Results
    #
    #                               Original VCF    Refined VCF     Artifactual VCF
    # Mutations Count (%)               value           value           value
    # Cosine Similarity Score           value           value           value
    # Residual Sum of Squares (RSS)     value           value           value
    #
    df.mutalisk.results <- data.frame(
        list(" " = c("Mutations Count (%)","Cosine similarity score",
                     "Residual sum of squares (RSS)"),
            "Original VCF" = c(paste0(
                format(x = data$raw.muts.mutalisk.results$num.point.mutations, 
                    big.mark = ","), " (100%)"),
                format(x = data$raw.muts.mutalisk.results$cos.sim.score,
                    digits = 3),
                format(x = data$raw.muts.mutalisk.results$rss, digits = 3)),
            "Refined VCF" = c(paste0(
                format(x = data$refined.muts.mutalisk.results$num.point.mutations,
                    big.mark = ","), " (",
                format(x = refined.muts.proportion * 100, digits = 4), "%)"),
                format(x = data$refined.muts.mutalisk.results$cos.sim.score,
                    digits = 3),
                format(x = data$refined.muts.mutalisk.results$rss, digits = 3)),
            "Artifactual VCF" = c(paste0(
                format(x = artifactual.muts.count, big.mark = ","), " (", 
                format(x = artifactual.muts.proportion * 100, digits = 4), "%)"),
                format(x = data$artifactual.muts.mutalisk.results$cos.sim.score,
                    digits = 3),
                format(x = data$artifactual.muts.mutalisk.results$rss,
                    digits = 3))),
        stringsAsFactors = F,
        check.names = F)

    return(list(df.parameters = df.parameters,
                df.filter.cutoffs = df.filter.cutoffs,
                df.results = df.results,
                df.mutalisk.results = df.mutalisk.results))
}


#' @title PrepareReportPlotsfromConfig
#' @description Prepares report plots
#'
#' @param x.solution.decimal A numeric vector
#' @param data A list of main data from \code{\link{RunFIREVAT}}
#'
#' @return A list with following items:
#' \itemize{
#'  \item{optimization.iter.plot}{A plot for the trend of objective value during iteraction}
#'  \item{refined.muts.seq.art.iter.plot}{A plot for showing sequencing artifact weights in refined variants}
#'  \item{artifactual.muts.seq.art.iter.plot}{A plot for showing sequencing artifact types in artifactual variants}
#'  \item{raw.muts.mutalisk.plots}{A plot for mutalisk results of original variants}
#'  \item{refined.muts.mutalisk.plots}{A plot for mutalisk results of refined variants}
#'  \item{artifactual.muts.mutalisk.plots}{A plot for mutalisk results of artifactual variants}
#'  \item{original.vcf.stats.plots}{A histogram for original VCF stats}
#'  \item{refined.vcf.stats.plots}{A histogram for refined VCF stats}
#'  \item{artifact.vcf.stats.plots}{A histogram for artifactual VCF stats}
#'  \item{vcf.stats.boxplots}{A plot for showing comparison results among VCF stats of original, refined and artifactual VCFs}
#' }
#' @importFrom graphics hist
#' @export
PrepareReportPlotsfromConfig <- function(x.solution.decimal, data)  {
    # Prepares report plots
    #
    # Args:
    #   x.solution.decimal = numeric vector
    #   data = list with the following elements
    #
    # Returns:
    #
    #
    
    # 1. Refinement Optimization> Section
    # 1.1 Objective Function Optimization
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)
    df.optimization.logs.temp <- df.optimization.logs[,c(
        "iteration",
        "refined.muts.sequencing.artifact.signatures.weights.sum",
        "refined.muts.cosine.similarity.score",
        "refined.muts.proportion",
        "objective.value"
    )]
    
    colnames(df.optimization.logs.temp) <- c(
        "iteration",
        "Refined Mutations Artifact Signatures Weights Sum",
        "Refined Mutations Cosine Similarity Score",
        "Refined Mutations Proportion",
        "Objective Value"
    )
    
    columns.to.plot <- colnames(df.optimization.logs.temp)[2:length(
        colnames(df.optimization.logs.temp)
    )]
    
    optimization.iter.plot <- PlotOptimizationIterations(
    
        df = df.optimization.logs.temp,
        columns.to.plot = columns.to.plot,
        x.axis.var = "iteration",
        x.axis.title = "Iteration",
        y.axis.title = "Value",
        title = "Objective Function Optimization",
        x.max = max(df.optimization.logs$iteration),
        legend.ncol = 2,
        save.file = paste0(data$output.dir, 
            data$vcf.file.basename,
            "_FIREVAT_Optimization_Iterations.png"),
        connect.dots = T
    
    )

    # 2. Plot sequencing artifact weight (refined mutations)
    # Enumerate sequencing artifact signature weights at each iteration
    # Get all sequencing artifact signatures
    refined.muts.seq.artifact.sigs <- c()
    for (i in 1:nrow(df.optimization.logs))  {
        curr.row <- df.optimization.logs[i,]
        curr.row.seq.art.sigs <- as.character(curr.row$refined.muts.sequencing.artifact.signatures)
        if (nchar(curr.row.seq.art.sigs) > 0)  {
            curr.row.seq.art.sigs <- strsplit(curr.row.seq.art.sigs, "\\,")[[1]]
            
            refined.muts.seq.artifact.sigs <- c(
                refined.muts.seq.artifact.sigs, curr.row.seq.art.sigs
            )
        }
    }
    refined.muts.seq.artifact.sigs <- unique(refined.muts.seq.artifact.sigs)
    # Create a dataframe with columns 'iteration', 'SBS<NUM>', 'SBS<NUM>', ...
    
    df.plot.data <- data.frame(
        matrix(0, ncol = 1 + length(refined.muts.seq.artifact.sigs),nrow = 0)
    )
    
    colnames(df.plot.data) <- c("iteration", refined.muts.seq.artifact.sigs)
    
    # Fill in data for each iteration
    for (i in 1:nrow(df.optimization.logs))  {
        curr.row <- df.optimization.logs[i,]
        curr.row.iter <- as.integer(curr.row$iteration)
        curr.row.seq.art.sigs <- as.character(
            curr.row$refined.muts.sequencing.artifact.signatures)
        curr.row.seq.art.sigs.probs <- as.character(
            curr.row$refined.muts.sequencing.artifact.signatures.weights)
        
        df.plot.data.temp <- data.frame(
            matrix(0, ncol = 1 + length(refined.muts.seq.artifact.sigs), nrow = 1)
        )
        
        colnames(df.plot.data.temp) <- c(
            "iteration", refined.muts.seq.artifact.sigs)

        df.plot.data.temp['iteration'] <- curr.row.iter
        
        if (nchar(curr.row.seq.art.sigs) > 0)  {
            
            curr.row.seq.art.sigs <-as.character(
                strsplit(curr.row.seq.art.sigs, "\\,")[[1]]
            )
            
            curr.row.seq.art.sigs.probs <- as.numeric(
                strsplit(curr.row.seq.art.sigs.probs, "\\,")[[1]]
            )
            
            for (j in 1:length(curr.row.seq.art.sigs))  {
                curr.sig <- curr.row.seq.art.sigs[j]
                curr.sig.prob <- curr.row.seq.art.sigs.probs[j]
                df.plot.data.temp[curr.sig] <- curr.sig.prob
            }
        }
        df.plot.data <- rbind(df.plot.data, df.plot.data.temp)
    }
    refined.muts.seq.art.iter.plot <- PlotOptimizationIterations(
        df = df.plot.data,
        columns.to.plot = c(refined.muts.seq.artifact.sigs),
        x.axis.var = "iteration",
        x.axis.title = "Iteration",
        x.max = max(df.plot.data$iteration),
        y.max = max(df.plot.data) * 1.3,
        y.axis.title = "Weight",
        title = "Refined Mutations Sequencing Artifact Weights at Each Iteration",
        legend.ncol = 5,
        save.file = paste0(
            data$output.dir, 
            data$vcf.file.basename,
            "_FIREVAT_Optimization_Iterations_Refined_Muts","_Seq_Art_Sigs_Weights.png"
            ),
        connect.dots = T)

    # 3. Plot sequencing artifact weight (artifactual mutations)
    # Enumerate sequencing artifact signature weights at each iteration
    # Get all sequencing artifact signatures
    artifactual.muts.seq.artifact.sigs <- c()
    for (i in 1:nrow(df.optimization.logs))  {
        curr.row <- df.optimization.logs[i,]
        curr.row.seq.art.sigs <- as.character(
            curr.row$artifactual.muts.sequencing.artifact.signatures)
        if (nchar(curr.row.seq.art.sigs) > 0)  {
            curr.row.seq.art.sigs <- strsplit(curr.row.seq.art.sigs, "\\,")[[1]]
            artifactual.muts.seq.artifact.sigs <- c(
                artifactual.muts.seq.artifact.sigs, curr.row.seq.art.sigs
            )
        }
    }
    artifactual.muts.seq.artifact.sigs <- unique(
        artifactual.muts.seq.artifact.sigs
    )
    # Create a dataframe with columns 'iteration', 'SBS<NUM>', 'SBS<NUM>', ...
    df.plot.data <- data.frame(
        matrix(0, ncol = 1 + length(artifactual.muts.seq.artifact.sigs),
               nrow = 0)
    )
    colnames(df.plot.data) <- c("iteration", artifactual.muts.seq.artifact.sigs)
    # Fill in data for each iteration
    for (i in 1:nrow(df.optimization.logs))  {
        curr.row <- df.optimization.logs[i,]
        curr.row.iter <- as.integer(curr.row$iteration)
        curr.row.seq.art.sigs <- as.character(
            curr.row$artifactual.muts.sequencing.artifact.signatures)
        curr.row.seq.art.sigs.probs <- as.character(
            curr.row$artifactual.muts.sequencing.artifact.signatures.weights)
        df.plot.data.temp <- data.frame(
            matrix(0, ncol = 1 + length(artifactual.muts.seq.artifact.sigs),
                   nrow = 1)
        )
        colnames(df.plot.data.temp) <- c("iteration", 
            artifactual.muts.seq.artifact.sigs)
        df.plot.data.temp['iteration'] <- curr.row.iter

        if (nchar(curr.row.seq.art.sigs) > 0)  {
            curr.row.seq.art.sigs <-as.character(
                strsplit(curr.row.seq.art.sigs, "\\,")[[1]])
            curr.row.seq.art.sigs.probs <- as.numeric(
                strsplit(curr.row.seq.art.sigs.probs, "\\,")[[1]])
            for (j in 1:length(curr.row.seq.art.sigs))  {
                curr.sig <- curr.row.seq.art.sigs[j]
                curr.sig.prob <- curr.row.seq.art.sigs.probs[j]
                df.plot.data.temp[curr.sig] <- curr.sig.prob
            }
        }
        df.plot.data <- rbind(df.plot.data, df.plot.data.temp)
    }
    artifactual.muts.seq.art.iter.plot <- PlotOptimizationIterations(
        df = df.plot.data,
        columns.to.plot = c(artifactual.muts.seq.artifact.sigs),
            x.axis.var = "iteration",
            x.axis.title = "Iteration",
            x.max = max(df.plot.data$iteration),
            y.max = max(df.plot.data[,artifactual.muts.seq.artifact.sigs]) * 1.3,
            title = paste0(
                "Artifactual Mutations Sequencing Artifact",
                "Weights at Each Iteration"
                ),
            legend.ncol = 7,
            save.file = paste0(
                data$output.dir,
                data$vcf.file.basename,
                "_FIREVAT_Optimization_Iterations_Artifactual_Muts",
                "_Seq_Art_Sigs_Weights.png"
                ),
            connect.dots = T)
    
        
    # 3. Plot Mutalisk results
    # Update MuTect2 filter
    vcf.filter <- UpdateFilterFromConfig(vcf.filter = data$vcf.filter, 
                                         param.values = x.solution.decimal)
    # Filter vcf.data based on the updated vcf.filter
    vcf.objs <- FilterVCFFromConfig(vcf.obj = data$vcf.obj, 
                                    vcf.filter = vcf.filter,
                                    config.obj = data$config.obj, 
                                    verbose = T)

    # Plot Mutalisk results
    trinuc.max.y <- max(c(data$raw.muts.mutalisk.results$sub.types.spectrum,
                          data$raw.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$sub.types.spectrum,
                          data$refined.muts.mutalisk.results$residuals,
                          data$artifactual.muts.mutalisk.results$sub.types.spectrum,
                          data$artifactual.muts.mutalisk.results$residuals)) * 1.3
    trinuc.min.y <- min(c(data$raw.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$residuals))
    mut.type.max.y <- max(c(
        EnumerateTriNucCounts(
            data$raw.muts.mutalisk.results$identified.mut.sigs.spectrum),
        EnumerateTriNucCounts(
            data$refined.muts.mutalisk.results$identified.mut.sigs.spectrum),
        EnumerateTriNucCounts(
            data$artifactual.muts.mutalisk.results$identified.mut.sigs.spectrum)
        )) * 1.1
    all.sigs <- c(data$raw.muts.mutalisk.results$identified.mut.sigs,
                  data$refined.muts.mutalisk.results$identified.mut.sigs,
                  data$artifactual.muts.mutalisk.results$identified.mut.sigs)
    
    raw.muts.mutalisk.plots <- PlotMutaliskResults(
        mutalisk.results = data$raw.muts.mutalisk.results,
        signatures = all.sigs,
        trinuc.max.y = trinuc.max.y,
        trinuc.min.y = trinuc.min.y,
        mut.type.max.y = mut.type.max.y,
        title = "Original VCF"
    )
    
    refined.muts.mutalisk.plots <- PlotMutaliskResults(
        mutalisk.results = data$refined.muts.mutalisk.results,
        signatures = all.sigs,
        trinuc.max.y = trinuc.max.y,
        trinuc.min.y = trinuc.min.y,
        mut.type.max.y = mut.type.max.y,
        title = "Refined VCF"
    )
    
    artifactual.muts.mutalisk.plots <- PlotMutaliskResults(
        mutalisk.results = data$artifactual.muts.mutalisk.results,
        signatures = all.sigs,
        trinuc.max.y = trinuc.max.y,
        trinuc.min.y = trinuc.min.y,
        mut.type.max.y = mut.type.max.y,
        title = "Artifactual VCF"
    )
    
    # 4. VCF Statistics
    # Get VCF statistics

    original.vcf.stats <- GetVCFValues(vcf.obj = data$vcf.obj,
                                       vcf.filter = data$vcf.filter)
    refined.vcf.stats <- GetVCFValues(vcf.obj = vcf.objs$vcf.obj.filtered,
                                      vcf.filter = data$vcf.filter)
    artifact.vcf.stats <- GetVCFValues(vcf.obj = vcf.objs$vcf.obj.artifact,
                                       vcf.filter = data$vcf.filter)

    x.axis.labels <- names(data$vcf.filter)
    # Get the maximum frequency value for each VCF statistic variable
    stat.y.max.vals <- c()
    stat.x.max.vals <- c()
    
    for (param in x.axis.labels)  {
        
        stat.x.max.vals <- c(stat.x.max.vals, max(max(original.vcf.stats[[param]]),
                                                  max(refined.vcf.stats[[param]]),
                                                  max(artifact.vcf.stats[[param]])))
        orig.hist.data <- hist(original.vcf.stats[[param]],
                               plot = F,
                               breaks = seq(0, max(original.vcf.stats[[param]]) + 1, by = 1))$counts
        refi.hist.data <- hist(refined.vcf.stats[[param]],
                               plot = F,
                               breaks = seq(0, max(refined.vcf.stats[[param]]) + 1, by = 1))$counts
        arti.hist.data <- hist(artifact.vcf.stats[[param]],
                               plot = F,
                               breaks = seq(0, max(artifact.vcf.stats[[param]]) + 1, by = 1))$counts
        stat.y.max.vals <- c(stat.y.max.vals,
                             max(c(max(orig.hist.data),
                                   max(refi.hist.data),
                                   max(arti.hist.data))))
    }
    # Plot VCF statistics

    df.optimization.logs <- ReadOptimizationIterationReport(data = data)

    cutoff.values <- as.vector(unlist(sapply(
        names(vcf.filter), function(x) return(tail(df.optimization.logs[x], 1))
    )))

    original.vcf.stats.plots <- PlotVCFStatsHistograms(
        plot.values = original.vcf.stats,
        x.axis.labels = x.axis.labels,
        sample.id = data$vcf.file.basename,
        title = "Original VCF",
        cutoff.values = cutoff.values,
        plot.cutoff.value.lines = T,
        stat.y.max.vals = stat.y.max.vals,
        stat.x.max.vals = stat.x.max.vals,
        save.file = paste0(
            data$output.dir, 
            data$vcf.file.basename,
            "_FIREVAT_Original_VCF_Stats_Histograms.png"
        )
    )

    refined.vcf.stats.plots <- PlotVCFStatsHistograms(
        plot.values = refined.vcf.stats,
        x.axis.labels = x.axis.labels,
        title = "Refined VCF",
        plot.cutoff.line.color = "gainsboro",
        cutoff.values = cutoff.values,
        plot.cutoff.value.lines = T,
        stat.y.max.vals = stat.y.max.vals,
        stat.x.max.vals = stat.x.max.vals,
        sample.id = data$vcf.file.basename,
        save.file = paste0(
            data$output.dir, 
            data$vcf.file.basename,
            "_FIREVAT_Refined_VCF_Stats_Histograms.png"
        )
    )
    
    artifact.vcf.stats.plots <- PlotVCFStatsHistograms(
        plot.values = artifact.vcf.stats,
        x.axis.labels = x.axis.labels,
        title = "Artifactual VCF",
        plot.cutoff.line.color = "gainsboro",
        cutoff.values = cutoff.values,
        plot.cutoff.value.lines = T,
        stat.y.max.vals = stat.y.max.vals,
        stat.x.max.vals = stat.x.max.vals,
        sample.id = data$vcf.file.basename,
        save.file = paste0(
            data$output.dir, 
            data$vcf.file.basename,
            "_FIREVAT_Artifact_VCF_Stats_Histograms.png"
        )
    )
    
    vcf.stats.boxplots <- list()
    
    for (i in 1:length(x.axis.labels))  {
        
        param <- x.axis.labels[i]
        
        curr.var.orig.vcf.stats <- original.vcf.stats[[param]]
        curr.var.refi.vcf.stats <- refined.vcf.stats[[param]]
        curr.var.arti.vcf.stats <- artifact.vcf.stats[[param]]
        
        g <- PlotVCFStatsBoxPlots(
                original.vcf.stat.values = curr.var.orig.vcf.stats,
                refined.vcf.stat.values = curr.var.refi.vcf.stats,
                artifact.vcf.stat.values = curr.var.arti.vcf.stats,
                xlab = x.axis.labels[i]
            )
        
        vcf.stats.boxplots[[i]] <- g
    }
    
    return(list(optimization.iter.plot = optimization.iter.plot,
                refined.muts.seq.art.iter.plot = refined.muts.seq.art.iter.plot,
                artifactual.muts.seq.art.iter.plot = artifactual.muts.seq.art.iter.plot,
                raw.muts.mutalisk.plots = raw.muts.mutalisk.plots,
                refined.muts.mutalisk.plots = refined.muts.mutalisk.plots,
                artifactual.muts.mutalisk.plots = artifactual.muts.mutalisk.plots,
                original.vcf.stats.plots = original.vcf.stats.plots,
                refined.vcf.stats.plots = refined.vcf.stats.plots,
                artifact.vcf.stats.plots = artifact.vcf.stats.plots,
                vcf.stats.boxplots = vcf.stats.boxplots))
}


#' @title ReportFIREVATResultsfromConfig
#' @description Reports FIREVAT results with in html format generated from Rmd
#'
#' @param x.solution.decimal A numeric vector
#' @param data A list of main data from \code{\link{RunFIREVAT}}
#' @param verbose If true, provides process detail
#' @param report.format The format of FIREVAT report. In current version, we only support "html"
#'
#' @importFrom rmarkdown render
#' @export
ReportFIREVATResultsfromConfig <- function(x.solution.decimal,
                                           data, verbose = T,
                                           report.format = "html")  {
    
    if (verbose) print("Started generating FIREVAT report")
    
    # Run Mutalisk
    # Update filter
    vcf.filter <- UpdateFilterFromConfig(vcf.filter = data$vcf.filter, 
                                      param.values = x.solution.decimal)
    # Filter vcf.data based on the updated vcf.filter
    vcf.objs <- FilterVCFFromConfig(vcf.obj = data$vcf.obj,
                                    config.obj = data$config.obj, 
                                    vcf.filter = vcf.filter,
                                    verbose = T)
    df.ref.mut.sigs <- GetPCAWGMutSigs()
    target.mut.sigs <- GetPCAWGMutSigsNames()
    data$raw.muts.mutalisk.results <- RunMutalisk(
        vcf.obj = data$vcf.obj,
        df.ref.mut.sigs = df.ref.mut.sigs,
        target.mut.sigs = target.mut.sigs
    )
    
    data$refined.muts.mutalisk.results <- RunMutalisk(
        vcf.obj = vcf.objs$vcf.obj.filtered,
        df.ref.mut.sigs = df.ref.mut.sigs,
        target.mut.sigs = target.mut.sigs
    )
    
    data$artifactual.muts.mutalisk.results <- RunMutalisk(
        vcf.obj = vcf.objs$vcf.obj.artifact,
        df.ref.mut.sigs = df.ref.mut.sigs,
        target.mut.sigs = target.mut.sigs
    )
     
    # Prepare data to run report.Rmd
    report.plots <- PrepareReportPlotsfromConfig(
        x.solution.decimal = x.solution.decimal, data = data
    )
    
    report.data <- PrepareReportDatafromConfig(
        x.solution.decimal = x.solution.decimal, data = data
    )

    # Generate html report
    if (report.format == "html") {
        rmarkdown::render(input = "../inst/rmd_template/report.Rmd", 
                          params = list(data = data,
                                        report.data = report.data,
                                        report.plots = report.plots),
                          output_dir = data$output.dir,
                          output_file = paste0(data$vcf.file.basename,
                                               "_FIREVAT_Report.html"))
    }

    if (verbose) print("Finished generating FIREVAT report")
}


#' @title GetVCFValues
#' @description Get values of filtering parameters from vcf.obj
#'
#' @param vcf.obj A list of vcf data
#' @param vcf.filter  A list of vcf filtering information
#'
#' @return A list with filtering parameter values
#' @keywords internal
GetVCFValues <- function(vcf.obj, vcf.filter)  {
    params <- names(vcf.filter)
    output.list <- vcf.obj$data[,..params]
    return(output.list)
}

