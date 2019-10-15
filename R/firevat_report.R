# FIREVAT Report Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title ReadOptimizationIterationReport
#' @description Read optimization iteration report
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame of FIREVAT optimization logs
#'
#' @export
ReadOptimizationIterationReport <- function(data) {
    log.file <- paste0(data$output.dir,
                       data$vcf.file.basename, "_FIREVAT_Optimization_Logs.tsv")
    df <- read.table(log.file, sep = "\t",
                     header = TRUE,
                     na.strings = "",
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
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
#'
#' @export
GetOptimizedSignatures <- function(data,
                                   mutations.type = "refined",
                                   signatures = "all") {

    if (mutations.type == "refined") {
        target.signatures.str <- "refined.muts.target.signatures"
        target.signatures.weights.str <- "refined.muts.target.signatures.weights"
        artifact.signatures.str <- "refined.muts.sequencing.artifact.signatures"
        artifact.signatures.weights.str <- "refined.muts.sequencing.artifact.signatures.weights"
    } else if (mutations.type == "artifact") {
        target.signatures.str <- "artifactual.muts.target.signatures"
        target.signatures.weights.str <- "artifactual.muts.target.signatures.weights"
        artifact.signatures.str <- "artifactual.muts.sequencing.artifact.signatures"
        artifact.signatures.weights.str <- "artifactual.muts.sequencing.artifact.signatures.weights"
    } else {
        PrintLog(paste0("ERROR - unknown parameter passed.",
                        "'mutations.type' must be one of the following: ",
                        "'refined' or 'artifact'"), type = "ERROR")
        return(NULL)
    }

    # Get the signatures from the last optimization iteration from Mutational Patterns
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)

    target.signatures <- tail(df.optimization.logs[,target.signatures.str], 1)
    if (is.na(target.signatures)) {
        target.signatures <- c()
        target.signatures.weights <- c()
    } else {
        target.signatures <- strsplit(target.signatures, ",")[[1]]
        target.signatures.weights <- tail(df.optimization.logs[,target.signatures.weights.str], 1)
        target.signatures.weights <- strsplit(target.signatures.weights, ",")[[1]]
    }

    artifact.signatures <- tail(df.optimization.logs[,artifact.signatures.str], 1)
    if (is.na(artifact.signatures)) {
        artifact.signatures <- c()
        artifact.signatures.weights <- c()
    } else {
        artifact.signatures <- strsplit(artifact.signatures, ",")[[1]]
        artifact.signatures.weights <- tail(df.optimization.logs[,artifact.signatures.weights.str], 1)
        artifact.signatures.weights <- strsplit(artifact.signatures.weights, ",")[[1]]
    }

    if (signatures == "all") {
        df <- data.frame(list(signature = c(target.signatures, artifact.signatures),
                              weight = c(target.signatures.weights, artifact.signatures.weights)),
                         stringsAsFactors = F,
                         check.names = F)
    } else if (signatures == "target") {
        df <- data.frame(list(signature = target.signatures,
                              weight = target.signatures.weights),
                         stringsAsFactors = F,
                         check.names = F)
    } else if (signatures == "artifact") {
        df <- data.frame(list(signature = artifact.signatures,
                              weight = artifact.signatures.weights),
                         stringsAsFactors = F,
                         check.names = F)
    } else {
        PrintLog(paste0("ERROR - unknown parameter passed.",
                        "'signatures' must be one of the following: ",
                        "'all', 'target' or 'artifact'"), type = "ERROR")
        return(NULL)
    }
    return(df)
}


#' @title PrepareFilterCutoffsTable
#' @description Prepares filter cutoffs table for reporting
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareFilterCutoffsTable <- function(data) { # COMMON_1
    filter.param.values <- data$x.solution.decimal
    filter.direction <- sapply(
        names(data$vcf.filter), function(x) return(data$config.obj[[x]]["direction"])
    )
    filter.direction <- replace(filter.direction, filter.direction=="POS", ">=")
    filter.direction <- replace(filter.direction, filter.direction=="NEG", "<=")

    # Filter Variable       Filter Direction    Optimized Cutoff
    # Variable.1                 > or <             value
    df.filter.cutoffs <- data.frame(
        list(
            "Filter Variable" = names(data$vcf.filter),
            "Filter Direction" = as.vector(unlist(filter.direction)),
            "Optimized Cutoff" = as.vector(unlist(filter.param.values))
        ),
        stringsAsFactors = F,
        check.names = F)
    return(df.filter.cutoffs)
}


#' @title PrepareIdentifiedSignaturesPlotHelper
#' @description PrepareIdentifiedSignaturesPlot helper function
#'
#' @param mutalisk.results A list of elements returned from \code{\link{RunMutalisk}}
#' @param signatures A character vector of signature names
#' @param title String value
#' @param df.ref.mut.sigs.groups.colors A data.frame with signature and color hex value
#'
#' @return A ggplot object
#'
#' @keywords internal
PrepareIdentifiedSignaturesPlotHelper <- function(mutalisk.results,
                                                  signatures,
                                                  title,
                                                  df.ref.mut.sigs.groups.colors) {
    # Prepare plot data
    df.identified.mut.sigs <- data.frame(
        list(signature = mutalisk.results$identified.mut.sigs,
             weight = mutalisk.results$identified.mut.sigs.probs),
        stringsAsFactors = F,
        check.names = F
    )
    remaining.sigs <- setdiff(signatures, as.character(df.identified.mut.sigs$signature))
    df.remaining.mut.sigs <- data.frame(list(signature = remaining.sigs,
                                             weight = rep(0, length(remaining.sigs))),
                                        stringsAsFactors = F,
                                        check.names = F)
    df.identified.mut.sigs <- rbind(df.identified.mut.sigs,
                                    df.remaining.mut.sigs)
    df.identified.mut.sigs <- df.identified.mut.sigs[order(df.identified.mut.sigs$signature),]

    # Plot contribution probabilities
    fig <- PlotSignaturesContProbs(df.identified.mut.sigs,
                                   title = title,
                                   df.ref.sigs.groups.colors = df.ref.mut.sigs.groups.colors)
    return(fig)

}


#' @title PrepareIdentifiedSignaturesPlot
#' @description Prepares identified signatures plot for reporting
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggarrange object
#'
#' @export
PrepareIdentifiedSignaturesPlot <- function(data) { # COMMON_2
    # Identify all signatures across original vcf, refined vcf, and artifact vcf
    all.signatures <- c(data$raw.muts.mutalisk.results$identified.mut.sigs,
                        data$refined.muts.mutalisk.results$identified.mut.sigs,
                        data$artifactual.muts.mutalisk.results$identified.mut.sigs)

    # Generate plots
    raw.vcf.mutalisk.plot <- PrepareIdentifiedSignaturesPlotHelper(
        mutalisk.results = data$raw.muts.mutalisk.results,
        signatures = all.signatures,
        title = "Original VCF",
        df.ref.mut.sigs.groups.colors = data$df.ref.mut.sigs.groups.colors
    )
    refined.vcf.mutalisk.plot <- PrepareIdentifiedSignaturesPlotHelper(
        mutalisk.results = data$refined.muts.mutalisk.results,
        signatures = all.signatures,
        title = "Refined VCF",
        df.ref.mut.sigs.groups.colors = data$df.ref.mut.sigs.groups.colors
    )
    artifact.vcf.mutalisk.plot <- PrepareIdentifiedSignaturesPlotHelper(
        mutalisk.results = data$artifactual.muts.mutalisk.results,
        signatures = all.signatures,
        title = "Artifactual VCF",
        df.ref.mut.sigs.groups.colors = data$df.ref.mut.sigs.groups.colors
    )

    # Merge the three plots into one
    g <- ggarrange(raw.vcf.mutalisk.plot,
                   refined.vcf.mutalisk.plot,
                   artifact.vcf.mutalisk.plot,
                   ncol = 3, nrow = 1,
                   common.legend = TRUE, legend="bottom")

    return(g)
}


#' @title PrepareTrinucleotideSpectrumsTable
#' @description Prepares trinucleotide spectrums table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareTrinucleotideSpectrumsTable <- function(data) { # COMMON_3
    #                               Original VCF    Refined VCF     Artifactual VCF
    # Mutations Count (%)               value           value           value
    # Cosine Similarity Score           value           value           value
    # Residual Sum of Squares (RSS)     value           value           value
    original.muts.count <- nrow(data$vcf.obj$data)
    refined.muts.count <- nrow(data$refined.vcf.obj$data)
    refined.muts.proportion <- (refined.muts.count / original.muts.count)
    artifactual.muts.count <- nrow(data$artifactual.vcf.obj$data)
    artifactual.muts.proportion <- (artifactual.muts.count / original.muts.count)

    df <- data.frame(
        list(" " = c("Mutations Count (%)","Cosine Similarity Score",
                     "Residual Sum of Squares (RSS)"),
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
    return(df)
}


#' @title PrepareObservedSpectrumsPlot
#' @description Prepares observed spectrums plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggarrange object
#'
#' @export
PrepareObservedSpectrumsPlot <- function(data) { # COMMON_4
    # Compute the maximum y-axis value
    trinuc.max.y <- max(c(data$raw.muts.mutalisk.results$sub.types.spectrum,
                          data$raw.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$sub.types.spectrum,
                          data$refined.muts.mutalisk.results$residuals,
                          data$artifactual.muts.mutalisk.results$sub.types.spectrum,
                          data$artifactual.muts.mutalisk.results$residuals)) * 1.3

    # Plot trinucleotide spectrums (96 substitution subtypes) - original vcf
    g1 <- PlotTriNucSpectrum(sub.types = data$raw.muts.mutalisk.results$sub.types,
                             spectrum = data$raw.muts.mutalisk.results$sub.types.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Obs Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Original VCF")
    # Plot trinucleotide spectrums (96 substitution subtypes) - refined vcf
    g2 <- PlotTriNucSpectrum(sub.types = data$refined.muts.mutalisk.results$sub.types,
                             spectrum = data$refined.muts.mutalisk.results$sub.types.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Obs Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Refined VCF")
    # Plot trinucleotide spectrums (96 substitution subtypes) - artifactual vcf
    g3 <- PlotTriNucSpectrum(sub.types = data$artifactual.muts.mutalisk.results$sub.types,
                             spectrum = data$artifactual.muts.mutalisk.results$sub.types.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Obs Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Artifactual VCF")

    fig <- ggarrange(g1, g2, g3,
                     ncol = 1, nrow = 3)
    return(fig)
}


#' @title PrepareMLEReconstructedSpectrumsPlot
#' @description Prepares MLE reconstructed spectrums plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggarrange object
#'
#' @export
PrepareMLEReconstructedSpectrumsPlot <- function(data) { # COMMON_5
    # Compute the maximum y-axis value
    trinuc.max.y <- max(c(data$raw.muts.mutalisk.results$sub.types.spectrum,
                          data$raw.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$sub.types.spectrum,
                          data$refined.muts.mutalisk.results$residuals,
                          data$artifactual.muts.mutalisk.results$sub.types.spectrum,
                          data$artifactual.muts.mutalisk.results$residuals)) * 1.3

    # Plot trinucleotide spectrums (96 substitution subtypes) - original vcf
    g1 <- PlotTriNucSpectrum(sub.types = data$raw.muts.mutalisk.results$sub.types,
                             spectrum = data$raw.muts.mutalisk.results$identified.mut.sigs.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "MLE Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Original VCF")
    # Plot trinucleotide spectrums (96 substitution subtypes) - refined vcf
    g2 <- PlotTriNucSpectrum(sub.types = data$refined.muts.mutalisk.results$sub.types,
                             spectrum = data$refined.muts.mutalisk.results$identified.mut.sigs.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "MLE Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Refined VCF")
    # Plot trinucleotide spectrums (96 substitution subtypes) - artifactual vcf
    g3 <- PlotTriNucSpectrum(sub.types = data$artifactual.muts.mutalisk.results$sub.types,
                             spectrum = data$artifactual.muts.mutalisk.results$identified.mut.sigs.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "MLE Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Artifactual VCF")

    fig <- ggarrange(g1, g2, g3,
                     ncol = 1, nrow = 3)
    return(fig)
}


#' @title PrepareResidualSpectrumsPlot
#' @description Prepares residual spectrums plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggarrange object
#'
#' @export
PrepareResidualSpectrumsPlot <- function(data) { # COMMON_6
    # Compute the maximum y-axis value
    trinuc.max.y <- max(c(data$raw.muts.mutalisk.results$sub.types.spectrum,
                          data$raw.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$sub.types.spectrum,
                          data$refined.muts.mutalisk.results$residuals,
                          data$artifactual.muts.mutalisk.results$sub.types.spectrum,
                          data$artifactual.muts.mutalisk.results$residuals)) * 1.3
    trinuc.min.y <- min(c(data$raw.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$residuals,
                          data$refined.muts.mutalisk.results$residuals))

    # Plot trinucleotide spectrums (96 substitution subtypes) - original vcf
    g1 <- PlotTriNucSpectrum(sub.types = data$raw.muts.mutalisk.results$sub.types,
                             spectrum = data$raw.muts.mutalisk.results$residuals,
                             max.y.val = trinuc.max.y,
                             min.y.val = trinuc.min.y,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Res Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Original VCF")
    # Plot trinucleotide spectrums (96 substitution subtypes) - refined vcf
    g2 <- PlotTriNucSpectrum(sub.types = data$refined.muts.mutalisk.results$sub.types,
                             spectrum = data$refined.muts.mutalisk.results$residuals,
                             max.y.val = trinuc.max.y,
                             min.y.val = trinuc.min.y,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Res Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Refined VCF")
    # Plot trinucleotide spectrums (96 substitution subtypes) - artifactual vcf
    g3 <- PlotTriNucSpectrum(sub.types = data$artifactual.muts.mutalisk.results$sub.types,
                             spectrum = data$artifactual.muts.mutalisk.results$residuals,
                             max.y.val = trinuc.max.y,
                             min.y.val = trinuc.min.y,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Res Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = "Artifactual VCF")

    fig <- ggarrange(g1, g2, g3,
                     ncol = 1, nrow = 3)
    return(fig)
}


#' @title PrepareNucleotideSubstitutionTypesPlot
#' @description Prepares nucleotide substitution types plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggarrange object
#'
#' @export
PrepareNucleotideSubstitutionTypesPlot <- function(data) { # COMMON_7
    # Compute mutation type max y.axis value
    mut.type.max.y <- max(c(
        EnumerateTriNucCounts(
            data$raw.muts.mutalisk.results$sub.types.spectrum),
        EnumerateTriNucCounts(
            data$refined.muts.mutalisk.results$sub.types.spectrum),
        EnumerateTriNucCounts(
            data$artifactual.muts.mutalisk.results$sub.types.spectrum)
    )) * 1.1

    # Plot mutation types - original vcf
    g1 <- PlotMutationTypes(mutation.types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                            mutation.types.values = EnumerateTriNucCounts(
                                data$raw.muts.mutalisk.results$sub.types.spectrum),
                            mutation.types.colors = TriNuc.Mutation.Type.Hex.Colors,
                            max.y.val = mut.type.max.y,
                            convert.to.percentage = T,
                            show.legend = F,
                            title = "Original VCF")
    # Plot mutation types - refined vcf
    g2 <- PlotMutationTypes(mutation.types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                            mutation.types.values = EnumerateTriNucCounts(
                                data$refined.muts.mutalisk.results$sub.types.spectrum),
                            mutation.types.colors = TriNuc.Mutation.Type.Hex.Colors,
                            max.y.val = mut.type.max.y,
                            convert.to.percentage = T,
                            show.legend = F,
                            title = "Refined VCF")
    # Plot mutation types - artifactual vcf
    g3 <- PlotMutationTypes(mutation.types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                            mutation.types.values = EnumerateTriNucCounts(
                                data$artifactual.muts.mutalisk.results$sub.types.spectrum),
                            mutation.types.colors = TriNuc.Mutation.Type.Hex.Colors,
                            max.y.val = mut.type.max.y,
                            convert.to.percentage = T,
                            show.legend = F,
                            title = "Artifactual VCF")

    fig <- ggarrange(g1, g2, g3,
                     ncol = 3, nrow = 1)
    return(fig)
}


#' @title PrepareOptimizedVCFStatisticsPlot
#' @description Prepares optimized VCF statistics plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggarrange object
#'
#' @export
PrepareOptimizedVCFStatisticsPlot <- function(data) { # COMMON_8
    # Get VCF statistics
    original.vcf.stats <- GetVCFValues(vcf.obj = data$vcf.obj,
                                       vcf.filter = data$vcf.filter)
    refined.vcf.stats <- GetVCFValues(vcf.obj = data$refined.vcf.obj,
                                      vcf.filter = data$vcf.filter)
    artifact.vcf.stats <- GetVCFValues(vcf.obj = data$artifactual.vcf.obj,
                                       vcf.filter = data$vcf.filter)

    x.axis.labels <- names(data$vcf.filter)

    # Get the maximum frequency value for each VCF statistic variable
    stat.y.max.vals <- c()
    stat.x.max.vals <- c()
    for (param in x.axis.labels) {
        stat.x.max.vals <- c(stat.x.max.vals, max(max(original.vcf.stats[[param]][is.finite(original.vcf.stats[[param]])]),
                                                  max(refined.vcf.stats[[param]][is.finite(refined.vcf.stats[[param]])]),
                                                  max(artifact.vcf.stats[[param]][is.finite(artifact.vcf.stats[[param]])])))
        orig.hist.data <- hist(original.vcf.stats[[param]][is.finite(original.vcf.stats[[param]])],
                               plot = F,
                               breaks = seq(0, max(original.vcf.stats[[param]][is.finite(original.vcf.stats[[param]])]) + 1, by = 1))$counts
        refi.hist.data <- hist(refined.vcf.stats[[param]][is.finite(refined.vcf.stats[[param]])],
                               plot = F,
                               breaks = seq(0, max(refined.vcf.stats[[param]][is.finite(refined.vcf.stats[[param]])]) + 1, by = 1))$counts
        arti.hist.data <- hist(artifact.vcf.stats[[param]][is.finite(artifact.vcf.stats[[param]])],
                               plot = F,
                               breaks = seq(0, max(artifact.vcf.stats[[param]][is.finite(artifact.vcf.stats[[param]])]) + 1, by = 1))$counts
        stat.y.max.vals <- c(stat.y.max.vals,
                             max(c(max(orig.hist.data),
                                   max(refi.hist.data),
                                   max(arti.hist.data))))
    }

    # Plot VCF statistics - original vcf
    original.vcf.stats.plots <- PlotVCFStatsHistograms(
        plot.values = original.vcf.stats,
        x.axis.labels = x.axis.labels,
        sample.id = data$vcf.file.basename,
        title = "Original VCF",
        cutoff.values = data$x.solution.decimal,
        plot.cutoff.value.lines = T,
        stat.y.max.vals = stat.y.max.vals,
        stat.x.max.vals = stat.x.max.vals,
        save.file = paste0(
            data$output.dir,
            data$vcf.file.basename,
            "_FIREVAT_Original_VCF_Stats_Histograms.png"
        )
    )
    # Plot VCF statistics - refined vcf
    refined.vcf.stats.plots <- PlotVCFStatsHistograms(
        plot.values = refined.vcf.stats,
        x.axis.labels = x.axis.labels,
        title = "Refined VCF",
        plot.cutoff.line.color = "gainsboro",
        cutoff.values = data$x.solution.decimal,
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
    # Plot VCF statistics - artifactual vcf
    artifact.vcf.stats.plots <- PlotVCFStatsHistograms(
        plot.values = artifact.vcf.stats,
        x.axis.labels = x.axis.labels,
        title = "Artifactual VCF",
        plot.cutoff.line.color = "gainsboro",
        cutoff.values = data$x.solution.decimal,
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
    # Merge the plots
    vcf.stats.boxplots <- list()
    for (i in 1:length(x.axis.labels)) {
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

    # Boxplots
    vcf.stat.plots <- list()
    list.idx <- 1
    for (i in 1:length(original.vcf.stats.plots$graphs))  {
        vcf.stat.plots[[list.idx]] <- original.vcf.stats.plots$graphs[[i]]
        list.idx <- list.idx + 1

        vcf.stat.plots[[list.idx]] <- refined.vcf.stats.plots$graphs[[i]]
        list.idx <- list.idx + 1

        vcf.stat.plots[[list.idx]] <- artifact.vcf.stats.plots$graphs[[i]]
        list.idx <- list.idx + 1

        vcf.stat.plots[[list.idx]] <- vcf.stats.boxplots[[i]]
        list.idx <- list.idx + 1
    }

    fig <- ggarrange(plotlist = vcf.stat.plots,
                     ncol = 4,
                     nrow = length(original.vcf.stats.plots$graphs),
                     widths = c(2,2,2,1.5),
                     align = "h")
    return(list(fig = fig,
                params.count = length(original.vcf.stats.plots$graphs)))
}


#' @title PrepareGeneticAlgorithmParametersTable
#' @description Prepares Genetic Algorithm parameters table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareGeneticAlgorithmParametersTable <- function(data) { # GA_1
    # FIREVAT Genetic Algorithm (GA) Parameters
    # GA Population Size                                value
    # GA Maximum Iteration                              value
    # GA Run                                            value
    # GA Mutation Probability                           value
    #
    df <- data.frame(
        list("FIREVAT Genetic Algorithm (GA) Parameters" = c(
            "GA Population Size", "GA Maximum Iteration", "GA Run",
            "GA Mutation Probability"),
            " " = c(format(x = data$ga.pop.size, nsmall = 0, big.mark = ","),
                    format(x = data$ga.max.iter, nsmall = 0, big.mark = ","),
                    format(x = data$ga.run, nsmall = 0, big.mark = ","),
                    format(x = data$ga.pmutation, nsmall = 3))),
        stringsAsFactors = F,
        check.names = F)
    return(df)
}


#' @title PrepareOptimizationResultsTable
#' @description Prepares optimization results table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareOptimizationResultsTable <- function(data) { # Optional_1
    # Get optimization results
    df <- ReadOptimizationIterationReport(data = data)
    refined.muts.cos.sim <- tail(df$refined.muts.cosine.similarity.score, 1)
    refined.muts.count <- tail(df$refined.mutations.count, 1)
    refined.muts.proportion <- tail(df$refined.muts.proportion, 1)
    refined.muts.seq.art.sigs.weights.sum <- tail(df$refined.muts.sequencing.artifact.signatures.weights.sum, 1)
    artifactual.muts.cos.sim <- tail(df$artifactual.muts.cosine.similarity.score, 1)
    artifactual.muts.count <- tail(df$artifact.mutations.count, 1)
    artifactual.muts.proportion <- tail(df$artifactual.muts.proportion, 1)
    artifactual.muts.seq.art.sigs.weights.sum <- tail(df$artifactual.muts.sequencing.artifact.signatures.weights.sum, 1)
    objective.value <- tail(df$objective.value, 1)

    # Objective Value   C.refined   W.refined   C.artifact  W.artifacdt
    #       value         value       value       value       value
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
    return(df.results)
}


#' @title PrepareRefinedMutsOptimizationIterationsPlot
#' @description Prepares refined mutations optimization iterations plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggplot object
#'
#' @export
PrepareRefinedMutsOptimizationIterationsPlot <- function(data) { # Optional_2
    # Get optimization results
    df <- ReadOptimizationIterationReport(data = data)
    df.temp <- df[,c(
        "iteration",
        "refined.muts.sequencing.artifact.signatures.weights.sum",
        "refined.muts.cosine.similarity.score",
        "refined.muts.proportion",
        "objective.value"
    )]
    max.val <- max(df.temp$objective.value,
                   df.temp$refined.muts.sequencing.artifact.signatures.weights.sum,
                   df.temp$refined.muts.cosine.similarity.score,
                   df.temp$refined.muts.proportion)
    colnames(df.temp) <- c(
        "iteration",
        "Refined Mutations Artifact Signatures Weights Sum",
        "Refined Mutations Cosine Similarity Score",
        "Refined Mutations Proportion",
        "Objective Value"
    )

    # Plot refined mutations optimization iterations
    columns.to.plot <- colnames(df.temp)[2:length(colnames(df.temp))]
    optimization.iter.plot <- PlotOptimizationIterations(
        df = df.temp,
        columns.to.plot = columns.to.plot,
        x.axis.var = "iteration",
        x.axis.title = "Iteration",
        y.axis.title = "Value",
        y.max = 1.3 * max.val,
        title = "Objective Function Optimization",
        x.min = min(df.temp$iteration),
        x.max = max(df.temp$iteration),
        legend.ncol = 2,
        save.file = paste0(data$output.dir,
                           data$vcf.file.basename,
                           "_FIREVAT_Optimization_Iterations.png"),
        connect.dots = T
    )
    return(optimization.iter.plot)
}


#' @title PrepareArtifactualMutsOptimizationIterationsPlot
#' @description Prepares artifactual mutations optimization iterations plot
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A ggplot object
#'
#' @export
PrepareArtifactualMutsOptimizationIterationsPlot <- function(data) { # Optional_3
    # Get optimization results
    df.optimization.logs <- ReadOptimizationIterationReport(data = data)

    # Enumerate sequencing artifact signature weights at each iteration
    # Get all sequencing artifact signatures
    artifactual.muts.seq.artifact.sigs <- c()
    for (i in 1:nrow(df.optimization.logs)) {
        curr.row <- df.optimization.logs[i,]
        curr.row.seq.art.sigs <- as.character(curr.row$artifactual.muts.sequencing.artifact.signatures)
        if (is.na(curr.row.seq.art.sigs) == FALSE) {
            curr.row.seq.art.sigs <- strsplit(curr.row.seq.art.sigs, "\\,")[[1]]
            artifactual.muts.seq.artifact.sigs <- c(artifactual.muts.seq.artifact.sigs,
                                                    curr.row.seq.art.sigs)
        }
    }

    if (length(artifactual.muts.seq.artifact.sigs) == 0) {
        return(NULL)
    }

    artifactual.muts.seq.artifact.sigs <- unique(artifactual.muts.seq.artifact.sigs)
    # Create a dataframe with columns 'iteration', 'SBS<NUM>', 'SBS<NUM>', ...
    df.plot.data <- data.frame(
        matrix(0, ncol = 1 + length(artifactual.muts.seq.artifact.sigs), nrow = 0)
    )
    colnames(df.plot.data) <- c("iteration", artifactual.muts.seq.artifact.sigs)
    # Fill in data for each iteration
    for (i in 1:nrow(df.optimization.logs)) {
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

        if (is.na(curr.row.seq.art.sigs) == FALSE) {
            curr.row.seq.art.sigs <-as.character(
                strsplit(curr.row.seq.art.sigs, "\\,")[[1]])
            curr.row.seq.art.sigs.probs <- as.numeric(
                strsplit(curr.row.seq.art.sigs.probs, "\\,")[[1]])
            for (j in 1:length(curr.row.seq.art.sigs)) {
                curr.sig <- curr.row.seq.art.sigs[j]
                curr.sig.prob <- curr.row.seq.art.sigs.probs[j]
                df.plot.data.temp[curr.sig] <- curr.sig.prob
            }
        }
        df.plot.data <- rbind(df.plot.data, df.plot.data.temp)
    }

    # Plot sequencing artifact weight (artifactual mutations)
    artifactual.muts.seq.art.iter.plot <- PlotOptimizationIterations(
        df = df.plot.data,
        columns.to.plot = c(artifactual.muts.seq.artifact.sigs),
        x.axis.var = "iteration",
        x.axis.title = "Iteration",
        y.axis.title = "Value",
        x.min = min(df.plot.data$iteration),
        x.max = max(df.plot.data$iteration),
        y.max = max(df.plot.data[,artifactual.muts.seq.artifact.sigs]) * 1.3,
        title = paste0(
            "Artifactual Mutations Sequencing Artifact ",
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
    return(artifactual.muts.seq.art.iter.plot)
}


#' @title PrepareRefinedAnnotationTable
#' @description Prepares refined mutations annotation (filtered, queried) table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareRefinedAnnotationTable <- function(data) { # Optional_4
    if (data$annotate == FALSE) {
        return(data.frame())
    }

    # CHROM POS REF ALT <data$annotated.columns.to.display>
    cols.to.display <- unique(c(c("CHROM", "POS", "REF", "ALT"),
                                data$annotated.columns.to.display))
    df <- as.data.frame(data$refined.vcf.obj.annotated.queried$data)
    df <- df[, cols.to.display]
    df$POS <- format(x = df$POS, big.mark = ",")
    return(df)
}


#' @title PrepareArtifactAnnotationTable
#' @description Prepares artifactual mutations annotation (filtered, queried) table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareArtifactAnnotationTable <- function(data) { # Optional_5
    if (data$annotate == FALSE) {
        return(data.frame())
    }

    # CHROM POS REF ALT <data$annotated.columns.to.display>
    cols.to.display <- unique(c(c("CHROM", "POS", "REF", "ALT"),
                                data$annotated.columns.to.display))
    df <- as.data.frame(data$artifactual.vcf.obj.annotated.queried$data)
    df <- df[, cols.to.display]
    df$POS <- format(x = df$POS, big.mark = ",")
    return(df)
}


#' @title PrepareRefinedStrandBiasTable
#' @description Prepares refined mutations strand biased variants table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareRefinedStrandBiasTable <- function(data) {
    if (data$perform.strand.bias.analysis == FALSE) {
        return(data.frame())
    }

    if (data$strand.bias.perform.fdr.correction == TRUE) {
        include.array <- data$refined.vcf.obj$data$StrandBiasQValue < data$filter.by.strand.bias.analysis.cutoff
        vcf.objs <- FilterVCF(vcf.obj = data$refined.vcf.obj,
                              include.array = include.array,
                              force.include = TRUE,
                              verbose = FALSE)
        vcf.obj <- vcf.objs$vcf.obj.filtered
        stat.sig.value.col <- "StrandBiasQValue"
    } else {
        include.array <- data$refined.vcf.obj$data$StrandBiasPValue < data$filter.by.strand.bias.analysis.cutoff
        vcf.objs <- FilterVCF(vcf.obj = data$refined.vcf.obj,
                              include.array = include.array,
                              force.include = TRUE,
                              verbose = FALSE)
        vcf.obj <- vcf.objs$vcf.obj.filtered
        stat.sig.value.col <- "StrandBiasPValue"
    }

    cols.to.display <- c("CHROM", "POS", "REF", "ALT",
                         data$ref.forward.strand.var,
                         data$ref.reverse.strand.var,
                         data$alt.forward.strand.var,
                         data$alt.reverse.strand.var,
                         stat.sig.value.col)
    df <- as.data.frame(vcf.obj$data)
    df <- df[, cols.to.display]
    df$POS <- format(x = df$POS, big.mark = ",")
    df <- df[order(df[stat.sig.value.col]),]
    return(df)
}


#' @title PrepareArtifactStrandBiasTable
#' @description Prepares artifactual mutations strand biased variants table
#'
#' @param data A list of elements returned from \code{\link{RunFIREVAT}}
#'
#' @return A data.frame
#'
#' @export
PrepareArtifactStrandBiasTable <- function(data) {
    if (data$perform.strand.bias.analysis == FALSE) {
        return(data.frame())
    }
    if (data$strand.bias.perform.fdr.correction == TRUE) {
        include.array <- data$artifactual.vcf.obj$data$StrandBiasQValue < 0.05
        vcf.objs <- FilterVCF(vcf.obj = data$artifactual.vcf.obj,
                              include.array = include.array,
                              force.include = TRUE,
                              verbose = FALSE)
        vcf.obj <- vcf.objs$vcf.obj.filtered
        stat.sig.value.col <- "StrandBiasQValue"
    } else {
        include.array <- data$artifactual.vcf.obj$data$StrandBiasPValue < 0.05
        vcf.objs <- FilterVCF(vcf.obj = data$artifactual.vcf.obj,
                              include.array = include.array,
                              force.include = TRUE,
                              verbose = FALSE)
        vcf.obj <- vcf.objs$vcf.obj.filtered
        stat.sig.value.col <- "StrandBiasPValue"
    }

    cols.to.display <- c("CHROM", "POS", "REF", "ALT",
                         data$ref.forward.strand.var,
                         data$ref.reverse.strand.var,
                         data$alt.forward.strand.var,
                         data$alt.reverse.strand.var,
                         stat.sig.value.col)
    df <- as.data.frame(vcf.obj$data)
    df <- df[, cols.to.display]
    df$POS <- format(x = df$POS, big.mark = ",")
    df <- df[order(df[stat.sig.value.col]),]
    return(df)
}


#' @title ReportFIREVATResults
#' @description Reports FIREVAT results in html format (generated from Rmd)
#'
#' @param data A list of main data from \code{\link{RunFIREVAT}}
#'
#' @return An updated data list
#'
#' @importFrom rmarkdown render
#' @export
ReportFIREVATResults <- function(data) {
    if (data$verbose == TRUE) {
        PrintLog("* Started generating FIREVAT report")
    }

    # Generate html report
    if (data$report.format == "html") {
        if (data$mode == "ga") {
            report.items <- list(df.genetic.algo.params = PrepareGeneticAlgorithmParametersTable(data = data),
                                 df.filter.cutoffs = PrepareFilterCutoffsTable(data = data),
                                 df.optimization.results = PrepareOptimizationResultsTable(data = data),
                                 refined.muts.optimization.iter.plot = PrepareRefinedMutsOptimizationIterationsPlot(data = data),
                                 artifactual.muts.optimization.iter.plot = PrepareArtifactualMutsOptimizationIterationsPlot(data = data),
                                 identified.signatures.plot = PrepareIdentifiedSignaturesPlot(data = data),
                                 df.trinucleotide.spectrums = PrepareTrinucleotideSpectrumsTable(data = data),
                                 observed.spectrums.plot = PrepareObservedSpectrumsPlot(data = data),
                                 mle.reconstructed.spectrums.plot = PrepareMLEReconstructedSpectrumsPlot(data = data),
                                 residual.spectrums.plot = PrepareResidualSpectrumsPlot(data = data),
                                 nucleotide.substitution.types.plot = PrepareNucleotideSubstitutionTypesPlot(data = data),
                                 vcf.stats.plot = PrepareOptimizedVCFStatisticsPlot(data = data),
                                 df.refined.vcf.annotated = PrepareRefinedAnnotationTable(data = data),
                                 df.artifact.vcf.annotated = PrepareArtifactAnnotationTable(data = data),
                                 df.refined.vcf.strand.bias = PrepareRefinedStrandBiasTable(data = data),
                                 df.artifact.vcf.strand.bias = PrepareArtifactStrandBiasTable(data = data))
            rmarkdown.file <- system.file("rmd_template", "report_ga.Rmd", package = "FIREVAT")
        } else if (data$mode == "manual") {
            report.items <- list(df.filter.cutoffs = PrepareFilterCutoffsTable(data = data),
                                 identified.signatures.plot = PrepareIdentifiedSignaturesPlot(data = data),
                                 df.trinucleotide.spectrums = PrepareTrinucleotideSpectrumsTable(data = data),
                                 observed.spectrums.plot = PrepareObservedSpectrumsPlot(data = data),
                                 mle.reconstructed.spectrums.plot = PrepareMLEReconstructedSpectrumsPlot(data = data),
                                 residual.spectrums.plot = PrepareResidualSpectrumsPlot(data = data),
                                 nucleotide.substitution.types.plot = PrepareNucleotideSubstitutionTypesPlot(data = data),
                                 vcf.stats.plot = PrepareOptimizedVCFStatisticsPlot(data = data),
                                 df.refined.vcf.annotated = PrepareRefinedAnnotationTable(data = data),
                                 df.artifact.vcf.annotated = PrepareArtifactAnnotationTable(data = data),
                                 df.refined.vcf.strand.bias = PrepareRefinedStrandBiasTable(data = data),
                                 df.artifact.vcf.strand.bias = PrepareArtifactStrandBiasTable(data = data))
            rmarkdown.file <- system.file("rmd_template", "report_manual.Rmd", package = "FIREVAT")
        }

        rmarkdown::render(input = rmarkdown.file,
                          params = list(data = data,
                                        report.items = report.items),
                          output_dir = data$output.dir,
                          intermediates_dir = data$output.dir,
                          output_file = paste0(data$vcf.file.basename,
                                               "_FIREVAT_Report.html"))
    }

    if (data$verbose == TRUE) {
        PrintLog("* Finished generating FIREVAT report")
    }

    data$report.items <- report.items
    return(data)
}
