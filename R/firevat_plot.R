# FIREVAT Plot Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title PlotVCFStatsHistograms
#' @description Plots multiple VCF stats histograms into one figure
#'
#' @param plot.values A list of multiple numeric vectors
#' @param x.axis.labels A character vector of x axis labels
#' @param stat.y.max.vals A numeric vector of max y-axis values
#' @param stat.x.max.vals A numeric vector of max x-axis values
#' @param sample.id A string value of sample ID
#' @param save.file A string value of file to which the resulting plot will be saved
#' @param title A string value of plot title
#' @param cutoff.values A numeric vector of cutoff values
#' @param plot.boxplot A boolean value (default = False)
#' @param plot.cutoff.line.color A hex string value (default = "#D4012E")
#' @param plot.cutoff.value.lines A boolean value (default = False)
#' @param bin.width An integer value (default = 1; histogram bin width)
#' @param ncol An integer value (default = 4; ggarrange ncol)
#' @param nrow An integer value (default = 3; ggarrange nrow)
#' @param font.size.med An integer value (default = 10)
#' @param font.size.large An integer value (default = 12)
#' @param plot.margin A list (default = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
#'
#' @return A list with the following elements
#' \itemize{
#'  \item{f}{ = A ggarrange object}
#'  \item{graphs}{ = A list of length 3; each element is a ggplot histogram}
#' }
#'
#' @import ggplot2
#' @import ggpubr
#' @export
PlotVCFStatsHistograms <- function(plot.values,
                                   x.axis.labels,
                                   stat.y.max.vals,
                                   stat.x.max.vals,
                                   sample.id,
                                   save.file,
                                   title,
                                   cutoff.values,
                                   plot.boxplot = F,
                                   plot.cutoff.line.color = "#D4012E",
                                   plot.cutoff.value.lines = F,
                                   bin.width = 1,
                                   ncol = 4,
                                   nrow = 3,
                                   font.size.med = 10,
                                   font.size.large = 12,
                                   plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) {
    PrintLog("** Started plotting vcf stats histograms")
    graphs <- list()
    i <- 1
    for (curr.plot.values in plot.values) {
        df <- data.frame('x' = curr.plot.values,
                         stringsAsFactors = F)
        g <- ggplot(df, aes_string(x = 'x')) +
            ggtitle(paste0("\n",title)) +
            geom_histogram(binwidth = bin.width,
                           color="black", fill="white") +
            xlab("") +
            ylab("Count\n") +
            theme_minimal() + theme_classic() +
            theme(plot.title = element_text(size = font.size.large,
                                            hjust = 0.5),
                  axis.text.x = element_text(size = font.size.med),
                  axis.text.y = element_text(size = font.size.med),
                  axis.title.x = element_text(size = font.size.large),
                  axis.title.y = element_text(size = font.size.med)) +
            scale_y_continuous(limits = c(0, stat.y.max.vals[i])) +
            scale_x_continuous(limits = c(-0.5, stat.x.max.vals[i] + 1),
                               expand = c(0.1, 0.1))

        if (plot.cutoff.value.lines) {
            g <- g + geom_vline(xintercept = cutoff.values[i], colour = plot.cutoff.line.color)
        }

        if (plot.boxplot) {
            p <- ggplot(aes_string(x = "", y = 'x'), data = df) +
                stat_boxplot(geom ='errorbar', width = 0.4) +
                geom_boxplot(outlier.alpha = 0) +
                coord_flip() +
                theme_minimal() + theme_classic() +
                theme(axis.text.x = element_text(size = font.size.med),
                      axis.text.y = element_text(size = font.size.med),
                      axis.title.y = element_text(size = font.size.med),
                      axis.title.x = element_text(size = font.size.med),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank()) +
                scale_y_continuous(limits = c(-0.5, stat.x.max.vals[i]),
                                   expand = c(0.1, 0.1)) +
                stat_summary(fun.y = mean, colour = "purple", geom = "point",
                             shape = 18, size = 3) +
                labs(x = "", y = xlab(paste0(x.axis.labels[i],"\n")))

            g <- ggarrange(g, p, heights = c(3,1), nrow = 2, align = "v")
        } else {
            g <- g + xlab(paste0(x.axis.labels[i],"\n"))
        }
        graphs[[i]] <- g
        i <- i + 1
    }

    f <- ggarrange(plotlist = graphs, ncol = ncol, nrow = nrow)
    title <- paste0(sample.id, " (n=", length(plot.values[[1]]), ")")
    annotate_figure(f, top = text_grob(title,
                                       color = "black",
                                       size = font.size.large))
    ggsave(save.file, width = 16, height = 9)

    PrintLog("** Finished plotting vcf stats histograms")
    return(list(f = f, graphs = graphs))
}


#' @title PlotVCFStatsBoxPlots
#' @description Plots multiple (original, refined, artifact vcf) boxplots for single filter parameter
#'
#' @param original.vcf.stat.values A numeric vector corresponding to the original vcf.obj values of single filter parameter
#' @param refined.vcf.stat.values A numeric vector corresponding to the refined vcf.obj values of single filter parameter
#' @param artifact.vcf.stat.values A numeric vector corresponding to the artifact vcf.obj values of single filter parameter
#' @param xlab A string value (x-axis label)
#' @param axis.font.size An integer value (axis font size)
#' @param label.font.size An integer value (label font size)
#' @param title.font.size An integer value (title font size)
#'
#' @return A ggboxplot
#'
#' @import ggplot2
#' @export
PlotVCFStatsBoxPlots <- function(original.vcf.stat.values,
                                 refined.vcf.stat.values,
                                 artifact.vcf.stat.values,
                                 xlab,
                                 axis.font.size = 10,
                                 label.font.size = 10,
                                 title.font.size = 12) {
    PrintLog("** Started plotting vcf stats boxplots")

    # Prepare plot data
    df1 <- data.frame(type = rep('Original', length(original.vcf.stat.values)),
                      value = original.vcf.stat.values,
                      stringsAsFactors = F,
                      check.names = F)
    df2 <- data.frame(type = rep('Refined', length(refined.vcf.stat.values)),
                      value = refined.vcf.stat.values,
                      stringsAsFactors = F,
                      check.names = F)
    df3 <- data.frame(type = rep('Artifact', length(artifact.vcf.stat.values)),
                      value = artifact.vcf.stat.values,
                      stringsAsFactors = F,
                      check.names = F)
    df <- rbindlist(list(df1, df2, df3))

    # Get y-axis min and max
    original.quantiles <- quantile(original.vcf.stat.values)
    original.iqr <- original.quantiles[4] - original.quantiles[2]

    refined.quantiles <- quantile(refined.vcf.stat.values)
    refined.iqr <- refined.quantiles[4] - refined.quantiles[2]

    artifact.quantiles <- quantile(artifact.vcf.stat.values)
    artifactd.iqr <- artifact.quantiles[4] - artifact.quantiles[2]

    # y minimum = max((Q1 - 1.5*IQR) * 0.9, globa.y.min)
    global.y.min <- min(original.vcf.stat.values,
                        refined.vcf.stat.values,
                        artifact.vcf.stat.values)
    original.y.min <- original.quantiles[2] - (1.5 * original.iqr)
    refined.y.min <- refined.quantiles[2] - (1.5 * refined.iqr)
    artifact.y.min <- artifact.quantiles[2] - (1.5 * artifactd.iqr)
    y.min <- max(min(c(original.y.min,
                       refined.y.min,
                       artifact.y.min)) * 0.9,
                 global.y.min)

    # y maxmimum = max(Q3 + 1.5*IQR)
    global.y.max <- max(original.vcf.stat.values,
                        refined.vcf.stat.values,
                        artifact.vcf.stat.values)
    original.y.max <- original.quantiles[4] + (1.5 * original.iqr)
    refined.y.max <- refined.quantiles[4] + (1.5 * refined.iqr)
    artifact.y.max <- artifact.quantiles[4] + (1.5 * artifactd.iqr)
    y.max <- max(original.y.max,
                 refined.y.max,
                 artifact.y.max,
                 1)
    y.top.buffer <- ((y.max - y.min) * 1.3)

    # Plot
    comparisons <- list(c("Original", "Refined"),
                        c("Refined", "Artifact"),
                        c("Original", "Artifact"))
    g <- ggboxplot(df, x = "type", y = "value", outlier.shape = NA,
                   title = "\nComparisons",
                   palette = c("#000000", "#000000", "#000000")) +
        font("x.text", size = axis.font.size) +
        font("y.text", size = axis.font.size) +
        font("title", color = "black", size = title.font.size, hjust = 0.5) +
        font("ylab", size = axis.font.size) +
        font("xlab", size = title.font.size) +
        ylim(c(y.min, y.max + y.top.buffer)) + xlab(paste0(xlab,"\n")) +
        stat_compare_means(comparisons = comparisons,
                           label.y = c(y.max + y.top.buffer * 0.2,
                                       y.max + y.top.buffer * 0.4,
                                       y.max + y.top.buffer * 0.6),
                           label = "p.signif",
                           tip.length = 0) +
        stat_compare_means(label.y = y.max + y.top.buffer * 0.8,
                           label = "p.format") +
        rremove("legend") + rremove("ylab")

    PrintLog("** Finished plotting vcf stats boxplots")
    return(g)
}


#' @title PlotOptimizationIterations
#' @description Plots multiple scatter plots into one figure
#'
#' @param df A data.frame (from reading "FIREVAT_Optimization_Logs.tsv")
#' @param columns.to.plot A character vector (of column names to plot)
#' @param x.axis.var x axis variable
#' @param x.axis.title x axis title
#' @param x.min x axis minimum value
#' @param x.max x axis maximum value
#' @param save.file Filename (including full path) to which the plot will be saved
#' @param title Plot title
#' @param y.axis.title y axis title; Default = ""
#' @param y.max y axis maximum value; Default = 1
#' @param point.size Point size; Default = 1
#' @param connect.dots If True draws dots for each iteration; Default = True
#' @param plot.legend If True write legend of plot; Default = T
#' @param legend.ncol legend.n Default = 1
#' @param font.size.med Medium font size; Default = 14
#' @param font.size.large Large font size; Default = 16
#' @param plot.margin Margin vector for plot; Default = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import ggpubr
#' @importFrom data.table melt
#' @export
PlotOptimizationIterations <- function(df,
                                       columns.to.plot,
                                       x.axis.var,
                                       x.axis.title,
                                       x.min,
                                       x.max,
                                       save.file,
                                       title,
                                       y.axis.title = "",
                                       y.max = 1,
                                       point.size = 1,
                                       connect.dots = T,
                                       plot.legend = T,
                                       legend.ncol = 1,
                                       font.size.med = 14,
                                       font.size.large = 16,
                                       plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) {
    PrintLog("** Started plotting optimization iterations")

    df.temp <- df[,c(x.axis.var, columns.to.plot)]
    df.temp[,x.axis.var] <- as.numeric(as.character(df.temp[, x.axis.var]))
    m <- melt(df.temp, id.vars = x.axis.var)
    g <- ggplot(m, aes_string(x = x.axis.var, y = "value", colour = "variable")) +
        geom_point(size = point.size) +
        geom_line(size = 1) +
        ggtitle(title) +
        xlab(paste0("\n", x.axis.title)) +
        ylab(y.axis.title) +
        theme_minimal() +
        theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
              axis.text.x = element_text(size = font.size.med),
              axis.text.y = element_text(size = font.size.med),
              axis.title.x = element_text(size = font.size.med),
              axis.title.y = element_text(size = font.size.med),
              plot.title = element_text(size = font.size.large, hjust = 0.5),
              axis.ticks.length = unit(.25, "cm"),
              plot.margin = plot.margin,
              axis.ticks = element_line(colour = 'black', size = 0.5)) +
        scale_y_continuous(limits = c(0, y.max), expand = c(0, 0)) +
        scale_x_continuous(limits = c(x.min, x.max), expand = c(0, 0))

    if (plot.legend == FALSE) {
        g <- g + guides(fill = FALSE) +
            theme(legend.position = "none")
    } else {
        g <- g +
            guides(color = guide_legend(ncol = legend.ncol,
                                        nrow = ceiling(length(unique(columns.to.plot)) / legend.ncol))) +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = font.size.med),
                  legend.title = element_text(size = 0),
                  legend.key.size = unit(2, 'lines'),
                  legend.key.width = unit(2, 'lines'),
                  legend.key.height = unit(1.5, 'lines'),
                  legend.spacing = unit(2, 'lines'))
    }

    # Save file
    ggsave(save.file, width = 16, height = 9)

    PrintLog("** Finished plotting optimization iterations")
    return(g)
}


#' @title PlotTriNucSpectrum
#' @description
#' Plots the spectrum of 96 trinucleotide distribution
#' (C>A, C>G, C>T, T>A, T>C, T>G)
#' Please note that this function assumes that both sub.types and spectrum
#' are sorted in the following order: C>A, C>G, C>T, T>A, T>C, T>G
#'
#' @param sub.types A character vector (types of 96 trinucleotide substitutions)
#' @param spectrum A numeric vector (96 elements)
#' @param max.y.val y axis maximum value
#' @param min.y.val y axis minimum value
#' @param y.axis.title y axis title
#' @param draw.top.strip If True then draws top strip; Default = T
#' @param draw.x.axis.labels If True then draws x axis labels; Default = T
#' @param draw.y.axis.labels If True then draws y axis labels; Default = T
#' @param draw.y.axis.title If True then draws y axis title; Default = T
#' @param font.size.small Small font size; Default = 8
#' @param font.size.med Medium font size; Default = 14
#' @param plot.margin.top Top margin; Default = 0.5
#' @param plot.margin.bottom Bottom margin; Default = 0.5
#' @param plot.margin.left Left margin; Default = 0.5
#' @param plot.margin.right Right margin; Default = 0.5
#' @param title Plot title
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import ggpubr
#' @import extrafont
#' @export
PlotTriNucSpectrum <- function(sub.types,
                               spectrum,
                               max.y.val,
                               min.y.val,
                               y.axis.title,
                               draw.top.strip = T,
                               draw.x.axis.labels = T,
                               draw.y.axis.labels = T,
                               draw.y.axis.title = T,
                               font.size.small = 8,
                               font.size.med = 14,
                               plot.margin.top = 0.5,
                               plot.margin.bottom = 0.5,
                               plot.margin.left = 0.5,
                               plot.margin.right = 0.5,
                               title) {
    # extrafont
    loadfonts()

    PrintLog("** Started plotting trinucleotide spectrum")

    # Prepare plot data
    plot.colors <- TriNuc.Mutation.Type.Hex.Colors
    mutation.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

    df.plot.data <- list()
    for (i in 1:6) {
        start.idx <- ((i - 1) * 16) + 1
        end.idx <- i * 16
        df.temp <- data.frame(
            list("mutation.subtype" = sub.types[start.idx:end.idx],
                 "value" = spectrum[start.idx:end.idx],
                 "title" = rep(mutation.types[i], 16)),
            stringsAsFactors = F)
        df.plot.data[[i]] <- df.temp
    }

    # Plot
    if (draw.x.axis.labels) {
        axis.text.x <- element_text(size = font.size.small,
                                    family = "Arial",
                                    colour = "#888888",
                                    angle = 90,
                                    vjust = 0.5,
                                    hjust = 0.5)
    } else {
        axis.text.x <- element_blank()
    }

    if (draw.y.axis.labels) {
        axis.text.y <- element_text(size = font.size.med - 4,
                                    family = "Arial",
                                    colour = "#888888")
    } else {
        axis.text.y <- element_blank()
    }

    strip.text.x.colors <- c("white", "white", "white", "white", "white", "white")

    plots <- list()
    for (i in 1:6) {
        p <- ggplot(df.plot.data[[i]], aes_string(x = 'mutation.subtype', y = 'value')) +
            geom_bar(stat = "identity", width = 0.9,
                     fill = plot.colors[i]) + ylab(NULL) + xlab("") +
            scale_y_continuous(limits = c(min.y.val, max.y.val),
                           expand = c(0,0)) +
            theme(axis.text.x = axis.text.x,
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid.major.y = element_line(colour = "#DCDCDC",
                                                    linetype = "solid"),
                  panel.grid.major.x = element_blank(),
                  panel.background = element_blank())

        if (draw.top.strip) { # draw.top.strip
            p <- p +
                theme(strip.background = element_rect(fill = plot.colors[i],
                                                      colour = plot.colors[i]),
                      strip.text.x = element_text(size = font.size.med - 2,
                                                  family = "Arial",
                                                  face = "bold",
                                                  colour = strip.text.x.colors[i])) +
                facet_grid(~title)
        }

        if (i == 1) {
            p <- p + theme(axis.text.y = axis.text.y,
                           plot.margin = unit(c(plot.margin.top, 0,
                                                plot.margin.bottom,
                                                plot.margin.left), "cm"))
        } else if(i == 6) {
            p <- p + theme(axis.text.y = element_text(size = font.size.med - 4,
                                                      family = "Arial", colour = "#FFFFFF00"),
                           plot.margin = unit(c(plot.margin.top,
                                                plot.margin.right,
                                                plot.margin.bottom,
                                                0), "cm"))
        } else {
            p <- p + theme(axis.text.y = element_text(size = font.size.med - 4,
                                                      family = "Arial", colour = "#FFFFFF00"),
                           plot.margin = unit(c(plot.margin.top, 0,
                                                plot.margin.bottom, 0), "cm"))
        }
        plots[[i]] <- p
    }

    figure <- ggarrange(plots[[1]], plots[[2]], plots[[3]],
                        plots[[4]], plots[[5]], plots[[6]],
                        nrow = 1, ncol = 6, heights = c(rep(1, 6)),
                        align = "h")

    if (draw.y.axis.title) {
        figure <- annotate_figure(figure,
                                  top = text_grob(title,
                                                  family = "Arial",
                                                  color = "black",
                                                  face = "bold",
                                                  size = font.size.med),
                                  left = text_grob(y.axis.title,
                                                   vjust = 0.5,
                                                   size = font.size.med - 2,
                                                   color = "black", rot = 90))
    }

    PrintLog("** Finished plotting trinucleotide spectrum")
    return(figure)
}


#' @title PlotSignaturesContProbs
#' @description
#' Plots a horizontal barplot of identified mutational signatures
#'
#' @param df.identified.mut.sigs A data.frame of identified mutational signatures
#' @param df.ref.sigs.groups.colors A data.frame with 'signature', 'group', and 'color' columns
#' @param title Plot title
#' @param convert.to.percentage If true, convert y values to percentage (x 100); Default = T,
#' @param font.size.small Small font size; Default = 8,
#' @param font.size.med Medium font size; Default = 14,
#' @param plot.margin Margin vector for drawing plot; Default = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @examples
#' \dontrun{
#'  g <- PlotSignaturesContProbs(sigs = c(mutalisk.results$identified.mut.sigs),
#'  sigs.probs = c(mutalisk.results$identified.mut.sigs.probs),
#'  df.ref.sigs.groups.colors = GetPCAWGMutSigsEtiologiesColors())
#'  print(g)
#' }
#' @export
PlotSignaturesContProbs <- function(df.identified.mut.sigs,
                                    df.ref.sigs.groups.colors,
                                    title,
                                    convert.to.percentage = T,
                                    font.size.small = 8,
                                    font.size.med = 14,
                                    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) {
    PrintLog("** Started plotting signatures probabilities")

    # Prepare plot data
    x <- as.character(df.identified.mut.sigs$signature)
    y <- as.numeric(df.identified.mut.sigs$weight)
    groups <-c ()
    groups.colors <- c()
    for (curr.sig in x) {
        matched <- df.ref.sigs.groups.colors[df.ref.sigs.groups.colors$signature == curr.sig,]
        groups <- c(groups, as.character(matched$group))
        groups.colors <- c(groups.colors, as.character(matched$color))
    }

    if (convert.to.percentage) {
        y <- y * 100
        y <- round(y, 1)
        y.lim <- 100
        y.labels <- paste0(y, "%")
    } else {
        y <- round(y, 2)
        y.lim <- 1
        y.labels <- y
    }

    df.plot.group <- data.frame("x" = x,
                                "y" = y,
                                "label" = y.labels,
                                "group" = groups,
                                "group.color" = groups.colors,
                                stringsAsFactors = F)

    names(groups.colors) <- groups

    # Plot
    g <- ggplot(data = df.plot.group, aes_string(x = 'x', y = 'y', fill = 'group')) +
        geom_bar(stat="identity", width = 0.5,
                 aes_string(fill = 'group')) +
        ggtitle(title) + ylab("Weight (%)") +
        geom_text(aes_string(label = 'label'),
                  vjust = 0.5,
                  hjust = eval(parse(text = "ifelse((y / y.lim) < 0.5, -0.3, 1.3)")),
                  color = "black",
                  size = font.size.med / 3) +
        scale_y_continuous(limits = c(0, y.lim),
                           expand = c(0,0)) +
        theme(plot.margin = plot.margin,
              plot.title = element_text(size = font.size.med,
                                        hjust = 0.5),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = font.size.med),
              axis.text.y = element_text(size = font.size.med,
                                         hjust = 0,
                                         margin = unit(c(0, 5, 0, 0), "mm")),
              panel.grid.major.x = element_line(colour = "#DCDCDC", size = 0.5),
              panel.grid.major.y = element_blank(),
              panel.background = element_blank(),
              legend.text = element_text(size = font.size.med),
              legend.title = element_text(size = font.size.med),
              legend.key.size = unit(1, 'lines'),
              legend.key.width = unit(1, 'lines'),
              legend.key.height = unit(1, 'lines'),
              legend.position = "bottom") +
        coord_flip() + xlab("") + ylab("") +
        scale_fill_manual(values = groups.colors) +
        guides(fill = guide_legend(title = "Aetiologies",
                                   title.theme = element_text(size = font.size.med,
                                                              hjust = 0.5),
                                   ncol = 2,
                                   nrow = ceiling(length(unique(groups)) / 2),
                                   title.position = "top",
                                   title.hjust = 0))

    PrintLog("** Finished plotting signatures probabilities")
    return(g)
}


#' @title PlotMutationTypes
#' @description
#' Plots a horizontal barplot of mutation types
#'
#' @param mutation.types Mutation types; Default = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
#' @param mutation.types.values Mutation count for each mutation type
#' @param mutation.types.colors A color vector for indicating mutation types
#' @param max.y.val y axis maximum value
#' @param title Plot title
#' @param convert.to.percentage if True convert y values to percentage (x 100); Default = T
#' @param show.legend If True, show legend; Default = T
#' @param font.size.small Small font size; Default = 8
#' @param font.size.med Medium font size; Default = 14
#' @param plot.margin Margin vector for drawing plot; Default = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' p <- PlotMutationTypes(mutation.types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
#'                   mutation.types.values = c(0.3, 0.3, 0.1, 0.1, 0.1, 0.1),
#'                   mutation.types.colors = TriNuc.Mutation.Type.Hex.Colors,
#'                   max.y.val = 0.5,
#'                   convert.to.percentage = T,
#'                   show.legend = T,
#'                   font.size.small = 8,
#'                   font.size.med = 14,
#'                   plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
#' print(p)
#' }
#' @import ggplot2
#' @export
PlotMutationTypes <- function(mutation.types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              mutation.types.values,
                              mutation.types.colors,
                              max.y.val,
                              title,
                              convert.to.percentage = T,
                              show.legend = T,
                              font.size.small = 8,
                              font.size.med = 14,
                              plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) {
    PrintLog("** Started plotting mutation types")

    # Prepare plot data
    x <- mutation.types
    y <- mutation.types.values
    if (convert.to.percentage) {
        y <- y * 100
        y <- round(y, 1)
        y.labels <- paste0(y, "%")
        max.y.val <- max.y.val * 100
    } else {
        y <- round(y, 2)
        y.labels <- y
        y.lim <- 1
    }
    names(mutation.types.colors) <- x

    df.plot <- data.frame(list(x = x,
                               y = y,
                               label = y.labels),
                          stringsAsFactors = F)

    # Plot
    g <- ggplot(data = df.plot, aes(x = x, y = y, fill = x)) +
        geom_bar(stat="identity", width = 0.5) +
        ggtitle(title) +
        xlab("") + ylab("Contribution (%)") +
        scale_y_continuous(limits = c(0, max.y.val),
                           expand = c(0,0)) +
        theme(axis.ticks = element_blank(),
              plot.title = element_text(size = font.size.med,
                                        hjust = 0.5),
              axis.text.x = element_text(size = font.size.small),
              axis.text.y = element_text(size = font.size.med),
              axis.title.y = element_text(size = font.size.med),
              panel.grid.major.y = element_line(colour = "#DCDCDC", size = 0.5),
              panel.grid.major.x = element_blank(),
              panel.background = element_blank(),
              plot.margin = plot.margin) +
        scale_fill_manual(values = mutation.types.colors)

    if (show.legend) {
        g <- g + guides(fill = guide_legend(title = "Mutation Type")) +
            theme(legend.title = element_text(size = font.size.med),
                  legend.text = element_text(size = font.size.med))
    } else {
        g <- g + guides(fill = FALSE)
    }

    PrintLog("** Finished plotting signatures probabilities")
    return(g)
}


#' @title PlotMutaliskResults
#' @description
#' Plots Mutalisk results
#'
#' @param mutalisk.results A list obtained from \code{\link{RunMutalisk}}
#' @param signatures A character vector of mutational signatures names
#' @param df.ref.sigs.groups.colors A data.frame with signature groups and colors
#' @param trinuc.max.y A numeric value (maximum y-axis value)
#' @param trinuc.min.y A numeric value (minimum y-axis value)
#' @param mut.type.max.y A numeric value
#' @param title A string value
#' @param font.size.small A numeric value
#' @param font.size.med A numeric value
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#'   df.ref.mut.sigs <- GetPCAWGMutSigs()
#'   target.mut.sigs <- GetPCAWGMutSigsNames()
#'   vcf.obj <- ReadVCF(vcf.file = "../data/sample/P-233-CT.final.vcf")
#'   mutalisk.results <- RunMutalisk(vcf.obj = vcf.obj,
#'                                   df.ref.mut.sigs = df.ref.mut.sigs,
#'                                   target.mut.sigs = target.mut.sigs)
#'   p <- PlotMutaliskResults(mutalisk.results = mutalisk.results)
#'   print(p)
#' }
#' @import ggplot2
#' @export
PlotMutaliskResults <- function(mutalisk.results,
                                signatures,
                                df.ref.sigs.groups.colors,
                                trinuc.max.y,
                                trinuc.min.y,
                                mut.type.max.y,
                                title,
                                font.size.small = 8,
                                font.size.med = 14) {
    PrintLog("** Started plotting Mutalisk results")

    # 1. Plot contribution probabilities of identified signatures
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
    f1 <- PlotSignaturesContProbs(df.identified.mut.sigs,
                                  title = title,
                                  df.ref.sigs.groups.colors = df.ref.sigs.groups.colors,
                                  font.size.small = font.size.small,
                                  font.size.med = font.size.med)

    # 2. Plot trinucleotide spectrums (96 substitution subtypes)
    g1 <- PlotTriNucSpectrum(sub.types = mutalisk.results$sub.types,
                             spectrum = mutalisk.results$sub.types.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Obs Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = title,
                             font.size.small = font.size.small,
                             font.size.med = font.size.med)
    g2 <- PlotTriNucSpectrum(sub.types = mutalisk.results$sub.types,
                             spectrum = mutalisk.results$identified.mut.sigs.spectrum,
                             max.y.val = trinuc.max.y,
                             min.y.val = 0,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "MLE Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = title,
                             font.size.small = font.size.small,
                             font.size.med = font.size.med)
    g3 <- PlotTriNucSpectrum(sub.types = mutalisk.results$sub.types,
                             spectrum = mutalisk.results$residuals,
                             max.y.val = trinuc.max.y,
                             min.y.val = trinuc.min.y,
                             draw.top.strip = T,
                             draw.x.axis.labels = F,
                             y.axis.title = "Res Frac",
                             plot.margin.top = 0.05,
                             plot.margin.right = 0.5,
                             plot.margin.bottom = 0.05,
                             plot.margin.left = 0.5,
                             title = title,
                             font.size.small = font.size.small,
                             font.size.med = font.size.med)

    # 3. Plot mutation types (6 substitution types)
    f3 <- PlotMutationTypes(mutation.types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                            mutation.types.values = EnumerateTriNucCounts(
                                mutalisk.results$identified.mut.sigs.spectrum),
                            mutation.types.colors = TriNuc.Mutation.Type.Hex.Colors,
                            max.y.val = mut.type.max.y,
                            convert.to.percentage = T,
                            show.legend = F,
                            title = title,
                            font.size.small = font.size.small,
                            font.size.med = font.size.med)

    PrintLog("** Finished plotting Mutalisk results")
    return(list(f1 = f1,
                f2.1 = g1,
                f2.2 = g2,
                f2.3 = g3,
                f3 = f3))
}


#' @title PlotTable
#' @description
#' Plots basic statistics table
#'
#' @param df = A data.frame where the first column is header and the second column is data value
#' @param padding Padding size; Default = 20
#' @param font.size Font size; Default = 14
#'
#' @return A plot
#'
#' @importFrom gtable gtable_add_grob
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom grid segmentsGrob grid.newpage gpar
#' @export
PlotTable <- function(df,
                      padding = 20,
                      font.size = 14) {
    PrintLog("** Started plotting table")

    # Add padding
    padding.char <- rep(" ", padding)
    padding.char <- paste0(padding.char, collapse = "")
    df[,2] <- paste0(padding.char, df[,2])

    padding.char <- rep(" ", nchar(colnames(df)[1]))
    padding.char <- paste0(padding.char, collapse = "")
    colnames(df)[1] <- paste0(colnames(df)[1], padding.char)

    # Render tableGrob
    g <- tableGrob(df[1:nrow(df),],
                   rows = NULL,
                   theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0)),
                                          colhead = list(fg_params=list(hjust = 0, x = 0))))

    # Align second column to right
    id <- which(grepl("core-fg", g$layout$name) & g$layout$l == 2)
    for (i in id) {
        g$grobs[[i]]$x <- unit(1, "npc")
        g$grobs[[i]]$hjust <- 1
    }

    # Set font size
    lapply(1:length(g$grobs),function(i){g$grobs[[i]]$gp$fontsize <<- font.size})

    # Add a horizontal line in each row
    for (i in 1:(nrow(df))) {
        lwd <- 0.5
        if (i == 1) {
            lwd <- 1
        }
        g <- gtable_add_grob(g,
                             grobs = segmentsGrob( # line across the bottom
                                 x0 = unit(0,"npc"),
                                 y0 = unit(0,"npc"),
                                 x1 = unit(1,"npc"),
                                 y1 = unit(0,"npc"),
                                 gp = eval(parse(text = "gpar(lwd = lwd)"))),
                             t = i, b = 1, l = 1, r = 2)
    }
    grid.newpage()

    PrintLog("** Finished plotting table")
    return(g)
}
