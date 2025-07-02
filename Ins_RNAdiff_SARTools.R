#-------------------------------------------------------------------------------#
# Mylan Ansel, Institut Pasteur, 2024

# Modified version of template_script_DESeq2.r from SARTools (https://github.com/PF2-pasteur-fr/SARTools/tree/master)
# containing modified functions for aesthetic purpose (Volcano plot + PCA)
# Original functions are also conserved to generate the report.

# It also contains DESeq2 tests for the different confounding factors.

#-------------------------------------------------------------------------------#
#---------------------- Import library for modified functions ------------------#

library(devtools)
#devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")
library(SARTools)
library(gridExtra)
library(grid)
library(ggdendro)
library(scales)
library(viridis)
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(readr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(tibble)
library(ggalt) 
library("RColorBrewer") 
library(apeglm)
library(genefilter)
library(tidyr)
library(stringr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(limma)
#-------------------------------------------------------------------------------#
#------------------------- Modified functions ----------------------------------#
#-------------------------------------------------------------------------------#

# Modify exploreCount function to use new PCA and Volcano plot functions for aesthetic purpose
# and to save svg plots
exploreCounts_1 <- function(object, group, genotype = NULL, condition = NULL,
                            typeTrans = "VST", gene.selection = "pairwise",
                            col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"),
                            shapes = c(16, 17, 18, 15),
                            ggplot_theme = theme_light(),
                            legend_labels = NULL,
                            genotype_labels = NULL,
                            condition_labels = NULL,
                            save_pca_svg = FALSE,
                            pca_svg_filename = "pca_plot.svg",
                            do_batch_correction = FALSE,
                            pca_components = c("PC1", "PC2")) { # New parameter for selecting PCA components
  
  # Input validation for pca_components (delegated to PCAPlot_1, but good to have here too if needed)
  
  if (class(object) == "DESeqDataSet") {
    vsd_or_rld_obj <- NULL 
    
    if (typeTrans == "VST") {
      vsd_or_rld_obj <- DESeq2::varianceStabilizingTransformation(object)
      counts.trans <- assay(vsd_or_rld_obj)
    } else {
      vsd_or_rld_obj <- DESeq2::rlogTransformation(object)
      counts.trans <- assay(vsd_or_rld_obj)
    }
    
    # PCA Plot - Pass all parameters, including new ones for batch correction and PCA components
    pca_plot <- PCAPlot_1(counts.trans = counts.trans,
                          group = group,
                          genotype = genotype,
                          condition = condition,
                          col = col,
                          shapes = shapes,
                          ggplot_theme = ggplot_theme,
                          legend_labels = legend_labels,
                          genotype_labels = genotype_labels,
                          condition_labels = condition_labels,
                          do_batch_correction = do_batch_correction,
                          vsd_obj = vsd_or_rld_obj,
                          pca_components = pca_components) # Pass the new parameter
    print(pca_plot)
    
    # Save PCA plot as SVG if requested
    if (save_pca_svg) {
      # It's good practice to reflect the plotted components in the filename if it's dynamic
      # For now, keeping original filename logic, but consider enhancing it.
      # Example: pca_svg_filename = paste0("PCA_", pca_components[1], "-", pca_components[2], ".svg")
      cowplot::save_plot(filename = pca_svg_filename, plot = pca_plot,
                         base_width = 10, base_height = 8)
      message(paste("PCA plot saved as:", pca_svg_filename))
    }
    
    # Clustering Plot (keeping original functionality)
    clusterPlot(counts.trans = counts.trans, group = group, ggplot_theme = ggplot_theme)
    
  } else if (class(object) == "DGEList") {
    
    # MDS Plot for edgeR DGEList object
    MDSPlot(dge = object, group = group, col = col, gene.selection = gene.selection, ggplot_theme = ggplot_theme)
    # Clustering Plot for edgeR DGEList object
    clusterPlot(counts.trans = edgeR::cpm(object, prior.count = 2, log = TRUE), group = group, ggplot_theme = ggplot_theme)
    
  } else {
    stop("The object is not a DESeqDataSet nor a DGEList")
  }
}

PCAPlot_1 <- function(counts.trans, group, genotype = NULL, condition = NULL,
                      n = min(500, nrow(counts.trans)),
                      col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"),
                      shapes = c(16, 17, 18, 15),
                      ggplot_theme = theme_minimal(),
                      legend_labels = NULL,
                      genotype_labels = NULL,
                      condition_labels = NULL,
                      do_batch_correction = FALSE,
                      vsd_obj = NULL,
                      pca_components = c("PC1", "PC2")) { # New parameter for selecting PCA components
  
  # Input validation for pca_components
  if (!is.character(pca_components) || length(pca_components) != 2 || 
      !all(grepl("^PC[0-9]+$", pca_components))) {
    stop("pca_components must be a character vector of length 2, e.g., c('PC1', 'PC2').")
  }
  
  # Extract the numeric part of the component names (e.g., 1 from "PC1")
  pc_x_idx <- as.numeric(gsub("PC", "", pca_components[1]))
  pc_y_idx <- as.numeric(gsub("PC", "", pca_components[2]))
  
  # Check if batch correction is requested and vsd_obj is provided
  if (do_batch_correction) {
    if (is.null(vsd_obj)) {
      stop("To perform batch correction, 'vsd_obj' must be provided with the 'vsd$batch' information.")
    }
    if (!inherits(vsd_obj, "DESeqTransform") && !inherits(vsd_obj, "RangedSummarizedExperiment")) {
      warning("vsd_obj is not a DESeqTransform or RangedSummarizedExperiment object. Ensure it contains 'assay' and '$batch' for batch correction.")
    }
    
    # Ensure that vsd_obj$batch exists and is correctly defined
    if (!"batch" %in% names(colData(vsd_obj))) {
      stop("The 'vsd_obj' does not contain a 'batch' column in its colData. Batch correction cannot be performed.")
    }
    
    message("Performing batch correction using removeBatchEffect...")
    counts.trans.corrected <- limma::removeBatchEffect(assay(vsd_obj), batch = vsd_obj$batch)
    
    # Replace the original counts.trans with the corrected data for PCA
    counts.trans <- counts.trans.corrected
  }
  
  # PCA on the most variable features
  rv = apply(counts.trans, 1, var, na.rm = TRUE)
  pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), ][1:n,]))
  
  # Calculate proportion of variance for all PCs, and ensure we have enough PCs
  prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
  
  if (max(pc_x_idx, pc_y_idx) > length(prp)) {
    stop(paste0("Cannot plot ", paste(pca_components, collapse = " and "), 
                ". Only ", length(prp), " principal components were found."))
  }
  
  # Extract proportions for the selected components
  prp_x <- round(prp[pc_x_idx], 2)
  prp_y <- round(prp[pc_y_idx], 2)
  
  # Data preparation - use selected PCA components
  d <- data.frame(x = pca$x[, pc_x_idx], y = pca$x[, pc_y_idx],
                  sample = factor(rownames(pca$x), levels = rownames(pca$x)))
  
  # Check if separate genotype and condition are provided
  if (!is.null(genotype) && !is.null(condition)) {
    # Check lengths match
    if (length(genotype) != nrow(d)) {
      stop(paste("Length of genotype vector (", length(genotype), 
                 ") does not match number of samples (", nrow(d), ")", sep = ""))
    }
    if (length(condition) != nrow(d)) {
      stop(paste("Length of condition vector (", length(condition), 
                 ") does not match number of samples (", nrow(d), ")", sep = ""))
    }
    
    # Use separate genotype and condition variables
    d$genotype <- genotype
    d$condition <- condition
    d$group <- group # Keep original group for backward compatibility
    
    # PCA scatter plot with color for genotype and shape for condition
    pca_plot <- ggplot(data = d, aes(x = x, y = y, color = genotype, shape = condition)) +
      geom_point(show.legend = TRUE, size = 3) +
      scale_colour_manual(values = col, 
                          labels = if (!is.null(genotype_labels)) genotype_labels else waiver(),
                          name = "Genotype") +
      scale_shape_manual(values = shapes,
                         labels = if (!is.null(condition_labels)) condition_labels else waiver(),
                         name = "Condition") +
      xlab(paste0(pca_components[1], " (", prp_x, "%)")) + # Use selected component name
      ylab(paste0(pca_components[2], " (", prp_y, "%)")) + # Use selected component name
      theme_void() +
      theme(panel.border = element_rect(fill = NA, color = "black"),
            axis.title = element_text(size = 12))
    
    # For density plots, use genotype for coloring (consistent with main plot)
    d$group_combined <- d$genotype
    
  } else {
    # Use original group variable
    d$group <- group
    
    # Original PCA scatter plot
    pca_plot <- ggplot(data = d, aes(x = x, y = y, color = group)) +
      geom_point(show.legend = TRUE, size = 3) +
      labs(color = "") +
      scale_colour_manual(values = col, labels = if (!is.null(legend_labels)) legend_labels else waiver()) +
      xlab(paste0(pca_components[1], " (", prp_x, "%)")) + # Use selected component name
      ylab(paste0(pca_components[2], " (", prp_y, "%)")) + # Use selected component name
      theme_void() +
      theme(panel.border = element_rect(fill = NA, color = "black"),
            axis.title = element_text(size = 12))
    
    d$group_combined <- d$group
  }
  
  # Density plots - use only genotype colors to match the number of colors provided
  density_x <- ggplot(d, aes(x = x, fill = group_combined)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = col[1:length(unique(d$group_combined))]) +
    theme_void() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.margin = margin(0, 5.5, 0, 5.5)) +
    xlab(NULL)
  
  density_y <- ggplot(d, aes(x = y, fill = group_combined)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = col[1:length(unique(d$group_combined))]) +
    coord_flip() +
    theme_void() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.margin = margin(5.5, 0, 5.5, 0)) +
    ylab(NULL)
  
  # Combine plots
  combined <- cowplot::insert_xaxis_grob(
    cowplot::insert_yaxis_grob(pca_plot, density_y, grid::unit(0.2, "null"), position = "right"),
    density_x, grid::unit(0.2, "null"), position = "top"
  )
  
  # Draw the plot
  cowplot::ggdraw(combined)
}

volcanoPlot2 <- function(complete, alpha=0.05, outfile=FALSE, padjlim=NULL, ggplot_theme=theme_gray()) {
  ncol <- min(2, length(complete))
  nrow <- ceiling(length(complete) / ncol)
  
  for (name in names(complete)) {
    complete.name <- complete[[name]]
    complete.name$padj[which(complete.name$padj == 0)] <- .Machine$double.xmin
    complete.name <- complete.name[which(!is.na(complete.name$padj)), ]
    complete.name$DE <- factor(ifelse(complete.name$padj <= alpha, "yes", "no"), levels=c("no", "yes"))
    
    reverselog_trans <- function(base = exp(1)) {
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      trans_new(paste0("reverselog-", format(base)), trans, inv,
                log_breaks(base=base),
                domain=c(.Machine$double.xmin, Inf))
    }
    
    plot <- ggplot(data=complete.name, 
                   aes(x=.data$log2FoldChange, y=.data$padj, color=.data$DE)) +
      geom_point(show.legend=FALSE, alpha=0.5) +
      geom_hline(yintercept=alpha, color="black", linetype="dashed", linewidth=0.5) +  # Horizontal line at alpha threshold
      geom_vline(xintercept=0, color="black", linetype="dashed", linewidth=0.5) +  # Horizontal line at alpha threshold
      scale_y_continuous(trans=reverselog_trans(10),
                         breaks=trans_breaks("log10", function(x) 10^x),
                         labels=trans_format("log10", math_format(~10^.x))) +
      scale_colour_manual(values=c("no"="#999999", "yes"="darkblue"), drop=FALSE) +
      xlab(expression(log[2]~fold~change)) +
      ylab("Adjusted P-value") +
      ggtitle(paste0("Volcano plot - ", gsub("_", " ", name))) +
      ggplot_theme +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.title = element_blank(),
            axis.text.x = element_text(size = 10, angle = 0, vjust = 0, hjust = 0.4, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title =  element_text(size = 13, color = "black"),
            axis.line = element_line(color="black"),
            legend.position = "none") 
    
    # Print the plots directly to the current graphical device
    print(plot)
    #ggsave(paste0("Volcano_", name, ".png"), plot = plot)
    ggsave(filename = paste0("Volcano_", name, ".svg"), plot = plot, width = 4, height = 6)
  }
}

summarizeResults.DESeq2_1 <- function(out.DESeq2, group, independentFiltering=TRUE, cooksCutoff=TRUE,
                                    alpha=0.05, col=c("lightblue","orange","MediumVioletRed","SpringGreen"),
                                    log2FClim=NULL, padjlim=NULL, ggplot_theme=theme_light()) {
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # Diagnostic of the size factors
  diagSizeFactorsPlots(dds=dds, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # Boxplots before and after normalisation
  countsBoxplots(dds, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # Dispersions plot
  dispersionsPlot(dds=dds, ggplot_theme=ggplot_theme)
  
  # Results of the independent filtering
  if (independentFiltering) {
    tabIndepFiltering <- tabIndepFiltering(results)
    cat("Number of features discarded by the independent filtering:\n")
    print(tabIndepFiltering, quote=FALSE)
  } else {
    tabIndepFiltering <- NULL
  }
  
  # Exporting results of the differential analysis
  complete <- exportResults.DESeq2_1(out.DESeq2, group=group, alpha=alpha)
  
  # Small table with number of differentially expressed features
  nDiffTotal <- nDiffTotal(complete=complete, alpha=alpha)
  cat("\nNumber of features down/up and total:\n")
  print(nDiffTotal, quote=FALSE)
  
  # Histograms of raw p-values
  rawpHist(complete=complete, ggplot_theme=ggplot_theme)
  
  # MA-plots
  MAPlot(complete=complete, alpha=alpha, log2FClim=log2FClim, ggplot_theme=ggplot_theme)
  
  # Volcano plots (now directly displayed)
  volcanoPlot2(complete=complete, alpha=alpha, padjlim=padjlim, ggplot_theme=ggplot_theme)
  
  return(list(complete=complete, tabIndepFiltering=tabIndepFiltering, nDiffTotal=nDiffTotal))
}

# Change export to FALSE to not save the files
exportResults.DESeq2_1 <- function(out.DESeq2, group, alpha=0.05, export=FALSE){
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # comptages bruts et normalis?s
  counts <- data.frame(Id=rownames(counts(dds)), counts(dds), round(counts(dds, normalized=TRUE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
  # baseMean avec identifiant
  bm <- data.frame(Id=rownames(results[[1]]),baseMean=round(results[[1]][,"baseMean"],2))
  # merge des info, comptages et baseMean selon l'Id
  base <- merge(counts, bm, by="Id", all=TRUE)
  tmp <- base[,paste("norm", colnames(counts(dds)), sep=".")]
  for (cond in levels(group)){
    base[,cond] <- round(apply(as.data.frame(tmp[,group==cond]),1,mean),0)
  }
  
  complete <- list()
  for (name in names(results)){
    complete.name <- base
    
    # ajout d'elements depuis results
    res.name <- data.frame(Id=rownames(results[[name]]),
                           FoldChange=round(2^(results[[name]][,"log2FoldChange"]), 3),
                           log2FoldChange=round(results[[name]][,"log2FoldChange"], 3),
                           stat=round(results[[name]][,"stat"], 3),
                           pvalue=results[[name]][,"pvalue"],
                           padj=results[[name]][,"padj"])
    complete.name <- merge(complete.name, res.name, by="Id", all=TRUE)
    # ajout d'elements depuis mcols(dds)
    mcols.add <- data.frame(Id=rownames(counts(dds)),dispGeneEst=round(mcols(dds)$dispGeneEst,4),
                            dispFit=round(mcols(dds)$dispFit,4),dispMAP=round(mcols(dds)$dispMAP,4),
                            dispersion=round(mcols(dds)$dispersion,4),betaConv=mcols(dds)$betaConv,
                            maxCooks=round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by="Id", all=TRUE)
    complete[[name]] <- complete.name
    
    if (export){
      # s?lection des up et down
      up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0),]
      up.name <- up.name[order(up.name$padj),]
      down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange<=0),]
      down.name <- down.name[order(down.name$padj),]
      
      # exports
      name <- gsub("_","",name)
      write.table(complete.name, file=paste0("tables/",name,".complete.txt"), sep="\t", row.names=FALSE, dec=".", quote=FALSE)
      write.table(up.name, file=paste0("tables/", name,".up.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
      write.table(down.name, file=paste0("tables/", name,".down.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
    }
  }
  
  return(complete)
}


# ------------------------ Original template ----------------------------------#

################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### March 23rd, 2022
### designed to be executed with SARTools 1.8.1
################################################################################


################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
#rm(list=ls())                                     # remove all the objects from the R session

workDir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses"      # working directory for the R session

projectName <- "Chantal2025_KO_insertion"           # name of the project
author <- "Ansel Mylan"                                 # author of the statistical analysis/report

targetFile <- "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/design.txt"       # path to the design/target file
rawDir <- "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/all_features_count"                  # path to the directory containing raw counts files
featuresToRemove <- NULL # NULL if no feature to remove

varInt <- "genotype"                                 # factor of interest
condRef <- "WT"                                      # reference biological condition
batch <- "batch"                                        # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                       # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "shorth"                                 # "median" (default) or "shorth" to estimate the size factors

colors <- c("blue","#86c5da","#86dac5",
            "red","orange","#e632c2")


forceCairoGraph <- FALSE

################################################################################
###                             running script                               ###
################################################################################

#setwd(workDir)
library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots (fail because more than 12 samples)
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# Save dds results for downstream analysis
saveRDS(out.DESeq2$dds, file = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Ins_out_deseq2_dds.rds")

# ------------------------ Test of potential confondant factors ----------------------------------#

# Test if puromycin explain some variability
dds_full <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = target,
                                   design = ~ puro + bacteria + batch + genotype)
dds_full <- DESeq(dds_full, test = "LRT", reduced = ~ bacteria + batch +genotype)
res_batch <- results(dds_full)
summary(res_batch)

# Test if puromycin explain some variability (without including bacteria as well)
dds_full <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = target,
                                   design = ~ puro + batch + genotype)
dds_full <- DESeq(dds_full, test = "LRT", reduced = ~ batch +genotype)
res_batch <- results(dds_full)
summary(res_batch)

# Test if chain or rosette (bacteria) explain some variability
dds_full <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = target,
                                   design = ~ bacteria + puro + batch + genotype)
dds_full <- DESeq(dds_full, test = "LRT", reduced = ~ puro + batch + genotype)
res_batch <- results(dds_full)
summary(res_batch)

# Test if chain or rosette (bacteria) explain some variability (without including puromycin as well)
dds_full <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = target,
                                   design = ~ bacteria + batch + genotype)
dds_full <- DESeq(dds_full, test = "LRT", reduced = ~ batch + genotype)
res_batch <- results(dds_full)
summary(res_batch)

# Test if batch (date of RNA extraction) explain some variability
dds_full <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = target,
                                   design = ~ cassette + bacteria + puro + batch + genotype)
dds_full <- DESeq(dds_full, test = "LRT", reduced = ~ bacteria + puro + genotype)
res_batch <- results(dds_full)
summary(res_batch)

# PCA and Heatmap to see the contribution of clone, bacteria and puromycin (puro)
vsd <- vst(out.DESeq2$dds)
plotPCA(vsd, intgroup=c("genotype", "clone"))
plotPCA(vsd, intgroup=c("puro"))
plotPCA(vsd, intgroup=c("bacteria"))
plotPCA(vsd, intgroup=c("batch"))

pcaData <- plotPCA(vsd, intgroup = c("genotype", "batch"), returnData = TRUE)

# 2. Calculate percent variance for axis labels
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 3. Custom ggplot: color by genotype, shape by batch
ggplot(pcaData, aes(x = PC1, y = PC2, color = batch, shape = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# ------------------------ Downstream analysis ----------------------------------#
# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt],
                                          col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, 
                   workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)

# ------------------------ PCA and volcano plots with modified scripts ----------------------------------#

# PCA + clustering (use modified version for asthetic purpose)
# Get sample names from your DESeq2 object
sample_names <- colnames(out.DESeq2$dds)

# Remove everything after the second underscore to get the base name
base_names <- sub("_[^_]*$", "", sample_names)  # Remove last _xxx part

# Determine condition based on presence of "a" or "algo" after first underscore
condition_vector <- ifelse(grepl("_a_|_algo_", base_names), 
                           "A.machipongonensis", "E.pacifica")

# Create genotype vector by extracting the prefix before first underscore
genotype_prefixes <- sub("_.*", "", base_names)  # Get everything before first underscore

# Apply genotype transformations
genotype_vector <- genotype_prefixes
genotype_vector <- ifelse(genotype_vector == "Sr", "WT", genotype_vector)
genotype_vector <- ifelse(startsWith(genotype_vector, "H"), "Hippo-KO", genotype_vector)
genotype_vector <- ifelse(startsWith(genotype_vector, "Y"), "Yorkie-KO", genotype_vector)
genotype_vector <- ifelse(startsWith(genotype_vector, "W1"), "Warts-KO", genotype_vector)
genotype_vector <- ifelse(genotype_vector == "Nssp", "Neg2", genotype_vector)
genotype_vector <- ifelse(genotype_vector == "Nap", "Neg2 puro", genotype_vector)

desired_genotype_order <- c("WT", "Neg2","Neg2 puro", "Warts-KO","Yorkie-KO","Hippo-KO")

# Convert genotype_vector to a factor with the desired order
genotype_vector <- factor(genotype_vector, levels = desired_genotype_order)

# Define your desired display order for condition labels
# Make sure these labels exactly match the values in your condition_vector
desired_condition_order <- c("A.machipongonensis", "E.pacifica")

# Convert condition_vector to a factor with the desired order
condition_vector <- factor(condition_vector, levels = desired_condition_order)

summary_table <- data.frame(
  sample = sample_names,
  base_name = base_names,
  genotype = genotype_vector,
  condition = condition_vector,
  stringsAsFactors = FALSE
)

# Call the function with the new parameters
exploreCounts_1(object = out.DESeq2$dds,
                group = target[, varInt],  # Keep for backward compatibility
                genotype = genotype_vector,
                condition = condition_vector,
                typeTrans = typeTrans,
                col = colors,
                shapes = c(16, 17),  # Shapes for the 2 conditions (circle, triangle)
                genotype_labels = desired_genotype_order,
                condition_labels = desired_condition_order,
                save_pca_svg = TRUE,
                pca_svg_filename = "W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/PCA1-2.svg",
                do_batch_correction = TRUE,
                pca_components = c("PC1", "PC2"))

# modified version to just plot adapted volcano plot
summaryResults_1 <- summarizeResults.DESeq2_1(out.DESeq2, group=target[,varInt],
                                          col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)


# save image of the R session
#save.image(file=paste0(projectName, ".RData"))
