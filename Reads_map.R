#-------------------------------------------------------------------------------#
# Mylan Ansel, Institut Pasteur, 2025
# This use Gviz to vizualize reads mapping from RNAseq
# 
#-------------- Install and load libraries -------------------------------------#

#BiocManager::install(c("Gviz", "GenomicAlignments", "Rsamtools"))
library(Gviz)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)

#-------------- Load GFF and Reads -------------------------------------------------------#
gff <- rtracklayer::import("W:/mansel/RNAseq/For_publication/Salpingoeca_rosetta_gca_gff3.gff3", format = "gff3")
txdb <- makeTxDbFromGFF("W:/mansel/RNAseq/For_publication/Salpingoeca_rosetta_gca_gff3.gff3", format = "gff3")

# RNAseq of KO per insertion
bam_files_WT_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Sr_1_S6.sorted.bam", 
                      "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Sr_2_S20.sorted.bam", 
                      "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Sr_4_S47.sorted.bam")

bam_files_Warts_KO_1_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/W1C3_1_S10.sorted.bam", 
                              "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/W1C3_2_S24.sorted.bam", 
                              "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/W1C3_3_S70.sorted.bam")

bam_files_Warts_KO_2_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/W1F9_2_S25.sorted.bam", 
                              "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/W1F9_3_S71.sorted.bam")

bam_files_Yorkie_KO_1_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Y2B4_1_S2.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Y2B4_2_S16.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Y2B4_3_S30.sorted.bam")

bam_files_Yorkie_KO_2_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Y2E2_1_S1.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Y2E2_2_S15.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/Y2E2_3_S29.sorted.bam")

bam_files_Hippo_KO_1_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/H1A12_1_S7.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/H1A12_2_S21.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/H1A12_3_S67.sorted.bam")

bam_files_Hippo_KO_2_ins <- c("W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/H1C5_1_S13.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/H1C5_2_S27.sorted.bam", 
                               "W:/mansel/RNAseq/For_publication/KO_Insertion/input_for_R_analysis/BAM/H1C5_3_S33.sorted.bam")

# RNAseq of KO per deletion
bam_files_WT_del <- c("W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/Sr_1_S1.sorted.bam", 
                  "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/Sr_2_S2.sorted.bam", 
                  "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/Sr_3_S3.sorted.bam", 
                  "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/Sr_4_S4.sorted.bam")

bam_files_Warts_KO_1_del <- c("W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/WB7_1_S15.sorted.bam", 
                  "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/WB7_2_S16.sorted.bam", 
                  "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/WB7_3_S17.sorted.bam")

bam_files_Warts_KO_2_del <- c("W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/WBA_1_S12.sorted.bam", 
                          "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/WBA_2_S13.sorted.bam", 
                          "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/WBA_3_S14.sorted.bam")

bam_files_Yorkie_KO_2_del <- c("W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/YC7_1_a_S8.sorted.bam", 
                        "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/YC7_2_S10.sorted.bam", 
                        "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/YC7_3_S11.sorted.bam")

bam_files_Yorkie_KO_1_del <- c("W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/YC5_1_S5.sorted.bam", 
                           "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/YC5_2_S6.sorted.bam", 
                           "W:/mansel/RNAseq/For_publication/KO_Deletion/input_for_R_analysis/BAM/YC5_3_S7.sorted.bam")


#-------------- Function to make Tracks -------------------------------------------------------#

getCoverage <- function(bam, region) {
  aln <- readGAlignments(bam, param = ScanBamParam(which = region))
  cov_list <- coverage(aln)
  chr <- as.character(seqnames(region))
  if (!(chr %in% names(cov_list))) {
    warning(paste("No reads in", bam, "on", chr))
    return(Rle(rep(0, width(region))))
  }
  cov_chr <- cov_list[[chr]]
  return(Views(cov_chr, start(region), end(region))[[1]])
}

makeRawTrack <- function(bam, label, color, region, axis_color = "black") {
  aln <- readGAlignments(bam, param = ScanBamParam(which = region))
  cov_list <- coverage(aln)
  chr <- as.character(seqnames(region))
  
  if (!(chr %in% names(cov_list))) {
    warning(paste("No reads in", bam, "on", chr))
    cov <- Rle(rep(0, width(region)))
  } else {
    cov_chr <- cov_list[[chr]]
    cov <- Views(cov_chr, start(region), end(region))[[1]]
  }
  
  raw_cov <- as.numeric(cov)
  max_val <- max(raw_cov, na.rm = TRUE)
  
  pos_vector <- start(region):(end(region))
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos_vector, width = 1))
  
  DataTrack(
    range = gr,
    data = raw_cov,
    genome = "yourGenome",
    name = sprintf(label),
    type = "histogram",
    col.histogram = color,
    fill.histogram = color,
    yaxis = TRUE,
    col.axis="black",
    background.title = "white",
    background.panel = "white",     
    col.frame = "black",
    col.title="black",
    showTitle = TRUE,
    titleWidth = 0.2   
  )
}

#---------------------------------------------------------------------------------#
#-------------- For KO per insertion ----------------------------------------------#
#---------------------------------------------------------------------------------#

#-------------- For Warts -------------------------------------------------------#

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_04961", gff$ID) | grepl("PTSG_04961", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT <- mapply(
  makeRawTrack, bam_files_WT_ins, paste0("WT", 1:3),
  MoreArgs = list(color = "blue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Warts_KO_1_ins, paste0("WartsKO1_", 1:3),
  MoreArgs = list(color = "red", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Warts_KO_2_ins, paste0("WartsKO2_", 1:2),
  MoreArgs = list(color = "red", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)


gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_04961",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)

# Define CRISPR guide coordinates
crispr_guides_df <- data.frame(
  chr = "GL832965",
  start = c(1158337),
  end = c(1158337 + 20 - 1), # 20 bases long
  id = c("CRISPR_guide_1")
)

crispr_guides_gr <- with(crispr_guides_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
crispr_track <- AnnotationTrack(
  range = crispr_guides_gr,
  name = "CRISPR Guides",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "darkgreen",
  col = "darkgreen",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)

# Define first alternative ATG coordinates
ATG_df <- data.frame(
  chr = "GL832965",
  start = c(1158359),
  end = c(1158359 + 3 - 1), # 3 bases long
  id = c("alternative ATG")
)

ATG_gr <- with(ATG_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
ATG_track <- AnnotationTrack(
  range = ATG_gr,
  name = "Alternative ATG",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "brown",
  col = "brown",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)


svg("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Reads/Inser_Warts-PTSG_04961_read_coverage.svg", width = 10, height = 8)
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track,crispr_track,ATG_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  chromosome = chrom
)
dev.off()

## ------------ WartsKo on Yorkie gene -------------------##

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_06057", gff$ID) | grepl("PTSG_06057", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT <- mapply(
  makeRawTrack, bam_files_WT_ins, paste0("WT", 1:3),
  MoreArgs = list(color = "blue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Warts_KO_1_ins, paste0("WartsKO1_", 1:3),
  MoreArgs = list(color = "red", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Warts_KO_2_ins, paste0("WartsKO2_", 1:3),
  MoreArgs = list(color = "red", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)


gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_06057",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)


svg("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Reads/Inser_Warts-PTSG_06057_read_coverage.svg", width = 10, height = 8) 
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  reverseStrand = (strand == "-"),  
  chromosome = chrom
)
dev.off()

#-------------- For Yorkie -------------------------------------------------------#

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_06057", gff$ID) | grepl("PTSG_06057", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT <- mapply(
  makeRawTrack, bam_files_WT_ins, paste0("WT", 1:3),
  MoreArgs = list(color = "blue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Yorkie_KO_1_ins, paste0("YorkieKO1_", 1:3),
  MoreArgs = list(color = "orange", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Yorkie_KO_2_ins, paste0("YorkieKO2_", 1:3),
  MoreArgs = list(color = "orange", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_06057",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)

# Define CRISPR guide coordinates
crispr_guides_df <- data.frame(
  chr = "GL832969",
  start = c(1249173),
  end = c(1249173 + 20 - 1), # 20 bases long
  id = c("CRISPR_guide_1")
)

crispr_guides_gr <- with(crispr_guides_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
crispr_track <- AnnotationTrack(
  range = crispr_guides_gr,
  name = "CRISPR Guides",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "darkgreen",
  col = "darkgreen",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)

# Define first alternative ATG coordinates
ATG_df <- data.frame(
  chr = "GL832969",
  start = c(1247217),
  end = c(1247217 + 3 - 1), # 3 bases long
  id = c("alternative ATG")
)

ATG_gr <- with(ATG_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
ATG_track <- AnnotationTrack(
  range = ATG_gr,
  name = "Alternative ATG",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "brown",
  col = "brown",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)


svg("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Reads/Insert_Yorkie-PTSG_06057_read_coverage.svg", width = 10, height = 8)  # set size in inches
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track,crispr_track,ATG_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  reverseStrand = (strand == "-"),  
  chromosome = chrom
)
dev.off()

#-------------- For Hippo -------------------------------------------------------#

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_10780", gff$ID) | grepl("PTSG_10780", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT <- mapply(
  makeRawTrack, bam_files_WT_ins, paste0("WT", 1:3),
  MoreArgs = list(color = "blue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Hippo_KO_1_ins, paste0("HippoKO1_", 1:3),
  MoreArgs = list(color = "#b61696", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Hippo_KO_2_ins, paste0("HippoKO2_", 1:3),
  MoreArgs = list(color = "#b61696", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_10780",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)

# Define CRISPR guide coordinates
crispr_guides_df <- data.frame(
  chr = "GL832988",
  start = c(550497),
  end = c(550497 + 20 - 1), # 20 bases long
  id = c("CRISPR_guide_1", "CRISPR_guide_2")
)

crispr_guides_gr <- with(crispr_guides_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
crispr_track <- AnnotationTrack(
  range = crispr_guides_gr,
  name = "CRISPR Guides",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "darkgreen",
  col = "darkgreen",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)

# Define first alternative ATG coordinates
ATG_df <- data.frame(
  chr = "GL832988",
  start = c(550164),
  end = c(550164 + 3 - 1), # 3 bases long
  id = c("alternative ATG")
)

ATG_gr <- with(ATG_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
ATG_track <- AnnotationTrack(
  range = ATG_gr,
  name = "Alternative ATG",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "brown",
  col = "brown",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)

svg("W:/mansel/RNAseq/For_publication/KO_Insertion/output_R_analyses/Reads/Insert_Hippo-PTSG_10780_read_coverage.svg", width = 10, height = 8)  # set size in inches
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track,crispr_track,ATG_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  reverseStrand = (strand == "-"),  
  chromosome = chrom
)
dev.off()

#---------------------------------------------------------------------------------#
#-------------- For KO per deletion ----------------------------------------------#
#---------------------------------------------------------------------------------#

#-------------- For Warts -------------------------------------------------------#

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_04961", gff$ID) | grepl("PTSG_04961", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT <- mapply(
  makeRawTrack, bam_files_WT_del, paste0("WT", 1:4),
  MoreArgs = list(color = "darkblue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Warts_KO_1_del, paste0("KO1_", 1:3),
  MoreArgs = list(color = "darkred", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Warts_KO_2_del, paste0("KO2_", 1:3),
  MoreArgs = list(color = "darkred", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)


gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_04961",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)

# Define CRISPR guide coordinates
crispr_guides_df <- data.frame(
  chr = "GL832965",
  start = c(1158337, 1160398),
  end = c(1158337 + 20 - 1, 1160398 + 20 - 1), # 20 bases long
  id = c("CRISPR_guide_1", "CRISPR_guide_2")
)

crispr_guides_gr <- with(crispr_guides_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
crispr_track <- AnnotationTrack(
  range = crispr_guides_gr,
  name = "CRISPR Guides",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "darkgreen",
  col = "darkgreen",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)

svg("W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Reads/Warts-PTSG_04961_read_coverage.svg", width = 10, height = 8)
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track,crispr_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  chromosome = chrom
)
dev.off()

## ------------ WartsKo on Yorkie gene -------------------##

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_06057", gff$ID) | grepl("PTSG_06057", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT_del <- mapply(
  makeRawTrack, bam_files_WT_del, paste0("WT", 1:4),
  MoreArgs = list(color = "darkblue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Warts_KO_1_del, paste0("KO1_", 1:3),
  MoreArgs = list(color = "darkred", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Warts_KO_2_del, paste0("KO2_", 1:3),
  MoreArgs = list(color = "darkred", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)


gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_06057",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)

svg("W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Reads/Warts-PTSG_06057_read_coverage.svg", width = 10, height = 8) 
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  reverseStrand = (strand == "-"),  
  chromosome = chrom
)
dev.off()

#-------------- For Yorkie -------------------------------------------------------#

# Filter for gene feature with ID or Name matching "PTSG_04961"
ptsg_gene <- gff[grepl("PTSG_06057", gff$ID) | grepl("PTSG_06057", gff$Name), ]

# Get basic coordinates
chrom <- as.character(seqnames(ptsg_gene)[1])
start_pos <- min(start(ptsg_gene))
end_pos <- max(end(ptsg_gene))
strand <- as.character(strand(ptsg_gene)[1])

options(ucscChromosomeNames = FALSE)

region <- GRanges(chrom, IRanges(start_pos - 500, end_pos + 500))

tracks_WT <- mapply(
  makeRawTrack, bam_files_WT_del, paste0("WT", 1:4),
  MoreArgs = list(color = "darkblue", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_1 <- mapply(
  makeRawTrack, bam_files_Yorkie_KO_1_del, paste0("KO1_", 1:3),
  MoreArgs = list(color = "darkorange", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

tracks_KO_2 <- mapply(
  makeRawTrack, bam_files_Yorkie_KO_2_del, paste0("KO2_", 1:3),
  MoreArgs = list(color = "darkorange", region = region, axis_color = "black"),
  SIMPLIFY = FALSE
)

gene_track <- GeneRegionTrack(
  txdb,
  genome = "S. rosetta",
  chromosome = chrom,
  start = start_pos - 500,
  end = end_pos + 500,
  name = "PTSG_06057",
  background.title = "white",     # white title background
  background.panel = "white",     # white plotting area
  col.title = "black",            # title text
  col.axis = "black",             # axis ticks
  col = "black",                  # gene model lines
  fill = "black",
  col.line="black",
  arrowHeadWidth=40,
  showTitle = TRUE,
  titleWidth = 0.2 # Adjust as needed
)

# Define CRISPR guide coordinates
crispr_guides_df <- data.frame(
  chr = "GL832969",
  start = c(1249173, 1245385),
  end = c(1249173 + 20 - 1, 1245385 + 20 - 1), # 20 bases long
  id = c("CRISPR_guide_1", "CRISPR_guide_2")
)

crispr_guides_gr <- with(crispr_guides_df, GRanges(chr, IRanges(start, end), id = id))

# Create an AnnotationTrack for CRISPR guides
crispr_track <- AnnotationTrack(
  range = crispr_guides_gr,
  name = "CRISPR Guides",
  genome = "yourGenome", # Use the same genome as your DataTrack (e.g., "S. rosetta")
  fill = "darkgreen",
  col = "darkgreen",
  stacking = "dense",
  showTitle = TRUE,
  titleWidth = 0.2,
  background.title = "white",     # white title background
  background.panel = "white", 
  col.title="black",
  fontcolor.item = "black", # color of the labels on the track
  cex.item = 1,
  showID=TRUE, # show the IDs of the guides
  just.group = "right"
)

svg("W:/mansel/RNAseq/For_publication/KO_Deletion/output_R_analyses/Reads/Yorkie-PTSG_06057_read_coverage.svg", width = 10, height = 8)  # set size in inches
plotTracks(
  c(list(GenomeAxisTrack(  col.title = "black",            # title text
                           col.axis = "black",             # axis ticks
                           col = "black",                  # gene model lines
                           fill = "black",
                           col.line="black",
                           col.frame="black",
                           fontcolor="black"), gene_track,crispr_track), tracks_WT, tracks_KO_1,tracks_KO_2),
  from = start(region),
  to = end(region),
  reverseStrand = (strand == "-"),  
  chromosome = chrom
)
dev.off()


