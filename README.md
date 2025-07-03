# RNA-seq Analyses of *A Selection-Based Knockout Method for a Choanoflagellate Reveals Regulation of Multicellular Development by Hippo Signaling*

This repository contains the code used to generate data for the paper:  
*A selection-based knockout method for a choanoflagellate reveals regulation of multicellular development by Hippo signaling.*

## Data Availability

**- Raw sequencing data (fastq and BAM files) are available in SRA (coming soon).**
**- The files required to run the scripts in this repository, intermediate results, Sequana and SARTools reports, and final plots are publicly available on [Figshare](https://figshare.com/articles/dataset/A_selection-based_knockout_method_for_a_choanoflagellate_reveals_regulation_of_multicellular_development_by_Hippo_signaling/29401835). Input and output files for the analysis of pacI mutants (insertion) and pacΔ mutants (deletion) are located in the directories `KO_Insertion_DEG_analysis` and `KO_Deletion_DEG_analysis`, respectively.**

---

## Installation and Dependencies

The RNA-seq analyses were performed on a SLURM cluster using the [Sequana RNA-seq pipeline](https://github.com/sequana/sequana_rnaseq).  
Please refer to the official tutorial for installation and usage: https://sequana.readthedocs.io/en/main/tutorial.html#rna-seq

Differential expression analyses were performed using [SARTools](https://github.com/PF2-pasteur-fr/SARTools/), which wraps around DESeq2 in R.

---

## Pipeline Overview

### 1️/ Sequana RNA-seq Analysis

Steps included:
- Adapter trimming: `Cutadapt 2.7`
- Alignment to *Salpingoeca rosetta* genome (GCA_000188695.1): `STAR 2.7.3a`
- Read counting: `FeatureCounts 1.6.4` (strand-specific, NCBI annotation v92)
- Quality control: `MultiQC 1.6`

> Based on QC, one WT, one `hippopacI1`, and one `yorkiepacI` sample were excluded from downstream analysis of the insertion dataset.

**Input:**
- FASTQ files must be downloaded from SRA

**Output:**
- BAM files can be founded in SRA
- Feature_count files: `input_for_R_analysis/all_feature_counts/`
- Sequana reports: `output_sequana_summary`

---

### 2️/ DEG Analysis Using SARTools (DESeq2)

Performed separately for insertion and deletion experiments using DESeq2 (v1.42.1) with Benjamini-Hochberg correction (FDR < 0.05). 
Original R scripts templates were modified to include Likelihood ratio tests (LRT) and modified functions to adapt generated Volcano plots and PCA plots for esthetic purpose, notably removal of batch effect via limma (version 3.58.1) (Ritchie et al., 2015). 
Significant covariates: Only RNA extraction date (included as batch effect).

Modifications from original SARTools templates:
- Likelihood Ratio Tests (LRT)
- Custom PCA (with Batch correction using `limma` (v3.58.1)) and Volcano plots

**Scripts:**
- Insertion: `Ins_RNAdiff_SARTools.R`
- Deletion: `Del_RNAdiff_SARTools.R`

**Input:**
- `input_for_R_analysis/`

**Output:**
- `output_for_R_analysis/`

---

### 3️/ Read Mapping

Reads were visualized using [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html).  
Alternative start codons identified using Geneious Prime® 2024.0.4 based on relaxed Kozak consensus (NNNRNNATGN).

**Script:**
- `Reads_map.R`

**Input:**
- BAM files must be downloaded from SRA ``

**Output:**
- Plots in `output_for_R_analysis/Reads/`

---

### 4️/ Gene Filtering

To remove confounding effects of puromycin resistance and other unindentified factors:
- Genes differentially expressed between WT and `epi-pac` or inconsistently identified between WT and edited strains were filtered.

To get high confident DEGs for Warts-KO
- “Warts-KO consistent DEGs” = DEGs shared between Warts-KO insertion and deletion conditions.

**Scripts:**
- Insertion: `Ins_Filter_DEG_Heatmap.R`
- Deletion: `Del_Filter_DEG_Heatmap.R`

**Input:**
- Pre-filter DEG tables: `output_for_R_analysis/tables/`

**Output:**
- Filtered DEG lists and heatmap plots: `output_for_R_analysis/Heatmap/filtered_deg_lists/`

---

### 5️/ Heatmaps

Gene annotations were retrieved via UniProt API, except:
- `PTSG_07358` (named "couscous" from previous study) and `PTSG_01648` (unannotated) but with structural annotation via [AlphaFold](https://alphafold.ebi.ac.uk) and Foldseek using AFDB50 and PDB databases

**Scripts:**
- Same as in the filtering step

**Input:**
- Pre-/Post-filter DEG tables:  
  `output_for_R_analysis/tables/` and  
  `output_for_R_analysis/Heatmap/filtered_deg_lists/`

**Output:**
- Heatmaps in `output_for_R_analysis/Heatmap/`

---

### 6️/ Volcano Plots

**Scripts:**
- `Del_VolcanoPlot.R`: significant and filtered DEGs
- `Del_VolcanoPlot_YORKIE.R`: highlights Yorkie DEGs overlapping with Warts-KO DEGs

**Input:**
- DEG tables from:  
  `output_for_R_analysis/tables/` and  
  `output_for_R_analysis/Heatmap/filtered_deg_lists/`

**Output:**
- Plots in `output_for_R_analysis/VolcanoPlot/`


### 7/ Log2FoldChange Plot

**Scripts:**
- `Del_VolcanoPlot_YORKIE.R`: Compare log2FoldChange of Yorkie-KO genes which are also consistent Warts-KO DEGs with all other Yorkie-KO genes

**Input:**
- DEG tables from:  
  `output_for_R_analysis/tables/` and  
  `output_for_R_analysis/Heatmap/filtered_deg_lists/`

**Output:**
- Plot: `output_for_R_analysis/Dotplot_YorkieKO-all_vs_YorkieKO-Warts_log2FC.svg`
---

## Citation

If you use this repository, please cite our paper:  
*A selection-based knockout method for a choanoflagellate reveals regulation of multicellular development by Hippo signaling*  
_(full citation coming soon)_

---

## Contact

For questions, please contact:  
**Mylan Ansel**  
[mylan.ansel@pasteur.fr]

---


## License
This work is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
