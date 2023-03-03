# Copertino et al.: The role of HODHbt for HIV CTL response in kick-and-kill approaches

The repo here contains the code that was used for the data processing of human PBMC single-cell RNA-seq data.

## 5' scRNA-seq Library Preparation and Sequencing

Following in vitro treatment with IL-15 + DMSO (control), IL-15 + HODHBt (both from R&D Systems, or DMSO,
PBMC from 3 ARV-suppressed PWH following in vitro treatment with were resuspended at a density of 1000 cells/ul
in PBS plus 0.04% bovine serum albumin on ice and loaded into the 10x Genomics Chromium Controller (10x Genomics, USA) 
with a target capture of ca. 5,000 cells per condition/donor using the single cell immune profiling 5' chip and reagent/gel
bead kits according to the manufacturer’s protocol. Barcoded sample libraries were quantified
and pooled using Qubit fluorometric quantification (Thermo Fisher Scientific, USA) and Bioanalyzer (Agilent, USA).
Libraries were sequenced on an Illumina Novaseq in a 26 X 8 X 91 bp configuration.

## 5' scRNA-seq data analysis

FASTQ files were processed using Cellranger 6.1.1 and mapped to a custom combined reference with HXB2 HIV reference genome 
added to human GRCh38 reference FASTA and GTF files.
Following the workflow described by Amezquita et al., single-cell RNA-seq analysis was carried out with R/Bioconductor packages (47).
Quality control was carried out for each sample separately with functions from the scuttle package v.1.4.0. (48), cells with low
gene content (below 10e2.5 to 10e2.75, depending on the sample) and high mitochondrial gene content (>3 median absolute deviations; 
default scuttle setting) were removed from further analyses. Genes with zero expression across all cells were removed from the matrix. 
Cell cycle scores were calculated with the scuttle::cyclone() function using human cell cycle marker genes provided by scuttle.
The different count matrices across all samples were then scaled for sequencing depth differences and log-transformed using multiBatchNorm
from the batchelor v.1.10.0 (49, 50) package. Cell types were annotated with SingleR v. 1.8.1. using celldex::HumanPrimaryCellAtlasData() (51).
Differentially expressed genes (FDR<0.05) were detected using the pseudoBulkDGE function from the scran package (52), which uses the
quasi-likelihood method implemented by edgeR. Gene Ontology (GO) analysis was performed using the enrichGO function from
clusterProfiler v.3.10.1 (54).

For more details, see the sessionInfo.txt file here.

## References

47. Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith ML, Huber W, Morgan M, Gottardo R, Hicks SC. Publisher Correction: Orchestrating single-cell analysis with Bioconductor. Nat Methods. 2020;17(2):242. Epub 2019/12/13. doi: 10.1038/s41592-019-0700-8. PubMed PMID:
851 31827272.
48. McCarthy DJ, Campbell KR, Lun AT, Wills QF. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics.  2017;33(8):1179-86. Epub 2017/01/16. doi: 10.1093/bioinformatics/btw777. PubMed PMID: 28088763; PMCID: PMC5408845.
49. Lun AT, Bach K, Marioni JC. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biol. 2016;17:75. Epub 2016/04/29. doi:10.1186/s13059-016-0947-7. PubMed PMID: 27122128; PMCID: PMC4848819.
50. Lun AT. Further MNN algorithm development. https://MarioniLabgithubio/FurtherMNN2018/theory/descriptionhtml. 2018.
51. Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, Chak S, Naikawadi RP, Wolters PJ, Abate AR, Butte AJ, Bhattacharya M. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat Immunol. 2019;20(2):163-72. Epub 2019/01/16. doi: 10.1038/s41590-018-0276-y. PubMed PMID: 30643263; PMCID: PMC6340744.
52. Lun AT, McCarthy DJ, Marioni JC. A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Res. 2016;5:2122. Epub 2016/12/06. doi:10.12688/f1000research.9501.2. PubMed PMID: 27909575; PMCID: PMC5112579.
53. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-40. Epub 2009/11/17. doi: 10.1093/bioinformatics/btp616. PubMed PMID: 19910308; PMCID: PMC2796818.
54. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-7. Epub 2012/03/30. doi:10.1089/omi.2011.0118. PubMed PMID: 22455463; PMCID: PMC3339379.
