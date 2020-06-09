# Cindrova-Davies_Burton_Menstrual

RNASeq  analysis for Endometrial organoids compare to Menstrual  with 7 individuals using DESeq2 package in R (v3.6.2).


Code Release to accompany paper: [![DOI](https://zenodo.org/badge/xx.svg)](https://zenodo.org/badge/latestdoi/xx)


## Step 1: Get the samples information
  ### [CTR_gjb2_0010-Control_Menstrual_SampleTable.txt](Figures_Tables/CTR_gjb2_0010-Control_Menstrual_SampleTable.txt)

## Step 2: DESeq2 Analysis

                              Design formula ~ individual + condition


| Files | Name   |
| ----------------------------- | --- |
|CTR_gjb2_0010-PCAplot_rld.menstrual_TopN2000.pdf | [[PDF](Figures_Tables/CTR_gjb2_0010-PCAplot_rld.menstrual_TopN2000.pdf)] |
|CTR_gjb2_0010-Menstrual_Volcano_Plot.pdf |  [[PDF](Figures_Tables/CTR_gjb2_0010-Menstrual_Volcano_Plot.pdf)]|
|CTR_gjb2_0010-Heatmap_Control_Menstrual_padj0.05_lfc1_top101.pdf |[[PDF](Figures_Tables/CTR_gjb2_0010-Heatmap_Control_Menstrual_padj0.05_lfc1_top101.pdf)]|
|CTR_gjb2_0010-Menstrual_Control_sigDEGs_N101_List.csv|[[CSV](Figures_Tables/CTR_gjb2_0010-Menstrual_Control_sigDEGs_N101_List.csv)]|
|CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stress_N61.pdf|[[PDF](Figures_Tables/CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stress_N61.pdf)]|
|CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stress_N61_B70_B75.pdf|[[PDF](Figures_Tables/CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stress_N61_B70_B75.pdf)]|
|CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stemcell_N38_Cat7_B70_B75.pdf|[[PDF](Figures_Tables/CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stemcell_N38_Cat7_B70_B75.pdf)]|
|CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stress_N61_B70_B75.pdf|[[PDF](Figures_Tables/CTR_gjb2_0010-Heatmap_endometrial_Menstrual_stress_N61_B70_B75.pdf)]|
|Menstrual_StemCell_markers.xls|[[XLS](Figures_Tables/Menstrual_StemCell_markers.xls)]|




## Software Versions & Methods

````
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.4

Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /storage/Software/packages/R-3.6.2/lib/libRblas.so
LAPACK: /storage/Software/packages/R-3.6.2/lib/libRlapack.so

Random number generation:
 RNG:     Mersenne-Twister
 Normal:  Inversion
 Sample:  Rounding

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ComplexHeatmap_2.5.1        apeglm_1.8.0               
 [3] limma_3.42.2                ggalt_0.4.0                
 [5] dplyr_0.8.5                 plyr_1.8.6                 
 [7] biomaRt_2.42.1              reshape2_1.4.4             
 [9] ggrepel_0.8.2               pheatmap_1.0.12            
[11] cowplot_1.0.0               RColorBrewer_1.1-2         
[13] ggplot2_3.3.0               DESeq2_1.26.0              
[15] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
[17] BiocParallel_1.20.1         matrixStats_0.56.0         
[19] Biobase_2.46.0              GenomicRanges_1.38.0       
[21] GenomeInfoDb_1.22.1         IRanges_2.20.2             
[23] S4Vectors_0.24.4            BiocGenerics_0.32.0               
````

## Contact

Contact Xiaohui Zhao (xz289 -at- cam.ac.uk)
