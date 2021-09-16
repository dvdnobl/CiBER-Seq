# CiBER-Seq
Code to reproduce data analysis for identifying guides which cause significant differential expression in CiBER-Seq.

 - Data from the CiBER-Seq screen is put in the "data" directory. 

 - Metadata for the data from the CiBER-Seq screen needed for DESeq2 analysis is put in the "meta" directory. 

 - Plots generated from EDA, QC, and DESeq2 are put in the "plots" directory.

 - Code for EDA, QC, and running DESeq2 is saved as the "DESeq2_his4.R" file.

 - Code used for meta-analysis is in "meta-analysis.R"

R version 4.0.5 is required.

Require dependencies (latest version for all):
 - tidyverse
 - metaRNAseq
 - metap
 - DESeq2
 - pheatmap
 - metaRNASeq

Paper outlining the theory and protocol for CiBER-Seq can be found here: https://www.science.org/doi/10.1126/science.abb9662

![370_abb9662_fa](https://user-images.githubusercontent.com/49042863/133697211-6f69f274-5f5c-468c-8998-69ae6ef37ebe.jpeg)
