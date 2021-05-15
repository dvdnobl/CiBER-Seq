library(tidyverse)
library(metaRNASeq)
library(metap)

res_table_geno <- read.delim('results/deseq2_results_geno.txt', sep='\t')
guide.good.targets <- read.delim('data/guide.good.targets_plus_empty.txt')[,1:4]
mpralm_his4 <- read.delim('data/his4_pooled_sum_mpralm.txt')
guide_scores <- drop_na(read.delim('data/rescore.txt'))

## Adding guide scores to guide.good.targets
guide.good.targets.scores <- merge(guide.good.targets, guide_scores, by='TargetLoc')

## Looking at significant results
sig_barcodes <- res_table_geno[res_table_geno$padj <= .05,]

res_valid_guides <- drop_na(res_table_geno, guide)
res_valid_guides <- res_valid_guides[order(res_valid_guides$guide),]

dim(res_valid_guides)
length(unique(res_valid_guides$guide))

## Adding yorf info onto result table
res_yorf <- merge(res_valid_guides, guide.good.targets[,1:3], by.x = 'guide',
                  by.y = 'Guide')
res_yorf <- drop_na(res_yorf, Yorf1)
colnames(res_yorf)[9] <- 'gene'

sig_genes <- res_yorf[res_yorf$padj <= .05,]
length(unique(sig_genes$Yorf1))

## Creating a new dataframe with combined p-values using Fisher's method
res_combined <- data.frame(gene = unique(res_yorf$gene))          

combined_p_vals <- c()
avg_log2fc <- c()
for (orf in res_combined$gene) {
  v <- res_yorf[res_yorf$gene == orf,]
  f_method <- fishercomb(v$padj)$adjpval
  combined_p_vals <- c(combined_p_vals, f_method)
  avg_log2fc <- c(avg_log2fc, mean(v[,4]))
}
#fishercomb(res_yorf[res_yorf$gene == 'YCR012W', "padj"])

res_combined$avg_log2fc <- avg_log2fc
res_combined$comb_adjpval <- combined_p_vals

res_combined_sig <- res_combined[res_combined$comb_adjpval <= .05,]

write.table(res_combined, 'results/res_combined_p.txt', sep='\t')
write.table(res_combined_sig, 'results/res_combined_p_sig.txt', sep='\t')


## Comparing result tables to see if genes were found to be DE with F's method

res_uncombined <- data.frame(gene = unique(res_yorf$gene))
n = c()
min_pval <- c()
for (orf in res_combined$gene){
  v <- res_yorf[res_yorf$gene == orf, "padj"]
  min_pval <- c(min_pval, min(v))
  n = c(n, length(v))
}


res_uncombined$min_padj <- min_pval
res_uncombined_sig <- res_uncombined[res_uncombined$min_padj <= .05,]
res_uncombined_sig[order(res_uncombined_sig$min_padj),]

setdiff(res_combined_sig$gene, res_uncombined_sig$gene)
## Fisher's method did not show any DE genes that wouldn't have been found
## by looking at just padj values

## Comparing with Ryan's analysis

mpralm_his4_yorf <- na.omit(data.frame(yorf=unique(mpralm_his4$Yorf1)))
min_padjs <- c()
for (gene in mpralm_his4_yorf$yorf) {
  v <- mpralm_his4[mpralm_his4$Yorf1 == gene, "adj.P.Val"]
  min_padjs <- c(min_padjs, min(v))
}
mpralm_his4_yorf$min_padj_val <- min_padjs

mpralm_his4_yorf_sig <- mpralm_his4_yorf[mpralm_his4_yorf$min_padj_val <= .05,]

mpralm_unique_sig <- na.omit(mpralm_his4[mpralm_his4$adj.P.Val <= .05, c(1, 9)])
length(unique(mpralm_unique_sig$Guide))

mpralm_sig_genes <- unique(mpralm_unique_sig$Yorf1)
deseq2_sig_genes <- unique(sig_genes$Yorf1)

1 - (length(setdiff(deseq2_sig_genes, mpralm_sig_genes)) / length(deseq2_sig_genes))

## Adding guide info to results
res_guide <- merge(res_table_geno, guide.good.targets.scores, by.x='guide',
                   by.y='Guide')
yorf_stouffer <- data.frame(gene=unique(res_guide$Yorf1))

stouffer_p <- c()
for (orf in yorf_stouffer$gene){
  guides <- res_guide[res_guide$Yorf1 == orf,]
  if (nrow(guides) == 1){
    stouffer_p <- c(stouffer_p, guides$padj[1])
  }
  else {
    stouffer_p <- c(stouffer_p, sumz(guides$padj, weights=guides$response)$p[1,1])
  }
}

yorf_stouffer$combined_p <- stouffer_p
yorf_stouffer_sig <- yorf_stouffer[yorf_stouffer$combined_p <= .05,]

setdiff(yorf_stouffer_sig$gene, res_uncombined_sig$gene)
     
