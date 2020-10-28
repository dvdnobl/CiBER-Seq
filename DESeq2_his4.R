library(tidyverse)
library(DESeq2)
library(pheatmap)

library(metaRNASeq)

## Reading in count tables, extracting RNA count columns
his4_raw_counts <- read.delim("data/all_his4_moreseq_counts.txt")[,c(1,8,9,6,7)]

pgk1_raw_counts <- read.delim("data/all_pgk1_moreseq_counts.txt")[,c(1,8,9,6,7)]

## Changing column names
colnames(his4_raw_counts)
colnames(his4_raw_counts) <- c("barcode", "PreL", "PreR", "PostL", "PostR")

colnames(pgk1_raw_counts)
colnames(pgk1_raw_counts) <- c("barcode", "PreL", "PreR", "PostL", "PostR")

head(his4_raw_counts)
head(pgk1_raw_counts)

## Creating meta tables
his4_meta <- data.frame("sampletype" = colnames(his4_raw_counts)[2:5],
                        "guide_induction" = c("pre", "pre", "post", "post"),
                        "turb" = c("left", "right", "left", "right"))
row.names(his4_meta) <- his4_meta$sampletype
his4_meta <- his4_meta[,2:3]

pgk1_meta <- data.frame("sampletype" = colnames(pgk1_raw_counts)[2:5],
                        "guide_induction" = c("pre", "pre", "post", "post"),
                        "turb" = c("left", "right", "left", "right"))
row.names(pgk1_meta) <- pgk1_meta$sampletype
pgk1_meta <- pgk1_meta[,2:3]


write.table(his4_meta, "meta/his4_meta.txt", sep="\t")
write.table(pgk1_meta, "meta/pgk1_meta.txt", sep="\t")


## Visualizing distribution of expression counts
ggplot(his4_raw_counts) +
  geom_histogram(aes(x=PostL), stat="bin", bins=200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

## Visualizing mean vs. variance in data
his4_mean_counts <- apply(his4_raw_counts[,4:5], 1, mean)
his4_var_counts <- apply(his4_raw_counts[,4:5], 1, var)

his4_mean_vs_var <- data.frame(his4_mean_counts, his4_var_counts)

pgk1_mean_counts <- apply(pgk1_raw_counts[,4:5], 1, mean)
pgk1_var_counts <- apply(pgk1_raw_counts[,4:5], 1, var)

pgk1_mean_vs_var <- data.frame(pgk1_mean_counts, pgk1_var_counts)

ggplot(his4_mean_vs_var) +
  geom_point(aes(x=his4_mean_counts, y=his4_var_counts)) +
  geom_line(aes(x=his4_mean_counts, y=his4_mean_counts), color="blue") +
  scale_y_log10() +
  scale_x_log10()

ggplot(pgk1_mean_vs_var) +
  geom_point(aes(x=pgk1_mean_counts, y=pgk1_var_counts)) +
  geom_line(aes(x=pgk1_mean_counts, y=pgk1_mean_counts), color="green") +
  scale_y_log10() +
  scale_x_log10()


## Creating count matrices
his4_matrix <- his4_raw_counts[,2:5]
pgk1_matrix <- pgk1_raw_counts[,2:5]

## Checking the meta tables
all(colnames(his4_matrix) %in% rownames(his4_meta))
all(colnames(his4_matrix) == row.names(his4_meta))

all(colnames(pgk1_matrix) %in% rownames(pgk1_meta))
all(colnames(pgk1_matrix) == row.names(pgk1_meta))

## Creating DESeq objects
dds_his4 <- DESeqDataSetFromMatrix(countData = his4_matrix, colData = his4_meta, 
                                   design = ~ guide_induction)
dds_pgk1 <- DESeqDataSetFromMatrix(countData = pgk1_matrix, colData = pgk1_meta, 
                                   design = ~ guide_induction)

## Taking log2 transformation of counts for visualization
rld_his4 <- vst(dds_his4, blind=TRUE)
rld_pgk1 <- vst(dds_pgk1, blind=TRUE)

## Plotting PCA 1 and 2 for both matrices
plotPCA(rld_his4, intgroup="guide_induction")
plotPCA(rld_his4, intgroup="turb")

plotPCA(rld_pgk1, intgroup="guide_induction")
plotPCA(rld_pgk1, intgroup="turb")

## Further PCA analysis
rld_mat_his4 <- assay(rld_his4)
rld_mat_pgk1 <- assay(rld_pgk1)

pca_his4 <- prcomp(t(rld_mat_his4))
pca_pgk1 <- prcomp(t(rld_mat_pgk1))

summary(pca_his4)
summary(pca_pgk1)

## Finding correlations and creating heatmaps
rld_cor_his4 <- cor(rld_mat_his4)
rld_cor_pgk1 <- cor(rld_mat_pgk1)

head(rld_cor_his4)
head(rld_cor_pgk1)

pheatmap(rld_cor_his4)
pheatmap(rld_cor_pgk1)


## Uploading barcode to guide matching data, number No_gRNA entries
barcode_guide <- read.delim('data/grna-assign-barcode-grna-good.txt')
barcode_guide$guide[barcode_guide$guide == "No_gRNA"] <- paste("No_gRNA", seq(1:788), sep=" ")


## Attaching guide info into count data sets
## Filtering out rows without guide, file size too large
his4_guide_counts <- inner_join(his4_raw_counts, barcode_guide, by=c("barcode"))
head(his4_guide_counts)
dim(his4_guide_counts)

pgk1_guide_counts <- inner_join(pgk1_raw_counts, barcode_guide, by=c("barcode"))
head(pgk1_guide_counts)
dim(pgk1_guide_counts)

## Aggregating counts by guide
his4_agg <- aggregate(his4_guide_counts[,c(2,3,4,5)], by=list(his4_guide_counts$guide), FUN=sum)
names(his4_agg)[1] <- "guide"
head(his4_agg)
dim(his4_agg)

pgk1_agg <- aggregate(pgk1_guide_counts[,c(2,3,4,5)], by=list(pgk1_guide_counts$guide), FUN=sum)
names(pgk1_agg)[1] <- "guide"
head(pgk1_agg)
dim(pgk1_agg)


## Creating count matrix with entries from his4 and pgk1 matched by guide
his4_pgk1 <- inner_join(his4_agg, pgk1_agg, by=c("guide"), suffix=c(".his4", ".pgk1"))
head(his4_pgk1)
dim(his4_pgk1)


his4_pgk1_meta <- data.frame("sampletype" = colnames(his4_pgk1)[2:9], 
                             "guide_induction" = factor(c("Pre", "Pre", "Post", "Post", 
                                                   "Pre", "Pre", "Post", "Post"),
                                                   levels = c("Pre", "Post")),
                             "turb" = factor(c("L", "R", "L", "R", "L", "R", "L", "R"),
                                              levels = c("L", "R")),
                             "geno" = factor(c("his4", "his4", "his4", "his4", "pgk1",
                                        "pgk1","pgk1","pgk1"), levels=c('pgk1', 'his4')))
rownames(his4_pgk1_meta) <- his4_pgk1_meta$sampletype
his4_pgk1_meta <- his4_pgk1_meta[,c(2,3,4)]

write.table(his4_pgk1_meta, "meta/his4_pgk1_meta.txt", sep="\t")

## Creating count matrix
his4_pgk1_mtx = his4_pgk1[,2:9]

## Creating DESeq object
dds_his4_pgk1 <- DESeqDataSetFromMatrix(countData = his4_pgk1_mtx, colData = his4_pgk1_meta, 
                                        design = ~ guide_induction * geno + turb)

## Taking log2 transformation of counts for visualization
rld_his4_pgk1 <- vst(dds_his4_pgk1, blind=TRUE)

## Plotting PCA 
plotPCA(rld_his4_pgk1, intgroup="guide_induction")
plotPCA(rld_his4_pgk1, intgroup="turb")
plotPCA(rld_his4_pgk1, intgroup="geno")

## Further PCA analysis
rld_mat_his4_pgk1 <- assay(rld_his4_pgk1)

pca_his4_pgk1 <- prcomp(t(rld_mat_his4_pgk1))
summary(pca_his4_pgk1)

pca_his4_pgk1_df <- cbind(his4_pgk1_meta, pca_his4_pgk1$x)
ggplot(pca_his4_pgk1_df) + geom_point(aes(x=PC1, y=PC2, color = geno))



## Trying DESeq2 on counts by barcode, not guide
his4_pgk1_b <- full_join(his4_raw_counts, pgk1_raw_counts, by='barcode', 
                         suffix = c('.his4', '.pgk1'))
his4_pgk1_b_filt <- dplyr::filter(his4_pgk1_b,
                              PreL.his4 > 31, PreR.his4 > 31, PreL.pgk1 > 31, PreR.pgk1 > 31)

his4_pgk1_b_filt[is.na(his4_pgk1_b_filt)] <- 0
dim(his4_pgk1_b_filt)

write.table(his4_pgk1_b_filt, "data/his4_pgk1_barcode_filt.txt", sep="\t")

his4_pgk1_b_mtx <- his4_pgk1_b_filt[,2:9]

dds_his4_pgk1_b <- DESeqDataSetFromMatrix(countData = his4_pgk1_b_mtx,
                                          colData = his4_pgk1_meta,
                                          design = ~ guide_induction * geno + turb)
dds_his4_pgk1_b <- DESeq(dds_his4_pgk1_b)

## Plotting dispersion estimates
plotDispEsts(dds_his4_pgk1_b)

## Creating results dataframe for pre vs post guide induction
resultsNames(dds_his4_pgk1_b)
results_atc <- results(dds_his4_pgk1_b, name = 'guide_induction_Post_vs_Pre',
                        alpha = 0.05)

results_geno <- results(dds_his4_pgk1_b, name = 'guide_inductionPost.genohis4')

res_table_atc <- as.data.frame(results_atc)

res_table_geno <- as.data.frame(results_geno)

res_table_atc$barcode <- his4_pgk1_b_filt$barcode
res_table_atc$guide <- grna.assign.barcode.grna.good[match(res_table_atc$barcode,
                                                           grna.assign.barcode.grna.good$barcode),
                                                     "guide"]

res_table_geno$barcode <- his4_pgk1_b_filt$barcode
res_table_geno$guide <- grna.assign.barcode.grna.good[match(res_table_geno$barcode,
                                                           grna.assign.barcode.grna.good$barcode),
                                                     "guide"]
res_table_atc <- res_table_atc[,c(8,7,1,2,3,4,5,6)]

res_table_geno <- res_table_geno[,c(8,7,1,2,3,4,5,6)]

write.table(res_table_atc, "results/deseq2_results_atc.txt",
            sep="\t")

write.table(res_table_geno, "results/deseq2_results_geno.txt",
            sep="\t")

plotMA(results_atc_unshrunken)

plotMA(results_geno, ylim=c(-4,4))


## Looking at significant results
sig_barcodes <- res_table_geno[res_table_geno$padj <= .1,]
dim(sig_barcodes)  
dim(res_table_atc)

res_valid_guides <- drop_na(res_table_geno, guide)
res_valid_guides <- res_valid_guides[order(res_valid_guides$guide),]

write.table(res_valid_guides, "results/deseq2_results_geno_guides.txt",
            sep="\t")

dim(res_valid_guides)
length(unique(res_valid_guides$guide))

## Adding yorf info onto result table
res_yorf <- merge(res_valid_guides, guide.good.targets[,c(1,3)], by.x = 'guide',
                  by.y = 'Guide')
res_yorf <- drop_na(res_yorf, Yorf1)
colnames(res_yorf)[9] <- 'gene'

## Creating a new dataframe with combined p-values using Fisher's method
res_combined <- data.frame(gene = unique(res_yorf$gene))          

combined_p_vals <- c()
avg_log2fc <- c()
for (orf in res_combined$gene) {
  v <- res_yorf[res_yorf$gene == orf,]
  f_method <- fishercomb(v[,8])$adjpval
  combined_p_vals <- c(combined_p_vals, f_method)
  avg_log2fc <- c(avg_log2fc, mean(v[,4]))
}

res_combined$avg_log2fc <- avg_log2fc
res_combined$comb_adjpval <- combined_p_vals

res_combined_sig <- res_combined[res_combined$comb_adjpval <= .05,]

write.table(res_combined, 'results/res_combined_p.txt', sep='\t')
write.table(res_combined_sig, 'results/res_combined_p_sig.txt', sep='\t')


## Comparing result tables to see if genes were found to be DE with F's method

res_uncombined <- data.frame(gene = unique(res_yorf$gene))
min_pval <- c()
for (orf in res_combined$gene)

