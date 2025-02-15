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


#write.table(his4_meta, "meta/his4_meta.txt", sep="\t")
#write.table(pgk1_meta, "meta/pgk1_meta.txt", sep="\t")


## Visualizing distribution of expression counts
ggplot(his4_raw_counts) +
  geom_histogram(aes(x=PostL), stat="bin", bins=200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

sum(his4_raw_counts$PostL == 0) / length(his4_raw_counts$PostL)



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

length(unique(barcode_guide$guide))
length(unique(guide.good.targets$Guide))
length(unique(his4_raw_counts$barcode))

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
                             "replicate" = factor(c("L", "R", "L", "R", "L", "R", "L", "R"),
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
                                        design = ~ guide_induction * geno + replicate)

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

## Recreating histogram of distribution of counts after filtering
ggplot(his4_pgk1_b_filt) +
  geom_histogram(aes(x=PostL.his4), stat="bin", bins=200) +
  xlab("Raw expression counts") +
  ylab("Number of guides")

## Recreating mean vs. variance plots after filtering

his4_filt_mean <- apply(his4_pgk1_b_filt[,2:5], 1, mean)
his4_filt_var <- apply(his4_pgk1_b_filt[,2:5], 1, var)

pgk1_filt_mean <- apply(his4_pgk1_b_filt[,6:9], 1, mean)
pgk1_filt_var <- apply(his4_pgk1_b_filt[,6:9], 1, var)

his4_filt_mean_vs_var <- data.frame(his4_filt_mean, his4_filt_var)
pgk1_filt_mean_vs_var <- data.frame(pgk1_filt_mean, pgk1_filt_var)

ggplot(data=his4_filt_mean_vs_var, aes(x=his4_filt_mean, y=his4_filt_var)) +
  geom_point(size=0.2) +
  geom_smooth(method='lm', formula=y~x, color='red') +
  geom_line(aes(x=his4_filt_mean, y=his4_filt_mean), color="blue") + 
  scale_y_log10() +
  scale_x_log10()

ggplot(data=pgk1_filt_mean_vs_var, aes(x=pgk1_filt_mean, y=pgk1_filt_var)) +
  geom_point(size=0.2) +
  geom_smooth(method='lm', formula=y~x, color="red") +
  geom_line(aes(x=pgk1_filt_mean, y=pgk1_filt_mean), color="blue") + 
  scale_y_log10() +
  scale_x_log10() 


## Creating DESeq object

his4_pgk1_b_filt[is.na(his4_pgk1_b_filt)] <- 0
dim(his4_pgk1_b_filt)

write.table(his4_pgk1_b_filt, "data/his4_pgk1_barcode_filt.txt", sep="\t")

his4_pgk1_b_mtx <- his4_pgk1_b_filt[,2:9]

dds_his4_pgk1_b <- DESeqDataSetFromMatrix(countData = his4_pgk1_b_mtx,
                                          colData = his4_pgk1_meta,
                                          design = ~ geno * guide_induction + replicate)
dds_his4_pgk1_b <- DESeq(dds_his4_pgk1_b)

## Plotting dispersion estimates
plotDispEsts(dds_his4_pgk1_b)




## Creating results dataframe for pre vs post guide induction
resultsNames(dds_his4_pgk1_b)
results_atc <- results(dds_his4_pgk1_b, name = 'guide_induction_Post_vs_Pre',
                        alpha = 0.05)

results_geno <- results(dds_his4_pgk1_b, name = 'genohis4.guide_inductionPost',
                        alpha = 0.05)

res_table_atc <- as.data.frame(results_atc)

res_table_geno <- as.data.frame(results_geno)

res_table_atc$barcode <- his4_pgk1_b_filt$barcode
res_table_atc$guide <- grna.assign.barcode.grna.good[match(res_table_atc$barcode,
                                                           grna.assign.barcode.grna.good$barcode),
                                                     "guide"]

res_table_geno$barcode <- his4_pgk1_b_filt$barcode
res_table_geno$guide <- barcode_guide[match(res_table_geno$barcode,
                                                           barcode_guide$barcode),
                                                     "guide"]
res_table_atc <- res_table_atc[,c(8,7,1,2,3,4,5,6)]

res_table_geno <- res_table_geno[,c(8,7,1,2,3,4,5,6)]

#write.table(res_table_atc, "results/deseq2_results_atc.txt",
#            sep="\t")

write.table(res_table_geno, "results/deseq2_results_geno.txt",
            sep="\t")


plotMA(results_geno, ylim=c(-8,8))

res_table_geno_sig <- res_table_geno[res_table_geno$padj <= .05,]


# Results for shrunken LFC
coefs = coef(dds_his4_pgk1_b)
shrunk_res = lfcShrink(dds_his4_pgk1_b, coef='genohis4.guide_inductionPost', res=results_geno)

plotMA(shrunk_res)



