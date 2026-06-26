#!/bin/env -S Rscript --vanilla
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH --exclude cc011        # unreliable node
#SBATCH --mail-type END,FAIL   # get an email when this finishes or fails

library(rtracklayer)
library(ggplot2)
library(ggplot2)
bw <- import(format = "BigWig", con = "ENCFF444JRQ.bigWig")

bins_50kb <- tileGenome(seqlengths(bw), tilewidth = 50000, cut.last.tile.in.chrom=TRUE)

bin_result <- binnedAverage(bins_50kb, numvar = coverage(bw, weight = "score"), varname = "avg_score")

data_ranges_df <- bin_result %>% as.data.frame()

df <- read.csv("HCT116.50000.E0-E128.trans.eigvecs.pq.no_empty_fields.csv")

clust <- read.table("HCT116.50000.E1-E128.kmeans_sm.tsv", header=TRUE)
k9 <- clust[, c("chrom", "start", "end", "kmeans_sm9") ]
combo <- merge(df, k9)

small <- combo[,c("chrom","start","end","GC","arm","armlen","centel","centel_abs","is_bad","kmeans_sm9","E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10")]
library(tidyr)
longer <- small %>% pivot_longer(cols=c("E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"), values_to = "value")
bin_result
head(data_ranges_df)
longer[1:20, ]
longer[1:40, ]
longer[30:40, ]
longer[3010:3040, ]
longer2 <- longer
longer2 <- longer %>% mutate(start = start +1 )
longer2
big_combo <- merge(longer2, data_ranges_df)
ls()
head(small)
small2 <- small %>% mutate(start = start +1)
big_combo <- merge(small2, data_ranges_df)
head(big_combo)
head(data_ranges_df)
data_ranges_df2 <- data_ranges_df %>% mutate(chrom = seqnames)

# Did we want all?
big_combo <- merge(small2, data_ranges_df2, all = TRUE)
big_combo <- merge(small2, data_ranges_df2, all = FALSE)

plot_H3K9 <- ggplot(big_combo, aes(x=centel, y=1, color = avg_score)) +
                geom_tile() +
                scale_color_viridis_c() +
                facet_grid(.~kmeans_sm9)
plot_H3K9
