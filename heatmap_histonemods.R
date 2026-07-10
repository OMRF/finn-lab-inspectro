#!/bin/env -S Rscript --vanilla
#SBATCH --cpus-per-task 24
#SBATCH --mem-per-cpu 4G
#SBATCH --exclude cc011        # unreliable node
#SBATCH --mail-type END,FAIL   # get an email when this finishes or fails

library(rtracklayer)
library(tidyverse) #load dplyr, ggplot2, 

#Read in arguments from the command line, starting after this script's name
args <- commandArgs(trailingOnly = TRUE)

cluster_file <- args[1]
eigvecs_file <- args[2]
bigwig_file <- args [3]
mod_name <- args[4]

clusterbasename <- basename(cluster_file)
base_out_name <- strsplit(clusterbasename,"\\.") [[1]][1]

print("DEBUG: we made it to 21")

# Use rtracklayer's import function to read in a BigWig file
bw <- import(format = "BigWig", con = bigwig_file)


# Create 50kb bins
bins_50kb <- tileGenome(seqlengths(bw), tilewidth = 50000, cut.last.tile.in.chrom=TRUE)

# squish the results into 50kb bins
bin_result <- binnedAverage(bins_50kb, numvar = coverage(bw, weight = "score"), varname = "avg_score")

print("DEBUG: we made it to 28!!")

# convert binned result into a data frame
data_ranges_df <- bin_result %>% as.data.frame()

df <- read.csv(eigvecs_file)

# Get kmeans_sm9 cluster labels
clust <- read.table(cluster_file, header=TRUE)
k9 <- clust[, c("chrom", "start", "end", "kmeans_sm9") ]

# Add kmeans_sm9 cluster labels to eigenvecs info
combo <- merge(df, k9)

# Leave out most of the Eigenvector columns
small <- combo[,c("chrom","start","end","GC","arm","armlen","centel","centel_abs","is_bad","kmeans_sm9","E0","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10")]

# Convert coordinates to 1-based (instead of 0-based)
small2 <- small %>% mutate(start = start +1)

# add column named "chrom" that replicates the seqnames column
data_ranges_df2 <- data_ranges_df %>% mutate(chrom = seqnames)

print("DEBUG: we made it to 57!!")

#big_combo <- merge(small2, data_ranges_df2, all = TRUE) # Did we want all?
big_combo <- merge(small2, data_ranges_df2, all = FALSE)

plot_Histone <- ggplot(big_combo, aes(x=centel, y=1, color = avg_score)) +
                geom_tile() +
                scale_color_viridis_c() +
                facet_grid(.~kmeans_sm9)

print("DEBUG: we made it to 66!!")

pdf_filename <-paste0(base_out_name,".", mod_name,".kmeanssm9.pdf")
pdf(pdf_filename, width =18, height =2)
print(plot_Histone)
dev.off()

print("DEBUG: we finished!!")
