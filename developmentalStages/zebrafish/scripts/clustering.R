# Cluster the genes with respect to their expression pattern (slope) along time series

library(tscR)

# load raw expression data

# load raw expression data
table.all <- read.table('../output/dfsm_tpm_picked.txt', sep = '\t', header = T)
table.all[is.na(table.all)] <- 0

df.all <- table.all[table.all$Gene_ID != "", ]

# do clustering
df.cluster <- df.all[, 5:9]
## replace 0 with 1e-6
df.cluster[df.cluster == 0] <- 1e6
df.cluster[df.cluster == 1e6] <- min(df.cluster)/2
df.cluster <- log2(df.cluster)

time.points <- colnames(df.cluster)
times <- seq(1, 5)

sDist <- slopeDist(df.cluster, times)
sclust <- getClusters(sDist, k = 16)

for (i in seq(0, 3)) {
  pdf(paste("../output/cluster_group", as.character(i + 1), ".pdf", sep = ""))
  plotCluster(data = df.cluster, clust = sclust, ncluster = c(seq(1, 4) + 4 * i))
  dev.off()
}

gene.id <- df.all$Gene_ID
gene.name <- df.all$Gene_Name
category <- df.all$category

for (i in seq(1, 16)) {
  index.cluster <- which(sclust$clustering == i)
  gene.id.cluster <- gene.id[index.cluster]
  gene.name.cluster <- gene.name[index.cluster]
  category.cluster <- category[index.cluster]
  df.output <- data.frame("gene.id" = gene.id.cluster, "gene.name" = gene.name.cluster, "category" = category.cluster)
  
  write.table(df.output, paste("../output/cluster_", i, "_genes.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
}
