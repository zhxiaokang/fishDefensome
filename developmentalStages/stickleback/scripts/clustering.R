# Cluster the genes with respect to their expression pattern (slope) along time series

library(tscR)

# load raw expression data
table.early.morula <- read.table('../output/dfsm_early_morula_gene_norm.tsv')
table.late.morula <- read.table('../output/dfsm_late_morula_gene_norm.tsv')
table.mid.gastrula <- read.table('../output/dfsm_mid_gastrula_gene_norm.tsv')
table.early.organogenesis <- read.table('../output/dfsm_early_organogenesis_gene_norm.tsv')
table.hph24 <- read.table('../output/dfsm_24_hours_post_hatching_gene_norm.tsv')

# df of all stages
df.all.stages <- data.frame(gene.id = table.early.morula$V3, category = table.early.morula$V1, 
                            gene.name = table.early.morula$V2)
df.all.stages$early_morula <- apply(table.early.morula[, 4:6], 1, mean)
df.all.stages$late.morula <- apply(table.late.morula[, 4:6], 1, mean)
df.all.stages$mid.gastrula <- apply(table.mid.gastrula[, 4:6], 1, mean)
df.all.stages$early.organogenesis <- apply(table.early.organogenesis[, 4:6], 1, mean)
df.all.stages$hph24 <- apply(table.hph24[, 4:6], 1, mean)

# do clustering
df.cluster <- df.all.stages[, 4:8]
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

gene.id <- df.all.stages$gene.id
gene.name <- df.all.stages$gene.name
category <- df.all.stages$category

for (i in seq(1, 16)) {
  index.cluster <- which(sclust$clustering == i)
  gene.id.cluster <- gene.id[index.cluster]
  gene.name.cluster <- gene.name[index.cluster]
  category.cluster <- category[index.cluster]
  df.output <- data.frame("gene.id" = gene.id.cluster, "gene.name" = gene.name.cluster, "category" = category.cluster)
  
  write.table(df.output, paste("../output/cluster_", i, "_genes.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
}
