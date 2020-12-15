# plot expression change along time points

library(ggplot2)
library(reshape2)
library(cowplot)

rm(list=ls())

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


# plot
df.plot <- data.frame(t(df.all.stages[, 4:8]))
colnames(df.plot) <- df.all.stages$gene.id
df.plot$stage <- rownames(df.plot)
df.plot$stage <- factor(df.plot$stage, levels = df.plot$stage)

melted.df <- melt(df.plot, id.vars = 'stage')
melted.df$category <- rep(df.all.stages$category, each=5)


pdf(file='../output/expression_change_along_time.pdf', width = 8, height = 5)
print(ggplot(melted.df, aes(x=stage, y=value, group = variable, color = category)) + geom_line())
dev.off()

pdf(file='../output/expression_change_along_time_log2.pdf', width = 8, height = 5)
print(ggplot(melted.df, aes(x=stage, y=log2(value), group = variable, color = category)) + geom_line())
dev.off()

categories <- unique(melted.df$category)

plots <- list()
i <- 1
for (category in categories){
  melted.df.category <- melted.df[melted.df$category==category,]
  p <- ggplot(data=melted.df.category,
              aes(x=stage, y=value, group = variable, color=variable)) +
    geom_line() + theme(legend.position = "none")
  plots[[i]] <- p
  i = i+1
}

pdf(file='../output/expression_change_along_time_seperate_category.pdf', width = 18, height = 14)
print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], ncol = 3, nrow = 4, labels = categories))
dev.off()


plots <- list()
i <- 1
for (category in categories){
  melted.df.category <- melted.df[melted.df$category==category,]
  p <- ggplot(data=melted.df.category,
              aes(x=stage, y=log2(value), group = variable)) +
    geom_line(color='grey35') + theme(legend.position = "none")
  plots[[i]] <- p
  i = i+1
}

pdf(file='../output/expression_change_along_time_seperate_category_log2.pdf', width = 18, height = 14)
print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], ncol = 3, nrow = 4, labels = categories))
dev.off()


