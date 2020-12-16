# plot expression change along time points

library(ggplot2)
library(reshape2)
library(cowplot)

rm(list=ls())

# load raw expression data
table.all <- read.table('../output/dfsm_tpm_picked.txt', sep = '\t', header = T)
table.all[is.na(table.all)] <- 0

df.all <- table.all[table.all$Gene_ID != "", ]

# plot
df.plot <- data.frame(t(df.all[, 5:9]))
colnames(df.plot) <- df.all$Gene_ID
df.plot$stage <- rownames(df.plot)
df.plot$stage <- factor(df.plot$stage, levels = df.plot$stage)

melted.df <- melt(df.plot, id.vars = 'stage')
melted.df$category <- rep(df.all$category, each=5)

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
         geom_line(colour='grey35') + theme(legend.position = "none")
  plots[[i]] <- p
  i = i+1
}

pdf(file='../output/expression_change_along_time_seperate_category_log2.pdf', width = 18, height = 14)
print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], ncol = 3, nrow = 4, labels = categories))
dev.off()

