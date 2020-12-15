DEA <- function(count.table, norm.control, norm.treat, file.dea) {
  count.table <- read.table(count.table, header = TRUE, row.names = 1, sep = '\t')
  num.samples <- ncol(count.table)
  num.sample.control <- num.samples/2
  num.sample.treat <- num.samples/2
  count.control <- count.table[, 1:num.sample.control]
  count.treat <- count.table[, (num.sample.control+1):(num.sample.control+num.sample.treat)]
  
  # number of samples in control and treat groups (should be the same if it's a pair test)
  num.sample <- ncol(count.table)
  num.sample.control <- ncol(count.control)
  num.sample.treat <- ncol(count.treat)
  
  # samples of two groups
  sample.control <- colnames(count.control)
  sample.treat <- colnames(count.treat)
  
  # save gene list in gene.list for extracting gene names later
  gene.list <- rownames(count.table)
  
  # get the sample id
  samples <- colnames(count.table)
  
  # define the group
  group <- factor(c(rep('control', num.sample.control), rep('treat', num.sample.treat)))
  
  # The design matrix
  design <- model.matrix(~group)
  
  # normalize the two groups and save the normalized count table
  y.control <- DGEList(counts = count.control, genes = gene.list)
  y.treat <- DGEList(counts = count.treat, genes = gene.list)
  
  y.control <- calcNormFactors(y.control, method="TMM")
  count.table.control.norm <- cpm(y.control)
  write.table(count.table.control.norm, norm.control, quote = FALSE, sep = "\t")

  y.treat <- calcNormFactors(y.treat, method="TMM")
  count.table.treat.norm <- cpm(y.treat)
  write.table(count.table.treat.norm, norm.treat, quote = FALSE, sep = "\t")
  
  # Put the data into a DGEList object
  y <- DGEList(counts = count.table, genes = gene.list)
  
  # do DEA
  
  # Filtering
  countsPerMillion <- cpm(y)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) > 1)
  y <- y[keep, ]
  
  # Normalization
  y <- calcNormFactors(y, method="TMM")
  
  y$samples$group <- group
  
  rownames(design) <- colnames(y)
  
  # Estimating the dispersion
  
  # estimate the NB dispersion for the dataset
  y <- estimateDisp(y, design, robust = TRUE)
  
  # Differential expression
  
  # determine differentially expressed genes
  # fit genewise glms
  fit <- glmFit(y, design)
  
  # conduct likelihood ratio tests for tumour vs normal tissue differences and show the top genes
  lrt <- glmLRT(fit)
  
  # the DEA result for all the genes
  # dea <- lrt$table
  toptag <- topTags(lrt, n = nrow(y$genes), p.value = 1)
  dea <- toptag$table  # just to add one more column of FDR
  dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
  
  # differentially expressed genes
  toptag <- topTags(lrt, n = nrow(y$genes), p.value = 0.05)
  deg <- toptag$table
  if (!is.null(deg)) {
    deg <- deg[order(deg$FDR, -abs(deg$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC
  }
  
  # save the DEA result and DEGs to files
  write.table(dea, file.dea, row.names = F, quote = FALSE, sep = '\t')
}