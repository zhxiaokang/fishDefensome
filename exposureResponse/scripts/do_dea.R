# do DEA

library('edgeR')
source('DEA.R')

args <- commandArgs(TRUE)

count.table <- args[1]
norm.control <- args[2]
norm.treat <- args[3]
file.dea <- args[4]

DEA(count.table, norm.control, norm.treat, file.dea)