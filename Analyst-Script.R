Data <- read.table("gene_exp_2.diff", header = TRUE)
#data ordered by q_value ascending 
q_ordered_DEGs <- Data[order(Data$q_value),]

#top 10 differentially expressed genes
top_10 <- as.data.frame(q_ordered_DEGs[1:10,])
top_10 <- top_10[,c(2,3,8,9,10, 12,13)]

#top_10 log fold hist 
Top10_log2.fold_change <-top_10$log2.fold_change.
Top_10_DEGs_hist <-hist(Top10_log2.fold_change, breaks = 20, main = 'Top Ten DEGs')

#data frame of only significant genes 
significant_genes <- subset(Data, significant == 'yes')

#significant genes hist
significant_DEGs_log2.fold_change <-significant_genes$log2.fold_change.
plot_2 <-hist(significant_DEGs_log2.fold_change, breaks =20, main = 'Significant DEGs')

#DEGs having p <0.01 
pdegs <- subset(significant_genes, p_value <0.01)

#makinf subsets of up and down regulated genes based of significant genes that have p<0.01
up_pdegs <- subset(pdegs, log2.fold_change. >= 0)
down_pdegs <- subset(pdegs, log2.fold_change. < 0)

#making subsets of up and down regulated genes 
up_significant_DEGs <-subset(significant_genes, log2.fold_change. >= 0)
down_significant_DEGs <- subset(significant_genes, log2.fold_change. < 0)

up_DEG_names <-write.csv(up_significant_DEGs$gene, file = 'upregulated_genes', row.names = FALSE, quote = FALSE)
down_DEG_names <-write.csv(down_significant_DEGs$gene, file = 'downregulated_genes', row.names = FALSE, quote = FALSE)

up_DAVID <- read.csv("upregulated-cluster.csv", sep = "\t", header = F)
down_DAVID <- read.csv("downregulated-cluster.csv", sep = "\t", header = F)
