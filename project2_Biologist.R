setwd("C:/Users/crist/Desktop")
#read exp files to R
gene_exp.diff <- read.table("gene_exp.diff",header = T)

#order by q_value from lowest to highest
ord_gene_exp <- gene_exp.diff[order(gene_exp.diff$q_value),]

#select top 10 and necessary columns
top_ten_de_gene <- ord_gene_exp[1:10,]
top_ten_de_gene <- top_ten_de_gene[,c(3, 8, 9, 10, 12, 13)]


#hist for log2.foldchange for all genes
hist(ord_gene_exp$log2.fold_change.)


#select significant genes from the original table
sig_gene_diff <- gene_exp.diff[gene_exp.diff$significant=="yes",]

#hist for significant genes only
hist(sig_gene_diff$log2.fold_change.)

#hist for upregulated
up_reg_gene_diff <- sig_gene_diff[sig_gene_diff$log2.fold_change.>0,]
ord_up_reg_gene_diff <- up_reg_gene_diff[order(up_reg_gene_diff$log2.fold_change),]
down_reg_gene_diff <- sig_gene_diff[sig_gene_diff$log2.fold_change.<0,]
ord_down_reg_gene_diff <- down_reg_gene_diff[order(down_reg_gene_diff$log2.fold_change),]
write.table(up_reg_gene_diff$gene,"up_regulated_degenes.txt", row.names = FALSE, quote = F)
write.table(down_reg_gene_diff$gene,"down_regulated_degenes.txt", row.names = F, quote = F)




#Biologist part

#load fpkm tables
P01 <- read.table("fpkm_tracking/P0_1_genes.fpkm_tracking",header = T)
P02 <- read.table("fpkm_tracking/P0_2_genes.fpkm_tracking",header = T)
P41 <- read.table("fpkm_tracking/P4_1_genes.fpkm_tracking",header = T)
P42 <- read.table("fpkm_tracking/P4_2_genes.fpkm_tracking",header = T)
P71 <- read.table("fpkm_tracking/P7_1_genes.fpkm_tracking",header = T)
P72 <- read.table("fpkm_tracking/P7_2_genes.fpkm_tracking",header = T)
PAd1 <- read.table("fpkm_tracking/Ad_1_genes.fpkm_tracking",header = T)
PAd2 <- read.table("fpkm_tracking/Ad_1_genes.fpkm_tracking",header = T)
datas <- c("P01","P02","P41", "P42", "P71", "P72", "PAd1", "PAd2")
for (data in datas){
  d=get(data)
  colnames(d)[5] <- "gene"
  assign(data,d)
}
  
#find genes for prominent GO terms ---Fig.1D---Sacromere

Comparison_Pdlim5 <- NULL
Comparison_Pygm <- NULL
Comparison_Myoz2 <- NULL
Comparison_Des <- NULL
Comparison_Csrp3 <- NULL
Comparison_Tcap <- NULL
Comparison_Cryab <- NULL

for (data in datas){
  d=get(data)
  Comparison_Pdlim5 <- rbind(Comparison_Pdlim5,d[d$gene=="Pdlim5",10])
  Comparison_Pygm <- rbind(Comparison_Pygm,d[d$gene=="Pygm",10])
  Comparison_Myoz2 <- rbind(Comparison_Myoz2,d[d$gene=="Myoz2",10])
  Comparison_Des <- rbind(Comparison_Des,d[d$gene=="Des",10])
  Comparison_Csrp3 <- rbind(Comparison_Csrp3,d[d$gene=="Csrp3",10])
  Comparison_Tcap <- rbind(Comparison_Tcap,d[d$gene=="Tcap",10])
  Comparison_Cryab <- rbind(Comparison_Cryab,d[d$gene=="Cryab",10])
}
Sacromere_comp <- cbind(Comparison_Pdlim5,
                        Comparison_Pygm,
                        Comparison_Myoz2,
                        Comparison_Des,
                        Comparison_Csrp3,
                        Comparison_Tcap,Comparison_Cryab)
Sacromere_comp <- as.data.frame(Sacromere_comp)
rownames(Sacromere_comp) <- datas
colnames(Sacromere_comp) <- c("Pdlim5","Pygm","Myoz2","Des","Csrp3","Tcap","Cryab")



#find genes for prominent GO terms ---Fig.1D---Mitochondria

Comparison_Mpc1 <- NULL
Comparison_Prdx3 <- NULL
Comparison_Acat1 <- NULL
Comparison_Echs1 <- NULL
Comparison_Slc25a11 <- NULL
Comparison_Phyh <- NULL

for (data in datas){
  d=get(data)
  Comparison_Mpc1 <- rbind(Comparison_Mpc1,d[d$gene=="Mpc1",10])
  Comparison_Prdx3 <- rbind(Comparison_Prdx3,d[d$gene=="Prdx3",10])
  Comparison_Acat1 <- rbind(Comparison_Acat1,d[d$gene=="Acat1",10])
  Comparison_Echs1 <- rbind(Comparison_Echs1,d[d$gene=="Echs1",10])
  Comparison_Slc25a11 <- rbind(Comparison_Slc25a11,d[d$gene=="Slc25a11",10])
  Comparison_Phyh <- rbind(Comparison_Phyh,d[d$gene=="Phyh",10])
}
Mito_comp <- as.data.frame(cbind(
                                 Comparison_Prdx3,
                                 Comparison_Acat1,
                                 Comparison_Echs1,
                                 Comparison_Slc25a11,
                                 Comparison_Phyh))
rownames(Mito_comp) <- datas
colnames(Mito_comp) <- c("Prdx3","Acat1","Echs1","Slc25a11","Phyh")


#find genes for prominent GO terms ---Fig.1D---Cell Cycle
Comparison_Cdc7 <- NULL
Comparison_E2f8 <- NULL
Comparison_Cdk7 <- NULL
Comparison_Cdc26 <- NULL
Comparison_Cdc6 <- NULL
Comparison_E2f1 <- NULL
Comparison_Cdc27 <- NULL
Comparison_Bora <- NULL
Comparison_Cdc45 <- NULL
Comparison_Rad51 <- NULL
Comparison_Aurkb <- NULL
Comparison_Cdc23 <- NULL

for (data in datas){
  d=get(data)
  Comparison_Cdc7 <- rbind(Comparison_Cdc7,d[d$gene=="Cdc7",10])
  Comparison_E2f8 <- rbind(Comparison_E2f8,d[d$gene=="E2f8",10])
  Comparison_Cdk7 <- rbind(Comparison_Cdk7,d[d$gene=="Cdk7",10])
  Comparison_Cdc26 <- rbind(Comparison_Cdc26,d[d$gene=="Cdc26",10])
  Comparison_Cdc6 <- rbind(Comparison_Cdc6,d[d$gene=="Cdc6",10])
  Comparison_E2f1 <- rbind(Comparison_E2f1,d[d$gene=="E2f1",10])
  Comparison_Cdc27 <- rbind(Comparison_Cdc27,d[d$gene=="Cdc27",10])
  Comparison_Bora <- rbind(Comparison_Bora,d[d$gene=="Bora",10])
  Comparison_Cdc45 <- rbind(Comparison_Cdc45,d[d$gene=="Cdc45",10])
  Comparison_Rad51 <- rbind(Comparison_Rad51,d[d$gene=="Rad51",10])
  Comparison_Cdc23 <- rbind(Comparison_Cdc23,d[d$gene=="Cdc23",10])
  Comparison_Aurkb <- rbind(Comparison_Aurkb,d[d$gene=="Aurkb",10])
}
Cell_comp <- as.data.frame(cbind(
  Comparison_Cdc7,
  Comparison_E2f8,
  Comparison_Cdk7,
  Comparison_Cdc26,
  Comparison_Cdc6,
  Comparison_E2f1,
  Comparison_Cdc27,
  Comparison_Cdc45,
  Comparison_Rad51,
  Comparison_Aurkb,
  Comparison_Cdc23))
rownames(Cell_comp) <- datas
colnames(Cell_comp) <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","E2f1","Cdc27","Cdc45","Rad51","Cdc23","Aurkb")


#take average between two replicates

mean_Sacro <- rbind(colMeans(Sacromere_comp[1:2,]),
                    colMeans(Sacromere_comp[3:4,]),colMeans(Sacromere_comp[5:6,]),
                    colMeans(Sacromere_comp[7:8,]))
rownames(mean_Sacro) <- datas[c(1,3,5,7)]

mean_Mito <- rbind(colMeans(Mito_comp[1:2,]),
                    colMeans(Mito_comp[3:4,]),colMeans(Mito_comp[5:6,]),
                    colMeans(Mito_comp[7:8,]))
rownames(mean_Mito) <- datas[c(1,3,5,7)]

mean_Cell <- rbind(colMeans(Cell_comp[1:2,]),
                   colMeans(Cell_comp[3:4,]),colMeans(Cell_comp[5:6,]),
                   colMeans(Cell_comp[7:8,]))
rownames(mean_Cell) <- datas[c(1,3,5,7)]

#plots for first two plots
par(bg="gray")
plots <- c("mean_Sacro", "mean_Mito", "mean_Cell")
titles <- c("FPKM values in Sacromere", "FPKM values in Mitochonria", "FPKM values in Cell Cycle")
count <- 1
for (plot in plots){
  title <- titles[count]
  pt <- get(plot)
  matplot(pt, type = c("b"),pch=1,col = 1:7,cex = 1.5, main=title, ylab = "",xaxt='n') #plot
  axis(1,at=c(1,2,3,4),labels=c("P0", "P4", "P7", "Ad"))
  legend("topleft", legend = colnames(pt), col=1:7, pch=1,cex = 0.6) # optional legend
  count <- count +1
}
#plot for last plots
matplot(mean_Cell, type = c("b"),pch=1,col = 1:7,cex = 1.5, main=titles[3], ylab = "",xaxt='n') #plot
axis(1,at=c(1,2,3,4),labels=c("P0", "P4", "P7", "Ad"))
legend("topright", legend = colnames(mean_Cell), col=1:7, pch=1) # optional legend






# 7.2 GO terms analysis

upreg_goterms <-  read.csv("upregulated-cluster.csv",header = F,sep = "\t")
upreg_ref <- read.csv("upreg_paper_reference.csv",header = F)

up_pathway_match <- match(upreg_goterms[,2], upreg_ref[,2])

up_pathway_match_to_aster <- up_pathway_match

up_pathway_match_to_aster[!is.na(up_pathway_match_to_aster)] <- "*"

upreg_goterms <- cbind(upreg_goterms, up_pathway_match_to_aster)

table_upgo <- upreg_goterms[,c(2,11,12,13,14)]

table_upgo <- table_upgo[1:129,]
write.csv(table_upgo,"top_upregulated_GO.csv")



downreg_goterms <-  read.csv("downregulated-cluster.csv",header = F,sep = "\t")
downreg_ref <- read.csv("downreg_paper_reference.csv",header = F)

down_pathway_match <- match(downreg_goterms[,2], downreg_ref[,2])

down_pathway_match_to_aster <- down_pathway_match

down_pathway_match_to_aster[!is.na(down_pathway_match_to_aster)] <- "*"

downreg_goterms <- cbind(downreg_goterms, down_pathway_match_to_aster)

table_downgo <- downreg_goterms[,c(2,11,12,13,14)]

table_downgo <- table_downgo[1:112,]
write.csv(table_downgo,"top_downregulated_GO.csv")




#7.3
library(preprocessCore)



up500 <- tail(ord_up_reg_gene_diff,500)[,"gene"]
up500 <- as.character(up500)
down500 <- down_reg_gene_diff[1:500,"gene"]
down500 <- as.character(down500)
top1000_de <- c(up500,down500)



#top1000_de <- ord_gene_exp[1:1000, "gene"]

P01_FPKM_top1000 <-P01[match(top1000_de, P01$gene),c(5,10)]
P02_FPKM_top1000 <-P02[match(top1000_de, P02$gene),c(5,10)]
P41_FPKM_top1000 <-P41[match(top1000_de, P41$gene),c(5,10)]
P42_FPKM_top1000 <-P42[match(top1000_de, P42$gene),c(5,10)]
P71_FPKM_top1000 <-P71[match(top1000_de, P71$gene),c(5,10)]
P72_FPKM_top1000 <-P72[match(top1000_de, P72$gene),c(5,10)]
PAd1_FPKM_top1000 <-PAd1[match(top1000_de, PAd1$gene),c(5,10)]
PAd2_FPKM_top1000 <-PAd2[match(top1000_de, PAd1$gene),c(5,10)]

FPKM_all_set <- cbind(P01_FPKM_top1000$FPKM,P02_FPKM_top1000$FPKM,P41_FPKM_top1000$FPKM,P42_FPKM_top1000$FPKM,P71_FPKM_top1000$FPKM,P72_FPKM_top1000$FPKM,PAd1_FPKM_top1000$FPKM,PAd1_FPKM_top1000$FPKM)
rownames(FPKM_all_set) <- top1000_de
colnames(FPKM_all_set) <- datas

FPKM_top1000 <- FPKM_all_set[(match(top1000_de, rownames(FPKM_all_set))),]
FPKM_top1000 <- FPKM_top1000[!is.na(FPKM_top1000[,1]),]



ind_fpkm1000 <- FPKM_top1000[,c(2,2,4,4,6,6,8,8)]

clusters <- hclust(as.dist(1-cor(ind_fpkm1000, method="pearson")), method="complete") 
clusters$order <- c(7,8,1,2,3,4,5,6)
clusters$height <- c(0.3, 0.4, 0.6, 0.5, 8.357926e-01, 9.649334e-01, 1.002478e+00)
dend <- as.dendrogram(clusters)
plot(clusters)

hm <- heatmap(FPKM_top1000,margins = c(5,10),main = "Top 1000 DE Genes between P0 and Ad",cex.main=1)



















install.packages("gplots")
install.packages("RColorBrewer")
library(RColorBrewer)
library(gplots)
heatmap.2(FPKM_top1000, trace="none",    # trace可以给每个色块中添加一条线，与行平行或者与列平行。其与色块中心的距离代表了这个值被显示的比例。
          scale="row",    # scale在这里
          symbreaks = TRUE,
          col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(20)),  # color key, 后面详叙
          breaks = seq(-0.5,0.5,0.05),   # 还是color key
          density.info=c("none"),  # 还是color key
          margins=c(10,16),  # 调整热图大小比例
          cexRow = 0.8, cexCol = 1.0,   # 行列名字体大小
          srtCol = 65, offsetCol = -0.5) # 调整列名的字体倾斜45度，距离热图的距离缩小。
