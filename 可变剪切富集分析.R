setwd("C:/Users/yudonglin/Desktop/BBRC/可变剪切/new")
pro_SE<- read.table("./pro_beta_wt/SE.MATS.JC.txt",header = T)
dim(pro_SE)
pro_SE<-unique(pro_SE$geneSymbol)
dim(pro_SE)
length(pro_SE)
pro_A3SS <- read.table("./pro_beta_wt/A3SS.MATS.JC.txt",header = T)
dim(pro_A3SS)
pro_A3SS<-unique(pro_A3SS$geneSymbol)
dim(pro_A3SS)
length(pro_A3SS)
pro_A5SS <- read.table("./pro_beta_wt/A5SS.MATS.JC.txt",header = T)
dim(pro_A5SS)
pro_A5SS<-unique(pro_A5SS$geneSymbol)
length(pro_A5SS)
pro_MXE <- read.table("./pro_beta_wt/MXE.MATS.JC.txt",header = T)
dim(pro_MXE)
pro_MXE<-unique(pro_MXE$geneSymbol)
length(pro_MXE)
pro_RI <- read.table("./pro_beta_wt/RI.MATS.JC.txt",header = T)
dim(pro_RI)
pro_RI<-unique(pro_RI$geneSymbol)
length(pro_RI)
length(pro_SE)
dim(pro_A3SS)
length(pro_A3SS)
length(pro_A5SS)
length(pro_SE)
length(pro_A3SS)
length(pro_A5SS)
length(pro_MXE)
length(pro_RI)
baso_SE<- read.table("./baso_beta_wt/SE.MATS.JC.txt",header = T)
dim(baso_SE)
baso_SE<-unique(baso_SE$geneSymbol)
length(baso_SE)
baso_A3SS <- read.table("./baso_beta_wt/A3SS.MATS.JC.txt",header = T)
dim(baso_A3SS)
baso_A3SS<-unique(baso_A3SS$geneSymbol)
length(baso_A3SS)
baso_A5SS <- read.table("./baso_beta_wt/A5SS.MATS.JC.txt",header = T)
dim(baso_A5SS)
baso_A5SS<-unique(baso_A5SS$geneSymbol)
length(baso_A5SS)
baso_MXE <- read.table("./baso_beta_wt/MXE.MATS.JC.txt",header = T)
dim(baso_MXE)
baso_MXE<-unique(baso_MXE$geneSymbol)
length(baso_MXE)
baso_RI <- read.table("./baso_beta_wt/RI.MATS.JC.txt",header = T)
dim(baso_RI)
baso_RI<-unique(baso_RI$geneSymbol)
length(baso_SE)
length(baso_A3SS)
length(baso_A5SS)
length(baso_MXE)
length(baso_RI)
poly_SE<- read.table("./poly_beta_wt/SE.MATS.JC.txt",header = T)
dim(poly_SE)
poly_SE<-unique(poly_SE$geneSymbol)
length(poly_SE)
poly_A3SS <- read.table("./poly_beta_wt/A3SS.MATS.JC.txt",header = T)
dim(poly_A3SS)
poly_A3SS<-unique(poly_A3SS$geneSymbol)
poly_SE<- read.table("./poly_beta_wt/SE.MATS.JC.txt",header = T)
dim(poly_SE)
poly_SE<-unique(poly_SE$geneSymbol)
length(poly_SE)
poly_A5SS <- read.table("./poly_beta_wt/A5SS.MATS.JC.txt",header = T)
dim(poly_A5SS)
poly_A5SS<-unique(poly_A5SS$geneSymbol)
length(poly_A5SS)
poly_MXE <- read.table("./poly_beta_wt/MXE.MATS.JC.txt",header = T)
dim(poly_MXE)
poly_MXE<-unique(poly_MXE$geneSymbol)
length(poly_MXE)
poly_RI <- read.table("./poly_beta_wt/RI.MATS.JC.txt",header = T)
dim(poly_RI)
poly_RI<-unique(poly_RI$geneSymbol)
length(poly_RI)
ortho_SE<- read.table("./ortho_beta_wt/SE.MATS.JC.txt",header = T)
dim(ortho_SE)
ortho_SE<-unique(ortho_SE$geneSymbol)
length(ortho_SE)
ortho_A3SS <- read.table("./ortho_beta_wt/A3SS.MATS.JC.txt",header = T)
dim(ortho_A3SS)
ortho_A3SS<-unique(ortho_A3SS$geneSymbol)
length(ortho_A3SS)
ortho_A5SS <- read.table("./ortho_beta_wt/A5SS.MATS.JC.txt",header = T)
dim(ortho_A5SS)
ortho_A5SS<-unique(ortho_A5SS$geneSymbol)
length(ortho_A5SS)
ortho_MXE <- read.table("./ortho_beta_wt/MXE.MATS.JC.txt",header = T)
dim(ortho_MXE)
ortho_MXE<-unique(ortho_MXE$geneSymbol)
length(ortho_MXE)
ortho_RI <- read.table("./ortho_beta_wt/RI.MATS.JC.txt",header = T)
dim(ortho_RI)
ortho_RI<-unique(ortho_RI$geneSymbol)


length(pro_SE)
length(pro_A3SS)
length(pro_A5SS)
length(pro_MXE)
length(pro_RI)

length(baso_SE)
length(baso_A3SS)
length(baso_A5SS)
length(baso_MXE)
length(baso_RI)



length(poly_SE)
length(poly_A3SS)
length(poly_A5SS)
length(poly_MXE)
length(poly_RI)

length(ortho_SE)
length(ortho_A3SS)
length(ortho_A5SS)
length(ortho_MXE)
length(ortho_RI)

pdf("可变剪切基因的韦恩图.pdf",width = 10,height = 10)
library (ggvenn)
x<-list("pro SE"=pro_SE,"baso SE"=baso_SE,"poly SE"=poly_SE,"ortho SE"=ortho_SE)
ggvenn(x)
x<-list("pro RI"=pro_RI,"baso RI"=baso_RI,"poly RI"=poly_RI,"ortho RI"=ortho_RI)
ggvenn(x)
x<-list("pro MXE"=pro_MXE,"baso MXE"=baso_MXE,"poly MXE"=poly_MXE,"ortho MXE"=ortho_MXE)
ggvenn(x)
x<-list("pro A5SS"=pro_A5SS,"baso A5SS"=baso_A5SS,"poly A5SS"=poly_A5SS,"ortho A5SS"=ortho_A5SS)
ggvenn(x)
x<-list("pro A3SS"=pro_A3SS,"baso A3SS"=baso_A3SS,"poly A3SS"=poly_A3SS,"ortho A3SS"=ortho_A3SS)
ggvenn(x)
dev.off()
View(x)



intersect<-intersect(merge_tf0,merge_tf1)
#pro特异基因
pro_SE_specific<-setdiff(pro_SE,baso_SE)
pro_SE_specific<-setdiff(pro_SE_specific,poly_SE)
pro_SE_specific<-setdiff(pro_SE_specific,ortho_SE)
head(pro_SE_specific)
write.csv(pro_SE_specific,"pro_SE_specific.csv")
pro_RI_specific<-setdiff(pro_RI,baso_RI)
pro_RI_specific<-setdiff(pro_RI_specific,poly_RI)
pro_RI_specific<-setdiff(pro_RI_specific,ortho_RI)
head(pro_RI_specific)
write.csv(pro_RI_specific,"pro_RI_specific.csv")
pro_A3SS_specific<-setdiff(pro_A3SS,baso_A3SS)
pro_A3SS_specific<-setdiff(pro_A3SS_specific,poly_A3SS)
pro_A3SS_specific<-setdiff(pro_A3SS_specific,ortho_A3SS)
head(pro_A3SS_specific)
write.csv(pro_A3SS_specific,"pro_A3SS_specific.csv")
pro_A5SS_specific<-setdiff(pro_A5SS,baso_A5SS)
pro_A5SS_specific<-setdiff(pro_A5SS_specific,poly_A5SS)
pro_A5SS_specific<-setdiff(pro_A5SS_specific,ortho_A5SS)
head(pro_A5SS_specific)
write.csv(pro_A5SS_specific,"pro_A5SS_specific.csv")
pro_MXE_specific<-setdiff(pro_MXE,baso_MXE)
pro_MXE_specific<-setdiff(pro_MXE_specific,poly_MXE)
pro_MXE_specific<-setdiff(pro_MXE_specific,ortho_MXE)
head(pro_MXE_specific)
write.csv(pro_MXE_specific,"pro_MXE_specific.csv")

#baso特异基因
baso_SE_specific<-setdiff(baso_SE,pro_SE)
baso_SE_specific<-setdiff(baso_SE_specific,poly_SE)
baso_SE_specific<-setdiff(baso_SE_specific,ortho_SE)
head(baso_SE_specific)
write.csv(baso_SE_specific,"baso_SE_specific.csv")
baso_RI_specific<-setdiff(baso_RI,pro_RI)
baso_RI_specific<-setdiff(baso_RI_specific,poly_RI)
baso_RI_specific<-setdiff(baso_RI_specific,ortho_RI)
head(baso_RI_specific)
write.csv(baso_RI_specific,"baso_RI_specific.csv")
baso_A3SS_specific<-setdiff(baso_A3SS,pro_A3SS)
baso_A3SS_specific<-setdiff(baso_A3SS_specific,poly_A3SS)
baso_A3SS_specific<-setdiff(baso_A3SS_specific,ortho_A3SS)
head(baso_A3SS_specific)
write.csv(baso_A3SS_specific,"baso_A3SS_specific.csv")
baso_A5SS_specific<-setdiff(baso_A5SS,pro_A5SS)
baso_A5SS_specific<-setdiff(baso_A5SS_specific,poly_A5SS)
baso_A5SS_specific<-setdiff(baso_A5SS_specific,ortho_A5SS)
head(baso_A5SS_specific)
write.csv(baso_A5SS_specific,"baso_A5SS_specific.csv")
baso_MXE_specific<-setdiff(baso_MXE,pro_MXE)
baso_MXE_specific<-setdiff(baso_MXE_specific,poly_MXE)
baso_MXE_specific<-setdiff(baso_MXE_specific,ortho_MXE)
head(baso_MXE_specific)
write.csv(baso_MXE_specific,"baso_MXE_specific.csv")

#poly特异基因
poly_SE_specific<-setdiff(poly_SE,baso_SE)
poly_SE_specific<-setdiff(poly_SE_specific,pro_SE)
poly_SE_specific<-setdiff(poly_SE_specific,ortho_SE)
head(poly_SE_specific)
write.csv(poly_SE_specific,"poly_SE_specific.csv")
poly_RI_specific<-setdiff(poly_RI,baso_RI)
poly_RI_specific<-setdiff(poly_RI_specific,pro_RI)
poly_RI_specific<-setdiff(poly_RI_specific,ortho_RI)
head(poly_RI_specific)
write.csv(poly_RI_specific,"poly_RI_specific.csv")
poly_A3SS_specific<-setdiff(poly_A3SS,baso_A3SS)
poly_A3SS_specific<-setdiff(poly_A3SS_specific,pro_A3SS)
poly_A3SS_specific<-setdiff(poly_A3SS_specific,ortho_A3SS)
head(poly_A3SS_specific)
write.csv(poly_A3SS_specific,"poly_A3SS_specific.csv")
poly_A5SS_specific<-setdiff(poly_A5SS,baso_A5SS)
poly_A5SS_specific<-setdiff(poly_A5SS_specific,pro_A5SS)
poly_A5SS_specific<-setdiff(poly_A5SS_specific,ortho_A5SS)
head(poly_A5SS_specific)
write.csv(poly_A5SS_specific,"poly_A5SS_specific.csv")
poly_MXE_specific<-setdiff(poly_MXE,baso_MXE)
poly_MXE_specific<-setdiff(poly_MXE_specific,pro_MXE)
poly_MXE_specific<-setdiff(poly_MXE_specific,ortho_MXE)
head(poly_MXE_specific)
write.csv(poly_MXE_specific,"poly_MXE_specific.csv")

#ortho特异基因
ortho_SE_specific<-setdiff(ortho_SE,baso_SE)
ortho_SE_specific<-setdiff(ortho_SE_specific,pro_SE)
ortho_SE_specific<-setdiff(ortho_SE_specific,poly_SE)
head(ortho_SE_specific)
write.csv(ortho_SE_specific,"ortho_SE_specific.csv")
ortho_RI_specific<-setdiff(ortho_RI,baso_RI)
ortho_RI_specific<-setdiff(ortho_RI_specific,pro_RI)
ortho_RI_specific<-setdiff(ortho_RI_specific,poly_RI)
head(ortho_RI_specific)
write.csv(ortho_RI_specific,"ortho_RI_specific.csv")
ortho_A3SS_specific<-setdiff(ortho_A3SS,baso_A3SS)
ortho_A3SS_specific<-setdiff(ortho_A3SS_specific,poly_A3SS)
ortho_A3SS_specific<-setdiff(ortho_A3SS_specific,pro_A3SS)
head(ortho_A3SS_specific)
write.csv(ortho_A3SS_specific,"ortho_A3SS_specific.csv")
ortho_A5SS_specific<-setdiff(ortho_A5SS,baso_A5SS)
ortho_A5SS_specific<-setdiff(ortho_A5SS_specific,poly_A5SS)
ortho_A5SS_specific<-setdiff(ortho_A5SS_specific,pro_A5SS)
head(ortho_A5SS_specific)
write.csv(ortho_A5SS_specific,"ortho_A5SS_specific.csv")
ortho_MXE_specific<-setdiff(ortho_MXE,baso_MXE)
ortho_MXE_specific<-setdiff(ortho_MXE_specific,poly_MXE)
ortho_MXE_specific<-setdiff(ortho_MXE_specific,pro_MXE)
head(ortho_MXE_specific)
write.csv(ortho_MXE_specific,"ortho_MXE_specific.csv")

dir.create("Batch_Enrichment")
setwd("./Batch_Enrichment")


a <- read.csv("pro_SE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_SE_specific.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_SE_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_SE_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_SE_specific.csv")
pdf("KEGG_pro_SE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("pro_A5SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_A5SS_specific.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_A5SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_A5SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_A5SS_specific.csv")
pdf("KEGG_pro_A5SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("pro_A3SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_A3SS_specific.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_A3SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_A3SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_A3SS_specific.csv")
pdf("KEGG_pro_A3SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("pro_MXE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_MXE.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_MXE.csv")
pdf("KEGG_pro_MXE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("pro_RI_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_RI_specific.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_RI_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_RI_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_RI_specific.csv")
pdf("KEGG_pro_RI_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()





a <- read.csv("baso_SE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_SE_specific.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_SE_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_SE_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_SE_specific.csv")
pdf("KEGG_baso_SE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("baso_A5SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_A5SS_specific.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_A5SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_A5SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_A5SS_specific.csv")
pdf("KEGG_baso_A5SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("baso_A3SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_A3SS_specific.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_A3SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_A3SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_A3SS_specific.csv")
pdf("KEGG_baso_A3SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("baso_MXE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_MXE.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_MXE.csv")
pdf("KEGG_baso_MXE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("baso_RI_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_RI_specific.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_RI_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_RI_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_RI_specific.csv")
pdf("KEGG_baso_RI_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()





a <- read.csv("poly_SE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_SE_specific.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_SE_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_SE_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_SE_specific.csv")
pdf("KEGG_poly_SE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("poly_A5SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_A5SS_specific.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_A5SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_A5SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_A5SS_specific.csv")
pdf("KEGG_poly_A5SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("poly_A3SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_A3SS_specific.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_A3SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_A3SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_A3SS_specific.csv")
pdf("KEGG_poly_A3SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("poly_MXE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_MXE.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_MXE.csv")
pdf("KEGG_poly_MXE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("poly_RI_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_RI_specific.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_RI_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_RI_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_RI_specific.csv")
pdf("KEGG_poly_RI_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()





a <- read.csv("ortho_SE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_SE_specific.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_SE_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_SE_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_SE_specific.csv")
pdf("KEGG_ortho_SE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_SE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("ortho_A5SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_A5SS_specific.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_A5SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_A5SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_A5SS_specific.csv")
pdf("KEGG_ortho_A5SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_A5SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("ortho_A3SS_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_A3SS_specific.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_A3SS_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_A3SS_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_A3SS_specific.csv")
pdf("KEGG_ortho_A3SS_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_A3SS_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("ortho_MXE_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_MXE.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_MXE.csv")
pdf("KEGG_ortho_MXE_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_MXE_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.csv("ortho_RI_specific.csv",header = T)
b<-unique(a$x)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_RI_specific.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_RI_specific.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_RI_specific.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_RI_specific.csv")
pdf("KEGG_ortho_RI_specific.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_RI_specific.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()








#加载包
library(reshape2)#融合数据
library(ggplot2)#绘图工具
#读入数据，此例中数据的格式为xlsx
mydata <- read.csv("number.csv")
mydata <- melt(mydata,id="phase")
colnames(mydata) <- c("phase","AS_type","value")#更改列名
mydata$phase<-factor(mydata$phase,levels = c("pro","baso","poly","ortho"))
ggplot(data = mydata,aes(x=phase,y=value,group = AS_type,color=AS_type,shape=AS_type))+
  geom_point()+
  geom_line()+
  xlab("phase")+#横坐标名称
  ylab("AS gene number")+#纵坐标名称
  theme_bw() #去掉背景灰色
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
#         text = element_text(family = "STXihei"),#设置中文字体的显示
#         legend.position = c(.075,.915),#更改图例的位置，放至图内部的左上角
#         legend.box.background = element_rect(color="black"))+#为图例田间边框线
#   scale_x_continuous(limits = c(2000,2018),breaks = seq(2000,2018,1))#更改横坐标刻度值
# #点击zoom查看大图
#https://www.freesion.com/article/10101402279/
ggsave("number.pdf",width = 5,height = 5)


#加载包
library(reshape2)#融合数据
library(ggplot2)#绘图工具
#读入数据，此例中数据的格式为xlsx
mydata <- read.csv("event.csv")
mydata <- melt(mydata,id="phase")
colnames(mydata) <- c("phase","AS_type","value")#更改列名
mydata$phase<-factor(mydata$phase,levels = c("pro","baso","poly","ortho"))
ggplot(data = mydata,aes(x=phase,y=value,group = AS_type,color=AS_type,shape=AS_type))+
  geom_point()+
  geom_line()+
  xlab("phase")+#横坐标名称
  ylab("AS event number")+#纵坐标名称
  theme_bw() #去掉背景灰色
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
#         text = element_text(family = "STXihei"),#设置中文字体的显示
#         legend.position = c(.075,.915),#更改图例的位置，放至图内部的左上角
#         legend.box.background = element_rect(color="black"))+#为图例田间边框线
#   scale_x_continuous(limits = c(2000,2018),breaks = seq(2000,2018,1))#更改横坐标刻度值
# #点击zoom查看大图
#https://www.freesion.com/article/10101402279/
ggsave("event.pdf",width = 5,height = 5)






a <- read.table("./pro_beta_wt/SE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_SE.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_SE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_SE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_SE.csv")
pdf("KEGG_pro_SE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./pro_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_A5SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_A5SS.csv")
pdf("KEGG_pro_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.table("./pro_beta_wt/A3SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_A3SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_A3SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_A3SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_A3SS.csv")
pdf("KEGG_pro_A3SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./pro_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_A5SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_A5SS.csv")
pdf("KEGG_pro_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./pro_beta_wt/MXE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_MXE.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_MXE.csv")
pdf("KEGG_pro_MXE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./pro_beta_wt/RI.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_pro_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_pro_RI.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_pro_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_pro_RI.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_pro_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_pro_RI.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_pro_RI.csv")
pdf("KEGG_pro_RI.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_pro_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()









a <- read.table("./baso_beta_wt/SE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_SE.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_SE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_SE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_SE.csv")
pdf("KEGG_baso_SE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./baso_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_A5SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_A5SS.csv")
pdf("KEGG_baso_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.table("./baso_beta_wt/A3SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_A3SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_A3SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_A3SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_A3SS.csv")
pdf("KEGG_baso_A3SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./baso_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_A5SS.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_A5SS.csv")
pdf("KEGG_baso_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./baso_beta_wt/MXE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_MXE.csv")
#生物过程(biologicalbasocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_MXE.csv")
pdf("KEGG_baso_MXE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./baso_beta_wt/RI.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_baso_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_baso_RI.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_baso_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_baso_RI.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_baso_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_baso_RI.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_baso_RI.csv")
pdf("KEGG_baso_RI.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_baso_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()






a <- read.table("./poly_beta_wt/SE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_SE.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_SE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_SE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_SE.csv")
pdf("KEGG_poly_SE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./poly_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_A5SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_A5SS.csv")
pdf("KEGG_poly_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.table("./poly_beta_wt/A3SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_A3SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_A3SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_A3SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_A3SS.csv")
pdf("KEGG_poly_A3SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./poly_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_A5SS.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_A5SS.csv")
pdf("KEGG_poly_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./poly_beta_wt/MXE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_MXE.csv")
#生物过程(biologicalpolycess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_MXE.csv")
pdf("KEGG_poly_MXE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./poly_beta_wt/RI.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_poly_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_poly_RI.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_poly_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_poly_RI.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_poly_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_poly_RI.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_poly_RI.csv")
pdf("KEGG_poly_RI.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_poly_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()








a <- read.table("./ortho_beta_wt/SE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_SE.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_SE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_SE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_SE.csv")
pdf("KEGG_ortho_SE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_SE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./ortho_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_A5SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_A5SS.csv")
pdf("KEGG_ortho_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

a <- read.table("./ortho_beta_wt/A3SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_A3SS.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_A3SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_A3SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_A3SS.csv")
pdf("KEGG_ortho_A3SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_A3SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./ortho_beta_wt/A5SS.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_A5SS.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_A5SS.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_A5SS.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_A5SS.csv")
pdf("KEGG_ortho_A5SS.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_A5SS.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./ortho_beta_wt/MXE.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_MXE.csv")
#生物过程(biologicalorthocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_MXE.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_MXE.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_MXE.csv")
pdf("KEGG_ortho_MXE.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_MXE.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
a <- read.table("./ortho_beta_wt/RI.MATS.JC.txt",header = T)
b<-unique(a$geneSymbol)
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_MF_ortho_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_MF_ortho_RI.csv")
#生物过程(biologicalprocess)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_BP_ortho_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_BP_ortho_RI.csv")
#细胞组成(cellularcomponent)
ego <- enrichGO(
  gene          = gene,
  keyType = "ENTREZID",
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
pdf("GO_CC_ortho_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
dev.off()
write.csv(ego,file="GO_CC_ortho_RI.csv")
ekegg <- enrichKEGG(
  gene          = gene,
  keyType     = "kegg",
  organism   = "mmu",
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
write.csv(ekegg,file="kegg_ortho_RI.csv")
pdf("KEGG_ortho_RI.pdf",width=20,height=20)
barplot(ekegg, showCategory =50)
dotplot(ekegg, showCategory =50)
dev.off()
ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
pdf("GO_ortho_RI.pdf",width=20,height=20)
barplot(ego, showCategory =50)
dotplot(ego, showCategory =50)
#柱状图,各自展示前8个term
barplot(ego, showCategory =50,drop = TRUE,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
#气泡图,各自展示前8个term
dotplot(ego, showCategory =50,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()



Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)


#BP，以热图展示
#合并SE对应的富集分析结果
marker0<-read.csv("GO_BP_pro_SE.csv")
marker0_data<-marker0[,c("ID","Description","GeneRatio")]

marker0_data$ID<-paste(marker0_data$ID,marker0_data$Description,sep="_")
marker0_data<-marker0_data[,c("ID","GeneRatio")]
marker0_data$GeneRatio<-sapply(marker0_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker0_data)
dim(marker0_data)

marker1<-read.csv("GO_BP_baso_SE.csv")
marker1_data<-marker1[,c("ID","Description","GeneRatio")]

marker1_data$ID<-paste(marker1_data$ID,marker1_data$Description,sep="_")
marker1_data<-marker1_data[,c("ID","GeneRatio")]
marker1_data$GeneRatio<-sapply(marker1_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker1_data)
dim(marker1_data)

marker2<-read.csv("GO_BP_poly_SE.csv")
marker2_data<-marker2[,c("ID","Description","GeneRatio")]

marker2_data$ID<-paste(marker2_data$ID,marker2_data$Description,sep="_")
marker2_data<-marker2_data[,c("ID","GeneRatio")]
marker2_data$GeneRatio<-sapply(marker2_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker2_data)
dim(marker2_data)

marker3<-read.csv("GO_BP_ortho_SE.csv")
marker3_data<-marker3[,c("ID","Description","GeneRatio")]

marker3_data$ID<-paste(marker3_data$ID,marker3_data$Description,sep="_")
marker3_data<-marker3_data[,c("ID","GeneRatio")]
marker3_data$GeneRatio<-sapply(marker3_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker3_data)
dim(marker3_data)


marker4<-read.csv("GO_BP_pro_RI.csv")
marker4_data<-marker4[,c("ID","Description","GeneRatio")]

marker4_data$ID<-paste(marker4_data$ID,marker4_data$Description,sep="_")
marker4_data<-marker4_data[,c("ID","GeneRatio")]
marker4_data$GeneRatio<-sapply(marker4_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker4_data)
dim(marker4_data)


marker5<-read.csv("GO_BP_baso_RI.csv")
marker5_data<-marker5[,c("ID","Description","GeneRatio")]

marker5_data$ID<-paste(marker5_data$ID,marker5_data$Description,sep="_")
marker5_data<-marker5_data[,c("ID","GeneRatio")]
marker5_data$GeneRatio<-sapply(marker5_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker5_data)
dim(marker5_data)


marker6<-read.csv("GO_BP_poly_RI.csv")
marker6_data<-marker6[,c("ID","Description","GeneRatio")]

marker6_data$ID<-paste(marker6_data$ID,marker6_data$Description,sep="_")
marker6_data<-marker6_data[,c("ID","GeneRatio")]
marker6_data$GeneRatio<-sapply(marker6_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker6_data)
dim(marker6_data)

marker7<-read.csv("GO_BP_ortho_RI.csv")
marker7_data<-marker7[,c("ID","Description","GeneRatio")]

marker7_data$ID<-paste(marker7_data$ID,marker7_data$Description,sep="_")
marker7_data<-marker7_data[,c("ID","GeneRatio")]
marker7_data$GeneRatio<-sapply(marker7_data$GeneRatio, function(x) eval(parse(text=x)))
head(marker7_data)
dim(marker7_data)





#合并不同的矩阵，列名弄成一样长，然后补上0值。

ID<-c(marker0_data$ID,marker1_data$ID,marker2_data$ID,marker3_data$ID,marker4_data$ID,marker5_data$ID,marker6_data$ID,marker7_data$ID)
name<-as.data.frame(ID)
name<-name[!duplicated(name$ID),] #删掉所有列上都重复的
name<-as.data.frame(name)
colnames(name)<-c("ID")
head(name)



DATA_cluster0<-dplyr::bind_rows(name,marker0_data)
DATA_cluster0<-DATA_cluster0[order(DATA_cluster0[,2],decreasing=T),]#以第一列降序排列
DATA_cluster0<-DATA_cluster0[!duplicated(DATA_cluster0$ID),] #删掉所有列上都重复的
DATA_cluster0[is.na(DATA_cluster0)] <- 0#把NA值全部替换为0
colnames(DATA_cluster0)<-c("ID","pro_SE")
head(DATA_cluster0)


DATA_cluster1<-dplyr::bind_rows(name,marker1_data)
DATA_cluster1<-DATA_cluster1[order(DATA_cluster1[,2],decreasing=T),]#以第一列降序排列
DATA_cluster1<-DATA_cluster1[!duplicated(DATA_cluster1$ID),] #删掉所有列上都重复的
DATA_cluster1[is.na(DATA_cluster1)] <- 0#把NA值全部替换为0
colnames(DATA_cluster1)<-c("ID","baso_SE")
head(DATA_cluster1)


DATA_cluster2<-dplyr::bind_rows(name,marker2_data)
DATA_cluster2<-DATA_cluster2[order(DATA_cluster2[,2],decreasing=T),]#以第一列降序排列
DATA_cluster2<-DATA_cluster2[!duplicated(DATA_cluster2$ID),] #删掉所有列上都重复的
DATA_cluster2[is.na(DATA_cluster2)] <- 0#把NA值全部替换为0
colnames(DATA_cluster2)<-c("ID","poly_SE")
head(DATA_cluster2)


DATA_cluster3<-dplyr::bind_rows(name,marker3_data)
DATA_cluster3<-DATA_cluster3[order(DATA_cluster3[,2],decreasing=T),]#以第一列降序排列
DATA_cluster3<-DATA_cluster3[!duplicated(DATA_cluster3$ID),] #删掉所有列上都重复的
DATA_cluster3[is.na(DATA_cluster3)] <- 0#把NA值全部替换为0
colnames(DATA_cluster3)<-c("ID","ortho_SE")
head(DATA_cluster3)


DATA_cluster4<-dplyr::bind_rows(name,marker4_data)
DATA_cluster4<-DATA_cluster4[order(DATA_cluster4[,2],decreasing=T),]#以第一列降序排列
DATA_cluster4<-DATA_cluster4[!duplicated(DATA_cluster4$ID),] #删掉所有列上都重复的
DATA_cluster4[is.na(DATA_cluster4)] <- 0#把NA值全部替换为0
colnames(DATA_cluster4)<-c("ID","pro_RI")
head(DATA_cluster4)


DATA_cluster5<-dplyr::bind_rows(name,marker5_data)
DATA_cluster5<-DATA_cluster5[order(DATA_cluster5[,2],decreasing=T),]#以第一列降序排列
DATA_cluster5<-DATA_cluster5[!duplicated(DATA_cluster5$ID),] #删掉所有列上都重复的
DATA_cluster5[is.na(DATA_cluster5)] <- 0#把NA值全部替换为0
colnames(DATA_cluster5)<-c("ID","baso_RI")
head(DATA_cluster5)


DATA_cluster6<-dplyr::bind_rows(name,marker6_data)
DATA_cluster6<-DATA_cluster6[order(DATA_cluster6[,2],decreasing=T),]#以第一列降序排列
DATA_cluster6<-DATA_cluster6[!duplicated(DATA_cluster6$ID),] #删掉所有列上都重复的
DATA_cluster6[is.na(DATA_cluster6)] <- 0#把NA值全部替换为0
colnames(DATA_cluster6)<-c("ID","poly_RI")
head(DATA_cluster6)


DATA_cluster7<-dplyr::bind_rows(name,marker7_data)
DATA_cluster7<-DATA_cluster7[order(DATA_cluster7[,2],decreasing=T),]#以第一列降序排列
DATA_cluster7<-DATA_cluster7[!duplicated(DATA_cluster7$ID),] #删掉所有列上都重复的
DATA_cluster7[is.na(DATA_cluster7)] <- 0#把NA值全部替换为0
colnames(DATA_cluster7)<-c("ID","ortho_RI")
head(DATA_cluster7)



cluster01<-merge(DATA_cluster0,DATA_cluster1,by=c("ID"))
cluster23<-merge(DATA_cluster2,DATA_cluster3,by=c("ID"))
cluster45<-merge(DATA_cluster4,DATA_cluster5,by=c("ID"))
cluster67<-merge(DATA_cluster6,DATA_cluster7,by=c("ID"))

cluster0123<-merge(cluster01,cluster23,by=c("ID"))
cluster4567<-merge(cluster45,cluster67,by=c("ID"))

DATA<-merge(cluster0123,cluster4567,by=c("ID"))
# DATA<-DATA[,c("ID","cluster0","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7")]
# 
# dim(DATA)
# head(DATA)
# DATA<-DATA[!duplicated(DATA$ID),] #删掉所有列上都重复的
rownames(DATA)<-DATA[,1]
DATA<-DATA[,-1]
DATA<-as.data.frame(DATA)
# https://blog.csdn.net/c1z2w3456789/article/details/79467095

head(DATA)

write.csv(DATA,"GO_SE_RI.csv")



pheatmap::pheatmap(DATA,
                   cluster_rows = F,
                   cluster_cols = F,show_rownames = T)


library(pheatmap)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(viridis)
library(cowplot)


pheatmap(
  mat   = DATA,
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = F,
  drop_levels   = TRUE,
  fontsize  = 8,
  main = "GO"
)

pheatmap(DATA,show_rownames = F,
         cluster_rows = T,border_color = NA,
         cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))



