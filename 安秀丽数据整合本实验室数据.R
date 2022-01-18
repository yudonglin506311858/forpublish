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
library(pheatmap)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
setwd("D:/anxiiuli_data/count_mouse_new")


#读入数据并合并
#安秀丽第一批wt数据：小鼠
WT_pro_1<-read.table("mm_proerythroblast_1.count",header=T)
WT_baso_1<-read.table("mm_basophilic_1.count",header=T)
WT_poly_1<-read.table("mm_polychromatic_1.count",header=T)
WT_ortho_1<-read.table("mm_orthochromatic_1.count",header=T)

#安秀丽第二批wt数据（新测序）
WT_pro_2<-read.table("mm_proerythroblast_2.count",header=T)
WT_baso_2<-read.table("mm_basophilic_2.count",header=T)
WT_poly_2<-read.table("mm_polychromatic_2.count",header=T)
WT_ortho_2<-read.table("mm_orthochromatic_2.count",header=T)

WT_pro_3<-read.table("mm_proerythroblast_3.count",header=T)
WT_baso_3<-read.table("mm_basophilic_3.count",header=T)
WT_poly_3<-read.table("mm_polychromatic_3.count",header=T)
WT_ortho_3<-read.table("mm_orthochromatic_3.count",header=T)



#合并pro的数据
pro<-merge(WT_pro_1,WT_pro_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro<-merge(pro,WT_pro_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并baso的数据
baso<-merge(WT_baso_1,WT_baso_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
baso<-merge(baso,WT_baso_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并poly的数据
poly<-merge(WT_poly_1,WT_poly_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly<-merge(poly,WT_poly_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并ortho的数据
ortho<-merge(WT_ortho_1,WT_ortho_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
ortho<-merge(ortho,WT_ortho_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
##合并五个时期的数据
pro_baso<-merge(pro,baso,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly_ortho<-merge(poly,ortho,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro_baso_poly_ortho<-merge(pro_baso,poly_ortho,by=c("Geneid","Chr","Start","End","Strand","Length"))


WT_count<-pro_baso_poly_ortho

head(WT_count)
dim(WT_count)
colnames(WT_count)


#删除ensmusg的版本号

#WT_count$Geneid<-gsub("\\.*","",WT_count$Geneid)
#WT_count$Geneid <- gsub("\\.[0-9]*$", "", WT_count$Geneid)
rownames(WT_count)<-WT_count$Geneid
#https://www.nhooo.com/note/qa02b6.html
#https://www.biostars.org/p/178726/
WT_count_all<-WT_count[,c(1,6:18)]
head(WT_count_all)
colnames(WT_count_all)<-c("Geneid","Length","PRO_1","PRO_2","PRO_3",
                          "BASO_1","BASO_2","BASO_3",
                          "POLY_1","POLY_2","POLY_3",
                          "ORTHO_1","ORTHO_2","ORTHO_3")

cts<-WT_count_all
head(cts)
dim(cts)

write.csv(cts,"ensembl_gene_expression_count_anxiuli.csv")



#基因ID转换
library('biomaRt')
library("curl")
library(ensembldb)
library(dplyr)
library(AnnotationHub)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
saveRDS(mart,"mart_mouse.rds")
mart<-readRDS("mart_mouse.rds")

gene<-read.csv("ensembl_gene_expression_count_anxiuli.csv")
gene<-as.matrix(gene$X)
head(gene)
colnames(gene)[1]<-"ensembl_gene_id"
#listAttributes(mart)
id_con<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = gene, mart = mart)
head(id_con)
write.csv(id_con,"mouse_gene_ensembl_transition.csv")






library(stringr)
#cts$Geneid<-str_sub(cts$Geneid,1,str_locate(cts$Geneid,"\\.")[1]-1)
id_con<-read.csv("mouse_gene_ensembl_transition.csv")
id_con<-id_con[,c(2,3)]
colnames(cts)[1]<-"ensembl_gene_id"
head(cts)
head(id_con)
cts<-merge(id_con,cts,by=c("ensembl_gene_id"))
dim(cts)
write.csv(cts,file = "ternimalE_count_genesymbol.csv",row.names = F)
cts<-cts[,-1]
cts$external_gene_name<-make.names(cts$external_gene_name, unique = TRUE)
rownames(cts)<-cts$external_gene_name
write.table(cts,"data.txt", sep = "\t")



#读取count文件：
write.csv(cts,file="ternimalE_count.csv")

#读取tpm文件：
head(cts)
kb <- cts$Length / 1000
head(kb)
countdata <- cts[,3:14]
rpk <- countdata / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm,file="ternimalE_tpm.csv")

#读取fpkm文件：
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm,file="ternimalE_fpkm.csv")


#PCA
condition<-factor(c("PRO","PRO","PRO",
                    "BASO","BASO","BASO",
                    "POLY","POLY","POLY",
                    "ORTHO","ORTHO","ORTHO"),
                  levels = c("PRO","BASO","POLY","ORTHO"))

head(cts)
tmp<-cts[,3:14]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:14]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:14],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")

#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


dir.create("./merge")
setwd("D:/anxiiuli_data/count_mouse_new/merge")


anxiuli<-cts
head(anxiuli)

#整合本实验室数据
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
setwd("D:/paper/bulkRNA-seq/count_own")


#读入数据并合并
#师姐第一批wt数据
WT_pro_1<-read.table("XH-Z2wt-pro-19-11-26_FKDL192547546-1a-33.count",header=T)
WT_baso_1<-read.table("XH-Z2wt-Baso-19-11-26_FKDL192547546-1a-34.count",header=T)
WT_poly_1<-read.table("XH-Z2wt-poly-19-11-26_FKDL192547546-1a-35.count",header=T)
WT_ortho_1<-read.table("XH-Z2wt-ortho-19-11-26_FKDL192547546-1a-36.count",header=T)
WT_pro_2<-read.table("XH-Z2wt-pro-19-11-27_FKDL192547546-1a-41.count",header=T)
WT_baso_2<-read.table("XH-Z2wt-Baso-19-11-27_FKDL192547546-1a-42.count",header=T)
WT_poly_2<-read.table("XH-Z2wt-poly-19-11-27_FKDL192547546-1a-43.count",header=T)
WT_ortho_2<-read.table("XH-Z2wt-ortho-19-11-27_FKDL192547546-1a-44.count",header=T)

#师姐第二批wt数据（新测序）
WT_pro_3<-read.table("wt-pro-19-12-24_FKDL202555457-1a-A1.count",header=T)
WT_baso_3<-read.table("wt-Baso-19-12-24_FKDL202555457-1a-A2.count",header=T)
WT_poly_3<-read.table("wt-poly-19-12-24_FKDL202555457-1a-A3.count",header=T)
WT_ortho_3<-read.table("wt-ortho-19-12-24_FKDL202555457-1a-A4.count",header=T)

#师兄的两个wt数据
WT_pro_4<-read.table("XH-PRO-1_HL3NLCCXY_L3.count",header=T)
WT_pro_5<-read.table("XH-PRO-2_HL3NLCCXY_L3.count",header=T)
WT_baso_4<-read.table("XH-BASO-1_HL3NLCCXY_L3.count",header=T)
WT_baso_5<-read.table("XH-BASO-2_HL3NLCCXY_L3.count",header=T)
WT_poly_4<-read.table("XH-POLY-1_HL3NLCCXY_L3.count",header=T)
WT_poly_5<-read.table("XH-POLY-2_HL3NLCCXY_L3.count",header=T)
WT_ortho_4<-read.table("XH-ORTHO-1_HL3NLCCXY_L3.count",header=T)
WT_ortho_5<-read.table("XH-ORTHO-2_HL3NLCCXY_L3.count",header=T)

#合并pro的数据
pro<-merge(WT_pro_1,WT_pro_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro<-merge(pro,WT_pro_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro<-merge(pro,WT_pro_4,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro<-merge(pro,WT_pro_5,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并baso的数据
baso<-merge(WT_baso_1,WT_baso_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
baso<-merge(baso,WT_baso_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
baso<-merge(baso,WT_baso_4,by=c("Geneid","Chr","Start","End","Strand","Length"))
baso<-merge(baso,WT_baso_5,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并poly的数据
poly<-merge(WT_poly_1,WT_poly_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly<-merge(poly,WT_poly_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly<-merge(poly,WT_poly_4,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly<-merge(poly,WT_poly_5,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并ortho的数据
ortho<-merge(WT_ortho_1,WT_ortho_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
ortho<-merge(ortho,WT_ortho_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
ortho<-merge(ortho,WT_ortho_4,by=c("Geneid","Chr","Start","End","Strand","Length"))
ortho<-merge(ortho,WT_ortho_5,by=c("Geneid","Chr","Start","End","Strand","Length"))
##合并四个时期的数据
pro_baso<-merge(pro,baso,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly_ortho<-merge(poly,ortho,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro_baso_poly_ortho<-merge(pro_baso,poly_ortho,by=c("Geneid","Chr","Start","End","Strand","Length"))


WT_count<-pro_baso_poly_ortho

head(WT_count)
dim(WT_count)
colnames(WT_count)

rownames(WT_count)<-WT_count$Geneid
WT_count_all<-WT_count[,c(1,6:26)]
head(WT_count_all)
colnames(WT_count_all)<-c("Geneid","Length","PRO_1","PRO_2","PRO_3","PRO_4","PRO_5",
                          "BASO_1","BASO_2","BASO_3","BASO_4","BASO_5",
                          "POLY_1","POLY_2","POLY_3","POLY_4","POLY_5",
                          "ORTHO_1","ORTHO_2","ORTHO_3","ORTHO_4","ORTHO_5")


cts<-WT_count_all
head(cts)
dim(cts)

write.csv(cts,"ensembl_gene_expression_count_ydl.csv")


#基因ID转换
library('biomaRt')
library("curl")
library(ensembldb)
library(dplyr)
library(AnnotationHub)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
saveRDS(mart,"mart_mouse.rds")
mart<-readRDS("mart_mouse.rds")

gene<-read.csv("ensembl_gene_expression_count_anxiuli.csv")
gene<-as.matrix(gene$X)
head(gene)
colnames(gene)[1]<-"ensembl_gene_id"
#listAttributes(mart)
id_con<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = gene, mart = mart)
head(id_con)
write.csv(id_con,"mouse_gene_ensembl_transition.csv")








library(stringr)
#cts$Geneid<-str_sub(cts$Geneid,1,str_locate(cts$Geneid,"\\.")[1]-1)
id_con<-read.csv("mouse_gene_ensembl_transition.csv")
id_con<-id_con[,c(2,3)]
colnames(cts)[1]<-"ensembl_gene_id"
head(cts)
head(id_con)
cts<-merge(id_con,cts,by=c("ensembl_gene_id"))
dim(cts)
write.csv(cts,file = "ternimalE_count_genesymbol.csv",row.names = F)
cts<-cts[,-1]
cts$external_gene_name<-make.names(cts$external_gene_name, unique = TRUE)
rownames(cts)<-cts$external_gene_name
write.table(cts,"data.txt", sep = "\t")



#读取count文件：
write.csv(cts,file="ternimalE_count.csv")

#读取tpm文件：
head(cts)
kb <- cts$Length / 1000
head(kb)
#countdata <- cts[,3:34]
countdata <- cts[,3:22]
rpk <- countdata / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm,file="ternimalE_tpm.csv")

#读取fpkm文件：
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm,file="ternimalE_fpkm.csv")


#PCA
condition<-factor(c("PRO","PRO","PRO","PRO","PRO",
                    "BASO","BASO","BASO","BASO","BASO",
                    "POLY","POLY","POLY","POLY","POLY",
                    "ORTHO","ORTHO","ORTHO","ORTHO","ORTHO"),
                  levels = c("PRO","BASO","POLY","ORTHO"))

head(cts)
tmp<-cts[,3:22]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:22]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:22],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")



#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



anxiuli<-anxiuli[,-1]

colnames(anxiuli)<-c("external_gene_name","Length","PRO_6","PRO_7","PRO_8",
                     "BASO_6","BASO_7","BASO_8",
                     "POLY_6","POLY_7","POLY_8",
                     "ORTHO_6","ORTHO_7","ORTHO_8")

head(anxiuli)
head(cts)

data<-cbind(cts,anxiuli)

head(data,3)
dim(data)

new_data<-data[c(1,2,3,4,5,6,7,25,26,27,8,9,10,11,12,28,29,30,13,14,15,16,17,31,32,33,18,19,20,21,22,34,35,36)]
head(new_data)
cts<-new_data




#PCA
condition<-factor(c("PRO","PRO","PRO","PRO","PRO","PRO","PRO","PRO",
                    "BASO","BASO","BASO","BASO","BASO","BASO","BASO","BASO",
                    "POLY","POLY","POLY","POLY","POLY","POLY","POLY","POLY",
                    "ORTHO","ORTHO","ORTHO","ORTHO","ORTHO","ORTHO","ORTHO","ORTHO"),
                  levels = c("PRO","BASO","POLY","ORTHO"))

head(cts)
tmp<-cts[,3:34]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:34]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:34],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")



#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)










top100<-read.table("top100.txt",header = T)

library (ggvenn)
x<-list("pro"=top100$pro,
        "baso"=top100$baso,
        "poly"=top100$poly,
        "ortho"=top100$ortho)
ggvenn(x)
pdf("GO term的交集和补集.pdf")
ggvenn(x)
dev.off()


intersect<-intersect(top100$baso,top100$poly)



library(clusterProfiler)
library(org.Mm.eg.db)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
## 加载R包
library("Mfuzz")
rm(list=ls())
getwd()


normal<-cts[,-c(1,2)]
normal$pro<-apply(normal[,1:8],1, mean, na.rm = T) 
normal$baso<-apply(normal[,9:16],1, mean, na.rm = T) 
normal$poly<-apply(normal[,17:24],1, mean, na.rm = T) 
normal$ortho<-apply(normal[,25:32],1, mean, na.rm = T) 

normal <-normal[,33:36]

write.csv(normal,"红系终末分化四个时期的平均值.csv")
normal <- data.matrix(normal)
eset <- new("ExpressionSet",exprs = normal)
## 过滤缺失超过25%的基因
gene.r <- filter.NA(eset, thres=0.25)
## mean填补缺失
gene.f <- fill.NA(gene.r,mode="mean")
## knn/wknn方法表现更好，但是计算起来比较复杂
gene.f <- fill.NA(gene.r,mode="knn")
gene.f <- fill.NA(gene.r,mode="wknn")
## 过滤标准差为0的基因
tmp <- filter.std(gene.f,min.std=0)
## 标准化
gene.s <- standardise(tmp)
## 聚类个数
c <- 9
## 计算最佳的m值
m <- mestimate(gene.s)
## 聚类
cl <- mfuzz(gene.s, c = c, m = m)
## 查看每类基因数目
cl$size
## 查看每类基因ID
cl$cluster[cl$cluster == 1]
## 输出基因ID
write.table(cl$cluster,"output.txt",quote=F,row.names=T,col.names=F,sep="\t")
## 绘制折线图
pdf("try.pdf")
mfuzz.plot(gene.s,cl,mfrow=c(3,3),new.window= FALSE)
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


data<-cts[,-c(1,2)]
#data<-read.table("data.txt", sep = "\t",check.names = FALSE)
data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况

YDL <- CreateSeuratObject(counts = data,project = "ternimalE")#min.cells=10,min.genes=200,

Idents(YDL)<-YDL@meta.data$orig.ident
YDL@meta.data$orig.ident<- factor(YDL@meta.data$orig.ident,levels=c("PRO","BASO","POLY","ORTHO"))



##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中，mit-开头的为线粒体基因
YDL[["percent.mt"]] <- PercentageFeatureSet(object = YDL, pattern = "^mt-")
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


head(YDL@meta.data)
summary(YDL@meta.data$nCount_RNA)
summary(YDL@meta.data$nFeature_RNA)
summary(YDL@meta.data$percent.mt)




#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 2,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL)#选择top20个PC
pcSelect=2
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution =1)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)
DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5)
# YDL@meta.data$celltype<-c("PRO","PRO","PRO","PRO","PRO",
#                           "BASO","BASO","BASO","BASO","BASO",
#                           "POLY","POLY","POLY","POLY","POLY",
#                           "ORTHO","ORTHO","ORTHO","ORTHO","ORTHO")
YDL@meta.data$celltype<-factor(c("PRO","PRO","PRO","PRO","PRO","PRO","PRO","PRO",
                                 "BASO","BASO","BASO","BASO","BASO","BASO","BASO","BASO",
                                 "POLY","POLY","POLY","POLY","POLY","POLY","POLY","POLY",
                                 "ORTHO","ORTHO","ORTHO","ORTHO","ORTHO","ORTHO","ORTHO","ORTHO"),
                               levels = c("PRO","BASO","POLY","ORTHO"))

rownames(YDL@meta.data)
DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5,group.by = "celltype")
DimPlot(YDL,reduction = "pca",label = F,pt.size = 1.5,group.by = "celltype")

DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5)
saveRDS(YDL,"mm_ternimalE.rds")


##细胞周期归类

YDL<- CellCycleScoring(object = YDL, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = YDL@meta.data)
pdf(file="CellCycle_proe.pdf",width=6.5,height=6)
DimPlot(YDL,reduction = "pca",label = TRUE,group.by="Phase",pt.size = 1.5)
DimPlot(YDL,reduction = "pca",label = TRUE,group.by="Phase",pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)
dev.off()


library(rSuperCT)
library(ggplot2)
pred_obj <- ImportData(YDL)
dir.create('./models', showWarnings = FALSE)
pred_obj <- PredCellTypes(pred_obj, species = 'mouse', model = 'generic_38celltypes',results.dir = './models')
table(pred_obj@meta.data$pred_types)# 举例来展示不同细胞的比例
g <- plotHist(pred_obj) + scale_fill_manual(values = rep('blue', 13))
write.csv(g$data, 'pred_types.hist.csv',row.names = F)
#选择特定类型的细胞
Idents(YDL) <- pred_obj@meta.data$pred_types
#YDL <- subset(YDL, idents = 'Redblood')
#save(YDL, file = 'YDL.Redblood.Rdata')
Idents(YDL) <- YDL@meta.data$seurat_clusters
#另一个可视化的方法
DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5)

Idents(YDL) <- YDL@meta.data$celltype
pdf("cell.pdf")
#绘制UMI的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("PC_1","PC_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = PC_1,y =PC_2,colour = log(nCount_RNA)))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggplot2)
mydata<- FetchData(YDL,vars = c("PC_1","PC_2","nCount_RNA"))
a <- ggplot(mydata,aes(x = PC_1,y =PC_2,colour = nCount_RNA))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#绘制基因的分布图
library(ggplot2)
mydata<- FetchData(YDL,vars = c("PC_1","PC_2","nFeature_RNA"))
a <- ggplot(mydata,aes(x = PC_1,y =PC_2,colour = nFeature_RNA))+geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))

a+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5)

dev.off()

pdf("干性基因_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Cd34","Meis1","Hoxa9","Gata2"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, slot = "counts",reduction = "pca",pt.size = 1.5,features = c("Cd34","Meis1","Hoxa9","Gata2"),cols = c("gray", "red"))#actin


FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Cd34","Meis1","Hoxa9","Gata2"),cols = c("gray", "red"))#actin
dev.off()


pdf("血红蛋白_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Hbb.bs","Hbb.bt","Hba.a1","Hba.a2"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2"),cols = c("gray", "red"))#actin
dev.off()

pdf("血红蛋白_YDL1 .pdf")
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bs","Gata1","Xpo7","Cd79a"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Hbb-bs","Gata1","Xpo7","Cd79a"),cols = c("gray", "red"))#actin
dev.off()


pdf("转录因子_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Gata1","Klf1","Tal1","Nfe2"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Gata1","Klf1","Tal1","Nfe2"),cols = c("gray", "red"))#actin
dev.off()

pdf("晚期基因_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Xpo7","Tent5c","Cd36","Trim58"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Xpo7","Tent5c","Cd36","Trim58"),cols = c("gray", "red"))#actin
dev.off()
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Xpo7"),cols = c("gray", "red"),split.by = "orig.ident",ncol = 2)#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Hbb-bs"),cols = c("gray", "red"),split.by = "orig.ident",ncol = 2)#actin
pdf("B细胞1_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Igkc","Ly6d","Cd79a","Cd79b"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Igkc","Ly6d","Cd79a","Cd79b"),cols = c("gray", "red"))#actin
dev.off()

pdf("中性粒细胞_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Elane","Mpo","Lyz2","Prtn3"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("S100a8","S100a9","Lgals1","Lgals3"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Elane","Mpo","Lyz2","Prtn3"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("S100a8","S100a9","Lgals1","Lgals3"),cols = c("gray", "red"))#actin
dev.off()


pdf("tsne_result_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Actb","Tmsb4x","Arpc1b","Tagln2","Vim","Cnn2","Gsn","Fmnl1","Hck"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Actb","Tmsb4x","Arpc1b","Tagln2","Vim","Cnn2","Gsn","Fmnl1","Hck"),cols = c("gray", "red"))#actin

FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Actb"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Tmsb4x"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Arpc1b"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Flnb"),cols = c("gray", "red"))#Filamin-alpha
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Tagln2"),cols = c("gray", "red"))#Trangelin2
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Vim"),cols = c("gray", "red"))#Vimentin
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Cnn2"),cols = c("gray", "red"))#Calponin2
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Gsn"),cols = c("gray", "red"))#Gelsolin
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Fmnl1"),cols = c("gray", "red"))#Formin-like1
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hck"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Zyx"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Ear1"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Ear2"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Ear6"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Elane"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("S100a8"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bs"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bt"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bh1"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-bh2"),cols = c("gray", "red"))
FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = c("Hbb-y"),cols = c("gray", "red"))
dev.off()

pdf("炎性红细胞_YDL .pdf")
FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Cd79a","S100a8","Cd44","Ptprc"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("Cd79a","S100a8","Cd44","Ptprc"),cols = c("gray", "red"))#actin
dev.off()


FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Ptp4a3"),cols = c("gray", "red"))#actin


FeaturePlot(object = YDL, reduction = "pca",pt.size = 1.5,features = c("Cd79a","S100a8","Cd44","Ptprc"),cols = c("gray", "red"))#actin

Idents(YDL) <- YDL@meta.data$celltype
YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 't')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker
table(YDL.markers$cluster)
write.csv(YDL.markers,file="allmarker_ternimalE.csv")
#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_all.pdf",width=9,height=12)
DoHeatmap(object = YDL, features = top10$gene)
DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
DoHeatmap(object = YDL, features = c(top10$gene,"Ptp4a3")) + NoLegend()

DoHeatmap(object = YDL, features = YDL.markers$gene,)+ NoLegend()
dev.off()




gene<-read.table("磷酸化酶.txt",header = T)
DoHeatmap(object = YDL, features = gene$gene) + NoLegend()

library(monocle)
#准备monocle分析需要的文件
monocle.matrix=as.matrix(YDL@assays$RNA@data)
monocle.matrix=cbind(id=row.names(monocle.matrix),monocle.matrix)
write.table(monocle.matrix,file="monocleMatrix.txt",quote=F,sep="\t",row.names=F)
monocle.sample=as.matrix(YDL@meta.data)
monocle.sample=cbind(id=row.names(monocle.sample),monocle.sample)
write.table(monocle.sample,file="monocleSample.txt",quote=F,sep="\t",row.names=F)
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.geneAnn=cbind(id=row.names(monocle.geneAnn),monocle.geneAnn)
write.table(monocle.geneAnn,file="monocleGene.txt",quote=F,sep="\t",row.names=F)
write.table(YDL.markers,file="monocleMarkers.txt",sep="\t",row.names=F,quote=F)


#设置工作目录
monocle.matrix=read.table("monocleMatrix.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.sample=read.table("monocleSample.txt",sep="\t",header=T,row.names=1,check.names=F)
monocle.geneAnn=read.table("monocleGene.txt",sep="\t",header=T,row.names=1,check.names=F)
marker=read.table("monocleMarkers.txt",sep="\t",header=T,check.names=F)

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表和基因注释表表
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)


#给其中一列数据重命名
names(pData(cds))[names(pData(cds))=="orig.ident"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
# saveRDS(cds,"cds.rds")
# rm(list = ls())  
# marker<-read.csv("allmarker.csv")
# cds<-readRDS("cds.rds")
#伪时间分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, marker$gene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,reduction_method = 'DDRTree')
cds <- orderCells(cds,root_state = 2)
pdf(file="cluster.trajectory_SINGLET.pdf",width=6.5,height=6)
cds$celltype<-factor(cds$celltype,
                     levels = c("PRO","BASO","POLY","ORTHO"))
plot_cell_trajectory(cds, color_by = "celltype")
plot_cell_trajectory(cds,color_by="orig.ident")+facet_wrap(~orig.ident,nrow=2,ncol = 5)
plot_cell_trajectory(cds,color_by="Cluster")+facet_wrap(~orig.ident,nrow=2,ncol = 3)

plot_cell_trajectory(cds, color_by = "YDL@active.ident")
plot_cell_trajectory(cds, color_by="Pseudotime", show_backbone=FALSE)
# 可以很明显看到细胞的发育轨迹 
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") +facet_wrap(~State, nrow = 1)
dev.off()

saveRDS(cds,"cds_SINGLET.rds")


#使用monocle做拟时序分析（单细胞谱系发育）
#https://mp.weixin.qq.com/s/bkYIAseELBRi8pQYS81mlQ
# 加入为了方便起见，直接挑选top2000的MAD基因。
ordering_genes=names(tail(sort(apply(cds@assayData$exprs,1,mad)),2000))

hm<-plot_pseudotime_heatmap(cds,num_clusters=4,show_rownames=F,return_heatmap = T)
save(hm,file="heatmap.Rdata")
c1 <- as.data.frame(cutree(hm$tree_row, k=4)) 
colnames(c1) <- "Cluster"
c1$Gene <- rownames(c1)
write.table(c1,"genes_in_heatmap_clusters.txt",sep="\t",row.names = F)




disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
dim(cds)

diff_test_res <-differentialGeneTest(cds, fullModelFormulaStr = '~Cluster')
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
head(ordering_genes)
ordering_genes<- row.names (subset(diff_test_res, use_for_ordering==TRUE))
p=plot_pseudotime_heatmap(cds[ordering_genes,],
                          num_clusters = 4,return_heatmap=T,
                          show_rownames = T)
dev.off()
clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

head(clustering)

#选择特定基因进行拟时序热图可视化
Time_genes<-c("Hbb-bs","Hbb-a1","Gypa","Xpo7","Gata1","Nfe2","Tmsb4x",
              "Tmsb10","Tyrobp","Apoe","Vcam1","Jun","Jund","Sox4","Cd44","Ptprc",
              "Lyz2","Lgals3","Elane","Mpo","S100a8","Hbb-bs")

A<-ordering_genes
A<-as.data.frame(A)
head(A)
B<-Time_genes
B<-as.data.frame(B)
head(B)
library (ggvenn)
x<-list("a"=A$A,"b"=B$B)
result<-ggvenn(x)
genes<-intersect(A$A,B$B)
p1 = plot_pseudotime_heatmap(cds[genes,], num_clusters=5, show_rownames=T, return_heatmap=T)
ggsave("Genes_plot1.pdf", plot = p1, width = 16, height = 8)
head(clustering)



data<-cts[,-c(1,2)]
newdata<-data[c(top10$gene),]
pheatmap::pheatmap(newdata,
                   cluster_rows = T,
                   cluster_cols = F)

pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,show_rownames = T,
                   cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,show_rownames = T,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

data<-AverageExpression(YDL,return.seurat = F)
data<-as.data.frame(data)
newdata<-data[c(YDL.markers$gene),]
colnames(newdata)<-c("PRO","BASO","POLY","ORTHO")
pheatmap::pheatmap(newdata,
                   cluster_rows = T,
                   cluster_cols = F)

pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,show_rownames = T,
                   cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,show_rownames = T,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = F,show_rownames = T,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))




normal<-cts[,-c(1,2)]
normal$pro<-apply(normal[,1:8],1, mean, na.rm = T) 
normal$baso<-apply(normal[,9:16],1, mean, na.rm = T) 
normal$poly<-apply(normal[,17:24],1, mean, na.rm = T) 
normal$ortho<-apply(normal[,25:32],1, mean, na.rm = T) 

normal <-normal[,33:36]

newdata<-normal[c(YDL.markers$gene),]
pheatmap::pheatmap(newdata,
                   cluster_rows = T,
                   cluster_cols = F)

pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,show_rownames = T,
                   cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,show_rownames = T,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = F,show_rownames = T,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

#低值为蓝色，高值为红色，中间值为白色：
