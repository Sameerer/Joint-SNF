
##rMKL生存分析
#label<-read.csv("D:/wenjian/rMKL_LPP/N9/lable.csv")
label<-read.csv("D:/label.csv")
rownames(label)<-label[,1]
label<-label[,-1]
library(survival)
library(survminer)
load("D:/shuju/SNF_analysis/KIRP/.RData")
KIRP_surv<-cbind(KIRP_clinical,label)
#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival,as.numeric(vital_status))~label,data=KIRP_surv)
p_value=1-pchisq(sdf$chisq,length(sdf$n) - 1)

fit_cox2<-coxph(Surv(overall_survival,as.numeric(vital_status))~factor(label),data=KIRP_surv)
fit_cox2
summary(fit_cox2)
fit<-survfit(Surv(overall_survival,as.numeric(vital_status))~label,data=KIRP_surv)
fit
ggsurvplot(fit,xlab='survivaltime(day)')

#diff=survdiff(Surv(futime, fustat) ~a,data = rt)
#pValue=1-pchisq(diff$chisq,df=1)
#pValue=round(pValue,5)  #保留五位小数点
#########SNF的生存分析
KIRP_surv<-cbind(KIRP_clinical,KIRP_SNF_lable[[3]])
sdf <-survdiff(Surv(overall_survival,as.numeric(vital_status))~KIRP_SNF_lable[[3]],data=KIRP_surv)
sdf

fit_cox2<-coxph(Surv(overall_survival,as.numeric(vital_status))~factor(KIRP_SNF_lable[[3]]),data=KIRP_surv)
fit_cox2
summary(fit_cox2)

fit<-survfit(Surv(overall_survival,as.numeric(vital_status))~KIRP_SNF_lable[[3]],data=KIRP_surv)
fit

ggsurvplot(fit,xlab='survivaltime(day)')
###################################################
#下游分析——基因差异表达分析
foldChange=1
padj=0.05

library(edgeR)
library(limma)
#group=c("normal","tumor","tumor","normal","tumor")
group=c(rep("normal",4),rep("tumor",178))                         #按照癌症和正常样品数目修改
design <-model.matrix(~0+factor(group))
KIRP_miRNA_T<-as.data.frame(t(KIRP_miRNA))

y <- DGEList(counts=KIRP_miRNA_T,design)


#y <- DGEList(counts=KIRP_miRNA_T,group=factor(label))    ###整合数据，使得包可识别
y <- calcNormFactors(y)     ##因子
y <- estimateCommonDisp(y)  ##估计正常、病人的方差、变异系数
y <- estimateTagwiseDisp(y)
et <- exactTest(y) 
#et <- exactTest(y,pair = c("normal","tumor"))   ##差异表达分析，内部差异、组间差异
topTags(et)   ##
ordered_tags <- topTags(et, n=437)   ##输出n个差异表达情况

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]   ##437x4矩阵
diff=allDiff
newData=y$pseudo.counts

write.table(diff,file="edgerOut.xls",sep="\t",quote=F)
#####输出显著差异表达基因

diffSig = diff[(diff$FDR < padj & (diff$logFC>foldChange | diff$logFC<(-foldChange))),]    ###输出显著差异表达基因
write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)
diffUp = diff[(diff$FDR < padj & (diff$logFC>foldChange)),]   ##3上调基因，FDR>0
write.table(diffUp, file="up.xls",sep="\t",quote=F)
diffDown = diff[(diff$FDR < padj & (diff$logFC<(-foldChange))),]   ##下调基因
write.table(diffDown, file="down.xls",sep="\t",quote=F)

normalizeExp=rbind(id=colnames(newData),newData)
write.table(normalizeExp,file="normalizeExp.txt",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）
diffExp=rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）

heatmapData <- newData[rownames(diffSig),]
hmExp=log10(heatmapData+0.001)
library('gplots')
hmMat=as.matrix(hmExp)
pdf(file="heatmap.pdf",width=60,height=90)
par(oma=c(10,3,3,7))
heatmap.2(hmMat,col='greenred',trace="none")



boxplot(data.frame (exprSet),co1="blue") ##画箱式图，比较数据分布情祝
exprSet[1:5,1:5]
group <- read.csv("D:/wenjian/rMKL_LPP/N9/lable.csv", header=TRUE,row.names=1, check.names=FALSE)
group <- group[,1] #定义比较组，按照癌症和正常样品数目修改#
design <-model.matrix(~0+factor(group))#把group设置成一-个mode1 matrix#
co1names(design)=1evel(factor(group))
rownames (design)=colnames (exprSet)

############################################
#############################################
#将FPKM转换为TPM
expMatrix <- a
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)

