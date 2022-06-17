
###生存分析
##rMKL生存分析

#label<-read.csv("D:/wenjian/rMKL_LPP/N9/lable.csv")
label<-read.csv("D:/shuju/lable_24.csv")   ###分类标签是字符
rownames(label)<-label[,1]
label<-label[,-1]

######################
library(survival)
library(survminer)
load("D:/shuju/SNF_analysis/KIRP/.RData")
KIRP_surv<-cbind(KIRP_clinical,label)   ###类别1为对照，2、3、4分别与之比较
#log rank test
sdf=NULL    
p_value=NULL
sdf=survdiff(Surv(overall_survival,as.numeric(vital_status))~label,data=KIRP_surv)
p_value=1-pchisq(sdf$chisq,length(sdf$n) - 1)
fit_cox2<-coxph(Surv(overall_survival,as.numeric(vital_status))~relevel(as.factor(label),ref = 3),data=KIRP_surv)
fit_cox2
summary(fit_cox2)
fit<-survfit(Surv(overall_survival,as.numeric(vital_status))~label,data=KIRP_surv)
fit
ggsurvplot(fit,xlab='survivaltime(day)')

###############试着指定分组，2为对照，3、4与之比较
KIRP_use1<-KIRP_surv[KIRP_surv$label=="type1",]
KIRP_use2<-KIRP_surv[KIRP_surv$label=="type2",]
KIRP_use3<-KIRP_surv[KIRP_surv$label=="type3",]
KIRP_use4<-KIRP_surv[KIRP_surv$label=="type4",]
##合并2，4
KIRP_cbind<-rbind(KIRP_use2,KIRP_use4)
KIRP_cbind$label<-"type2"
##合并后再
KIRP_use<-rbind(KIRP_cbind,KIRP_use3,KIRP_use1)

###1,4合并
KIRP_cbind1<-rbind(KIRP_use1,KIRP_use4)
KIRP_cbind1$label<-"type1"
KIRP_use<-rbind(KIRP_cbind1,KIRP_use2,KIRP_use3)
#KIRP_use<-rbind(KIRP_cbind1,KIRP_use2)


###########
KIRP_use<-rbind(KIRP_use3,KIRP_use1,KIRP_use2,KIRP_use4)

#log rank test
sdf=NULL
p_value=NULL
sdf=survdiff(Surv(overall_survival,as.numeric(vital_status))~label,data=KIRP_use)
p_value=1-pchisq(sdf$chisq,length(sdf$n) - 1)

fit_cox2<-coxph(Surv(overall_survival,as.numeric(vital_status))~factor(label),data=KIRP_use)
fit_cox2
summary(fit_cox2)
fit<-survfit(Surv(overall_survival,as.numeric(vital_status))~label,data=KIRP_use)
fit

ggsurvplot(fit,xlab='survivaltime(day)')
####################################################################
###############计算轮廓系数

library(cluster)
KIRP<-cbind(KIRP_mRNA,KIRP_miRNA,KIRP_methy)
#整合数据集
#样本顺序与组学数据不一样

label1<-read.csv("D:/label_924.csv")   #标签为数字
rownames(label1)<-label1[,1]
label1<-label1[,-1]

d<- data.frame(t(r))#转化行、列
d1<-dist(KIRP_W)#计算相异度矩阵，行为样本，列为变量
d1<-daisy(KIRP_methy, metric = "manhattan", stand = FALSE)#计算相异度矩阵
d1<-daisy(KIRP_methy, metric = "gower", stand = FALSE)#计算相异度矩阵
d1<-daisy(KIRP_mRNA, metric = "euclidean", stand = FALSE)#计算相异度矩阵


sil<-silhouette(label1,d1)#，Cluster为聚类结果（向量）；得到的是每个组别的silhouette值
summary(sil) #计算平均得分







##################FPKM转换TPM数据
#读取未处理的FPKM数据 
mRNA<-read.csv("FPKM_mRNA_KIRP.csv")
##按老师命令处理
.............
#lable the 0 individual in each danba,should check the statues of NA,here is 0
lengthcheck0=function(x){
  length(which(x==0))     #which(is.na(x))
}

mRNA_tumor_KIRP0=apply(mRNA_tumor_KIRP[,-c(1,2)],1,lengthcheck0)
lable_30_KIRP=which(mRNA_tumor_KIRP0>0.3*dim(mRNA_tumor_KIRP[,-c(1,2)])[2])
length(lable_30_KIRP)
if(length(lable_30_KIRP)==0){
  mRNA_tumor_KIRP30=mRNA_tumor_KIRP
}else {
  mRNA_tumor_KIRP30=mRNA_tumor_KIRP[-lable_30_KIRP,]#pay attention,if the lable_30_KIRP=NULL,it will delete all the features
  
}
###这里不进行log转换，直接对其进行KNN填补
is.na(mRNA_tumor_KIRP30[,-c(1,2)])=!mRNA_tumor_KIRP30[,-c(1,2)]#change 0 to NA
#here is already 0,we  need to change
mRNA_tumor_KIRP30_knn=mRNA_tumor_KIRP30
mRNA_tumor_KIRP30_knn[,-c(1,2)]=knnImputation(mRNA_tumor_KIRP30_knn[,-c(1,2)],k=10,scale=T,meth="weighAvg",distData = NULL)
##填补后的数据，进行行名修改
mRNA_tumor_KIRP30_knn_NEW<-mRNA_tumor_KIRP30_knn
rownames(mRNA_tumor_KIRP30_knn_NEW)<-mRNA_tumor_KIRP30_knn_NEW[,1]
mRNA_tumor_KIRP30_knn_NEW<-mRNA_tumor_KIRP30_knn_NEW[,-c(1:2)]
##############################################
##转化函数
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

tpms <- apply(mRNA_tumor_KIRP30_knn_NEW,2,fpkmToTpm)   ##转化后数据
tpms<-as.data.frame(tpms)
##输出TPM数据
save(tpms,file = "DESeq.Rdata")
write.csv(tpms,file = "mRNA_tpms.csv",row.names = F)
################和其他组学都有的样本取交集
tpms_mRNA<-tpms          
tpms_name<-substr(colnames(tpms_mRNA),1,12)
colnames(tpms_mRNA)<-tpms_name
tpms_mRNA<-tpms_mRNA[KIRP_share3]

##查看数据变异情况，看是否要log转化
boxplot(tpms_mRNA,outline=FALSE, notch=T,las=2)

###差异分析
group <- c(rep("normal",12),rep("type1",120), rep("type2",17)) 
ibrary(limma)
label<-read.csv("D:/label.csv")
rownames(label)<-label[,1]
label<-label[,-1] 
design <-model.matrix(~0+factor(label))#把group设置成一个mode1 matrix#
rownames(design)=colnames(tpms_mRNA) 



levels(design)

contrast.matrix <- makeContrasts(GY_7d - normal, 
                                 GY_1d - EC,
                                 GY_7d - GY_1d, 
                                 levels=design)     


group <- c(rep("EC",3),rep("GY_7d",3), rep("GY_1d",3)) 
head(group)
View(group)

group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(GY_7d - EC, 
                                 GY_1d - EC,
                                 GY_7d - GY_1d, 
                                 levels=design)


library(limma)
rownames(label)<-label[,1]
label<-label[,-1] 
design <-model.matrix(~0+factor(label))#把group设置成一个mode1 matrix#
rownames(design)=rownames(KIRP_clinical) 
colnames(design)=c('type1','type2','type3','type4')

##cont.matrix=makeContrasts('type1-type4','type2-type4','type3-type4',levels = design)

label[c('type1','type2','type3','type4')]



########################################################################
# step3  log2(RPM+1) transformation          #
########################################################################

mRNA_tumor_KIRP30_log2=mRNA_tumor_KIRP30
mRNA_tumor_KIRP30_log2[,-c(1,2)]=log2(mRNA_tumor_KIRP30_log2[,-c(1,2)]+1)


########################################################################
