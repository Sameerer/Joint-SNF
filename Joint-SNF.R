

library(SNFtool)
library(r.jive)

##set all the parameters:
K = 20; ##number of neighbors, usually (10~30)
alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
T = 20; ###Number of iterations, usually (10~50)
nperms=100;# number of permutations for rank estimation
alpha=0.05;# the quantile to use for rank estimation.

###Step1 Extraction of joint structures by JIVE
#A list of multiple data matrices measured on same set of samples.
#For each data matrix in the list, samples should be on rows and features should be on columns.
data_list<-list(data1,data2,data3) 
Results <- jive(data_list,nperms=100,alpha=0.05)
#summary(Results);Provides a summary of JIVE output
##### Visualization
#showVarExplained(Results);Display Variance Explained
#showHeatmaps(Results);Heatmaps for JIVE Decompositions
#####joint structures
J1<-Results[["joint"]][[1]] 
J2<-Results[["joint"]][[2]]
J3<-Results[["joint"]][[3]]
J1<-data.frame(t(J1))
J2<-data.frame(t(J2))
J3<-data.frame(t(J3))
###### individual structure 
#I1<-Results[["individual"]][[1]]
#I2<-Results[["individual"]][[2]]
#I3<-Results[["individual"]][[3]]
#I1<-data.frame(t(I1))
#I2<-data.frame(t(I2))
#I3<-data.frame(t(I3))

###Step2 Similarity Network Fusion
J1_sd=standardNormalization(J1)
J2_sd=standardNormalization(J2)
J3_sd=standardNormalization(J3)
# Calculate distance matrices
Joint_SNF_Dist1=(dist2(as.matrix(J1_sd),as.matrix(J1_sd)))^(1/2)
Joint_SNF_Dist2=(dist2(as.matrix(J2_sd),as.matrix(J2_sd)))^(1/2)
Joint_SNF_Dist3=(dist2(as.matrix(J3_sd),as.matrix(J3_sd)))^(1/2)
# construct similarity graphs
Joint_SNF_W1=affinityMatrix(Joint_SNF_Dist1,K,alpha)
Joint_SNF_W2=affinityMatrix(Joint_SNF_Dist2,K,alpha)
Joint_SNF_W3=affinityMatrix(Joint_SNF_Dist3,K,alpha)
# next, we fuse all the graphs
Joint_SNF_list<-list(Joint_SNF_W1,Joint_SNF_W2,Joint_SNF_W3)
Joint_SNF_W=SNF(Joint_SNF_list,K,T)
# get cluster labels for each data point by spectral clustering
# C=3; number of clusters
Joint_lable<-spectralClustering(Joint_SNF_W,C=3)


