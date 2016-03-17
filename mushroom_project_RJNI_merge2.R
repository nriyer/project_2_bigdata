library(MASS)
library(dplyr)
library(rvgtest)
library(Hmisc)
library(reshape2)
library(caret)
library(FactoMineR)
library(ggplot2)
library(GGally)

#set wd
setwd("~/Desktop/Spring2016/6907BigDataAnalysis/Project2")
getwd()
mr <- read.table('agaricus-lepiota.data',sep=",",header=FALSE)
str(mr)

#########################################
# DEFINE FUNCTIONS THAT WILL BE USED    #
#########################################
#define function for preprocessing - using 'range'

scale <- function(df){
  pre_range <- preProcess(df,method="range")
  processed <- predict(pre_range,df)
  return(data.frame(processed))
}

#define function for wss and bss plots

wss_and_bss <- function(df){
  #within sum of squares
  wss <- (nrow(df)-1)*sum(apply(df,2,var))
  for (i in 1:12) wss[i] <- sum(kmeans(df, 
                                       centers=i)$withinss)
  print(plot(1:12, wss, type="b", xlab="Number of Clusters",
             ylab="Within groups sum of squares"))
  
  #between sum of squares
  bss <- (nrow(df)-1)*sum(apply(df,2,var))
  for (i in 1:12) bss[i] <- sum(kmeans(df, 
                                       centers=i)$betweenss)
  print(plot(1:12, bss, type="b", xlab="Number of Clusters",
             ylab="Between groups sum of squares"))
}


#define function for final output
kmeans_iyer <- function(df,k,n=20){
  fit <-kmeans(df,k)
  #get cluster means:this illustrates amount of each characteristic in each cluster
  aggregate(df,by=list(fit$cluster), FUN=mean)
  #append cluster assignment
  cluster_assignment <- data.frame(df, fit$cluster)
  #table of how many variables fit in each cluster
  #print(table(fit$cluster))
  km.out<- kmeans(cluster_assignment,k,nstart=n)
  print(km.out$tot.withinss)
  print(km.out$betweenss)
  print(km.out)
}

################################################################################################

#set names
names(mr) <- c("class", "capshape", "capsurface", "capcolor","bruises", "odor","gillattachment",
                "gillspacing","gillsize","gillcolor","stalkshape","stalkroot","stalksurfaceabovering",
                "stalksurfacebelowring","stalkcolorabovering","stalkcolorbelowring","veiltype",
                "veilcolor","ringnumber","ringtype","sporeprintcolor","population","habitat")
head(mr)

#####################################
## CHANGE ALPHA VARS TO NUMERIC    ##
#####################################
numeric <- mr
numeric[numeric == '?'] <- NA ##Protect the null values
numeric <- sapply(numeric, as.numeric)
#this translates column from factor to number
#loop to translate all to factors
str(numeric)
summary(numeric)

#####################################################################################
## You will need to divide your data set into a training set and a test set.       ##
## Use samples of 50-50, 60-40, and 70-30 for the training-test ratios             ##
##This is further down, right before clustering but noted based on project outline ##
#####################################################################################


#####################################################################################
## Try plotting the data using several plotting functions to see what it looks like## 
## Use pairs (e.g., 2D plots) or 3 variables (3D plots) based on the packages.     ##
## Try to filter the data by selecting samples with only certain attribute values  ##
## and plotting them.##                                                            ##
#####################################################################################

#data exploration
#plot all variables vs. class and ftable for each against class
for (i in names(mr)){
  plot(class ~ mr[,i], data=mr,xlab = i)
  print(ftable(mr$class,mr[,i]))
}


#look at pairwise correlations (look through all)
library(GGally)
library(ggplot2)
ggpairs(numeric, columns = c(1, 2, 3)) #cap
ggpairs(numeric, columns = c(4,5)) #bruises and odor
ggpairs(numeric, columns = c(6, 7, 8, 9)) #gill
ggpairs(numeric, columns = c(10, 11, 12, 13, 14, 15)) #stalk 
ggpairs(numeric, columns = c(16, 17)) #veil
ggpairs(numeric, columns = c(18, 19, 20)) #ring and spore
ggpairs(numeric, columns = c(21, 22)) #pop, habitat
numeric2 <- subset(numeric, select = -c(class) )
cor(numeric2)

#Continue to look at the data prior to PCA and Clustering 
mush <- subset(numeric, select = -c(class) )
dim(mush)
str(mush)
head(mush)
table1 <- (xtabs(~capshape+capsurface+capcolor, data=mush))
table1
plot(table1)

table2 <- (xtabs(~bruises+odor, data=mush))
table2
plot(table2)

table3 <- (xtabs(~gillattachment+gillspacing+gillcolor+gillsize, data=mush))
table3
plot(table3)

table4 <- (xtabs(~stalkshape+stalkroot+stalksurfaceabovering+stalksurfacebelowring+stalkcolorabovering+stalkcolorbelowring, data=mush))
table4
plot(table4)

table5 <- (xtabs(~veiltype+veilcolor, data=mush))
table5
plot(table5)

table6 <- (xtabs(~ringnumber+ringtype+sporeprintcolor))
table6
plot(table6)

table7 <- (xtabs(~population+habitat))
table7
plot(table7)

str(numeric)
##melt and try ggplot to see each var against target
numeric$class <- mr$class
melt.mr <- melt(numeric)
#value represented by 'stat="identity)
ggplot(data=melt.mr,aes(x=variable,y=value)) + geom_bar(stat="identity")
qplot(variable, color=value, data=melt.mr, geom='density')



#############################################################################################
## Try data reduction to eliminate some attributes through Principal Components Analysis.  ##
#############################################################################################
#####################################################
# remove variables with zero or near-zero variance  #
#####################################################

nzv <- nearZeroVar(numeric, saveMetrics = TRUE)
## only keep those columns that significant variance
mr_nzv <- numeric[,c(rownames(nzv[nzv$nzv == FALSE,]))]



#################################
# Principal Component Analysis  #
#################################


mr_nzv <- subset( mr_nzv, select = -c(class) )
pr <- PCA( mr_nzv, ncp = dim(mr_nzv)[2])

## take the first components that have
## eigenvalues greater than 1
index <- max(which(pr$eig$eigenvalue > 1))
mr.t <- pr$ind$coord[,1:index]
mr.t <- as.data.frame(mr.t)

## plot that shows how nice this transformation is
library(car)
f <- as.factor(numeric[,1])
scatter3d(mr.t[,1],mr.t[,2],mr.t[,3], surface = FALSE, groups = f)

######################################
# INITIAL CLUSTERING - WITH CLASS    #
######################################

numeric <- mr
#this translates column from factor to number
#loop to translate all to factors
for (i in names(numeric)){
  numeric[,i] <- as.numeric(numeric[,i])
}
head(numeric)
dim(numeric)

initial <- numeric
set.seed(12345)
row <- nrow(initial)
col <- col(initial)
#train/test 70-30
trainindex <- sample(row, 0.7*row, replace=FALSE)
train_initial <- initial[trainindex,]
test_initial <- initial[-trainindex,]
dim(train_initial)
##plot wss and bss
wss_and_bss(train_initial)
## optimum number here is 4 clusters
testclust <- kmeans_iyer(train_initial,4)
sort(testclust$size)
#K-means clustering with 4 clusters of sizes 1536, 1988, 1332, 830
#Within cluster sum of squares by cluster:
#  [1] 52552.72 55329.35 22906.27 28094.94
#(between_SS / total_SS =  47.2 %)
kmeans_iyer(test_initial,4)
#Within cluster sum of squares by cluster:
#  [1] 25666.340  8925.573 10536.847 21587.822
#(between_SS / total_SS =  47.8 %)


##Now to run PCA and reduce dimensions, then try clustering again using
## K=3, K=5 and K=7 AND two other clustering methods - K-Nearest Neighbors
## and Hierarchical clustering


#############################################################################
## 2. This will involve some statistical analysis and some clustering.     ##
## Use the r packages and functions in the notes as well as the ones below ## 
## CLUSTERING - now on reduced dimensionality set, we will use mush3       ##
## except for on the optimal number of clusters, where we will try         ##
## all three PCA sets #                                                    ##
#############################################################################
## Include:(a)clustering of the samples into 2 (edible or poisonous)classes##
## using three different clustering methods, such kmeans k-nearest neighbor##
## one other b. a clustering of the samples into N = 3, 5, 7 classes using ##
## three different clustering methods, such kmeans, k-nearest neighbor, one##
## other with number of points per cluster                                 ##
#############################################################################


#take out target var and turn rest into numeric for clustering
cluster <- mr.t
#cluster <- mush1[,-1]
#no variance in vtype
cluster$vtype <- NULL

#split into training and test
row <- nrow(cluster)
col <- col(cluster)
set.seed(12345)

##train and test at 50-50 split
trainindex <- sample(row, 0.5*row, replace=FALSE)
train_50 <- cluster[trainindex,]
test_50 <- cluster[-trainindex,]
dim(train_50)


##train and test at 60-40 split
trainindex <- sample(row, 0.6*row, replace=FALSE)
train_60 <- cluster[trainindex,]
test_60 <- cluster[-trainindex,]
dim(train_60)


##train and test set at 70-30 split
trainindex <- sample(row, 0.7*row, replace=FALSE)
train_70 <- cluster[trainindex,]
test_70 <- cluster[-trainindex,]
dim(train_70)

#preProcess, use 'range' and then also try 
#to use center,scale for std. dev 1 and mean 0
#to get each factor between 0 and 1 for distance measures
#to be more accurate

####################
# 50-50 train/test #
####################

##range 
train1 <- scale(train_50)
test1 <- scale(test_50)

##turn both into matrix
train1 <- as.matrix(train1)
test1 <- as.matrix(test1)

wss_and_bss(train1)
#optimal looks like 4 clusters. Will use K=4, K=5 and K=7

set.seed(12345)
#K=3
kmeans_iyer(train1,3)

#K-means clustering with 4 clusters of sizes 1436, 863, 1763
# Within cluster sum of squares by cluster:
# [1] 8033.1858  473.2955 4839.9445
# (between_SS / total_SS =  41.3 %)

kmeans_iyer(test1,4)
# Within cluster sum of squares by cluster:
#   [1] 1804.6730 6733.1266  898.5788  712.2011
# (between_SS / total_SS =  58.4 %)

#K=5
kmeans_iyer(train1,5)
#K-means clustering with 5 clusters of sizes 1277, 715, 338, 775, 957

# Within cluster sum of squares by cluster:
#   [1] 3881.9698  749.3130  758.7224 1079.6333  630.6668
# (between_SS / total_SS =  75.3 %)


kmeans_iyer(test1,5)
# Within cluster sum of squares by cluster:
#   [1] 1804.7616 1187.0728 1046.9565  712.2011 2310.7929
# (between_SS / total_SS =  75.4 %)

#K=7
kmeans_iyer(train1,7)
#K-means clustering with 7 clusters of sizes 453, 532, 666, 325, 854, 857, 375
# Within cluster sum of squares by cluster:
#   [1] 611.5559 453.6008 630.8959 709.4998 442.2659 374.6532 284.5403
# (between_SS / total_SS =  90.4 %)

kmeans_iyer(test1,7)
# Within cluster sum of squares by cluster:
#   [1] 350.8962 498.7762 503.0216 469.8624 632.9195 417.7128 638.1423
# (between_SS / total_SS =  90.9 %)

###7 does the best but could be overfitting also may have something to do with 50-50 split


####################
# 60-40 train/test #
####################
##range 
train2 <- scale(train_60)
test2 <- scale(test_60)

##turn both into matrix
train2 <- as.matrix(train2)
test2 <- as.matrix(test2)

wss_and_bss(train2)
#optimal looks like 4 clusters. Will use K=4, K=5 and K=7

set.seed(12345)
#K=3
kmeans_iyer(train2,3)

#K-means clustering with 4 clusters of sizes 1618, 2228, 1028
# Within cluster sum of squares by cluster:
# [1] 9105.5893 6409.5309  549.7098
# (between_SS / total_SS =  42.9 %)

kmeans_iyer(test2,3)
# Within cluster sum of squares by cluster:
# [1] 3294.8793 7193.6698  426.3172
# (between_SS / total_SS =  39.1 %)

#K=5
kmeans_iyer(train2,5)
#K-means clustering with 5 clusters of sizes 869, 1217, 814, 1035, 939

# Within cluster sum of squares by cluster:
#   [1] 4944.5984  819.5048  870.6196  571.8846 1285.0629
# (between_SS / total_SS =  74.4 %)

kmeans_iyer(test2,5)
# Within cluster sum of squares by cluster:
#   [1]  305.3508  712.1916  625.1542  400.0601 3840.9398
# (between_SS / total_SS =  75.0 %)

#K=7
kmeans_iyer(train2,7)
#K-means clustering with 7 clusters of sizes 816, 990, 714, 382, 445, 341, 1186
# Within cluster sum of squares by cluster:
#   [1]  876.00502 1488.07158  671.03950  849.68993  368.34202   69.32215  704.29210
# (between_SS / total_SS =  86.5 %)

kmeans_iyer(test2,7)
#Within cluster sum of squares by cluster:
# [1] 352.3810 472.7462 501.5998 396.1190 529.1403 309.3507 274.8304
# (between_SS / total_SS =  90.5 %)
###7 does the best again..

####################
# 70-30 train/test #
####################

##range 
train3 <- scale(train_70)
test3 <- scale(test_70)

##turn both into matrix
train3 <- as.matrix(train3)
test3 <- as.matrix(test3)

wss_and_bss(train3)
#optimal looks like 3 clusters. Will use K=3, K=5 and K=7

set.seed(12345)
#K=3
kmeans_iyer(train3,3)

#K-means clustering with 3 clusters of sizes 957, 1203, 3526
# Within cluster sum of squares by cluster:
#   [1]  1012.8757   666.3672 17173.1493
# (between_SS / total_SS =  41.5 %)


kmeans_iyer(test3,3)
# Within cluster sum of squares by cluster:
#   [1] 7198.9749  328.5411  424.3353
# (between_SS / total_SS =  41.5 %)

#K=5
kmeans_iyer(train3,5)
#K-means clustering with 5 clusters of sizes 1209, 975, 2516, 553, 433

# Within cluster sum of squares by cluster:
#   [1]  685.2579 1150.2749 6862.6231  467.2397  978.9314
# (between_SS / total_SS =  75.0 %)

kmeans_iyer(test3,5)
# Within cluster sum of squares by cluster:
#   [1]  377.9383  319.0392  403.2486  186.4270 2828.9784
# (between_SS / total_SS =  75.2 %)


#K=7
kmeans_iyer(train3,7)
#K-means clustering with 7 clusters of sizes 426, 1200, 1229, 906, 805, 555, 565

#Within cluster sum of squares by cluster:
# [1] 954.5145 661.5366 545.6186 874.6575 679.6864 477.9926 721.8379
# (between_SS / total_SS =  90.0 %)

kmeans_iyer(test3,7)
#Within cluster sum of squares by cluster:
# [1] 242.2359 399.5511 347.7086 387.1052 170.7704 260.8635 298.3935
# (between_SS / total_SS =  90.9 %)

#7 clusters does the best. Try 7 clusters with all 3 PCA sets??:



####################################
#HIERARCHICAL                     ##
####################################
#hierarchical starts with all variables and forms clusters based on cutoff criteria 
#and distance measures given

#use training and test 70-30 split, normalized using 'range'

train3 <- as.matrix(train3)
test3<- as.matrix(test3)

hc.complete=hclust(dist(train3), method="complete")
hc.average=hclust(dist(train3), method="average")
hc.single=hclust(dist(train3), method="single")


#this is a little hard to make sense of - trying another method below this
par(mfrow=c(1,3))
plot(hc.complete,main="Complete Linkage", xlab="", sub="", cex=.9)
plot(hc.average, main="Average Linkage", xlab="", sub="", cex=.9)
plot(hc.single, main="Single Linkage", xlab="", sub="", cex=.9)
dev.off()

#Ward Hierarchical Clustering
#Ward's minimum variance criterion minimizes the total within-cluster variance. 
#To implement this method, at each step find the pair of clusters that leads to minimum increase 
#in total within-cluster variance after merging. This increase is a weighted squared distance between 
#cluster centers. At the initial step, all clusters are singletons (clusters containing a single point). 
#To apply a recursive algorithm under this objective function, 
#the initial distance between individual objects must be (proportional to) squared Euclidean distance.

d <- dist(train_processed, method = "euclidean") # distance matrix
fit_ward <- hclust(d, method="ward.D") 
plot(fit_ward) # display dendogram
groups <- cutree(fit_ward, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit_ward, k=5, border="red")


######################################################
# K-NN nearest neigbors to edible and poisonous     ##
######################################################
# KNN uses the Euclidian distance measure in order to find the k-nearest neighbours to your new, 
# unknown instance. Here, the k parameter is one that you set yourself. 
# As mentioned before, new instances are classified by looking at the majority vote or weighted vote. 
# In case of classification, the data point with the highest score wins the battle and the 
# unknown instance receives the label of that winning data point. 
# If there is an equal amount of winners, the classification happens randomly.


##reference: https://www.datacamp.com/community/tutorials/machine-learning-in-r
library(class)
library(gmodels)
str(numeric)
#use numeric dataframe. predicting class, this is the only one i do not change into a factor.
table(numeric$class)
#classes are not biased

#training and test for knn
ind <- sample(2, nrow(numeric), replace=TRUE, prob=c(.67,.33))
knn.train <- numeric[ind==1,2:23]
knn.test <- numeric[ind==2,2:23]
#train and test labels for target var
knn.trainLabels <- numeric[ind==1,1]
knn.testLabels <- numeric[ind==2,1]

#k nearest on training using 10 as k
knn_pred <- knn(train = knn.train, test = knn.test, cl = knn.trainLabels, k=10)
knn_pred

#cross table
CrossTable(x = knn.testLabels, y = knn_pred, prop.chisq=FALSE)
#true negative/true positive = model score
#1292/1311 = .9855


