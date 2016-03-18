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

###kmeans###
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


###hierarchical clustering###

#starts with all variables and forms clusters based on cutoff criteria 
#and distance measures given
initial <- numeric
set.seed(12345)
row <- nrow(initial)
col <- col(initial)
#train/test 70-30
trainindex <- sample(row, 0.7*row, replace=FALSE)
train_initial <- initial[trainindex,]
test_initial <- initial[-trainindex,]
train_initial <- as.matrix(train_initial)
test_initial<- as.matrix(test_initial)

hc0.complete=hclust(dist(train_initial), method="complete")
hc0.average=hclust(dist(train_initial), method="average")
hc0.single=hclust(dist(train_initial), method="single")
plot(hc0.complete)

#cutree will cut the height of the cluster results/dendrogram 
cutree(hc0.complete, k = 15)




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
#Will use K=3, K=5 and K=7

set.seed(12345)
#K=3
kmeans_iyer(train1,3)

#K-means clustering with 3 clusters of sizes 1251, 1769, 1042
# Within cluster sum of squares by cluster:
#[1] 182.5771 216.5760 161.1758
#(between_SS / total_SS =  85.9 %)

set.seed(12345)
kmeans_iyer(test1,3)

#Within cluster sum of squares by cluster:
#  [1]  64.16686  27.87945 424.57828
#(between_SS / total_SS =  83.6 %)

#K=5
set.seed(12345)
kmeans_iyer(train1,5)
#5 clusters of sizes 385, 1546, 341, 926, 864
#Within cluster sum of squares by cluster:
#  [1]  15.35515 216.60085  34.77083  24.07242  20.05762
#(between_SS / total_SS =  96.8 %)

set.seed(12345)
kmeans_iyer(test1,5)
#Within cluster sum of squares by cluster:
#  [1]  42.72332  33.22492 309.30968  22.45656  31.46635
#(between_SS / total_SS =  96.2 %)

#K=7
set.seed(12345)
kmeans_iyer(train1,7)
#K-means clustering with 7 clusters of sizes 1015, 306, 277, 718, 929, 519, 298
# Within cluster sum of squares by cluster:
#[1]  24.072417 214.198568  10.810926   8.947576   7.518888   5.491708  14.469189
#(between_SS / total_SS =  98.2 %)

set.seed(12345)
kmeans_iyer(test1,7)
#Within cluster sum of squares by cluster:
#  [1]   2.349587  24.304690 101.517284   2.640166  56.058236  29.129220  10.465214
#(between_SS / total_SS =  98.5 %)

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
#Will use K=3, K=5 and K=7

set.seed(12345)
#K=3
kmeans_iyer(train2,3)

#K-means clustering with 3 clusters of sizes 1906, 1258, 1710
# Within cluster sum of squares by cluster:
#[1]  30.65981  21.86613 619.22763
#(between_SS / total_SS =  85.4 %)

set.seed(12345)
kmeans_iyer(test2,3)
#Within cluster sum of squares by cluster:
#  [1]  28.58319 339.18316  47.75000
#(between_SS / total_SS =  86.6 %)

#K=5
set.seed(12345)
kmeans_iyer(train2,5)
#K-means clustering with 5 clusters of sizes 5 clusters of sizes 1217, 1170, 822, 1039, 626

#Within cluster sum of squares by cluster:
#  [1] 140.25096  26.81310  36.22674  24.59775  67.11366
#(between_SS / total_SS =  97.7 %)

set.seed(12345)
kmeans_iyer(test2,5)
#Within cluster sum of squares by cluster:
#  [1]  14.26809 143.20761  23.93566  18.42919  25.01050
#(between_SS / total_SS =  97.3 %)

#K=7
set.seed(12345)
kmeans_iyer(train2,7)
#K-means clustering with 7 clusters of sizes 816, 990, 714, 382, 445, 341, 1186
#Within cluster sum of squares by cluster:
#  [1] 25.061781 17.640326 63.765930 37.614394 14.279230  7.414986 35.484587
#(between_SS / total_SS =  98.9 %)

set.seed(12345)
kmeans_iyer(test2,7)
#Within cluster sum of squares by cluster:
#  [1]   8.428006 180.325270  18.381682  15.081955  14.268094   7.489344   8.300694
#(between_SS / total_SS =  98.4 %)
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

#K-means clustering with 3 clusters of sizes 1906, 2813, 967
# Within cluster sum of squares by cluster:
#[1] 313.01196 406.58541  40.68651
#(between_SS / total_SS =  87.3 %)

set.seed(12345)
kmeans_iyer(test3,3)
#Within cluster sum of squares by cluster:
#  [1]  17.58074  44.67666 253.86466
#(between_SS / total_SS =  83.3 %)

#K=5
set.seed(12345)
kmeans_iyer(train3,5)
#K-means clustering with 5 clusters of sizes 448, 546, 2414, 476, 1802

#Within cluster sum of squares by cluster:
#  [1]  45.64282  15.56042 250.85349  17.37074 269.37067
#(between_SS / total_SS =  94.9 %)

set.seed(12345)
kmeans_iyer(test3,5)
#Within cluster sum of squares by cluster:
#  [1] 31.10798 15.51408 19.82691 13.61199 71.30005
#(between_SS / total_SS =  97.3 %)


#K=7
set.seed(12345)
kmeans_iyer(train3,7)
#K-means clustering with 7 clusters of sizes 717, 1239, 639, 905, 561, 1197, 428

#Within cluster sum of squares by cluster:
#  [1] 26.06610 22.40491 34.47960 34.33005 23.17867 26.66523 42.34009
#(between_SS / total_SS =  99.1 %)

set.seed(12345)
kmeans_iyer(test3,7)
#Within cluster sum of squares by cluster:
#  [1]  31.574552  12.224544   5.966243   2.052997 147.159923   8.296618  19.107347
#(between_SS / total_SS =  98.1 %)

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
plot(hc.complete)

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


