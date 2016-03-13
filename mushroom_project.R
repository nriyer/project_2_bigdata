library(MASS)
library(dplyr)
library(rvgtest)
library(Hmisc)
library(reshape2)
library(caret)
#set wd
setwd("~/Desktop/GW_spring2016/bigdata")
mr <- read.table('agaricus-lepiota.data',sep=",",header=FALSE)
str(mr)
#set names
names(mr) <- c("class","cshape","csurface","ccolor","bruises","odor","gattach","gspace",
"gsize","gcolor","sshape","sroot","ssabove","ssbelow","scabove","scbelow","vtype","vcolor",
"rnumber","rtype","spcolor","popnum","habitat")
names(mr)
#data exploration
#plot all variables vs. class and ftable for each against class
for (i in names(mr)){
  plot(class ~ mr[,i], data=mr,xlab = i)
  print(ftable(mr$class,mr[,i]))
}
#look at pairwise correlations (look through all)
pairs(class ~ ccolor, data = mr)

##-------get ready for clustering:
numeric <- mr
#this translates column from factor to number
#sorted (?)
numeric$cshape <- sort(as.numeric(numeric$cshape))

numeric$cshape <- as.numeric(numeric$cshape)
#loop to translate all to factors
for (i in names(numeric)){
  numeric[,i] <- as.numeric(numeric[,i])
}

##melt and try ggplot to see each var against target
numeric$class <- mr$class
melt.mr <- melt(numeric)
#value represented by 'stat="identity)
ggplot(data=melt.mr,aes(x=variable,y=value)) + geom_bar(stat="identity")
qplot(variable, color=value, data=melt.mr, geom='density')

#take out target var and turn rest into numeric for clustering
cluster <- numeric
cluster <- numeric[,-1]
#no variance in vtype
cluster$vtype <- NULL

#split into training and test
row <- nrow(cluster)
col <- col(cluster)
set.seed(12345)
##train and test set
trainindex <- sample(row, 0.7*row, replace=FALSE)
train <- cluster[trainindex,]
test <- cluster[-trainindex,]
dim(train)

#preProcess, use 'range' and then also try 
#to use center,scale for std. dev 1 and mean 0
#to get each factor between 0 and 1 for distance measures
#to be more accurate

##range 
str(train)
pre_range <- preProcess(train,method="range")
train_processed <- predict(pre_range, train)
str(test)
post_range <- preProcess(test,method="range")
test_processed <- predict(post_range, train)
test_processed_t <- as.data.frame(test_processed)


##centered and scaled for normal distribution
pre_center_scale <- preProcess(train,method=c('center','scale'))
post_center_scale <- preProcess(test,method=c('center','scale'))
train_processed_cs <- predict(pre_center_scale, train)
test_processed_cs <- predict(post_center_scale,test)

######################
##Clustering         #
######################

####################################
#K MEANS:                         ##
####################################

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
  print(table(fit$cluster))
  km.out<- kmeans(cluster_assignment,k,nstart=n)
  print(km.out$tot.withinss)
  print(km.out$betweenss)
  print(km.out)
}

library(NbClust)
#First use preprocessed range 0-1 data set:

#TRAIN MODEL using TRAINING SCALED BY 'RANGE'
train_processed <- as.matrix(train_processed)
#function to find optimal number of clusters (takes too long)
km.numcluster <- NbClust(train_processed, method = 'kmeans', min.nc=2,max.nc=12)

##START WITH BEST NUMBER OF K BY USING WSS AND BSS PLOTS
#instead look at within sum squares (wss) and between sum squares (bss)
wss_and_bss(train_processed)

#Both plots show that around 5 clusters is the point where sum of squares starts to level out
#SO we will use K=5, 5 clusters
#you want larger sum of squares between and within since this will make more distinct clusters

set.seed(12345)
#Kmeans for 5 clusters 
kmeans_iyer(train_processed,5)
check_5clusters <- cluster_assignment[which(cluster_assignment$cshape == .4),c(1,2,22)]
#output of number of points per cluster
#1    2    3    4    5 
#1591 1191  890 1475  539 

#Within cluster sum of squares by cluster:
#  [1] 2406.0151  353.2646 1003.6869  274.5883  979.7177
#(between_SS / total_SS =  77.9 %)

###TEST
kmeans_iyer(test_processed,5)
#output of number of points per cluster
#1    2    3    4    5 
#647  539 2450 1219  831 
#Within cluster sum of squares by cluster:
#  [1]  353.2646 2406.0151  274.5883  979.7177 1003.6869
#(between_SS / total_SS =  79.8 %)

######################
## 3 and 7 clusters  #
######################

##now try 3 clusters
kmeans_iyer(train_processed, 3)
#output of number of points per cluster
# K-means clustering with 3 clusters of sizes 1775, 2231, 1680
#Within cluster sum of squares by cluster:
#  [1]  994.2487 6149.4506  373.2694
#(between_SS / total_SS =  53.5 %)
kmeans_iyer(test_processed, 3)
#K-means clustering with 3 clusters of sizes 1775, 2231, 1680
#Within cluster sum of squares by cluster:
#  [1] 2281.522 2432.911 2738.665
#(between_SS / total_SS =  53.3 %)

##7 clusters
kmeans_iyer(train_processed, 7)
#K-means clustering with 7 clusters of sizes 396, 466, 795, 562, 539, 1515, 1413
#Within cluster sum of squares by cluster:
#  [1]   70.19977  579.30406  159.42366  292.17378  274.58829 2359.89559  795.79288
#(between_SS / total_SS =  87.3 %)
kmeans_iyer(test_processed, 7)
#K-means clustering with 7 clusters of sizes 1475, 950, 751, 656, 786, 529, 539
#Within cluster sum of squares by cluster:
#  [1]  979.7177 1553.7875  535.3927  137.0492 1185.7640  129.7308  274.5883
#(between_SS / total_SS =  85.7 %)

###5 clusters does the best, test the other pre processed data set on this:

wss_and_bss(train_processed_cs)
#looks around 6 but lets see performance of 5
kmeans_iyer(train_processed_cs,5)
#Within cluster sum of squares by cluster:
#  [1]  4193.817 24242.098 10881.830 16372.751  8120.760
#(between_SS / total_SS =  49.5 %)
kmeans_iyer(train_processed_cs,6)
#Within cluster sum of squares by cluster:
#  [1] 12068.212  7041.161 17683.127 10814.880  1944.343  5894.806
#(between_SS / total_SS =  59.9 %)

##6 does better, but i think 'range' is a better scaling technique for this data.

####################################
#HIERARCHICAL                     ##
####################################
#hierarchical starts with all variables and forms clusters based on cutoff criteria 
#and distance measures given

#use training and test normalized using 'range'

train_processed <- as.matrix(train_processed)
test_processed <- as.matrix(test_processed)

hc.complete=hclust(dist(train_processed), method="complete")
hc.average=hclust(dist(train_processed), method="average")
hc.single=hclust(dist(train_processed), method="single")


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

###################
# SVM             #
###################


