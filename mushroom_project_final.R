library(MASS)
library(dplyr)
library(rvgtest)
library(Hmisc)
library(reshape2)
library(caret)
library(class)
library(FactoMineR)
library(ggplot2)
library(GGally)
library(gmodels)
library(dendextend)

options(warn=-1)

#########################################
# DEFINE FUNCTIONS THAT WILL BE USED    #
#########################################
#define function for preprocessing - using 'range'

scale2 <- function(df){
  pre_range <- preProcess(df,method="range")
  processed <- predict(pre_range,df)
  return(data.frame(processed))
}

#define function for wss and bss plots
#in our final submission we do not use this function
#however it was very useful for our understanding
#so we leave it in
wss_and_bss <- function(df){
  #within sum of squares
  wss <- (length(df)-1)*sum(apply(df,2,var))
  for (i in 1:12) wss[i] <- sum(kmeans(df, 
                                       centers=i)$withinss)
  print(plot(1:12, wss, type="b", xlab="Number of Clusters",
             ylab="Within groups sum of squares"))
  
  #between sum of squares
  bss <- (length(df)-1)*sum(apply(df,2,var))
  for (i in 1:12) bss[i] <- sum(kmeans(df, 
                                       centers=i)$betweenss)
  print(plot(1:12, bss, type="b", xlab="Number of Clusters",
             ylab="Between groups sum of squares"))
}


#define function for final output
#that returns the algorithm score
#see write-up for score choice
# FOR KMEANS
#train data on train dataset
#take centers from training and cluster test data around those centroids
#  report back the withinness.tot/toalss for test data
#  report back the number of elements per cluster
kmeans_metric <- function(train,test,k){
  library(clue)
  fit <- kmeans(train,k)
  
  test.predict <- cl_predict(fit, test)
  fitted.test <- fit$centers[test.predict,]
  resid.test <- test - fitted.test
  
  ss <- function(x) sum(scale(x, scale = FALSE)^2)
  
  betweenness <- ss(fitted.test) ##compute sum square of betweeness
  withinness <- sapply(split(test, test.predict), ss) ## sum os square for withinness for each variable
  withinness.tot <- sum(withinness)
  
  cat("# (between_SS / total_SS =  ")
  cat(betweenness/(withinness.tot+betweenness) * 100, " %\n")
  cat("# Sizes: ", sort(table(test.predict)))

}

#define function for final output
#that returns the algorithm score
#see write-up for score choice
# FOR KMEANS
#train has so many columns, but the first is the labels for that column
#test has the same number of columns as train
#test.labels are the "true" labels for the test dataset for our score to be calculated
knn_metric <- function(train,test,train.labels,test.labels,k) { 

  knn_pred <- knn(train = train, test = test, cl = train.labels, k=k)
  
  table <- CrossTable(x = test.labels, y = knn_pred, prop.chisq=FALSE)$t
  cat("# (accuracy = ")
  cat(sum(diag(table))/sum(table) * 100)
  cat("%)\n")
  cat("# Sizes: ", table(knn_pred))
}

#define function for final output
#that returns the algorithm score
#see write-up for score choice
# FOR HIERARCHICAL
# return correlation metric between two training and test dendograms
# return counts of elements per group
hierarchical_metric <- function(train, test, k) {
  
  d.train <- dist(train) # distance matrix
  h.train <- hclust(d.train, method="ward.D") 
  groups <- cutree(h.train, k=k)
  
  ## get centroids of the clusters that the 
  ## hierarchical clustering found with training data
  centroids <- aggregate(train, list(groups), mean)[,2:6]
  centroids <- as.matrix(centroids)
  
  ## find which clustering centroid the test data is closest to
  distanceMatrix <- matrix(NA, nrow=dim(test)[1], ncol=dim(centroids)[1])
  for(i in 1:nrow(centroids)) {
    distanceMatrix[,i] <- sqrt(rowSums(t(t(test)-centroids[i,])^2))
  }
  clusters <- apply(distanceMatrix, 1, which.min)
  
  ## find sum of squares information and report that 
  fitted.test <- centroids[clusters,]
  resid.test <- test - fitted.test
  
  ss <- function(x) sum(scale(x, scale = FALSE)^2)
  
  betweenness <- ss(fitted.test) ##compute sum square of betweeness
  withinness <- sapply(split(test, test.predict), ss) ## sum os square for withinness for each variable
  withinness.tot <- sum(withinness)
  
  
  cat("# (between_SS / total_SS =  ")
  cat(betweenness/(withinness.tot+betweenness) * 100, " %\n")
  cat("# Sizes: ", sort(table(clusters)))
  
}

#####################################################################################
## Load and clean data                                                             ##
#####################################################################################

#load data
mr <- read.table('agaricus-lepiota.data',sep=",",header=FALSE)

#set names
names(mr) <- c("class", "capshape", "capsurface", "capcolor","bruises", "odor","gillattachment",
                "gillspacing","gillsize","gillcolor","stalkshape","stalkroot","stalksurfaceabovering",
                "stalksurfacebelowring","stalkcolorabovering","stalkcolorbelowring","veiltype",
                "veilcolor","ringnumber","ringtype","sporeprintcolor","population","habitat")

#####################################
## change alpha vars to numeric    ##
#####################################
numeric <- mr
numeric[numeric == '?'] <- NA ##Protect the null values
numeric <- as.data.frame(sapply(numeric, as.numeric))

numeric <- numeric-1
numeric <- scale2(numeric)

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

####################
### summary statistics
####################
summary(numeric) #entire dataset

edible <- numeric[which(numeric[,1] == 0),]
poisonous <- numeric[which(numeric[,1] == 1),]
summary(edible) #just edible ones
summary(poisonous) # just poisonous ones

####################
### PLOTS
####################

# simple scatter plots
plot(numeric[,2:3], main="capshape vs. capsurface")
f <- as.factor(numeric[,1])
scatter3d(numeric[,2],numeric[,3], numeric[,4], surface = FALSE, groups = f)

#plot all variables vs. class and ftable for each against class
for (i in names(mr)){
  plot(class ~ mr[,i], data=mr,xlab = i)
  print(ftable(mr$class,mr[,i]))
}

#look at pairwise correlations for the entire dataset (look through all)
library(GGally)
library(ggplot2)

par(c(1,3))
ggpairs(numeric, columns = c(1, 2, 3), title = "Entire Sample") #cap
ggpairs(edible, columns = c(1, 2, 3), title = "Edible Sample") #cap
ggpairs(poisonous, columns = c(1, 2, 3), title = "Poisonous Sample") #cap

ggpairs(numeric, columns = c(4,5), title = "Entire Sample") #bruises and odor
ggpairs(edible, columns = c(4,5), title = "Edible Sample") #bruises and odor
ggpairs(poisonous, columns = c(4,5), title = "Poisonous Sample") #bruises and odor

ggpairs(numeric, columns = c(6, 7, 8, 9), title = "Entire Sample") #gill
ggpairs(edible, columns = c(6, 7, 8, 9), title = "Edible Sample") #gill
ggpairs(poisonous, columns = c(6, 7, 8, 9), title = "Poisonous Sample") #gill

ggpairs(numeric, columns = c(10, 11, 12, 13, 14, 15), title = "Entire Sample") #stalk 
ggpairs(edible, columns = c(10, 11, 12, 13, 14, 15), title = "Edible Sample") #stalk 
ggpairs(poisonous, columns = c(10, 11, 12, 13, 14, 15), title = "Poisonous Sample") #stalk 

ggpairs(numeric, columns = c(16, 17), title = "Entire Sample") #veil
ggpairs(edible, columns = c(16, 17), title = "Edible Sample") #veil
ggpairs(poisonous, columns = c(16, 17), title = "Poisonous Sample") #veil

ggpairs(numeric, columns = c(18, 19, 20), title = "Entire Sample") #ring and spore
ggpairs(edible, columns = c(18, 19, 20), title = "Edible Sample") #ring and spore
ggpairs(poisonous, columns = c(18, 19, 20), title = "Poisonous Sample") #ring and spore

ggpairs(numeric, columns = c(21, 22), title = "Entire Sample") #pop, habitat
ggpairs(edible, columns = c(21, 22), title = "Edible Sample") #pop, habitat
ggpairs(poisonous, columns = c(21, 22), title = "Poisonous Sample") #pop, habitat


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
mr.test <- pr$ind$coord[,1:index]
mr.test <- as.data.frame(mr.test)

## how much variance are we explaining with this choice of basis?
pr$eig$`cumulative percentage of variance`[index]

## plot that shows how nice this transformation is
library(car)
f <- as.factor(numeric[,1])
scatter3d(mr.test[,1],mr.test[,2],mr.test[,3], surface = FALSE, groups = f, title(main="First three eignvectors as data basis"))


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


row <- dim(mr)[1]

##train and test at 50-50 split
trainindex50 <- sample(row, 0.5*row, replace=FALSE)
train_50 <- mr.test[trainindex50,]
train_50.label <- numeric[trainindex50,1]
test_50 <- mr.test[-trainindex50,]
test_50.label <- numeric[-trainindex50,1]


##train and test at 60-40 split
trainindex60 <- sample(row, 0.6*row, replace=FALSE)
train_60 <- mr.test[trainindex60,]
train_60.label <- numeric[trainindex60,1]
test_40 <- mr.test[-trainindex60,]
test_40.label <- numeric[-trainindex60,1]


##train and test set at 70-30 split
trainindex70 <- sample(row, 0.7*row, replace=FALSE)
train_70 <- mr.test[trainindex70,]
train_70.label <- numeric[trainindex70,1]
test_30 <- mr.test[-trainindex70,]
test_30.label <- numeric[-trainindex70,1]


####################
# k-means          #
####################

## N=2 ##

#50-50#
kmeans_metric(train_50,test_50,2)
# (between_SS / total_SS =  17.26922 %)
# Sizes:  530 3532

#60-40#
kmeans_metric(train_60,test_40,2)
# (between_SS / total_SS =  17.42808 %)
# Sizes:  1512 1738

#70-30#
kmeans_metric(train_70,test_30,2)
# (between_SS / total_SS =  26.95782 %)
# Sizes:  917 1521


## N=3 ##

#50-50#
kmeans_metric(train_50,test_50,3)
# (between_SS / total_SS =  41.05853 %)
# Sizes:  980 1288 1794

#60-40#
kmeans_metric(train_60,test_40,3)
# (between_SS / total_SS =  38.88912 %)
# Sizes:  575 698 1977

#70-30#
kmeans_metric(train_70,test_30,3)
# (between_SS / total_SS =  35.22559 %)
# Sizes:  321 544 1573


## N=5 ##

#50-50#
kmeans_metric(train_50,test_50,5)
# (between_SS / total_SS =  71.78903 %)
# Sizes:  326 377 712 902 1745

#60-40#
kmeans_metric(train_60,test_40,5)
# (between_SS / total_SS =  70.64303 %)
# Sizes:  408 594 699 768 781

#70-30#
kmeans_metric(train_70,test_30,5)
# (between_SS / total_SS =  72.89514 %)
# Sizes:  311 424 537 567 599



## N=7 ##

#50-50#
kmeans_metric(train_50,test_50,7)
# (between_SS / total_SS =  84.4691 %)
# Sizes:  312 377 473 483 669 866 882

#60-40#
kmeans_metric(train_60,test_40,7)
# (between_SS / total_SS =  74.33949 %)
# Sizes:  358 358 367 405 527 548 687

#70-30#
kmeans_metric(train_70,test_30,7)
# (between_SS / total_SS =  84.76866 %)
# Sizes:  253 273 279 306 397 400 530






######################################################
# K-NN nearest neigbors to edible and poisonous     ##
# performed on train/test of 50-50, 60-40 and 70-30 ##
######################################################
# KNN uses the Euclidian distance measure in order to find the k-nearest neighbours to your new, 
# unknown instance. Here, the k parameter is one that you set yourself. 
# As mentioned before, new instances are classified by looking at the majority vote or weighted vote. 
# In case of classification, the data point with the highest score wins the battle and the 
# unknown instance receives the label of that winning data point. 
# If there is an equal amount of winners, the classification happens randomly.

## k=2 ##

#50-50#
knn_metric(train_50,test_50, train_50.label, test_50.label, 2)
# (accuracy = 99.63072%)
# Sizes:  2098 1964

#60-40#
knn_metric(train_60,test_40, train_60.label, test_40.label, 2)
# (accuracy = 99.75385%)
# Sizes:  1702 1548

#70-30#
knn_metric(train_70,test_30, train_70.label, test_30.label, 2)
# (accuracy = 99.7539%)
# Sizes:  1232 1206


## k=3 ##

#50-50#
knn_metric(train_50,test_50, train_50.label, test_50.label, 3)
# (accuracy = 99.67996%)
# Sizes:  2096 1966

#60-40#
knn_metric(train_60,test_40, train_60.label, test_40.label, 3)
# (accuracy = 99.75385%)
# Sizes:  1704 1546

#70-30#
knn_metric(train_70,test_30, train_70.label, test_30.label, 3)
# (accuracy = 99.67186%)
# Sizes:  1234 1204


## k=5 ##

#50-50#
knn_metric(train_50,test_50, train_50.label, test_50.label, 5)
# (accuracy = 99.48301%)
# Sizes:  2104 1958

#60-40#
knn_metric(train_60,test_40, train_60.label, test_40.label, 5)
# (accuracy = 99.66154%)
# Sizes:  1707 1543

#70-30#
knn_metric(train_70,test_30, train_70.label, test_30.label, 5)
# (accuracy = 99.54881%)
# Sizes:  1237 1201


## k=7 ##

#50-50#
knn_metric(train_50,test_50, train_50.label, test_50.label, 7)
# (accuracy = 99.40916%)
# Sizes:  2101 1961

#60-40#
knn_metric(train_60,test_40, train_60.label, test_40.label, 7)
# (accuracy = 99.50769%)
# Sizes:  1712 1538

#70-30#
knn_metric(train_70,test_30, train_70.label, test_30.label, 7)
# (accuracy = 99.46678%)
# Sizes:  1237 1201



####################################
#HIERARCHICAL                     ##
####################################

#Ward Hierarchical Clustering 
#Ward's minimum variance criterion minimizes the total within-cluster variance. 
#To implement this method, at each step find the pair of clusters that leads to minimum increase 
#in total within-cluster variance after merging. This increase is a weighted squared distance between 
#cluster centers. At the initial step, all clusters are singletons (clusters containing a single point). 
#To apply a recursive algorithm under this objective function, 
#the initial distance between individual objects must be (proportional to) squared Euclidean distance.


## N=2 ##

#50-50#
hierarchical_metric(train_50,test_50,2)
# (between_SS / total_SS =  20.91378  %
# Sizes:  872 3190

#60-40#
hierarchical_metric(train_60,test_40,2)
# (between_SS / total_SS =  21.96063  %
# Sizes:  679 2571

#70-30#
hierarchical_metric(train_70,test_30,2)
# (between_SS / total_SS =  21.76016  %
# Sizes:  561 1877


## N=3 ##

#50-50#
hierarchical_metric(train_50,test_50,3)
# (between_SS / total_SS =  30.78866  %
# Sizes:  869 948 2245

#60-40#
hierarchical_metric(train_60,test_40,3)
# (between_SS / total_SS =  32.52458  %
# Sizes:  674 775 1801

#70-30#
hierarchical_metric(train_70,test_30,3)
# (between_SS / total_SS =  31.42805  %
# Sizes:  558 558 1322


## N=5 ##

#50-50#
hierarchical_metric(train_50,test_50,5)
# (between_SS / total_SS =  43.28707  %
# Sizes:  386 713 880 928 1155

#60-40#
hierarchical_metric(train_60,test_40,5)
# (between_SS / total_SS =  45.85888  %
# Sizes:  318 577 684 777 894

#70-30#
hierarchical_metric(train_70,test_30,5)
# (between_SS / total_SS =  43.85703  %
# Sizes:  227 419 551 565 676



## N=7 ##

#50-50#
hierarchical_metric(train_50,test_50,7)
# (between_SS / total_SS =  47.10419  %
# Sizes:  310 375 433 507 671 874 892

#60-40#
hierarchical_metric(train_60,test_40,7)
# (between_SS / total_SS =  49.19209  %
# Sizes:  234 307 363 388 543 677 738

#70-30#
hierarchical_metric(train_70,test_30,7)
# (between_SS / total_SS =  47.21947  %
# Sizes:  169 219 287 291 391 522 559
