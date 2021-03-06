library(MASS)
library(car)
library(ggplot2)
library(reshape)
require(graphics)
library(graphics)
library(plyr)

setwd("~/Desktop/Spring2016/6907BigDataAnalysis/Project2")
getwd()

mush <- read.csv(file="agaricus-lepiota.data", header=FALSE, sep=",")

str(mush)
mush <- as.data.frame(mush) 
summary(mush)
dim(mush)
myvars <- c("class", "capshape", "capsurface", "capcolor","bruises", "odor","gillattachment","gillspacing","gillsize","gillcolor","stalkshape","stalkroot","stalksurfaceabovering","stalksurfacebelowring","stalkcolorabovering","stalkcolorbelowring","veiltype","veilcolor","ringnumber","ringtype","sporeprintcolor","population","habitat")
colnames(mush) <- myvars
mush <- mush[myvars]
summary(mush)
colnames(mush) #should be 23 vars but need to drop class 
mush <- subset(mush, select = -c(class) )
dim(mush)
str(mush)


mush$capshape <- as.numeric(mush$capshape)
table(mush$capshape)
table(mush$capsurface)
table(mush$capcolor)
table(mush$capshbruisesape)

for (i in names(mush)){
  mush[,i] <- as.numeric(mush[,i])
}

str(mush)

#Explore the categorical predictors. 
str(mush)
summary(mush$capshape)
attach(mush)
sort(table(mush$capshape))
#c=1, s=2, b=3, k=4, f=5, x=6
plot(mush$capshape)

table1 <- (xtabs(~capshape+capsurface+capcolor, data=mush))
table1

table2 <- (xtabs(~bruises+odor, data=mush))
table2

table3 <- (xtabs(~gillattachment+gillspacing+gillcolor+gillsize, data=mush))
table3

table4 <- (xtabs(~stalkshape+stalkroot+stalksurfaceabovering+stalksurfacebelowring+stalkcolorabovering+stalkcolorbelowring, data=mush))
table4

table5 <- (xtabs(~veiltype+veilcolor, data=mush))
table5

table6 <- (xtabs(~ringnumber+ringtype+sporeprintcolor))
table6

table7 <- (xtabs(~population+habitat))
table7

#PCA and other reduction methods 
install.packages("principal")
library(principal) #not available (for R version 3.2.3)
install.packages("princomp")#not available (for R version 3.2.3)
library(princomp)
library(psych)
install.packages("psych")

summary(mush)
head(mush)
#factanal()

head(mush)
factanal(mush)
as.data.frame(mush)

kc <- kmeans(mush, 11, nstart=5)
kc

library(FactoMineR)
PCA1 <- PCA(mush, scale.unit = TRUE, ncp = 22, graph = TRUE,quali.sup=1:13 )
PCA1
PCA1$eig
PCA1$var
summary(PCA1) #we want components 1 - 4 because they have a eigenvalue of 1 and above 
dimdesc(PCA1, axes = 1:4) # to see which variables make an impact on the dimensions 
var(mush)
dimdesc(PCA1, axes = 1:4)
myvars1 <- c("sporeprintcolor","stalkcolorabovering","stalkcolorbelowring","ringnumber", "gillsize", "habitat", "ringtype", "population", "bruises","veilcolor", "gillattachment", "stalkshape")
mush1 <- mush[myvars1]
head(mush)
head(mush1)
str(mush1)
plot(mush1)


PCA2 <- PCA(mush, scale.unit = TRUE, ncp = 22, graph = TRUE)
PCA2
PCA2$eig #will use components 1-7 
PCA2$var
summary(PCA2) # 
dimdesc(PCA2, axes = 1:7)
summary(PCA2) # 
summary(PCA2)
PCA2$var$contrib
PCA2$call
PCA2$ind
PCA2$var$contrib
PCA2$ind$contrib
myvars2 <- c("sporeprintcolor","gillsize", "odor", "stalkshape","bruises","stalkroot", "gillcolor", "ringtype", "stalksurfaceabovering", "stalksurfacebelowring", "gillspacing", "population", "gillattachment","veilcolor", "habitat", "ringnumber", "capshape")
mush2 <- mush[myvars2]
str(mush2)


PCA3 <- PCA(mush, scale.unit = TRUE, ncp = 22, graph = TRUE,quali.sup=(c(4,6,7,8,10,16)))
summary(PCA3)
PCA3$eig #will use components 1-6
dimdesc(PCA3, axes = 1:6)
myvars3 <- c("ringtype","gillcolor", "stalkroot", "bruises","odor", "gillsize", "sporeprintcolor","stalksurfaceabovering","stalksurfacebelowring","bruises","stalkcolorbelowring","stalkcolorabovering","gillspacing","capcolor","population", "stalkshape","habitat","ringnumber","capsurface","veilcolor","gillattachment")
mush3 <- mush[myvars3]
dim(mush3)
