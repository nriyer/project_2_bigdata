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

#Explore the categorical predictors. 
str(mush)
summary(mush$capshape)
attach(mush)
sort(table(mush$capshape))
#c=1, s=2, b=3, k=4, f=5, x=6
plot(mush$capshape)
summary(capshape)
revalue(mush$capshape, c("c" = 1)) -> mush$capshape
revalue(mush$capshape, c("s" = 2)) -> mush$capshape
revalue(mush$capshape, c("b" = 3)) -> mush$capshape
revalue(mush$capshape, c("k" = 4)) -> mush$capshape
revalue(mush$capshape, c("f" = 5)) -> mush$capshape
revalue(mush$capshape, c("x" = 6)) -> mush$capshape


sort(table(mush$capsurface))
#g=1, f=2, s=3, y=4
plot(table(mush$capsurface))
revalue(mush$capsurface, c("g" = 1)) -> mush$capsurface
revalue(mush$capsurface, c("f" = 2)) -> mush$capsurface
revalue(mush$capsurface, c("s" = 3)) -> mush$capsurface
revalue(mush$capsurface, c("y" = 4)) -> mush$capsurface
sort(table(mush$capsurface))


sort(table(mush$capcolor))
#r=1, u=2, c=3, p=4, b=5,w=6, y=7, e=8, g=9, n=10
plot(table(mush$capcolor))
revalue(mush$capcolor, c("r" = 1)) -> mush$capcolor
revalue(mush$capcolor, c("u" = 2)) -> mush$capcolor
revalue(mush$capcolor, c("c" = 3)) -> mush$capcolor
revalue(mush$capcolor, c("p" = 4)) -> mush$capcolor
revalue(mush$capcolor, c("b" = 5)) -> mush$capcolor
revalue(mush$capcolor, c("w" = 6)) -> mush$capcolor
revalue(mush$capcolor, c("y" = 7)) -> mush$capcolor
revalue(mush$capcolor, c("e" = 8)) -> mush$capcolor
revalue(mush$capcolor, c("g" = 9)) -> mush$capcolor
revalue(mush$capcolor, c("n" = 10)) -> mush$capcolor
plot(table(mush$capcolor))


sort(table(mush$bruises))
#t=1, f=2
plot(table(mush$bruises))
revalue(mush$bruises, c("t" = 1)) -> mush$bruises
revalue(mush$bruises, c("f" = 2)) -> mush$bruises

sort(table(mush$odor))
#m=1, c=2, p=3, a=4, l=5, s=6, y=7, f=8, n=9
plot(table(mush$odor))
revalue(mush$odor, c("m" = 1)) -> mush$odor
revalue(mush$odor, c("c" = 2)) -> mush$odor
revalue(mush$odor, c("p" = 3)) -> mush$odor
revalue(mush$odor, c("a" = 4)) -> mush$odor
revalue(mush$odor, c("l" = 5)) -> mush$odor
revalue(mush$odor, c("s" = 6)) -> mush$odor
revalue(mush$odor, c("y" = 7)) -> mush$odor
revalue(mush$odor, c("f" = 8)) -> mush$odor
revalue(mush$odor, c("n" = 9)) -> mush$odor
plot(table(mush$odor))
sort(table(mush$odor))

sort(table(mush$gillattachment))
#a=1, f=2
plot(mush$gillattachment)
plot(table(mush$gillattachment))
sort(table(mush$gillattachment))
revalue(mush$gillattachment, c("a" = 1)) -> mush$gillattachment
revalue(mush$gillattachment, c("f" = 2)) -> mush$gillattachment
plot(mush$gillattachment)


sort(table(mush$gillspacing))
#w=1, c=2
plot(mush$gillspacing)
plot(table(mush$gillattachment))
revalue(mush$gillspacing, c("w" = 1)) -> mush$gillspacing
revalue(mush$gillspacing, c("c" = 2)) -> mush$gillspacing
plot(mush$gillspacing)


sort(table(mush$gillsize))
#n=1, b=2
plot(mush$gillsize)
revalue(mush$gillsize, c("n" = 1)) -> mush$gillsize
revalue(mush$gillsize, c("b" = 2)) -> mush$gillsize

sort(table(mush$gillcolor))
#r=1, o=2, y=3, e=4, k=5, u=6, h=7, g=8, n=9, w=10, p=11, b=12 
plot(mush$gillcolor)
revalue(mush$gillcolor, c("r" = 1)) -> mush$gillcolor
revalue(mush$gillcolor, c("o" = 2)) -> mush$gillcolor
revalue(mush$gillcolor, c("y" = 3)) -> mush$gillcolor
revalue(mush$gillcolor, c("e" = 4)) -> mush$gillcolor
revalue(mush$gillcolor, c("k" = 5)) -> mush$gillcolor
revalue(mush$gillcolor, c("u" = 6)) -> mush$gillcolor
revalue(mush$gillcolor, c("h" = 7)) -> mush$gillcolor
revalue(mush$gillcolor, c("g" = 8)) -> mush$gillcolor
revalue(mush$gillcolor, c("n" = 9)) -> mush$gillcolor
revalue(mush$gillcolor, c("w" = 10)) -> mush$gillcolor
revalue(mush$gillcolor, c("p" = 11)) -> mush$gillcolor
revalue(mush$gillcolor, c("b" = 12)) -> mush$gillcolor
sort(table(mush$gillcolor))

sort(table(mush$gillsize))
#n=1, b=2
plot(mush$gillsize)
revalue(mush$gillsize, c("n" = 1)) -> mush$gillsize
revalue(mush$gillsize, c("b" = 2)) -> mush$gillsize

sort(table(mush$stalkshape))
#e=1, t=2
plot(mush$stalkshape) 
revalue(mush$stalkshape, c("e" = 1)) -> mush$stalkshape
revalue(mush$stalkshape, c("t" = 2)) -> mush$stalkshape

sort(table(mush$stalkroot))
#r=1, c=2, e=3, ?=4, b=5
#what does ? mean?  Is this missing? 
plot(mush$stalkroot)
revalue(mush$stalkroot, c("r" = 1)) -> mush$stalkroot
revalue(mush$stalkroot, c("c" = 2)) -> mush$stalkroot
revalue(mush$stalkroot, c("e" = 3)) -> mush$stalkroot
revalue(mush$stalkroot, c("?" = ".")) -> mush$stalkroot
revalue(mush$stalkroot, c("b" = 5)) -> mush$stalkroot


sort(table(mush$stalksurfaceabovering))
#y=1, f=2, k=3, s=4
plot(mush$stalksurfaceabovering)
revalue(mush$stalksurfaceabovering, c("y" = 1)) -> mush$stalksurfaceabovering
revalue(mush$stalksurfaceabovering, c("f" = 2)) -> mush$stalksurfaceabovering
revalue(mush$stalksurfaceabovering, c("k" = 3)) -> mush$stalksurfaceabovering
revalue(mush$stalksurfaceabovering, c("s" = 4)) -> mush$stalksurfaceabovering

sort(table(mush$stalksurfacebelowring))
#y=1, f=2, k=3, s=4
plot(mush$stalksurfacebelowring)
revalue(mush$stalksurfacebelowring, c("y" = 1)) -> mush$stalksurfacebelowring
revalue(mush$stalksurfacebelowring, c("f" = 2)) -> mush$stalksurfacebelowring
revalue(mush$stalksurfacebelowring, c("k" = 3)) -> mush$stalksurfacebelowring
revalue(mush$stalksurfacebelowring, c("s" = 4)) -> mush$stalksurfacebelowring


sort(table(mush$stalkcolorabovering))
#y=1, c=2, e=3, o=4, b=5, n=6, g=7, p=8, w=9
plot(mush$stalkcolorabovering)
#w has the largest frequency
revalue(mush$stalkcolorabovering, c("y" = 1)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("c" = 2)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("e" = 3)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("o" = 4)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("b" = 5)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("n" = 6)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("g" = 7)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("p" = 8)) -> mush$stalkcolorabovering
revalue(mush$stalkcolorabovering, c("w" = 9)) -> mush$stalkcolorabovering

sort(table(mush$stalkcolorbelowring))
#y=1, c=2, e=3, o=4, b=5, n=6, g=7, p=8, w=9
plot(mush$stalkcolorbelowring)
revalue(mush$stalkcolorbelowring, c("y" = 1)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("c" = 2)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("e" = 3)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("o" = 4)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("b" = 5)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("n" = 6)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("g" = 7)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("p" = 8)) -> mush$stalkcolorbelowring
revalue(mush$stalkcolorbelowring, c("w" = 9)) -> mush$stalkcolorbelowring


sort(table(mush$veiltype))
#p=1
revalue(mush$veiltype, c("p" = 1)) -> mush$veiltype


sort(table(mush$veilcolor))
#y=1, n=2, o=3, w=4
revalue(mush$veilcolor, c("y" = 1)) -> mush$veilcolor
revalue(mush$veilcolor, c("n" = 2)) -> mush$veilcolor
revalue(mush$veilcolor, c("o" = 3)) -> mush$veilcolor
revalue(mush$veilcolor, c("w" = 4)) -> mush$veilcolor

sort(table(mush$ringnumber))
#n=1, t=2, o=3
revalue(mush$ringnumber, c("n" = 1)) -> mush$ringnumber
revalue(mush$ringnumber, c("t" = 2)) -> mush$ringnumber
revalue(mush$ringnumber, c("o" = 3)) -> mush$ringnumber

sort(table(mush$ringtype))
#n=1, f=2, l=3, e=4, p=5
revalue(mush$ringtype, c("n" = 1)) -> mush$ringtype
revalue(mush$ringtype, c("f" = 2)) -> mush$ringtype
revalue(mush$ringtype, c("l" = 3)) -> mush$ringtype
revalue(mush$ringtype, c("e" = 4)) -> mush$ringtype
revalue(mush$ringtype, c("p" = 5)) -> mush$ringtype


sort(table(mush$sporeprintcolor))
#b=1, o=2, u=3, y=4, r=5, h=6, k=7, n=8, w=9
revalue(mush$sporeprintcolor, c("b" = 1)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("o" = 2)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("u" = 3)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("y" = 4)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("r" = 5)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("h" = 6)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("k" = 7)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("n" = 8)) -> mush$sporeprintcolor
revalue(mush$sporeprintcolor, c("w" = 9)) -> mush$sporeprintcolor

sort(table(mush$population))
#c=1, a=2, n=3, s=4, y=5, v=6
revalue(mush$population, c("c" = 1)) -> mush$population
revalue(mush$population, c("a" = 2)) -> mush$population
revalue(mush$population, c("n" = 3)) -> mush$population
revalue(mush$population, c("s" = 4)) -> mush$population
revalue(mush$population, c("y" = 5)) -> mush$population
revalue(mush$population, c("v" = 6)) -> mush$population

sort(table(mush$habitat))
#w=1, m=2, u=3, l=4, p=5, g=6, d=7
revalue(mush$habitat, c("w" = 1)) -> mush$habitat
revalue(mush$habitat, c("m" = 2)) -> mush$habitat
revalue(mush$habitat, c("u" = 3)) -> mush$habitat
revalue(mush$habitat, c("l" = 4)) -> mush$habitat
revalue(mush$habitat, c("p" = 5)) -> mush$habitat
revalue(mush$habitat, c("g" = 6)) -> mush$habitat
revalue(mush$habitat, c("d" = 7)) -> mush$habitat

summary(mush)
attach(mush)

#You will need to divide your data set into a training set and a test set. 
#Use samples of 50-50, 60-40, and 70-30 for the training-test ratios

row<-nrow(mush)
col<-col(mush)
set.seed(12345)

#50/50 split 
trainindex1 <- sample(row, 0.5*row, replace=FALSE)
train50 <- mush[trainindex1,]
test50 <- mush[-trainindex1,]
dim(train50)
dim(test50)

#60/40 split 
trainindex2 <- sample(row, 0.6*row, replace=FALSE)
train60 <- mush[trainindex2,]
test40 <- mush[-trainindex2,]
dim(train60)
dim(test40)

#70/30 split 
trainindex <- sample(row, 0.7*row, replace=FALSE)
train70 <- mush[trainindex,]
test30 <- mush[-trainindex,]
dim(train70)
dim(test30)

#Try plotting the data using several plotting functions to see what it looks like. 
#Use pairs (e.g., 2D plots) or 3 variables (3D plots) based on the packages. 
dev.off()
dim(mush)
summary(mush)
library(ggplot2)
ggpairs(mush)
ggpairs(mush, columns = c(1, 2, 3)) #cap
ggpairs(mush, columns = c(4,5)) #bruises and odor
ggpairs(mush, columns = c(6, 7, 8, 9)) #gill
ggpairs(mush, columns = c(10, 11, 12, 13, 14, 15)) #stalk 
ggpairs(mush, columns = c(16, 17)) #veil
ggpairs(mush, columns = c(18, 19, 20)) #ring and spore
ggpairs(mush, columns = c(21, 22)) #pop, habitat

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
