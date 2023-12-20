
#****************************************************************************************
# Load packages
#****************************************************************************************
#For Trees
library(mclust)
library(tree)
library(randomForest)
#For PCA
library(maps)
library(akima)
library(fields)
library(tidyverse)
#Function to prune tree
source("PruneTree.r")
#Load this only when plotting a randomforest tree because
#it masks necssary commands: 
#source("RandomForestPlot.R")
#Load Data - monthly. 2008 to 2015, generally months 4-11
data=read.csv("lipid_data_cyto2.csv") %>% 
        mutate(Total_N = NO3+NO2+NH4)
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variaMGDGe (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$MGDG
data <- as.data.frame(cbind(Y,Xvar))


#****************************************************************************************
# CART: Fit tree prune, and scatterplot of observed vs CART estimated.
#****************************************************************************************

myTree <- tree(Y ~ ., data = data, model = T)
#Perform CV on tree object
cvTree <- cv.tree(myTree)
optTree <- which.min(cvTree$dev)
bestTree <- cvTree$size[optTree]

#Prune tree based on CV results
pruneTree <- prune.tree(myTree, best = bestTree)

#Use tree to predict original data
treePred <- predict(pruneTree)
treeResid <- resid(pruneTree)
myRange <- range(treePred, data$Y)

#Plot Unpruned Tree
plot(myTree)
text(myTree, cex = .75)
title(main = c("Unpruned Tree for MGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for MGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for MGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="MGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <-sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  2881
#> R2
#[1]0.75


#****************************************************************************************
# Random Forest
#****************************************************************************************
#Fit random forest and plot
forest.T = randomForest(Y ~ ., data = data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for MGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = ("Pred vs True for Random Forest of MGDG (ng/L)"))
lines(data$Y, data$Y, col = "MGDGack")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <-sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

varImpPlot(forest.T, pch = 20, main = "Importance of VariaMGDGes")
#Most important variaMGDGe is PC2.

#Can also see importance of the variaMGDGes as taMGDGe
importance(forest.T)

#Tree visuallization - look at tree  number - load source above
#Can plot a given tree number in the random forest:
tree_num = 350
x=tree_func(forest.T,tree_num)

#Can use getTree to to also visualize tree number 350, since the Error VS
#trees plot had stabilized by then.
getTree(forest.T, k=tree_num, labelVar=TRUE)










#************************************************************************************************************************************
#
#*## All depths and mesocosms rerun random forest with selected variables only
data=read.csv("lipid_data_cyto2.csv") %>% 
        mutate(Total_N = NO3+NO2+NH4)
#drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "POP", "pH", "Diat", "POC", "Prym", "SiO2", "PON", "Synecho_all", "DGDG", "Pelago", "Dino", "Crypto", "NO3", "NH4", "FL4_group", "Other", "Bsi", "NP_Ratio", "NO2", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "Synecho_low", "Synecho_high", "Synechococcus", "MikroI", "MikroII", "Chloro")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$MGDG
data <- as.data.frame(cbind(Y,Xvar))


#****************************************************************************************
# CART: Fit tree prune, and scatterplot of observed vs CART estimated.
#****************************************************************************************

myTree <- tree(Y ~ ., data = data, model = T)
#Perform CV on tree object
cvTree <- cv.tree(myTree)
optTree <- which.min(cvTree$dev)
bestTree <- cvTree$size[optTree]

#Prune tree based on CV results
pruneTree <- prune.tree(myTree, best = bestTree)

#Use tree to predict original data
treePred <- predict(pruneTree)
treeResid <- resid(pruneTree)
myRange <- range(treePred, data$Y)

#Plot Unpruned Tree
plot(myTree)
text(myTree, cex = .75)
title(main = c("Unpruned Tree for MG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for MG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for MG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed MG (ng/L)", ylab="Modeled MG (ng/L)",
     main = c("Pred vs True for MG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="MG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <- sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  335
#> R2
#[1]0.79


#****************************************************************************************
# Random Forest
#****************************************************************************************
#Fit random forest and plot
forest.T = randomForest(Y ~ ., data = data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for MG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed MG (ng/L)", ylab="Modeled MG (ng/L)",
     main = ("Pred vs True for Random Forest of MG (ng/L)"))
lines(data$Y, data$Y, col = "black")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <- sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

varImpPlot(forest.T, pch = 20, main = "Importance of Variables")
#Most important variable is PC2.

#Can also see importance of the variables as table
importance(forest.T)

#Tree visuallization - look at tree  number - load source above
#Can plot a given tree number in the random forest:
tree_num = 350
x=tree_func(forest.T,tree_num)

#Can use getTree to to also visualize tree number 350, since the Error VS
#trees plot had stabilized by then.
getTree(forest.T, k=tree_num, labelVar=TRUE)






























############################################
#### Including pacific





data=randomtree_pacific
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variaMGDGe (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$MGDG
data <- as.data.frame(cbind(Y,Xvar))
data <- data[-c(33), ]

#****************************************************************************************
# CART: Fit tree prune, and scatterplot of observed vs CART estimated.
#****************************************************************************************

myTree <- tree(Y ~ ., data = data, model = T)
#Perform CV on tree object
cvTree <- cv.tree(myTree)
optTree <- which.min(cvTree$dev)
bestTree <- cvTree$size[optTree]

#Prune tree based on CV results
pruneTree <- prune.tree(myTree, best = bestTree)

#Use tree to predict original data
treePred <- predict(pruneTree)
treeResid <- resid(pruneTree)
myRange <- range(treePred, data$Y)

#Plot Unpruned Tree
plot(myTree)
text(myTree, cex = .75)
title(main = c("Unpruned Tree for MGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for MGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for MGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "Black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="MGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <-sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  2881
#> R2
#[1]0.75


#****************************************************************************************
# Random Forest
#****************************************************************************************
#Fit random forest and plot
forest.T = randomForest(Y ~ ., data = data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for MGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = ("Pred vs True for Random Forest of MGDG (ng/L)"))
lines(data$Y, data$Y, col = "Black")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <-sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

varImpPlot(forest.T, pch = 20, main = "Importance of VariaMGDGes")
#Most important variaMGDGe is PC2.

#Can also see importance of the variaMGDGes as taMGDGe
importance(forest.T)






############################################
#### Including pacific





data=randomtree_pacific
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variaMGDGe (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "POP", "SiO2", "DGDG", "PO4", "Chla", "pH", "Crypto", "NO3", "NH4", "Light", "FL4_group", "Other", "Bsi", "NP", "N_star", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "O2", "Synecho_low", "Synecho_high", "Synechococcus", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$MGDG
data <- as.data.frame(cbind(Y,Xvar))
data <- data[-c(33), ]

#****************************************************************************************
# CART: Fit tree prune, and scatterplot of observed vs CART estimated.
#****************************************************************************************

myTree <- tree(Y ~ ., data = data, model = T)
#Perform CV on tree object
cvTree <- cv.tree(myTree)
optTree <- which.min(cvTree$dev)
bestTree <- cvTree$size[optTree]

#Prune tree based on CV results
pruneTree <- prune.tree(myTree, best = bestTree)

#Use tree to predict original data
treePred <- predict(pruneTree)
treeResid <- resid(pruneTree)
myRange <- range(treePred, data$Y)

#Plot Unpruned Tree
plot(myTree)
text(myTree, cex = .75)
title(main = c("Unpruned Tree for MGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for MGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for MGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "Black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="MGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <-sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  2881
#> R2
#[1]0.75


#****************************************************************************************
# Random Forest
#****************************************************************************************
#Fit random forest and plot
forest.T = randomForest(Y ~ ., data = data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for MGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = ("Pred vs True for Random Forest of MGDG (ng/L)"))
lines(data$Y, data$Y, col = "Black")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <-sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

varImpPlot(forest.T, pch = 20, main = "Importance of VariaMGDGes")
#Most important variaMGDGe is PC2.

#Can also see importance of the variaMGDGes as taMGDGe
importance(forest.T)



































############################################
#### Relative Abundances





data=randomtree_pacific2
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variaMGDGe (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$MGDG
data <- as.data.frame(cbind(Y,Xvar))

#****************************************************************************************
# CART: Fit tree prune, and scatterplot of observed vs CART estimated.
#****************************************************************************************

myTree <- tree(Y ~ ., data = data, model = T)
#Perform CV on tree object
cvTree <- cv.tree(myTree)
optTree <- which.min(cvTree$dev)
bestTree <- cvTree$size[optTree]

#Prune tree based on CV results
pruneTree <- prune.tree(myTree, best = bestTree)

#Use tree to predict original data
treePred <- predict(pruneTree)
treeResid <- resid(pruneTree)
myRange <- range(treePred, data$Y)

#Plot Unpruned Tree
plot(myTree)
text(myTree, cex = .75)
title(main = c("Unpruned Tree for MGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for MGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for MGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "Black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="MGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <-sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  2881
#> R2
#[1]0.75


#****************************************************************************************
# Random Forest
#****************************************************************************************
#Fit random forest and plot
forest.T = randomForest(Y ~ ., data = data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for MGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = ("Pred vs True for Random Forest of MGDG (ng/L)"))
lines(data$Y, data$Y, col = "Black")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <-sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

varImpPlot(forest.T, pch = 20, main = "Importance of VariaMGDGes")
#Most important variaMGDGe is PC2.

#Can also see importance of the variaMGDGes as taMGDGe
importance(forest.T)






############################################
#### Including pacific





data=randomtree_pacific2
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variaMGDGe (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "POC", "PON", "Dino", "POP", "SiO2", "DGDG", "PO4", "Chla", "pH", "Crypto", "NO3", "NH4", "Light", "FL4_group", "Other", "Bsi", "NP", "N_star", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "O2", "Synecho_low", "Synecho_high", "Synechococcus", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$MGDG
data <- as.data.frame(cbind(Y,Xvar))

#****************************************************************************************
# CART: Fit tree prune, and scatterplot of observed vs CART estimated.
#****************************************************************************************

myTree <- tree(Y ~ ., data = data, model = T)
#Perform CV on tree object
cvTree <- cv.tree(myTree)
optTree <- which.min(cvTree$dev)
bestTree <- cvTree$size[optTree]

#Prune tree based on CV results
pruneTree <- prune.tree(myTree, best = bestTree)

#Use tree to predict original data
treePred <- predict(pruneTree)
treeResid <- resid(pruneTree)
myRange <- range(treePred, data$Y)

#Plot Unpruned Tree
plot(myTree)
text(myTree, cex = .75)
title(main = c("Unpruned Tree for MGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for MGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for MGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "Black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="MGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <-sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  2881
#> R2
#[1]0.75


#****************************************************************************************
# Random Forest
#****************************************************************************************
#Fit random forest and plot
forest.T = randomForest(Y ~ ., data = data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for MGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed MGDG (ng/L)", ylab="Modeled MGDG (ng/L)",
     main = ("Pred vs True for Random Forest of MGDG (ng/L)"))
lines(data$Y, data$Y, col = "Black")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <-sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

varImpPlot(forest.T, pch = 20, main = "Importance of VariaMGDGes")
#Most important variaMGDGe is PC2.

#Can also see importance of the variaMGDGes as taMGDGe
importance(forest.T)


