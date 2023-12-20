
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
#Function to prune tree
source("PruneTree.r")
#Load this only when plotting a randomforest tree because
#it masks necssary commands: 
#source("RandomForestPlot.R")
#Load Data - monthly. 2008 to 2015, generally months 4-11
data=read.csv("lipid_data_cyto2.csv") %>% 
        mutate(Total_N = NO3+NO2+NH4)
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$DGDG
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
title(main = c("Unpruned Tree for DGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for DGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for DGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = c("Pred vs True for PC"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="DGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <- sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1] 3176
#> R2
#[1]0.64


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
plot(forest.T, main="Random Forest for DGDG (mg/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed TOC Concentration (mg/L)", ylab="Modeled TOC Concentration (mg/L)",
     main = ("Pred vs True for Random Forest of TOC (mg/L)"))
lines(data$Y, data$Y, col = "black")

#How good are predictions?
# R2
R2 <- (cor(data$Y,forest.T.Pred))^2
R2
##RMSE
RMSE <- sqrt(mean((data$Y -forest.T.Pred)^2))
RMSE

#Importance of the different variables
#Dotchart of variable importance as measured by a Random Forest
#Shows total decrease in node impurities from splitting on the variable, averaged over
#all trees. For classification, the node impurity is measured by the Gini index. 
#For regression, it is measured by residual sum of squares. Variables are ordered
#by importance were the top most variable is the most important.
#You can use this chart to determine which variables to include in a PCA by looking
#for a large break between the variables, which will help you to decide how many important
#vairables to choose.
varImpPlot(forest.T, pch = 20, main = "Importance of Variables")
#Most important variable is PC2.

#Can also see importance of the variables as table
importance(forest.T)

#Tree visuallization - look at tree  number - load source above
#Can plot a given tree number in the random forest:
tree_num = 350
x=tree_func(forest.T,tree_num)
#Note: I found this function online. Not sure if it works 100% correctly.
#https://shiring.github.io/machine_learning/2017/03/16/rf_plot_ggraph

#Can use getTree to to also visualize tree number 350, since the Error VS
#trees plot had stabilized by then.
getTree(forest.T, k=tree_num, labelVar=TRUE)







#************************************************************************************************************************************
#
#*## All depths and mesocosms rerun random forest with selected variables only
data=read.csv("lipid_data_cyto2.csv")
# data <- data[-c(4,36), ]
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "Crypto", "Diat", "Synechococcus", "MGDG", "Light", "PON", "POC", "SiO2", "PO4", "POP", "NP_Ratio", "Prym", "Temp", "NH4", "DGDG", "NO3", "Other", "Bsi", "NO2", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$DGDG
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
title(main = c("Unpruned Tree for DGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for DGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for DGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = c("Pred vs True for DGDG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="DGDG (ng/L)")
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
plot(forest.T, main="Random Forest for DGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = ("Pred vs True for Random Forest of DGDG (ng/L)"))
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
































###########################################
data=randomtree_pacific
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$DGDG
data <- as.data.frame(cbind(Y,Xvar))
data <- data[-c(56,75), ]

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
title(main = c("Unpruned Tree for DGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for DGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for DGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = c("Pred vs True for PC"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="DGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <- sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1] 3176
#> R2
#[1]0.64


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
plot(forest.T, main="Random Forest for DGDG (mg/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed TOC Concentration (mg/L)", ylab="Modeled TOC Concentration (mg/L)",
     main = ("Pred vs True for Random Forest of TOC (mg/L)"))
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







#************************************************************************************************************************************
#
#*## All depths and mesocosms rerun random forest with selected variables only
data=randomtree_pacific
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "Diat", "SiO2", "PON", "POP", "N_star", "PO4", "NP", "Prym", "Temp", "NH4", "DGDG", "Crypto", "Light", "NO3", "Other", "Bsi", "NO2", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$DGDG
data <- as.data.frame(cbind(Y,Xvar))
data <- data[-c(56, 75), ]

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
title(main = c("Unpruned Tree for DGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for DGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for DGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = c("Pred vs True for DGDG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="DGDG (ng/L)")
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
plot(forest.T, main="Random Forest for DGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = ("Pred vs True for Random Forest of DGDG (ng/L)"))
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




























###########################################
data=randomtree_pacific2
# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$DGDG
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
title(main = c("Unpruned Tree for DGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for DGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for DGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = c("Pred vs True for PC"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="DGDG (ng/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <- sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1] 3176
#> R2
#[1]0.64


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
plot(forest.T, main="Random Forest for DGDG (mg/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed TOC Concentration (mg/L)", ylab="Modeled TOC Concentration (mg/L)",
     main = ("Pred vs True for Random Forest of TOC (mg/L)"))
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







#************************************************************************************************************************************
#
#*## All depths and mesocosms rerun random forest with selected variables only
data=randomtree_pacific2
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "Diat", "SiO2", "Synecho_all", "PON", "POP", "N_star", "PO4", "NP", "Prym", "Temp", "NH4", "DGDG", "Crypto", "Light", "NO3", "Other", "Bsi", "NO2", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$DGDG
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
title(main = c("Unpruned Tree for DGDG"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for DGDG"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for DGDG"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = c("Pred vs True for DGDG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="DGDG (ng/L)")
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
plot(forest.T, main="Random Forest for DGDG (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed DGDG (ng/L)", ylab="Modeled DGDG (ng/L)",
     main = ("Pred vs True for Random Forest of DGDG (ng/L)"))
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

