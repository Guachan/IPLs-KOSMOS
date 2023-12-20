#******************************************************************************
#
# Fort Collins - Regression Trees
# Variable Modeled: TOC
# CART and Random Forest
# Kelsey Reeves
#
#******************************************************************************

setwd("C:/Users/Group-SplitX2/Desktop/Research/Locations/Fort Collins/Trees - TOC")

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
        mutate(Total_N = NO3+NO2+NH4)#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$PE
data <- cbind(Y,Xvar)


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
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed Turbidity (NTU)", ylab="Modeled Turbidity (NTU)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="TOC (mg/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <- sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  0.890991
#> R2
#[1]0.83333


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
plot(forest.T, main="Random Forest for TOC (mg/L)")

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
#drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "Pelago", "Synechococcus", "pH", "POC", "POP", "Bsi", "Dino", "MGDG", "SiO2", "PON", "DGDG", "Chla", "NP_Ratio", "Prym", "Nano", "Crypto", "Peuks", "Chains", "NO3", "Other", "Diat", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high", "MikroI", "MikroII", "Chloro")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$PE
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
title(main = c("Unpruned Tree for PE"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for PE"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for PE"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed PE (ng/L)", ylab="Modeled PE (ng/L)",
     main = c("Pred vs True for PE"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="PE (ng/L)")
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
plot(forest.T, main="Random Forest for PE (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed PE (ng/L)", ylab="Modeled PE (ng/L)",
     main = ("Pred vs True for Random Forest of PE (ng/L)"))
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









##########################

data=randomtree_pacific
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$PE
data <- cbind(Y,Xvar)


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
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed Turbidity (NTU)", ylab="Modeled Turbidity (NTU)",
     main = c("Pred vs True for MGDG"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="TOC (mg/L)")
summary(treeFit)

#How good are predictions? -had to place treeFun code in here to do
#R2
R2 <- (cor(data$Y,treePred))^2
R2
#RMSE
RMSE <- sqrt(mean((data$Y -treePred)^2))
RMSE
#> RMSE
#[1]  0.890991
#> R2
#[1]0.83333


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
plot(forest.T, main="Random Forest for TOC (mg/L)")

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
data=randomtree_pacific
#drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "Pelago", "N_star", "pH", "POC", "Bsi", "Dino", "MGDG", "SiO2", "PON", "NH4", "DGDG", "Chla", "NP_Ratio", "Prym", "Nano", "Crypto", "Peuks", "Chains", "Light", "NO3", "Other", "Diat", "NO2", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high", "MikroI", "MikroII", "Chloro")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$PE
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
title(main = c("Unpruned Tree for PE"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for PE"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for PE"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed PE (ng/L)", ylab="Modeled PE (ng/L)",
     main = c("Pred vs True for PE"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="PE (ng/L)")
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
plot(forest.T, main="Random Forest for PE (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed PE (ng/L)", ylab="Modeled PE (ng/L)",
     main = ("Pred vs True for Random Forest of PE (ng/L)"))
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

