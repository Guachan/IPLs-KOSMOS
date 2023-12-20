
#For Trees
library(mclust)
library(tree)
library(randomForest)
#For PCA
library(maps)
library(akima)
library(fields)




## All depths and mesocosms
data=read.csv("lipid_data_cyto2.csv") %>% 
        mutate(Total_N = NO3+NO2+NH4) %>% 
        mutate(Si_P = SiO2/PO4)
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
#drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other", "FL4_group", "Synecho_all", "Chains", "Crypto_all", "Prym", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$BL
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
title(main = c("Unpruned Tree for BL"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for BL"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for BL"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = c("Pred vs True for BL"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="BL (ng/L)")
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
plot(forest.T, main="Random Forest for BL (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = ("Pred vs True for Random Forest of BL (ng/L)"))
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













#************************************************************************************************************************************
#
#*## All depths and mesocosms rerun random forest with selected variables only
data=read.csv("lipid_data_cyto2.csv") %>% 
        mutate(Total_N = NO3+NO2+NH4) %>% 
        mutate(Si_P = SiO2/PO4)
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Light", "NO3", "Other", "Bsi", "Diat", "NP_Ratio", "NO2", "NH4", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "O2", "PO4", "Chains", "FL4_group", "Synecho_all", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "PON", "SiO2", "Dino", "NP_Ratio", "Total_N", "POP", "POC", "SQ", "Chla", "MGDG", "PO4", "POP", "Prym", "DGDG", "Light", "NO3", "Other", "Bsi", "Diat", "O2", "Synecho_all", "NO2", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "Chains", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$BL
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
title(main = c("Unpruned Tree for BL"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for BL"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for BL"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = c("Pred vs True for BL"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="BL (ng/L)")
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
plot(forest.T, main="Random Forest for BL (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = ("Pred vs True for Random Forest of BL (ng/L)"))
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


## All depths and mesocosms
data=randomtree_pacific
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
#drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other", "FL4_group", "Synecho_all", "Chains", "Crypto_all", "Prym", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$BL
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
title(main = c("Unpruned Tree for BL"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for BL"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for BL"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = c("Pred vs True for BL"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="BL (ng/L)")
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
plot(forest.T, main="Random Forest for BL (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = ("Pred vs True for Random Forest of BL (ng/L)"))
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













#************************************************************************************************************************************
#
#*## All depths and mesocosms rerun random forest with selected variables only
data=randomtree_pacific
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Light", "NO3", "Other", "Bsi", "Diat", "NP_Ratio", "NO2", "NH4", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "O2", "PO4", "Chains", "FL4_group", "Synecho_all", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "Chla", "MGDG", "PO4", "POP", "Prym", "NP", "N_star", "DGDG", "Light", "NO3", "Other", "Bsi", "Diat", "Crypto", "O2", "Synecho_all", "NP_Ratio", "NO2", "NH4", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "Chains", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$BL
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
title(main = c("Unpruned Tree for BL"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for BL"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for BL"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = c("Pred vs True for BL"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="BL (ng/L)")
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
plot(forest.T, main="Random Forest for BL (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = ("Pred vs True for Random Forest of BL (ng/L)"))
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
#Relative abundances
############################################

data=randomtree_pacific2
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other")
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
#drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Other", "FL4_group", "Synecho_all", "Chains", "Crypto_all", "Prym", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$BL
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
title(main = c("Unpruned Tree for BL"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for BL"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for BL"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = c("Pred vs True for BL"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="BL (ng/L)")
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
plot(forest.T, main="Random Forest for BL (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = ("Pred vs True for Random Forest of BL (ng/L)"))
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
# drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "MGDG", "DGDG", "Light", "NO3", "Other", "Bsi", "Diat", "NP_Ratio", "NO2", "NH4", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "O2", "PO4", "Chains", "FL4_group", "Synecho_all", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Chloro", "Pelago")
drop<-c("Mesocosm", "Depth","Day", "BL", "PG", "PC", "PE", "SQ", "Chla", "Dino", "MGDG", "POC", "PON", "pH", "PO4", "POP", "Prym", "NP", "DGDG", "Light", "NO3", "Other", "Bsi", "Diat", "Crypto", "O2", "NP_Ratio", "NO2", "NH4", "OrgN_OrgP_Ratio", "OrgC_OrgN_Ratio", "Chains", "FL4_group", "Crypto_all", "Synecho_low", "Fl4_group", "Synecho_high","Synechococcus", "Nano", "MikroI", "MikroII", "Pelago")
Xvar=data[,!(names(data)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data$BL
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
title(main = c("Unpruned Tree for BL"))

#Plot CV
plot(cvTree$size, cvTree$dev, type = "b", 
     main = c("Cross Validation for BL"))

#Plot Prunned Tree
plot(pruneTree)
text(pruneTree, cex = .75)
title(main =c("Pruned Tree for BL"))

#Plot Observed vs Predicted
plot(data$Y, treePred, xlim = myRange, ylim = myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = c("Pred vs True for BL"))
lines(data$Y, data$Y, col = "black")



#Fit tree
tree.T = tree(Y ~ ., data =data)
summary(tree.T)
#Use treeFun to prune the tree, predict from the tree, and get plots
treeFit = treeFun(data, toPlot=T, title="BL (ng/L)")
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
plot(forest.T, main="Random Forest for BL (ng/L)")

#Plot observed VS predicted
plot(data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed BL (ng/L)", ylab="Modeled BL (ng/L)",
     main = ("Pred vs True for Random Forest of BL (ng/L)"))
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

