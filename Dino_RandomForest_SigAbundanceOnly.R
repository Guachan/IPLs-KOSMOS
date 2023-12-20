
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
library(ggpubr)
library(ggpmisc)
#Function to prune tree
source("PruneTree.r")

load_data = function(path){
        raw_lipid = read.csv(paste0(path, 'lipid_abundances_012021.csv')) %>%
                select(-X)
        # pivot on lipid class
        lipid_class = raw_lipid %>%
                select(-'Compound') %>%
                group_by(Depth, Day, Mesocosm, Class) %>%
                summarize(total_L = sum(total_L), .groups='drop_last') %>%
                pivot_wider(names_from='Class', values_from='total_L')
        # pivot on lipid compound
        lipid_compound = raw_lipid %>%
                select(-'Class') %>%
                pivot_wider(names_from='Compound', values_from='total_L')
        phyto = read.csv(paste0(path, 'phytoplankton_abundances.csv'))
        water = read.csv(paste0(path, 'water chemistry.csv'))
        # join data frames
        data = left_join(water, phyto, by=c('Depth', 'Day', 'Mesocosm')) %>%
                left_join(lipid_class, by=c('Depth', 'Day', 'Mesocosm')) %>%
                left_join(lipid_compound, by=c('Depth', 'Day', 'Mesocosm'))
        # add features
        data[is.na(data)] = 0
        data$Depth = as.factor(data$Depth)
        data$Mesocosm = as.factor(data$Mesocosm)
        # add n-to-p ratio
        data = data %>%
                mutate(NP_Ratio=(NO3+NO2+NH4)/PO4)
        return(data)
}


rel_ab <- read.csv('lipid_rel_abundances.csv') %>% 
    filter(percentage >.99) %>% 
    filter(Mesocosm <9) %>% 
    select(Compound) %>% 
    unique() %>% 
    rename(Lipid = Compound)
# read in data
path = "C:/Users/canta/MLRs and GAMs -- KOSMOS/"
data = load_data(path) %>% 
    rename("DG.Total" = DGDG,
           "MG.Total" = MGDG,
           "SQ.Total" = SQ,
           "PC.Total" = PC,
           "PE.Total" = PE,
           "PG.Total" = PG,
           "BL.Total" = BL)
data2 <- data %>% gather(Lipid, Concentration, 29:193) %>% 
    merge(rel_ab, by = "Lipid") %>% 
    spread(Lipid, Concentration)


#Load this only when p

#it masks necssary commands: 
#source("RandomForestPlot.R")
#Load Data - monthly. 2008 to 2015, generally months 4-11

# data <- data[1:39]
#Remove unnessary water quality parameters, the dependent variable (TOC), and dates
drop<-c("Mesocosm", "Depth","Day", "Diat", "Dino", "Synechococcus", "Pelago", "Prym", "Chloro", "Crypto", "Light", "Chla", "Other", "NH4", "NO3", "NO2", "NP_Ratio", "O2", "Temp", "pH", "PO4", "SiO2", "PI-DAG-38:5")
Xvar=data2[,!(names(data2)%in% drop)]
# data=as.data.frame(cbind(Y,Xvar))
Y=data2$Dino
data <- as.data.frame(cbind(Y,Xvar))
names(data) <- gsub(x = names(data), patter = "\\-", replacement = ".")
names(data) <- gsub(x = names(data), patter = "\\:", replacement = ".")
names(data) <- gsub(x = names(data), patter = "PC", replacement = "PC.DAG")
names(data) <- gsub(x = names(data), patter = "PC", replacement = "PC.DAG")
names(data) <- gsub(x = names(data), patter = "DAG.DAG", replacement = "DAG")
names(data) <- gsub(x = names(data), patter = "PC.DAG.Total", replacement = "PC.Total")


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
plot(forest.T, main="Random Forest for Dino (mg/L)")

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
Dino_lipids = importance(forest.T) %>% as.data.frame() %>% 
        rename(RMSE = "%IncMSE") %>% 
        arrange(desc(RMSE)) %>% 
        slice(1:20)
        
dino_names = rownames(Dino_lipids)

    

revised_data <- data[dino_names] %>% 
        cbind(Y)
        

forest.T = randomForest(Y ~ ., data = revised_data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for Dino (mg/L)")

#Plot observed VS predicted
png(filename = "Dino_R2.png", width = 500, height = 500, units = "px")
par(mar=c(5,5,2,2))
plot(revised_data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed Dino Concentration (mg/L)", ylab="Modeled Dino Concentration (mg/L)", cex.lab =1.75, cex.axis = 1.5, cex = 2)
lines(revised_data$Y, revised_data$Y, col = "black")
dev.off()

#How good are predictions?
# R2
R2 <- (cor(revised_data$Y,forest.T.Pred))^2
R2 ## .67
##RMSE
RMSE <- sqrt(mean((revised_data$Y -forest.T.Pred)^2))
RMSE ##3925

varImpPlot(forest.T, pch = 20, main = "Importance of Variables")



### Remove < 2.5% RMSE


Dino_lipids = importance(forest.T) %>% as.data.frame() %>% 
    rename(RMSE = "%IncMSE") %>% 
    filter(RMSE > 2.49) %>% 
    arrange(desc(RMSE))

dino_names = rownames(Dino_lipids)



revised_data <- data[dino_names] %>% 
    cbind(Y)


forest.T = randomForest(Y ~ ., data = revised_data, importance = TRUE)
summary(forest.T)
print(forest.T)
forest.T.Pred = predict(forest.T)
forest.T.res = data$Y - forest.T.Pred  # no residual command so use: res = y - yhat
forest.myRange <- range(forest.T.Pred, data$Y)

#Plot error VS trees
plot(forest.T, main="Random Forest for Dino (mg/L)")

#Plot observed VS predicted
plot(revised_data$Y,  forest.T.Pred, xlim = forest.myRange, ylim = forest.myRange, xlab="Observed TOC Concentration (mg/L)", ylab="Modeled TOC Concentration (mg/L)",
     main = ("Pred vs True for Random Forest of TOC (mg/L)"))
lines(revised_data$Y, revised_data$Y, col = "black")

#How good are predictions?
# R2
R2 <- (cor(revised_data$Y,forest.T.Pred))^2
R2 ## .67
##RMSE
RMSE <- sqrt(mean((revised_data$Y -forest.T.Pred)^2))
RMSE ##3925

varImpPlot(forest.T, pch = 20, main = "Importance of Variables")






revised_data2 <- scale(revised_data) %>% tbl_df()
dino_lm <- lm(Y ~., data = revised_data)

pred_dino <- predict(dino_lm, newdata = data, se.fit = TRUE) %>% as_tibble() %>% 
        bind_cols(data)
formula <-  y ~ x
ggplot(data = pred_dino, aes(fit, Y))+
        geom_point()+
        geom_smooth(method = "lm", se = FALSE)+
        stat_regline_equation(label.y = 40000, aes(label = ..rr.label..), size =5)+
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text',
                        aes(label = paste("p = ", signif(..p.value.., digits = 2), sep = "")),
                        label.y = 38000, size = 5)
