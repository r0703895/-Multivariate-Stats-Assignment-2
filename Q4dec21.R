####### QUESTION 4 #######

###### DATA SETUP ######
rm(list=ls())
dev.off()
load(file.choose())
library(ca)
library(ggplot2)
# a<-setwd("C:/Users/mbrus/Desktop/2021 Sem 1/MVS/Assignment 2/") 
# load("C:/Users/mbrus/Desktop/2021 Sem 1/MVS/Assignment 2/datacar(1)")

###### EXPLORATION ######
chi2<-chisq.test(datacar)

###### CORRESPONDANCE ANALYSIS #######
print(ca<-ca(datacar))
casum<-summary(ca) #only for two dimensions
# #Checks:
round(sum(ca$sv^2),7)==round(chi2$statistic/sum(datacar),7) #total inertia is sum of sv's
round(sum(ca$sv^2),7)==round(sum(ca$rowinertia),7) #total inertia is sum of row inertias
round(sum(ca$sv^2),7)==round(sum(ca$colinertia),7) #total inertia is sum of column inertias
#There are 9 possible dimensions, = min(10-1,15-1)
#The first two dimensions explain:
sum(ca$sv[1:2]^2)/sum(ca$sv^2) #95% of the total inertia in the dataset
#Scree plot: 
plot(ca$sv^2/sum(ca$sv^2), main="Scree plot", ylab="% of total inertia explained", xlab="dimension")
lines(ca$sv^2/sum(ca$sv^2)) #There's a very clear elbow
#Quality of representation in the two dimensions (~communality)
par(mfrow=c(1,1))
barplot(casum$rows$` qlt`/1000, names.arg = lapply(strsplit(rownames(datacar), " "), `[[`, 1), las=2, main= "Qualities of rows (car models) [%]") #all rows seem pretty well represented except that one
barplot(casum$columns$` qlt`/1000, names.arg=colnames(datacar), las=2, cex.names=0.8, main="Qualities of columns (attributes) [%]")
#Interpreting the dimensions using corrs (=how much of the inertia in a point is explained by a specific dimension)
par(mfrow=c(1,2))
rowcorrs<-as.data.frame(
          cbind(
            c(
              casum[["rows"]][[6]]/1000,
              casum[["rows"]][[9]]/1000
              ),
            rep(unlist(lapply(strsplit(rownames(datacar), " "), `[[`, 1), use.names=FALSE),2),
            c(
              rep(1, 10),
              rep(2, 10)
              )
            ))
ggplot(data = rowcorrs, aes(x=as.character(V3), y=V2,fill=as.numeric(V1))) + 
  geom_tile()+labs(y= "squared row correlation", x = "Dimension", fill="corr^2")
colcorrs<-as.data.frame(
  cbind(
    c(
      casum[["columns"]][[6]]/1000,
      casum[["columns"]][[9]]/1000
    ),
    rep(colnames(datacar),2),
    c(
      rep(1, 15),
      rep(2, 15)
    )
  ))
ggplot(data = colcorrs, aes(x=as.character(V3), y=V2,fill=as.numeric(V1))) + 
  geom_tile()+labs(y= "squared column correlation", x = "Dimension", fill="corr^2")

#Note that these should hold but sometimes don't and idk why
# casum$rows$` qlt` == casum[["rows"]][[6]] + casum[["rows"]][[9]]
# casum$columns$` qlt` == casum[["columns"]][[6]] + casum[["columns"]][[9]]
par(mfrow=c(1,1))
barplot(casum$rows$mass, main="Row (car model) masses",names.arg = lapply(strsplit(rownames(datacar), " "), `[[`, 1), las=2 )
barplot(casum$columns$mass, main="Column (attributes) masses",names.arg=colnames(datacar), las=2,cex.names=0.8 )
#We can safely interpret the contribution values. If one row/column was disproportionately represented in the dataset
#(had a high mass), then a dimension could be disproportionately determined by that row/column.
#But it's not so that's good. The contributions:  
rowcons<-as.data.frame(
  cbind(
    c(
      casum[["rows"]][[7]]/1000,
      casum[["rows"]][[10]]/1000
    ),
    rep(unlist(lapply(strsplit(rownames(datacar), " "), `[[`, 1), use.names=FALSE),2),
    c(
      rep(1, 10),
      rep(2, 10)
    )
  ))
ggplot(data = rowcons, aes(x=as.character(V3), y=V2,fill=as.numeric(V1))) + 
  geom_tile()+labs(y= "Row contributions", x = "Dimension", fill="cont")
colcons<-as.data.frame(
  cbind(
    c(
      casum[["columns"]][[7]]/1000,
      casum[["columns"]][[10]]/1000
    ),
    rep(colnames(datacar),2),
    c(
      rep(1, 15),
      rep(2, 15)
    )
  ))
ggplot(data = colcons, aes(x=as.character(V3), y=V2,fill=as.numeric(V1))) + 
  geom_tile()+labs(y= "Column contributions", x = "Dimension", fill="corr^2")


########## BIPLOT ###########
library(factoextra)
fviz_ca_biplot(ca, repel = T,arrows = c(T, F), map="colprincipal")
fviz_ca_biplot(ca, repel = T,arrows = c(F, T), map="rowprincipal")

