setwd("P:/Subjects/SpatialStats/Project")
library(sp)
library(BayesX)
library(R2BayesX)


map1<-readRDS("DEU_adm1.rds")
kar1 = read.table(file = "yield.csv", header = T,
                  sep = ",",dec = ".",na.strings = ".")

#*********************************************************************
#hardest way to draw a map
#*********************************************************************
#normalizing potato yield for shades of gray
kar1$col=(kar1$Kartoffeln-min(kar1$Kartoffeln,na.rm = TRUE))/
  (max(kar1$Kartoffeln,na.rm = TRUE)-min(kar1$Kartoffeln,na.rm = TRUE))
#simplifying map1 definitions
map1Reg=map1["CCA_1"]
#converting region codes from string to number
map1Reg$CCA_1=as.numeric(map1Reg$CCA_1)
#creating color codes for each region according to their order in map
regCol=rep("",length(map1Reg$CCA_1))
for(i in 1:length(map1Reg$CCA_1)){
  entry=kar1$col[which(kar1$StateNo==map1Reg$CCA_1[i])]
  if(is.na(entry)){
    regCol[i]="red"
  }else{
    regCol[i]=gray(entry)
  }
}
#plotting regional(state) map of germany with corresponding potato yield
#from black to white potato yield increases, no data for red areas
plot(map1Reg,col = regCol, border = 'darkgrey')


#****************************************************************
#a simpler way for drawing the same map
#****************************************************************

map1Bnd=sp2bnd(spObject = map1,
               regionNames = as.character(as.numeric(map1$CCA_1)))


drawmap(data=kar1[which(!is.na(kar1$Kartoffeln)),],
        map=map1Bnd,main="DEU Kartoffeln",
        regionvar=1,plotvar=3)

# normalize the yield values
kar1$zscore=(kar1$Kartoffeln-mean(kar1$Kartoffeln,na.rm = TRUE))/
  sd(kar1$Kartoffeln,na.rm = TRUE)

drawmap(data=kar1[which(!is.na(kar1$zscore)),],
        map=map1Bnd,main="DEU Kartoffeln",
        regionvar=1,plotvar=4)

#***************************************************************
# districts of Germany
#***************************************************************

library(MBA)


data("GermanyBnd")
data2 = read.table(file = "yield2.csv", header = T,
                  sep = ",",dec = ".",na.strings = c(".","-"))

data2=data2[-which(is.na(data2$kartoffeln)),]
data2=data2[-which(data2$kartoffeln<1),]
data2$stateName=factor(data2$stateName)

# normalize the yield values
data2$zscore=(data2$kartoffeln-mean(data2$kartoffeln))/
  sd(data2$kartoffeln)

R2BayesX::plotmap(map = GermanyBnd, x = data2, c.select = 4,
                  type = "mba",interp = FALSE, extrap = TRUE,
                  legend = TRUE, missing = TRUE, 
                  names = FALSE, values = TRUE,
                  cex.names = 0.5, cex.values = 0.5,
                  digits = 2,
                  main="German Potato Yields")

# (i) determine the precision matrix via neigborhood structure

GermanyGraph=bnd2gra(GermanyBnd)
GermanyGraph[1:10,1:10] # looks like a precision matrix

#check properties of pmat
#names?
countyCodes=colnames(GermanyGraph)

#Ruegen,the only district of Germany which consists solely of islands.
#Neighbor should be attained to this district
GermanyGraph[362,362] # equals to zero meaning that no neighbor
countyCodes[362] #county with code 13061 has no neighbors
which(countyCodes=="13057")
countyCodes[358]
GermanyGraph[362,362]=1
GermanyGraph[362,358]=-1
GermanyGraph[358,362]=-1
GermanyGraph[358,358]=GermanyGraph[358,358]+1

#Symetric?
sqrt(sum(GermanyGraph==t(GermanyGraph)))
439*439 # number of graphs

# rowsum = 0 ?
sum(rowSums(GermanyGraph))
sum(colSums(GermanyGraph))

# (ii)
# Load the data
#deleting the counties which are not represented by the map
for(i in nrow(data2):1){
  if(!any(data2$stateCode[i]==countyCodes)){
    data2=data2[-i,]
  }
}
data2$stateName=factor(data2$stateName)

n <- nrow(data2) # number of observed yields


# Compute design matrix
Z <- matrix(0,n,439)
for(j in 1:439) {
  Z[which(data2$stateCode==countyCodes[j]),j] <- 1
}
sum(rowSums(Z))
sum(colSums(Z))

# (iii) Model estimation
# penalizied LS-Estimation for several values of lambda

lambda <- 1
beta1 <- as.vector(solve(t(Z)%*%Z+lambda*GermanyGraph)%*%
                     t(Z)%*%data2$zscore)
lambda <- 3
beta3 <- as.vector(solve(t(Z)%*%Z+lambda*GermanyGraph)%*%
                     t(Z)%*%data2$zscore)
lambda <- 5
beta5 <- as.vector(solve(t(Z)%*%Z+lambda*GermanyGraph)%*%
                     t(Z)%*%data2$zscore)

# Plots
res <- data.frame(countyCodes,beta1,beta3,beta5)

R2BayesX::plotmap(map = GermanyBnd, x = res, c.select = 2,
                  legend = TRUE, missing = TRUE, 
                  names = FALSE, values = TRUE,
                  cex.names = 0.5, cex.values = 0.5,
                  digits = 2,
                  main="Yields smoothed at lambda=1")

R2BayesX::plotmap(map = GermanyBnd, x = res, c.select = 3,
                  legend = TRUE, missing = TRUE, 
                  names = FALSE, values = TRUE,
                  cex.names = 0.5, cex.values = 0.5,
                  digits = 2,
                  main="Yields smoothed at lambda=3")

R2BayesX::plotmap(map = GermanyBnd, x = res, c.select = 4,
                  legend = TRUE, missing = TRUE, 
                  names = FALSE, values = TRUE,
                  cex.names = 0.5, cex.values = 0.5,
                  digits = 2,
                  main="Yields smoothed at lambda=5")



x=data2$kartoffeln
x=data2$zscore
h<-hist(x, breaks=20, col="red",main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)
rm(x,h,xfit,yfit)















