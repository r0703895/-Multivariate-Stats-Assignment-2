### Part 0: Setup ###
rm(list=ls())
dev.off()
load(file=file.choose())
library(smacof)
library(stringr)
library(corrplot)

### Part 1: MDS for different measurement levels ###
m1<-smacofSym(delta=dissim, ndim=2, type="ratio", init="torgerson") #ratio
m2<-smacofSym(delta=dissim, ndim=2, type="interval", init="torgerson") #interval
m3<-smacofSym(delta=dissim, ndim=2, type="ordinal", init="torgerson") #ordinal
m4<-smacofSym(delta=dissim, ndim=2 ,type="mspline",spline.degree =4 , spline.intKnots = 4, init="torgerson") #spline

### Part 2: Selecting best analysis and judging fit and stability ###
#STRESS VALUES
round(c(m1$stress,m2$stress,m3$stress,m4$stress),3) #stress-1 values
mt<-smacofSym(delta=dissim, ndim=2 ,type="mspline",spline.degree =40 ,
              spline.intKnots = 40, init="torgerson") #Maybe tuning the spline parameters can improve the model
round(c(mt$stress),3)
#STRESS NORM TEST (Appendix I)
set.seed(1) 
rstress<-matrix(rep(0, len=8), nrow = 2)
colnames(rstress)<-c("ratio","interval","ordinal","mspline")
rownames(rstress)<-c("mean", "sd")
ratio <-randomstress(n=36,ndim=2,nrep=300,type="ratio")
interval <-randomstress(n=36,ndim=2,nrep=300,type="interval")
ordinal <-randomstress(n=36,ndim=2,nrep=300,type="ordinal")
mspline <-randomstress(n=36,ndim=2,nrep=300,type="mspline")
u<-seq(from = m3$stress, to = 0.5, by=0.001)
matplot(x = u,
        y = cbind(
          dnorm(x=u, mean=mean(ratio), sd=sd(ratio)),
          dnorm(x=u, mean=mean(interval), sd=sd(interval)),
          dnorm(x=u, mean=mean(ordinal), sd=sd(ordinal)),
          dnorm(x=u, mean=mean(mspline), sd=sd(mspline))
        ),
        type = "l",
        lty = 1,
        col = c("red", "blue", "green", "orange"),
        xlab = "stress values",
        ylab = "densities",
        main = "Stress Norm Test")
text(x=mean(ratio), y=50, labels="Ratio~N(0.47,0.01)", cex=0.8, col="red")
text(x=mean(interval)+0.03, y=50, labels="Interval~N(0.36,0.010)", cex=0.8, col="blue")
text(x=mean(ordinal)-0.03, y=80, labels="Ordinal~N(0.35,0.0040)", cex=0.8, col="green")
text(x=mean(mspline-0.03), y=100, labels="Spline~N(0.36,0.0036)", cex=0.8, col="orange")
abline(v=m1$stress, col="red")
text(x=m1$stress, y=10, labels=round(m1$stress,3), cex=0.8, col="red")
abline(v=m2$stress, col="blue")
text(x=m2$stress, y=10, labels=round(m2$stress,3), cex=0.8, col="blue")
abline(v=m3$stress, col="green")
text(x=m3$stress+0.01, y=10, labels=round(m3$stress,3), cex=0.8, col="green")
abline(v=m4$stress, col="orange")
text(x=m4$stress, y=15, labels=round(m4$stress,3), cex=0.8, col="orange")
#PERMUTATION TESTS
rstress.perm<-matrix(rep(0, len=8), nrow = 2) 
colnames(rstress.perm)<-c("ratio","interval","ordinal","mspline")
rownames(rstress.perm)<-c("mean", "sd")
ratio.perm<-permtest(m1,nrep=100)
interval.perm<-permtest(m2,nrep=100)
ordinal.perm<-permtest(m3,nrep=100)
spline.perm<-permtest(m4,nrep=100)
u<-seq(from = m3$stress, to = 0.5, by=0.001)
matplot(x = u,
        y = cbind(
          dnorm(x=u, mean=mean(ratio.perm$stressvec), sd=sd(ratio.perm$stressvec)),
          dnorm(x=u, mean=mean(interval.perm$stressvec), sd=sd(interval.perm$stressvec)),
          dnorm(x=u, mean=mean(ordinal.perm$stressvec), sd=sd(ordinal.perm$stressvec)),
          dnorm(x=u, mean=mean(spline.perm$stressvec), sd=sd(spline.perm$stressvec))
        ),
        type = "l",
        lty = 1,
        col = c("red", "blue", "green", "orange"),
        xlab = "stress values",
        ylab = "densities",
        main = "Permutation Test")
text(x=mean(ratio.perm$stressvec)+0.04, y=50, labels="Ratio~N(0.36,0.0036)", cex=0.8, col="red")
text(x=mean(interval.perm$stressvec)+0.04, y=70, labels="Interval~N(0.36,0.0035)", cex=0.8, col="blue")
text(x=mean(ordinal.perm$stressvec)-0.03, y=80, labels="Ordinal~N(0.35,0.0037)", cex=0.8, col="green")
text(x=mean(spline.perm$stressvec)+0.04, y=100, labels="Spline~N(0.36,0.0035)", cex=0.8, col="orange")
abline(v=m1$stress, col="red")
text(x=m1$stress, y=10, labels=round(m1$stress,3), cex=0.8, col="red")
abline(v=m2$stress, col="blue")
text(x=m2$stress, y=10, labels=round(m2$stress,3), cex=0.8, col="blue")
abline(v=m3$stress, col="green")
text(x=m3$stress+0.01, y=10, labels=round(m3$stress,3), cex=0.8, col="green")
abline(v=m4$stress, col="orange")
text(x=m4$stress, y=15, labels=round(m4$stress,3), cex=0.8, col="orange")
#RESIDUAL AND SHEPHARD PLOTS (for ordinal and spline)
par(mfrow=c(2,2))
plot(m3,plot.type="resplot",main="residual plot ordinal MDS")
plot(m3,plot.type="Shepard",main="Shepard diagram ordinal MDS")
plot(m4,plot.type="resplot",main="residual plot spline MDS")
plot(m4,plot.type="Shepard",main="Shepard diagram spline MDS") 
par(mfrow=c(1,1))
#JACKKNIFE PLOTS (for ordinal and spline)
jack<-cbind(jackmds(m1)$stab,jackmds(m2)$stab,jackmds(m3)$stab,jackmds(m4)$stab)
colnames(jack)<-c("Ratio", "Interval", "Ordinal", "Spline")
par(mfrow=c(1,2))
plot(jackmds(m3),xlim=c(-1.2,1.2),ylim=c(-1,1), main="Ordinal")
text(x=1,y=0.8,cex=1, labels=round(jack[3],5))
plot(jackmds(m4),xlim=c(-1.2,1.2),ylim=c(-1,1), main="Spline")
text(x=1,y=0.8,cex=1, labels=round(jack[4],5))
#SCREE PLOT
Scree<-matrix(rep(0, len=40), nrow = 10)
colnames(Scree)<-c("Ratio", "Interval", "Ordinal", "Spline (d=4)")
for (i in 1:10) {
  #ratio
  Scree[i,1]<-smacofSym(delta=dissim, ndim=i, type="ratio", init="torgerson")$stress
  #interval
  Scree[i,2]<-smacofSym(delta=dissim, ndim=i, type="interval", init="torgerson")$stress
  #ordinal
  Scree[i,3]<-smacofSym(delta=dissim, ndim=i, type="ordinal", init="torgerson")$stress
  #spline
  Scree[i,4]<-smacofSym(delta=dissim, ndim=i ,type="mspline",spline.degree =4 ,
                        spline.intKnots = 4, init="torgerson")$stress
}
plot(Scree[,"Ratio"], col="red", ylab="Stress-1", xlab="Dimensions", main="Scree plot for MDS")
points(Scree[,"Interval"], col = "blue")
points(Scree[,"Ordinal"], col = "green")
points(Scree[,"Spline"], col = "orange")
abline(a=0.05,b=0)
abline(a=0.10, b=0)
abline(a=0.20, b=0)
text(x=2, y=0.22, labels="poor")
text(x=2, y = 0.12, labels = "fair")
text(x=2, y = 0.06, labels = "good")
legend(x="topright",legend=c("Ratio", "Interval", "Ordinal", "Spline"), col=c("red", "blue", "green", "orange"), pch=15)

### Part 3: Anaysis of the solution and externall variables
#CONFIGURATION PLOT
plot(m3,plot.type="conf",main="Morse (ordinal)")
#EXTERNAL VARIABLES
print(beeps<-nchar(colnames(dissim))-str_count(colnames(dissim)," ")) #Amount of Beeps
print(first<-str_sub(colnames(dissim),1,1)) #The first symbol (which we translate into binary) Appendix III
    print(firstbin<-factor(first, levels=c("·","-"), labels=c(0, 1)))
    firstbin[9] <- "0"
print(last<-str_sub(colnames(dissim),-1,-1)) #The last symbol (again translated into binary) Appendix III
    print(lastbin<-factor(last, levels=c("·","-"), labels=c(0, 1)))
    lastbin[9] <- "0"
print(dashes<-str_count(colnames(dissim), pattern = "-")) #The amount of dashes in the signal
dots<-str_count(colnames(dissim), pattern = "·") #The amount of dots in the signal Appendix III
    dots[9] <- 2
print(pauses<-beeps-1) #The amount of pauses in the signal Appendix IV
print(length<-0.05*dots+0.15*dashes+0.05*pauses) #The length of the signal (in time)
print(dotprop<-dots/beeps) #The proportion of dots in the signal
print(dashprop<-dashes/beeps) #The proportion of dashes in the signal
#CORRELATION BETWEEN EXTERNAL VARIABLES
Y<-cbind(beeps, 
         as.numeric(levels(firstbin))[firstbin], 
         as.numeric(levels(lastbin))[lastbin], 
         dashes,
         dots,
         length,
         dotprop 
         #dashprop,
         #pauses
)
colnames(Y)[c(2,3)]<-c("first","last")
cormatY<-cor(Y)
corrplot(cormatY, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
#THE BIPLOT
biCon <- biplotmds(m3, Y)   
coef(biCon)
plot(biCon, main = "Biplot Vector Representation", vecscale = 0.8, 
     xlim = c(-1.5, 1.5), vec.conf = list(col = "brown"), pch = 20, cex = 0.5)
cormat<-t(cor(m3$conf,Y))

### APPENDICES ###
#Appendix I 
shapiro.test(ratio)
shapiro.test(interval)
shapiro.test(ordinal)
shapiro.test(mspline)
#Appendix II
shapiro.test(ratio.perm$stressvec)
shapiro.test(interval.perm$stressvec)
shapiro.test(ordinal.perm$stressvec)
shapiro.test(spline.perm$stressvec)
#Appendix III
    #Now notice that the both firstbin and lastbin have an NA value.
    #If you look at
print(colnames(dissim))
    #You'll see that the dot-symbols used are not regular dots ".",
    #but some weird floating dot "·". Entry number 9 for some 
    #reason does use regular dots ".", so the factor-function did not
    #recognize it. So we'll just change that value manually.
#Appendix IV
    #Note that I do not use the same computation for the spaces as for variable #1.
    #The data is a bit messy, some names have spaces between the symbols and some don't.
    #For variable #1, we just need to get rid of actual spaces between symbols.
    #For variable #6, we need to know the amount of real life pauses in the signal,
    #which will always be the amount of beeps minus one.
