##Design for two factors for Incidence and Mortality for Race and Ethnicity types
#data
library(lattice)
library(gridExtra)
library(MASS)
library(ggplot2)
incidence <- c("White"= 604668,"Hispanic"=25561,"Black"=3471,"Asian"=1884)
mortality <- c("White"= 90180,"Hispanic"=190,"Black"=1356,"Asian"=450)

#incidence <- c("Hispanic"=25561,"Black"=3471,"Asian"=1884)
#mortality <- c("Hispanic"=190,"Black"=1356,"Asian"=450)
# 2 dotplots side by side
options(scipen=10000)
first <- dotplot(incidence,
        origin=0,type=c("p","h"),
        xlab="per 100,000 population",
        ylab="Incidence rate",
        scales=list(x=list(tck=0, alternating=FALSE)),
        panel=function(...){
          panel.fill(rgb(.9,.9,.9))
          panel.grid(h=0,v=-1,col="white",lwd=1)
          panel.dotplot(col=rgb(0,.5,1),cex=1.1,...)
        }
)
second <- dotplot(mortality,
        origin=0,type=c("p","h"),
        xlab="per 100,000 population",
        ylab="Mortality rate",
        scales=list(x=list(tck=0, alternating=FALSE)),
        panel=function(...){
          panel.fill(rgb(.9,.9,.9))
          panel.grid(h=0,v=-1,col="white",lwd=1)
          panel.dotplot(col=rgb(0,.5,1),cex=1.1,...)
        }
)
quartz(width=6, height=3)
grid.arrange(first,second,ncol=2)

##--------------------------------------------------------------------------------------
#Principal Component

# 0. Setup
library(XLConnect)
library(MASS)
library(ggplot2)
library(gridExtra)
library(micromapST)
source("fileNameWithheld")

hwThemeSmallText <- hw +
  theme(
    axis.text=element_text(size=rel(.8)),
    axis.title=element_text(size=rel(.8)),
    plot.title=element_text(size=rel(.8))
  )

# 1. Read data and check for missing values
mfcancer <- readWorksheetFromFile(file='fileNameWithheld', 
                                       sheet = 10)

mfcancer <- data.frame(mfcancer)
any(is.na(mfcancer))

rate <- mfcancer[, "Deaths"]
male <- mfcancer[, "Sex"]=="Male"
mat <- cbind(rate[male], rate[!male])
rownames(mat) <- mfcancer$State[1:51]
colnames(mat) <- c("Male","Female")
mat
df <- as.data.frame(mat)
head(df)

# 2. Center the variables 

matCen <- scale(mat,scale=FALSE)
dfCen <- as.data.frame(matCen)
head(dfCen)

# Save the means for reference lines
maleM <- mean(df$Male)
femaleM <- mean(df$Female)

# produce plots

title <-paste("States Cancer Mortality Rates from 2001-2011", sep="")   
p <- ggplot(df,aes(x=Male,y=Female))+
  geom_hline(yint=femaleM,color='blue')+ 
  geom_vline(xint=maleM,color='blue')+
  geom_point(fill="red",shape=21)+
  labs( title=title)+hwThemeSmallText

pCen<- ggplot(dfCen,aes(x=Male,y=Female))+
  geom_hline(yint=0,color='blue')+ 
  geom_vline(xint=0,color='blue')+
  geom_point(fill="red",shape=21)+
  labs(x="Male - Mean(Male)",
       y="Female - mean(Female)",
       title=title)+hwThemeSmallText

quartz(width=6, height=3)
grid.arrange(p,pCen,ncol=2)

# 3. Principal components and rotation 

pc <- prcomp(mat,scale=FALSE)
pcDat <- pc$x # principal components
pcRotate <- pc$rotation  #rotation matrix

det(pcRotate)  # -1 means there is also a reflection 

big <- max(abs(pcDat[,1]))

dfPC <- data.frame(PC1=pcDat[,1],PC2 = pcDat[,2])
rownames(dfPC) = rownames(dfCen)

# Find some extreme points to label
# with state postal codes

id1 <- which.min(dfPC$PC1) #DC
id2 <- which.max(dfPC$PC1) #CA
id3 <- which.min(dfPC$PC2) #Florida
subs <- c(id1,id2,id3)

dfPC$State <- NA
dfPC$State[subs] <- row.names(dfPC)[subs] 

dfCen$State <- NA
dfCen$State[subs] <- row.names(dfPC)[subs]

# Store the two plots and then juxtapose
#
pCen <- ggplot(dfCen,aes(x=Male,y=Female))+
  geom_hline(yint=0,color='blue')+ 
  geom_vline(xint=0,color='blue')+
  geom_point(fill="red",shape=21)+
  ylim(-big,big)+
  xlim(-big,big)+
  geom_text(aes(y=Female-5,label=State),size=4,vjust=1)+
  labs(x="Male - Mean(Male)",
       y="Female - mean(Female)",
       title=title)+hwThemeSmallText

pRot <- ggplot(dfPC,aes(x=PC1,y=PC2))+
  geom_hline(yint=0,color='blue')+ 
  geom_vline(xint=0,color='blue')+
  geom_point(fill="red",shape=21)+
  ylim(-big,big)+
  xlim(-big,big)+
  geom_text(aes(y=-PC2+5,label=State),size=4,vjust=0)+
  labs(x="PC1 for Male and Female Mortality Rates",
       y="PC2 for Male and Female Mortality Rates",
       title=title)+hwThemeSmallText

quartz(width=6, height=3) 
grid.arrange(pCen,pRot,ncol=2)
##------------------------------------------------------------------------------------------
#Linear Regression

source("fileNameWithheld")
library(MASS)
library(ISLR)
library(lattice)
library(hexbin)
library(ggplot2)
library(car)
library(XLConnect)

#1.  Simple Linear Regression


## Run
cancer <- readWorksheetFromFile(file='fileNameWithheld',sheet = 11)
View(cancer) 

## End

#1.2 Fit a single predictor model

## Run

lm.fit <- lm(Male_White~Female_White,data=cancer)
lm.fit

## End

#1.3 Look at the model summary

## Run
summary(lm.fit) 
## End

#1.4 Extracting estimates
## Run  

# component names
names(lm.fit)

coef(lm.fit)
predict(lm.fit)
residuals(lm.fit)

## End

#1.5 Compute confidence intervals

## Run

confint(lm.fit) # 95% is the default
confint(lm.fit,level=.99)

## End

## Run

predict(lm.fit,data.frame(Female_White=(c(5,10,15))),
        interval="confidence")    

predict(lm.fit,data.frame(Female_White=(c(5,10,15))),
        interval="prediction")

## End

#1.6 Some quick base level R graphics

## Run

attach(cancer)

quartz()
plot(Female_White,Male_White)
abline(lm.fit,lwd=3,col="red")

plot(Female_White,Male_White,col="red")
plot(Female_White,Male_White,pch=20)
plot(Female_White,Male_White,pch="+")

## End 

#1.7 Better graphics

## Run

(ggplot(cancer,aes(x=Female_White,y=Male_White))
+ geom_point(shape=21,fill="red",
             color="black",size=2)
+ stat_smooth(method=lm,
              color="blue",fill="cyan")
+ labs(
  x="Mortality rate for Female White population",
  y="Mortality rate for Male White population",
  title="Melanoma Cancer Mortality Rates for White population based on gender") 
+ hw
)
##____________________________ Graphs for Correlation and Covariance ______________________
# 5. Biplot


## Run 

head(pc$x)

pc$rot

## End

quartz(width=6,height=6)
biplot(pc,scale=0,las=1,
       main= "Covariance PC Loadings Top and Right Axes")

quartz(width=6, height=6)
pcCor <- prcomp(mat,scale=TRUE)
biplot(pcCor,scale=1,las=1,
       main="Correlation PC Loading Top and Right Axes")

# 6. Percent of variability represented

pc$sdev

# The screeplot shows variances
screeplot(pc)

# Another choice is to look 
# cumulative percent of total
# variance 
dfPC$State <- NULL
vec <- diag(var(dfPC))
100*cumsum(vec)/sum(vec)  

# 6. Ordering States by the first principal
#component

ord <- with(dfPC,order(PC1))
round( dfPC[ord,1:2],1)

##-----------------------------------------------------------------------------------------
#Micromap to justify PCA
library(micromapST)
cancer_1data <- readWorksheetFromFile(file='fileNameWithheld', 
                                          sheet = 12)
#cancer_1data <- round( dfPC[ord,1:2],1)
subs <- match(cancer_1data$FIPS,stateNamesFips$fips)
row.names(cancer_1data) <- stateNamesFips$ab[subs]

### Sunlight vs UV
panelDesc1 <- data.frame(type = c("mapcum","id","dot"), lab1=c("","","Mortality Rates per 100,000 population")
                         ,lab3=c("","","per 100,000 population"),col1= c(NA,NA,2))
ExtTitle1 <- c("State-wise Mortality for Gender data 2001-2011")
pdf(file = "/Users/divya/PC.pdf", width = 7.5, height = 10)
micromapST(cancer_1data, panelDesc1,sortVar = 2, ascend = TRUE,title = ExtTitle1, 
           plotNames = "ab")
dev.off()

View(cancer_1data)