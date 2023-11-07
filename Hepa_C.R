setwd("C:/Users/48507/Desktop/MGR_data_science/SEM 2/AS_project")
set.seed(123)
library(ggplot2)
library(mice)
library(car)
library(grid)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(JointAI)
library(EnvStats)
library(mrfDepth)
library(DMwR2)
library(reshape2)
library(car)
library(stats)
library(stats)
library(caret)
library(ltm)
library(lmtest)
library(mlr3measures)

data<- read.csv("hcvdat0.csv")
data=data[,-1]

missing_values<-sum(is.na(data))
perc_na <- (missing_values/(dim(data)[1]*dim(data)[2]))*100
#female<-data[data$Sex=='f',] # 38.7%
#male<-data[data$Sex=='m',] # 61.3%
#os<-data[data$Category=="0s=suspect Blood Donor",] # 1,138211#
#hepatisis<-data[data$Category=="1=Hepatitis",] # 3.9%
#fibrosis<-data[data$Category=="2=Fibrosis",] # 3.4%
#cirhhosis<-data[data$Category=="3=Cirrhosis",] # 4.8%
``
#heat.dat = data
#heat.dat[!is.na(heat.dat)]<- 0
#heat.dat[is.na(heat.dat)] <- 1
#heat.dat$Sex=as.numeric(heat.dat$Sex)
#heatmap(t(as.matrix(heat.dat[,-1])),scale="none",Rowv=NA,Colv=NA,col=c("azure3","red"),xlab = "Observations", ylab = "Variable",main="Missing data")


data%>%summarise_all(funs(sum(is.na(.))))
md_pattern(data,rotate.names = TRUE,color=c("azure3","tomato"),ylab="Number of observations per pattern",
           xlab="Number of missing values",border="white")

miss_data<-data[rowSums(is.na(data))>0,]


na<-mice(data,m=6,defaultMethod="pmm",maxit=5)
new_data=complete(na)
new_data[,2]=data[,3];new_data[,3]=data[,2];
colnames(new_data)[2]="Sex"
colnames(new_data)[3]="Age"
new_data$Sex <- ifelse(new_data$Sex=='m','male','female')

##### BAR PLOT #####
new_data$Category = gsub(".*\\=","",data$Category)

theme_set(theme_classic())
ggplot(new_data,aes(x=Category,fill=Sex))+
  geom_bar(position="stack",color="white",width=0.8)+
  labs(title=paste("Barplot for each group of patients"),x="Category",y=paste("Quantity"))+
  theme(plot.title=element_text(size=18,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=11.5, angle=45, hjust=1),
        axis.text.y =element_text(size=11.5))+
  scale_fill_manual(values=c("female"="lightpink1","male"="cadetblue3"))

##### Merging groups ####
new_data[new_data$Category=="suspect Blood Donor",1]="Blood Donor"
new_data[new_data$Category=="Blood Donor",1]="Blood Donor"
new_data[new_data$Category=="Hepatitis",1]="Hepatitis"
new_data[new_data$Category=="Fibrosis",1]="Hepatitis"
new_data[new_data$Category=="Cirrhosis",1]="Hepatitis"

new_data$Category=as.factor(new_data$Category)

theme_set(theme_classic())
ggplot(new_data,aes(x=Category,fill=Sex))+
  geom_bar(position="dodge",color="white",width=0.8)+
  labs(title=paste("Barplot for each group of patients"),x="Category",y=paste("Quantity"))+
  theme(plot.title=element_text(size=18,face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=11.5, angle=45, hjust=1),
        axis.text.y =element_text(size=11.5))+
  scale_fill_manual(values=c("female"="lightpink1","male"="cadetblue3"))


#write.xlsx(new_data,"new_data.xlsx",row.names=TRUE,col.names=TRUE)
setwd("C:/Users/48507/Desktop/AS_project")
#save(new_data, file="H_data.RData")
load("H_data.RData")

names=colnames(new_data)
theme_set(theme_classic())

##### BOXPLOTS ####
for (i in 3:length(new_data)){
   
  pdf(paste("boxplot",i,".pdf",sep=""),width=11,height=6)

    ggplot(new_data,aes(y=new_data[,i],x=Category,group=Category,fill=Category,color=Category))+
    geom_boxplot(size=1)+
    #geom_jitter(aes(colour=Sex),position=position_jitter(0.25),alpha=0.7)+
    labs(title=paste("",names[i], ""),x="Category",y="Values")+
    guides(color=guide_legend(override.aes=list(size=3.5)))+
      theme(plot.title=element_text(size=18,face="bold",hjust=0.5),
            axis.title=element_text(size=12),
            legend.position="none")+
      scale_fill_manual(values=c("Blood Donor"="lightblue","Hepatitis"="lightcoral"))+
    scale_color_manual(values=c("Blood Donor"="darkblue","Hepatitis"="darkred"))
    
    #file_name<-paste("box",i,".png",sep="")
    #png(file_name,width=1000,height=650)
    print(g)
    dev.off()
  #dev.off()
}

##### HISTOGRAMS #####
for (i in 3:length(new_data)){
    
    
    dat.don=data.frame(var=new_data[which(new_data$Category=="Blood Donor"),i])
    p1=ggplot(dat.don, aes(x=var))+
      geom_histogram(bins=round(sqrt(540),0),fill="lightblue", #position="identity", lightblue
                     size=0.5,color="black")+
      labs(title="Blood Donor",subtitle="",y="Freguency",x="Value")+
      theme(plot.title=element_text(size=12,face="bold",hjust = 0.5,vjust=-5),legend.position="none")
     
    
    dat.hep=data.frame(var=new_data[which(new_data$Category=="Hepatitis"),i])
    p2=ggplot(dat.hep, aes(x=var))+
      geom_histogram(bins=round(sqrt(75),0),fill="lightcoral", #position="identity", lightcoral
                     size=0.5,color="black")+
      labs(title="Hepatisis",subtitle="",y="Freguency",x="Value")+
      theme(plot.title=element_text(size=12,face="bold",hjust = 0.5,vjust=-5),legend.position="none")
     
    
    #pdf(paste("histogram",i,".pdf",sep=""),width=11,height=6)
    
    figure=ggarrange(p1,p2,ncol=1)
    anot=annotate_figure(figure,fig.lab=paste("",names[i], ""),fig.lab.face = "bold",fig.lab.size = 16)
    print(anot)

    dev.off()
    
    rm(dat.don,dat.hep,p1,p2)
}

##### Descriptive statistics #####
stat=c();
for (i in 3:length(new_data)){
  
  sum=as.matrix(aggregate(new_data[,i],list(new_data$Category),summary))
  sd1=as.matrix(aggregate(new_data[,i],list(new_data$Category),sd))
  var1=as.matrix(aggregate(new_data[,i],list(new_data$Category),var))
  cols=cbind(sum,sd1[,2],var1[,2]); 
  stat=rbind(stat,cols)
  rm(cols)
}
colnames(stat)=c("group","min","q1","median","mean","q3","max","sd","variance")
write.xlsx(stat,"descriptive_stat.xlsx",row.names=TRUE,col.names=TRUE)

##### BOX COX ####

box.cox.alt<-boxcox(new_data$ALT,objective.name="Log-Likelihood")
box.fit.alt<-as.data.frame(t(rbind(box.cox.alt$lambda,box.cox.alt$objective)))
lambda.alt<-box.cox.alt$lambda[max(box.cox.alt$objective)==box.cox.alt$objective]
plot(box.cox.alt,main="Box-Cox Transformation Results",sub="Log-Likelihood vs Lambda")# choose lambda 0, it will be: log(x)
abline(v=0,col="red",lty=3)

#box.cox.ast<-boxcox(new_data$AST,objective.name="Log-Likelihood")
#box.fit.ast<-as.data.frame(t(rbind(box.cox.ast$lambda,box.cox.ast$objective)))
#lambda.ast<-box.cox.ast$lambda[max(box.cox.ast$objective)==box.cox.ast$objective]
#plot(box.cox.ast)# choose lambda -1 will be: x^-1

box.cox.bil<-boxcox(new_data$BIL,objective.name="Log-Likelihood")
box.fit.bil<-as.data.frame(t(rbind(box.cox.bil$lambda,box.cox.bil$objective)))
lambda.bil<-box.cox.bil$lambda[max(box.cox.bil$objective)==box.cox.bil$objective]
plot(box.cox.bil)# choose lambda -0.5 == 1/sqrt(x)
abline(v=-0.5,col="red",lty=3)

box.cox.ggt<-boxcox(new_data$GGT,objective.name="Log-Likelihood")
box.fit.ggt<-as.data.frame(t(rbind(box.cox.ggt$lambda,box.cox.ggt$objective)))
lambda.ggt<-box.cox.ggt$lambda[max(box.cox.ggt$objective)==box.cox.ggt$objective]
plot(box.cox.ggt)# choose lambda -0.5 == 1/sqrt(x)
abline(v=-0.5,col="red",lty=3)

new_data[,"ALT"]=log(new_data[,"ALT"])
#new_data[,"AST"]=(new_data[,"ALT"])^(-1)
new_data[,"BIL"]=1/(sqrt(new_data[,"BIL"]))
new_data[,"GGT"]=1/(sqrt(new_data[,"GGT"]))

##### OUTLIERS DETECTION HUBERTA'S METHOD ####
### MC = medcouple, a robust measure of skewness fo univariate data

MC.b=medcouple(new_data[which(new_data$Category=="Blood Donor"),4:13],do.reflect = NULL) # od 4:13 bo nie usuwamy outlierów w AGE
MC.h=medcouple(new_data[which(new_data$Category=="Hepatitis"),4:13],do.reflect = NULL) # od 4:13 bo nie usuwamy outlierów w AGE

for (i in 4:length(new_data)){
    q1.b=quantile(new_data[new_data$Category=="Blood Donor",i],probs=0.25)
    q3.b=quantile(new_data[new_data$Category=="Blood Donor",i],probs=0.75)
    q1.h=quantile(new_data[new_data$Category=="Hepatitis",i],probs=0.25)
    q3.h=quantile(new_data[new_data$Category=="Hepatitis",i],probs=0.75)
    
    iqr.b = q3.b - q1.b
    iqr.h = q3.h - q1.h
    
    if (MC.b[i-3]>=0){
      lower.b = q1.b - 1.5*exp(-4*MC.b[i-3])*iqr.b
      upper.b = q3.b + 1.5*exp(3*MC.b[i-3])*iqr.b
    }
    else{
      lower.b = q1.b - 1.5*exp(-3*MC.b[i-3])*iqr.b
      upper.b = q3.b + 1.5*exp(4*MC.b[i-3])*iqr.b
    }
    idx.b = which(new_data$Category=="Blood Donor"& new_data[,i]>upper.b| new_data$Category=="Blood Donor"& new_data[,i]<lower.b)

    if (MC.h[i-3]>=0){
      lower.h = q1.h - 1.5*exp(-4*MC.h[i-3])*iqr.h
      upper.h = q3.h + 1.5*exp(3*MC.h[i-3])*iqr.h
    }
    else{
      lower.h = q1.h - 1.5*exp(-3*MC.h[i-3])*iqr.h
      upper.h = q3.h + 1.5*exp(4*MC.h[i-3])*iqr.h
    }
    idx.h = which(new_data$Category=="Hepatitis"& new_data[,i]>upper.h | new_data$Category=="Hepatitis"& new_data[,i]<lower.h)

    new_data[idx.b,i]=NA
    new_data[idx.h,i]=NA
}

rm(q1.b,q3.b,iqr.b,lower.b,upper.b,q1.h,q3.h,iqr.h,lower.h,upper.h)

miss_data<-new_data[rowSums(is.na(new_data))>0,]
missin_values<-sum(is.na(new_data))


donor=new_data[new_data$Category=="Blood Donor",-c(1,2)]
hepa=new_data[new_data$Category=="Hepatitis",-c(1,2)]

donor%>%summarise_all(funs(sum(is.na(.))))
md_pattern(donor,rotate.names = TRUE,color=c("azure3","tomato"),
           print_xaxis=TRUE,
           ylab="Number of observations per pattern",
           xlab="Number of outliers",border="white",
           legend.position = "none")

hepa%>%summarise_all(funs(sum(is.na(.))))
md_pattern(hepa,rotate.names = TRUE,color=c("azure3","tomato"),
           print_xaxis=TRUE,
           ylab="Number of observations per pattern",
           xlab="Number of outliers",border="white",
           legend.position = "none")


# Outlier imputation with knn

imputed_data=round(knnImputation(new_data[,4:13], k=5),2)
new_data2=cbind(new_data[,1:3],imputed_data)

#save(new_data2, file="H2_data.RData")
#load("H2_data.RData")

##### NORMALITY TESTING  ####
###### Q-Q plots #####

for (i in 3:length(new_data2)){
  
  pdf(paste("qqplolt",i,".pdf",sep=""),width=12,height=5.5)
  #file_name<-paste("qqplot_donors",i,".png",sep="")
  #png(file_name,width=1000,height=650)
  par(mfrow=c(1,2))
  print(qqnorm(new_data2[which(new_data2$Category=="Blood Donor"),i],pch=3,lwd=1,col="#191970",frame=FALSE, 
               main=bquote(paste("Blood donors Q-Q Plot for",~italic(.(names[i])), " variable"))),
  qqline(new_data2[which(new_data2$Category=="Blood Donor"),i],col="red",lwd=2))

  print(qqnorm(new_data2[which(new_data2$Category=="Hepatitis"),i],pch=3,col="#191970",frame=FALSE, 
               main=bquote(paste("Hepatitis Q-Q Plot for",~italic(.(names[i])), " variable"))),
        qqline(new_data2[which(new_data2$Category=="Hepatitis"),i],col="red",lwd=2))
  #dev.off()
  
} 

###### Shapiro - Wilk test #####
shap.test=c()
for (i in 3:length(new_data2)){
  donor.shapiro= shapiro.test(new_data2[which(new_data2$Category=="Blood Donor"),i])
  donor.p=donor.shapiro$p.value
  donor.t=donor.shapiro$statistic
  
  hepa.shapiro= shapiro.test(new_data2[which(new_data2$Category=="Hepatitis"),i])
  hepa.p=hepa.shapiro$p.value
  hepa.t=hepa.shapiro$statistic
  
  
  res=rbind(donor.p,donor.t,hepa.p,hepa.t)
  shap.test=cbind(shap.test,res)
  rm(res)
}
colnames(shap.test)=colnames(new_data)[3:13]
rownames(shap.test)=c("Blood Donor.pval","Blood Donor.tstat","Hepatitis.pval","Hepatitis.tstat")
write.xlsx(shap.test,"shapiro_test_new.xlsx",row.names=TRUE,col.names=TRUE)

##### VARIANCE HOMOGENITY TESTING ; levene ####
# null hypotesis = variances are equal

levene.res=c()
for (i in 3:length(new_data2)){
  l.data=melt(new_data2[,c(1,i)])
  levene= leveneTest(value~Category,data=l.data)
  p.levene=levene$`Pr(>F)`
  F.levene=levene$`F value`
  stat=rbind(p.levene[1],F.levene[1])
  levene.res=cbind(levene.res,stat)
  rm(stat)
}
colnames(levene.res)<-names[3:13]
rownames(levene.res)<-c("p-value","test statistic")

write.xlsx(levene.res,"levene_res_new.xlsx",row.names=TRUE,col.names=TRUE)

# other method to perfom non-parametric homogenity test based on ranks is: flinger-killeen test
fligner.res=c()
for (i in 3:length(new_data2)){
  l.data=melt(new_data2[,c(1,i)])
  fligner= fligner.test(value~Category,data=l.data)
  p.fligner=fligner[["p.value"]]
  Chi.fligner=fligner[["statistic"]]
  stat=rbind(p.fligner[1],Chi.fligner[1])
  fligner.res=cbind(fligner.res,stat)
  rm(stat)
}
colnames(fligner.res)<-names[3:13]

##### mean testing and CI ####
###### CI for median (non-normally distributed features) #####
library(DescTools)
a1=MedianCI(new_data2[new_data2$Category=="Blood Donor","Age"],conf.level=0.95)
a2=MedianCI(new_data2[new_data2$Category=="Hepatitis","Age"],conf.level=0.95)

a3=MedianCI(new_data2[new_data2$Category=="Blood Donor","ALB"],conf.level=0.95)
a4=MedianCI(new_data2[new_data2$Category=="Hepatitis","ALB"],conf.level=0.95)

a5=MedianCI(new_data2[new_data2$Category=="Blood Donor","ALP"],conf.level=0.95)
a6=MedianCI(new_data2[new_data2$Category=="Hepatitis","ALP"],conf.level=0.95)

a7=MedianCI(new_data2[new_data2$Category=="Blood Donor","ALT"],conf.level=0.95)
a8=MedianCI(new_data2[new_data2$Category=="Hepatitis","ALT"],conf.level=0.95)

a9=MedianCI(new_data2[new_data2$Category=="Blood Donor","AST"],conf.level=0.95)
a10=MedianCI(new_data2[new_data2$Category=="Hepatitis","AST"],conf.level=0.95)

a11=MedianCI(new_data2[new_data2$Category=="Blood Donor","BIL"],conf.level=0.95)
a12=MedianCI(new_data2[new_data2$Category=="Hepatitis","BIL"],conf.level=0.95)

a13=MedianCI(new_data2[new_data2$Category=="Blood Donor","CHE"],conf.level=0.95)
a14=MedianCI(new_data2[new_data2$Category=="Hepatitis","CHE"],conf.level=0.95)

a15=MedianCI(new_data2[new_data2$Category=="Blood Donor","CHOL"],conf.level=0.95)
a16=MedianCI(new_data2[new_data2$Category=="Hepatitis","CHOL"],conf.level=0.95)

a17=MedianCI(new_data2[new_data2$Category=="Blood Donor","CREA"],conf.level=0.95)
a18=MedianCI(new_data2[new_data2$Category=="Hepatitis","CREA"],conf.level=0.95)

cimed=data.frame(Group=rep(c("Blood Donor","Hepatits"),9),Variable=c("Age","Age","ALB","ALB","ALP","ALP","ALT","ALT","AST","AST","BIL","BIL","CHE","CHE","CHOL","CHOL","CREA","CREA"),
                 y.med=c(a1[[1]],a2[[1]],a3[[1]],a4[[1]],a5[[1]],a6[[1]],a7[[1]],a8[[1]],a9[[1]],
                         a10[[1]],a11[[1]],a12[[1]],a13[[1]],a14[[1]],a15[[1]],a16[[1]],a17[[1]],a18[[1]]),
                 y.min=c(a1[[2]],a2[[2]],a3[[2]],a4[[2]],a5[[2]],a6[[2]],a7[[2]],a8[[2]],a9[[2]],
                         a10[[2]],a11[[2]],a12[[2]],a13[[2]],a14[[2]],a15[[2]],a16[[2]],a17[[2]],a18[[2]]),
                 y.max=c(a1[[3]],a2[[3]],a3[[3]],a4[[3]],a5[[3]],a6[[3]],a7[[3]],a8[[3]],a9[[3]],
                         a10[[3]],a11[[3]],a12[[3]],a13[[3]],a14[[3]],a15[[3]],a16[[3]],a17[[3]],a18[[3]]))
ggplot(cimed,aes(x=Variable,y=y.med,col=Group))+
  geom_point()+
  geom_errorbar(aes(ymin=y.min,ymax=y.max))
  

###### CI for mean (normally distributed features)#####

meanB=mean(new_data2[new_data2$Category=="Blood Donor","PROT"])
meanH=mean(new_data2[new_data2$Category=="Hepatitis","PROT"])
sdB=sd(new_data2[new_data2$Category=="Blood Donor","PROT"])
sdH=sd(new_data2[new_data2$Category=="Hepatitis","PROT"])
nb=length(new_data2[new_data2$Category=="Blood Donor","PROT"])
nh=length(new_data2[new_data2$Category=="Hepatitis","PROT"])

marginB <- qt(0.95,df=nb-1)*sdB/sqrt(nb)
marginH <- qt(0.95,df=nh-1)*sdH/sqrt(nh)
#calculate lower and upper bounds of confidence interval
low <- meanB- marginB
low
high <- meanB + marginB
high

low <- meanH- marginH
low
high <- meanH + marginH
high



library(misty)
p.values=c()

###### U Mann Whitney test fro median #####
for (i in 3:(length(new_data2)-2)){
  mann.w = wilcox.test(new_data2[,i]~Category,data=new_data2,paired=FALSE)
  p.m = mann.w[["p.value"]]
  t.m = mann.w[["statistic"]]
  man = rbind(p.m,t.m)
  p.values=cbind(p.values,man)
}

library(effsize)
library(ltm)
library(effectsize)

###### Welch's test for mean #######
welch = test.welch(GGT~Category,data=new_data2,effsize=TRUE)
p.welch = welch[["result"]][["pval"]][2]
t.welch = welch[["result"]][["t"]][2]
p=rbind(p.welch,t.welch)
p.values=cbind(p.values,p)

welch1 = test.welch(PROT~Category,data=new_data2,effsize=TRUE)
p.welch1 = welch1[["result"]][["pval"]][2]
t.welch1 = welch1[["result"]][["t"]][2]
p1=rbind(p.welch1,t.welch1)
p.values=cbind(p.values,p1)

colnames(p.values)<-names[c(3:13)]
write.xlsx(p.values,"mean_testing.xlsx",row.names=TRUE,col.names=TRUE)

##### correction of p-value #####
p.val=p.values[1,]
p.adjusted = p.adjust(p.val,method="bonferroni")

p.adjusted
write.xlsx(p.adjusted,"p_adjusted.xlsx",row.names=TRUE,col.names=TRUE)

##### effect size #####
###### Hedge's ####
library(effsize)
hedges=c()
for (i in 3:length(new_data2)){
  eff.size = cohen.d(new_data2[,i],new_data2$Category,hedges.correction=TRUE)
  h=eff.size[["estimate"]]
  hedges=cbind(hedges,h)
}
colnames(hedges)<-names[c(3:13)]
write.xlsx(hedges,"size_eff.xlsx",row.names=TRUE,col.names=TRUE)


###### Rank biserial correlation coefficient #######
# rank-biserial r_rb correlation and Cliff's delta effect sizes for non-parametric (rank sum) differences.
rank_biserial(new_data2[new_data2$Category=="Blood Donor","CREA"],
              new_data2[new_data2$Category=="Hepatitis","CREA"])


pairs(new_data2[new_data2$Category=="Hepatitis",3:8],pch=19,col=4,main="Pair plot",
      gap=0, #subplot disatnce,
      labels=colnames(new_data2[,3:8]),
      row1attop=FALSE #diagonal direction,
      )
pairs(new_data2[new_data2$Category=="Blood Donor",3:8],pch=19,col=4,main="Pair plot",
      gap=0, #subplot disatnce,
      labels=colnames(new_data2[,3:8]),
      row1attop=FALSE #diagonal direction,
)

library(GGally)
ggcorr(new_data2[new_data2$Category=="Hepatitis",3:13])



pairs=ggpairs(new_data2,
        columns=3:13,
        aes(color=Category,alpha=0.5),
        upper="blank")
auxplot=ggplot(new_data2,aes(x=Age,y=GGT,color=Category))+
  geom_point()
legend=grab_legend(auxplot)
pairs1=putPlot(pairs,legend,9,-1)
show(pairs1)

library(corrplot)
pearson3 <- cor(new_data2[,3:13],method="pearson")
spearman <- cor(new_data2[,3:13],method="spearman")

ggcorrplot::ggcorrplot(pearson3,lab=TRUE,
                       type="lower",lab_size=4,
                       method="square",colors=c("blue","yellow","red"),
                       title = "Pearson correlation")

ggcorrplot::ggcorrplot(spearman,lab=TRUE,
                       type="lower",lab_size=4,
                       method="square",colors=c("blue","yellow","red"),
                       title = "Spearman correlation")

##### Fisher exact test ####
# examining the dependency between two categorical variables research
hepa=new_data[new_data$Category=="Hepatitis",c(1,2)]
don=new_data[new_data$Category=="Blood Donor",c(1,2)]
temp=rbind(don,hepa)
fisher.test(table(temp))

##### point-biserial corr #####
#examining the dependency/correlation between continous variable and dychotomous variable
biserial.cor(x=new_data2$Age,y=new_data2$Category)

##### Binary classification #####

sample = createDataPartition(new_data2$Category, p = 0.8)
sample = unlist(sample)
train = new_data2[sample, ]
test = new_data2[-sample, ]
train$Category = as.factor(ifelse(train$Category=="Hepatitis",1,0))
test$Category = as.factor(ifelse(test$Category=="Hepatitis",1,0))

# Select features using backward method based on BIC criterion
full.model=glm(Category~.,data=train,family="binomial")
summary(full.model)
BIC(full.model)

backward.model=step(full.model,k=log(nrow(train)))
summary(backward.model)
BIC(backward.model)

# Perform LRT test 
lrtest(backward.model) # H0 - null model is better, Ha - our model is better
lrtest(backward.model,full.model) # H0 - our model is better, Ha - full model is better

# Error bar for backward model selection;
CI=exp(confint(backward.model,level=0.95))
OR=exp(coef(backward.model))
CI=as.data.frame(CI)
CI=cbind(CI,OR)
CI=CI[-1,]
theme_set(theme_minimal())
ggplot(CI,aes(x=rownames(CI),lower=CI[,1],upper=CI[,2]))+
  geom_errorbar(data=CI,mapping=aes(x=rownames(CI),ymin=`2.5 %`,ymax=`97.5 %`),width=0.2,color="blue")+
  geom_point(data=CI,mapping=aes(x=rownames(CI),y=CI[,3]))+
  geom_hline(yintercept=1,linetype="dashed",col="red")+
  labs(title="Error bar for Backward model selection",x="Coefficients",y="Confidence Interval")



pred.train = ifelse(predict(backward.model,train,type="response")>treshold,1,0)
pred.train = as.factor(pred.train)
BACC = bacc(train$Category,pred.train)
conf.mat = table(train$Category, pred.train)
sens = sensitivity(train$Category,pred.train,positive='1')
spec = specificity(train$Category,pred.train,positive='1')

library(pROC)

roc(train$Category, predict(backward.model,train,type="response"),plot=TRUE,legacy.axes=TRUE,col="blue",lwd=3,
    print.auc=TRUE,main="ROC curve with AUC for training set")


# Testing
pred.test = ifelse(predict(backward.model,test,type="response")>treshold,1,0)
pred.test = as.factor(pred.test)
BACC.1 = bacc(test$Category,pred.test)
conf.mat = table(test$Category, pred.test)
sens.1 = sensitivity(test$Category,pred.test,positive='1')
spec.1 = specificity(test$Category,pred.test,positive='1')

roc(test$Category, predict(backward.model,test,type="response"),plot=TRUE,legacy.axes=TRUE,col="blue",lwd=3,
    print.auc=TRUE,main="ROC curve with AUC for testing set")


