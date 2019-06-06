suppressPackageStartupMessages(require(lme4))
suppressPackageStartupMessages(require(xtable))
suppressPackageStartupMessages(require(lattice))
suppressPackageStartupMessages(require(Rmisc))
# could have used t.test(x)$conf.int or bootstrapped
civals <- function(x) {
   lb <- mean(x) - qt(.975,length(x)-1)*sd(x)/sqrt(length(x))
   mx <- mean(x)
   ub <- mean(x) + qt(.975,length(x)-1)*sd(x)/sqrt(length(x))
return(c(lb,mx,ub))}
@

<<makedata,tidy=FALSE,echo=FALSE>>=
makedata <- function(reps=1000,nschool=500,nstud=50000,
  dag=1,ratio=.8,name="sim",writeit=FALSE){
effvals <- matrix(nrow=reps,ncol = 1+2*5)
pvals <- matrix(nrow=reps,ncol = 1+2) 
dirvals <- matrix(nrow=reps,ncol = 1+2) 
#during revision so just for a couple of models
dirvals[,1] <- pvals[,1] <- effvals[,1] <- rep(c(0,5),each=reps/2)

for (x in 1:reps){
 StudGr <- rbinom(nstud,1,.5)
 Achieve <- rnorm(nstud,-2+StudGr*2,sd=20)
 pre <- Achieve + rnorm(nstud,sd=20)
 PreGr <- as.numeric(pre > median(pre))
 school <- vector(length=nstud)
 SchGr <- rep(0:1,each=nschool/2)
 ifelse(dag==1,
    StSchGr <- rbinom(nstud,1,ratio*StudGr+(1-ratio)*(1-StudGr)),
    StSchGr <- rbinom(nstud,1,ratio*PreGr+(1-ratio)*(1-PreGr)))
 school[StSchGr == 0] <- 
    sample(1:(nschool/2),sum(StSchGr==0),replace=TRUE)
 school[StSchGr == 1] <- 
    sample((nschool/2 + 1):nschool,sum(StSchGr==1),replace=TRUE)
 seff <- rnorm(nschool,sd=5) + effvals[x,1] * SchGr*1
 SchData <- cbind(1:nschool,seff,SchGr)
 colnames(SchData) <- c("school","TrueVA","SchGr")
 StudData <- cbind(StudGr,Achieve,pre,school,StSchGr)
 Sim1Data <- merge(StudData,SchData,by="school",all.x=TRUE)
 Sim1Data$post <- Sim1Data$Achieve + rnorm(nstud,sd=20) + 
    Sim1Data$TrueVA 
 AvePrex <- aggregate(Sim1Data$pre,by=list(Sim1Data$school),mean)
 colnames(AvePrex) <- c("school","AvePre")
 Sim1 <- merge(Sim1Data,AvePrex,by="school",all.x=TRUE)
 m1 <- unlist(ranef(lmer(post~1+(1|school),data=Sim1)))
 effvals[x,2:3] <- tapply(m1,SchGr,mean)
 m2 <- unlist(ranef(lmer(post~pre+(1|school),data=Sim1)))
 pvals[x,2] <- t.test(m2~SchGr)$p.value
 dirvals[x,2] <- t.test(m2~SchGr)$statistic > 0
 effvals[x,4:5] <- tapply(m2,SchGr,mean)
 m3 <- unlist(ranef(lmer(post~pre+AvePre + 
     (1|school),data=Sim1)))
 effvals[x,6:7] <- tapply(m3,SchGr,mean)
 m4 <- unlist(ranef(lmer(post~pre+StSchGr + 
     (1|school),data=Sim1)))
 effvals[x,8:9] <- tapply(m4,SchGr,mean)
 m5 <- unlist(ranef(lmer(post - pre ~ 0 + 
     (1|school),data=Sim1)))
 pvals[x,3] <- t.test(m5~SchGr)$p.value
 dirvals[x,3] <- t.test(m5~SchGr)$statistic > 0
 effvals[x,10:11] <- tapply(m5,SchGr,mean)
}
if(writeit) write.csv(effvals,paste0("ef",name,".csv"),row.names=FALSE)
if(writeit) write.csv(pvals,paste0("pv",name,".csv"),row.names=FALSE)
if(writeit) write.csv(dirvals,paste0("dir",name,".csv"),row.names=FALSE)
return(list(effvals=effvals,pvals=pvals,dirvals=dirvals))}
