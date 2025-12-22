library(rstan)
library(parallel)

data<-read.table("data/data.txt",header=F)

colnames(data)<-c("subj","expt","item","cond","pos","word","region","rt")

#RC1: Distance manipulaton between relative pronoun and the relative clause verb in Subject RCs and Object RCs. Replication of the Levy Russian RC paper.

#No. of conditions: 4
#Condition a: Subject RC, Long (Canonical word-order)
#Condition b: Subject RC, Short (Non-Canonical word-order)
#Condition c: Object RC, Long (Canonical word-order)
#Condition d: Object RC, Short (Non-Canonical word-order)

#Total items: 24

#Critical region: RCVerb

RC1<-subset(data,expt=="RC1")

## nested contrasts:
RC1$dist.SR  <- ifelse(RC1$cond=="a",1,
                       ifelse(RC1$cond=="b",-1,0))
RC1$dist.OR  <- ifelse(RC1$cond=="c",1,
                       ifelse(RC1$cond=="d",-1,0))
RC1$RCType <- ifelse(RC1$cond%in%c("a","b"),-1,1)

## anova contrasts:
RC1$dist  <- ifelse(RC1$cond%in%c("a","c"),1,-1)
## RCType will be one contrast used here.
RC1$int <- ifelse(RC1$cond%in%c("a","d"),-1,1)

RC1<-subset(RC1,region!="1" & region!="0" & region!="-")
library(gdata)
RC1$region<-drop.levels(RC1$region)
## we will use log rt:
RC1$lrt<- log(RC1$rt)

RC1<-subset(RC1,region=="RCVerb")

RC1$region<-drop.levels(RC1$region)

write.table(RC1,file="expt1critdata.txt")

library(lme4)
m1<-lmer(lrt~dist+RCType+int+(1+dist+RCType+int|subj)+(1|item),RC1)
summary(m1)

## Stan/JAGS ready data:
dat <- list(mu_prior=c(0,0,0,0),
             subj=sort(as.integer(factor(RC1$subj))),
             lrt = RC1$lrt,
             distance = RC1$dist,
             rctype = RC1$RCType,
             interaction = RC1$int,
             N = nrow(RC1),
             I = length(unique(RC1$subj))
            )  

## to be used later: crossed subject and item:
dat2 <- list(mu_prior=c(0,0,0,0),
            subj=sort(as.integer(factor(RC1$subj))),
            item=sort(as.integer(factor(RC1$item))),
            lrt = RC1$lrt,
            distance = RC1$dist,
            rctype = RC1$RCType,
            interaction = RC1$int,
            N = nrow(RC1),
            I = length(unique(RC1$subj)),
            K = length(unique(RC1$item))
)  


#################################################
#e1.sm <- stan_model("expt1subj.stan", model_name = "e1subj")
#sflist <- mclapply(1:4, mc.cores = detectCores(),
#                   function(i) sampling(e1.sm, data = dat,
OA#                                        chains = 1, chain_id = i, seed = 12345))

#e1.sf <- sflist2stanfit(sflist)

#print(e1.sf,digits=4)
#                mean    sd      
#beta[1]         6.5546  0.0046
#beta[2]        -0.0289  0.0002
#beta[3]         0.0345  0.0002
#beta[4]         0.0004  0.0001

#pnorm(0,-0.0289,0.0002,lower.tail=F)
#pnorm(0,0.0345,0.0002,lower.tail=T)

#res<-as.matrix(e1.sf)

#head(res)

## model recovers all effects:
## effect of dist:
#hist(res[,2])
## effect of rc type:
#hist(res[,3])
# no interaction:
#hist(res[,4])

#plot(e1.sf)

## expt 4:
completiondata<-read.table("sentence_completion_stats.txt",header=T)
##probabilities:
probs<-c(0.5,0.3,
         1,0.2,
         1,0.3,
         0.9,0,
         0.9,0,
         1,0.4,
         0.8,0.3,
         0.2,0,
         1,0.7,
         1,0,
         1,0,
         1,0.1,
         1,0,
         0.5,0,
         1,0,
         1,0.5)

probabilities<-data.frame(item=rep(1:16,each=4),cond=rep(letters[1:4],16),probs=rep(probs,each=2))

meanprobs<-with(probabilities,tapply(probs,cond,mean))

#Experiment 4 (Complex predicate)
##
CP1<-subset(data,expt=="CP1")

## Data preparation for plotting:

## crit region: CPLightVerb/MainVerb:

a.regions<-c("PreNounPred","CPNounPred","PreLV","CPLightVerb", "PreCoordMV","CoordMainVerb")
b.regions<-c("PreNounPred","CPNounPred","PreLV","CPLightVerb", "PreCoordMV","CoordMainVerb")
c.regions<-c("PreMVObj","MVObj","PreMV","MainVerb", "PreCoordMV","CoordMainVerb")
d.regions<-c("PreMVObj","MVObj","PreMV","MainVerb", "PreCoordMV","CoordMainVerb")

regions<-c(a.regions,b.regions,
           c.regions,d.regions)
region.id<-rep(1:6,4)
cond.id<-rep(letters[1:4],each=6)

region.df<-data.frame(cond=factor(cond.id),
                      region.id=region.id,
                      region=factor(regions))

CP1<-subset(CP1,region!="1" & region!="0" & region!="-")

CP1$region<-drop.levels(CP1$region)
CP1$cond<-drop.levels(CP1$cond)

#update labels in the long conditions (PreLV1, PreMV1)
CP1$region[CP1$region == 'PreMV1'] <- 'PreMV'
CP1$region[CP1$region == 'PreLV1'] <- 'PreLV'

#Labeling error in the data; item==15, cond==a,b. 'PreMV' should be 'PreLV'
CP1$region[CP1$region == 'PreMV' & CP1$item == 15 & CP1$cond == 'a'] <- 'PreLV'
CP1$region[CP1$region == 'PreMV' & CP1$item == 15 & CP1$cond == 'b'] <- 'PreLV'

library(plyr)

#merge multiple "PreNounPred", "PreMVObj", "PreLV", "PreMV", "PreCoordMV" regions
CP1.uniq.reg<-ddply(CP1, .(subj, expt, item, cond, region), summarize, rt = sum(rt))

CP1.uniq.reg.merged<-merge(CP1.uniq.reg,region.df, by.x=c("cond","region"))

CP1.uniq.reg.merged$cond<-drop.levels(CP1.uniq.reg.merged$cond)


#head(CP1.uniq.reg.merged)

#with(CP1.uniq.reg.merged,tapply(rt,IND=list(cond,region.id),mean))

#head(CP1.uniq.reg.merged)

#unique(CP1.uniq.reg.merged$region)
#CP1: Expectation build-up for the verb due to nominal host of Complex Predicate. Effect of suprisal/locality on Normal vs Complex predicate.

#No. of conditions: 4
#Condition a: Expectation, Long
#Condition b: Expectation, Short
#Condition c: No Expectation, Long
#Condition d: No Expectation, Short

#Total items: 16

CP1<-CP1.uniq.reg.merged

#Critical region: CPLightVerb
## these are all region.id==4:
CP1.ab.crit<-subset(CP1,region=="CPLightVerb")
CP1.cd.crit<-subset(CP1,region=="MainVerb")
CP1.crit<-rbind(CP1.ab.crit,CP1.cd.crit)

CP1.crit$region2<-factor("crit")

CP1.postcrit<-subset(CP1,region.id==5)
CP1.postcrit$region2<-factor("postcrit")

CP1.postcrit2<-subset(CP1,region.id==6)
CP1.postcrit2$region2<-factor("postcrit2")

CP1.critdata<-rbind(CP1.crit,CP1.postcrit,CP1.postcrit2)

#completiondata<-read.table("sentence_completion_stats.txt",header=T)

#dim(completiondata)
## print out critical words:
#completiondata[31:32,c(1:2,68:69)]
#completiondata[31:32,]

## odds: (abandoned)
#c(0.5/0.5,0.3/0.7,
#  Inf,0.2/0.8,
#  Inf,0.3,
#  0.9/0.1,0,
#  0.9/0.1,0,
#  Inf,0.4/0.6,
#  0.8/0.2,0.3/0.7,
#  0.2/0.8,0,
#  Inf,0.7/0.3,
#  Inf,0,
#  Inf,0)


temp<-merge(CP1.critdata,probabilities,by.x=c("item","cond"),by.y=c("item","cond"))

CP1.critdata<-temp

## nested contrasts:
CP1.critdata$dist.exp<-ifelse(CP1.critdata$cond=="a",-1,ifelse(CP1.critdata$cond=="b",1,0))
CP1.critdata$dist.noexp<-ifelse(CP1.critdata$cond=="c",-1,ifelse(CP1.critdata$cond=="d",1,0))
CP1.critdata$exp<-ifelse(CP1.critdata$cond%in%c("a","b"),1,-1)

## anova contrasts:
CP1.critdata$dist<-ifelse(CP1.critdata$cond%in%c("a","c"),-1,1)
## exp is above:
CP1.critdata$int<-ifelse(CP1.critdata$cond%in%c("a","d"),-1,1)

CP1.critdata$lrt<-log(CP1.critdata$rt)

e2critdata<-subset(CP1.critdata,region2=="crit")

e2critdata$region2<-drop.levels(e2critdata$region2)

summary(e2critdata)

write.table(e2critdata,file="expt2critdata.txt")

#Remove items 1 and 8 because of the sentence completion study:

#CP1.critdata.reduced<-subset(CP1.critdata,item!=1)
#CP1.critdata.reduced<-subset(CP1.critdata.reduced,item!=8)

#sort(unique(CP1.critdata.reduced$item))

#CP1.critdata<-CP1.critdata.reduced
#head(CP1.critdata)

pred<-ifelse(CP1.critdata$probs<median(CP1.critdata$probs),"low","high")

CP1.critdata$pred<-factor(pred)

means.pred<-with(subset(CP1.critdata,region2=="crit"),tapply(rt,IND=list(factor(dist),pred),mean))

CP1.crit<-subset(CP1.critdata,region2=="crit")

dat2e2 <- list(mu_prior=c(0,0,0,0),
             subj=sort(as.integer(factor(CP1.crit$subj))),
             item=sort(as.integer(factor(CP1.crit$item))),
             lrt = log(CP1.crit$rt),
             distance = CP1.crit$dist,
             expectation = CP1.crit$exp,
             interaction = CP1.crit$int,
             N = nrow(CP1.crit),
             I = length(unique(CP1.crit$subj)),
             K = length(unique(CP1.crit$item))
)  

sink("e2results.txt")

set_cppo('fast')

e2.sm <- stan_model("expt2subjitem.stan", model_name = "e2subjitem")
sflist <- mclapply(1:4, mc.cores = detectCores(),
                   function(i) sampling(e2.sm, data = dat2e2,
                                        chains = 1, chain_id = i, seed = 12345))

e2.sf <- sflist2stanfit(sflist)

print(e2.sf,digits=4)
sink()