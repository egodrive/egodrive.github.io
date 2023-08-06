

##### P(A and B mate) ####  
#####____ Random mating simulation #####
library(data.table)
# How does N and SR affect probability of mating under random mating?
test2 = (
  lapply(10:50, function(N){
    lapply(seq(.1,.9,length.out = 20), function(x){
      CJ(Female = paste("F",c(1:(round(N*x))),sep = "_"), 
         Males = paste("M",c(1:(round(N*(1-x)))),sep = "_"))[,SR:=x][,N:=N]
    })  
  })
)
test2 = rbindlist(lapply(test2, rbindlist))
test2[,SR:=uniqueN(Males)/uniqueN(Female), c("SR", "N")]
test2[,nDyads:=.N, c("SR","N")]
test2 = test2[SR<1.5]
with(test2, plot(nDyads~SR))
test2[,Realised:=c(1, rep(0,(.N-1))), c("Female","N", "SR")]
glm(Realised~SR+N, data = test2, family = "binomial")
library(ggplot2);library(cowplot)
qplot(y = nDyads, x = SR, group = N, col = (N), data = test2) + stat_smooth(method = "lm", formula = y ~ x + I(x^2))

summary(glm(Realised~SR+N, 
            data = test2, family = "binomial"))
summary(glm(Realised~SR*N, 
            data = test2, family = "binomial"))

### Female preferences
Pref = function(trait, Old=FALSE, AssignMating = T, 
                Rate = c(0,6), MidPoint = c(0,1), linear = F){
  # AssignMating: assign mating to the one with highest prob, i.e. all females must mate
  # linear: usually assume a sigmod preference curve, but if 'linear = T', the rate parameter is the slope of a linear preference curve going through 
  # origo
  # trait = rnorm(n = 10)
  
  K = 1
  if(linear){
    prefs = sapply(1:length(Rate), function(x){
      prefs =MidPoint[x] + Rate[x]*trait
      prefs[prefs<0]<-0
      prefs[prefs>1]<-1
      prefs
      # if(uniqueN(prefs)==1){
      #   prefs <- rep(.5, length(prefs))
      # }else{
      #   prefs <- (prefs-min(prefs))/(max(prefs) - min(prefs))} # scale to be between 0->1
      # 
    })
  }else{
    prefs = sapply(1:length(Rate), function(x){
      prefs = K/(1 + exp(-(Rate[x])*(trait-(MidPoint[x]))))
    })
  }
  
  
  if(AssignMating){
    idx = apply(prefs, 2, which.max)
    for(x in 1:ncol(prefs)){
      prefs[idx[x],x]<-1
    }}
  attr(prefs, "Rate")<-Rate
  attr(prefs, "MidPoint")<-MidPoint
  prefs
}

library(vegan)
library(data.table)

# can try to use the morisita index also, which one can compare with a test statistics to see if it deviates from random mating at 1
# Tsuji & Tsuji 1998. 
Morisita = function(matings,...){
  n = length(matings)
  n*((sum(matings^2, na.rm = T) - sum(matings, na.rm = T))/((sum(matings, na.rm = T)^2) - sum(matings, na.rm = T)))
}

MatingDistribution = function(runs = 25, 
                              N = seq(20,50,length.out = 7), 
                              SR = seq(.3,1.5,length.out = 7),
                              MaleTraitsDistri = "uniform",
                              #FemaleAge = c(0,1),
                              Mids = c(0,0),
                              Rates = c(0,6),
                              IncludePopStructure=FALSE){
  rbindlist(lapply(1:runs, function(run){
    DT = rbindlist(lapply(SR, function(SR){
      rbindlist(lapply(N, function(N){
        # SR = .5;N = 40
        Nfemales = round(N/(SR + 1))
        Nmales = N - Nfemales
        if(MaleTraitsDistri=="uniform"){
          MaleTraits = sample(c(-10:10)/10, 
                              Nmales, replace = T)
        }else{
          MaleTraits = rnorm(n = Nmales)
        }
        
        
        MatingMatrix = Pref(MaleTraits, AssignMating = F,  
                            MidPoint = Mids, Rate = Rates)
        
        # Female encounter males in a random sequence
        MS = t(sapply(1:Nfemales, function(female){
          # set.seed(round(female*N*SR*run)) # so different combination have different random seeds
          #      print(female)
          FemaleAge = rbinom(n = 1, size = c(1,2), prob = c(.5,.5))+1 # The female is randomly assigned age
          Mated=0
          nMalesSampled = 0
          while(Mated==0){
            # She samples males before deciding
            nMalesSampled = nMalesSampled + 1
            Male = sample(1:Nmales, 1)
            #print(Male)
            Mated = sample(c(0,1), 1, replace = T,
                           prob = c(1-MatingMatrix[Male, FemaleAge],
                                    MatingMatrix[Male, FemaleAge]))
            
            
          }
          c(female, FemaleAge, Male, nMalesSampled)
        }))
        MS = data.table(MS)
        names(MS) = c("FemaleID", "FemaleAge", "MaleID", "nMalesSampled")
        MorI = dispindmorisita(sapply(1:Nmales, function(x){
          length(which(MS$MaleID %in% x))
        }))[[1]]
        
        
        
        if(IncludePopStructure){
          DT = data.table(MaleID = 1:Nmales, MaleTraits
          )[CJ(FemaleID=1:Nfemales, MaleID=1:Nmales
          ), on = "MaleID"]
          MS = MS[DT, on = c("FemaleID", "MaleID")
          ][,Realised:=ifelse(is.na(FemaleAge),0,1)
          ]
          
          MS[,FemaleAge:=imputeTS::na_locf(FemaleAge), "FemaleID"
          ][,nMalesSampled:=imputeTS::na_locf(nMalesSampled), "FemaleID"
          ]#[,MorI:=na.omit(MorI)]
          
        }
        MS$MorI = MorI
        MS$N = N
        MS$SR = SR
        MS$Nfemales = Nfemales
        MS$Nmales = Nmales
        MS$meanMaleTrait = mean(MaleTraits)
        MS$sdMaleTrait = sd(MaleTraits)
        MS = data.table(MS)
        MS
      }))
    }))
    DT[,Run := run]
    DT[,Mid:=ifelse(FemaleAge==2 & length(Mids)==2, Mids[2], Mids[1])]
    DT[,Rate:=ifelse(FemaleAge==2 & length(Rates)==2, Rates[2], Rates[1])]
    DT}))
}
Checks = MatingDistribution(runs = 10, 
                            IncludePopStructure = T)

Checks[,meanFemaleAge:=mean(FemaleAge), c("Run", "N", "SR")]
p = .2
Checks[,nOffspring:=ifelse(FemaleAge==2, 
                           sample(c(1,2),uniqueN(FemaleAge==2), prob = c(p,1-p)),
                           sample(c(1,2),uniqueN(FemaleAge==1), prob = c(1-p,p))
)]
with(Checks, boxplot(nOffspring~FemaleAge))

Checks[,nUnmated:=Nmales-uniqueN(.SD[Realised==1]$MaleID), c("Run", "N", "SR")]
Checks[,Run := factor(Run)]

gMod = glm(Realised~MaleTraits*SR*N + meanMaleTrait + (Run), data = Checks, family = "binomial", na.action = na.fail)
library(MuMIn)
dMods = dredge(gMod, fixed = "Run")
head(dMods)


with(Checks[,.SD[1], c("Run", "N", "SR")], cor.test(meanFemaleAge, meanMaleTrait))
library(MuMIn)
df = Checks[,.SD[1], c("Run", "N", "SR")]
df[,LogMorI:=log(MorI)][,CVAgeMale:=(sdMaleTrait^2)/(meanMaleTrait^2)]
pairs(df[MorI >0 & CVAgeMale <1000][,c("LogMorI", "meanFemaleAge", "meanMaleTrait", "CVAgeMale")])
lm(log(MorI)~meanFemaleAge+meanMaleTrait + sdMaleTrait, data = df[MorI>0])


set.seed(1)
windows()
par(mfrow = c(2,3), bty = "L")
dd = Checks[,.SD[1], c("Mid", "Rate")]
curve(Pref(trait = x, Rate = 0, Mid = 0), from =-2, to =2, col = "green",xlim = c(-1.75,1.75), ylim = c(0,1),
      xlab = "Male trait", ylab = "P(mating)", main = "Female preference curve")
curve(Pref(trait = x, Rate = 6, Mid = 0), from =-2, to =2, col = "red", add = T)
legend("topleft", bty = "n",
       lwd = 1,
       legend = c("Young females", "Old females"),
       col = c("green", "red")
)



with(Checks[,.SD[1], c("Run", "N", "SR")], boxplot((1-(nUnmated/Nmales))~cut(meanMaleTrait, seq(-.5,.5,.1)), xlab = "mean maletrait", ylab = "Proportion mated males"))
with(Checks[,.SD[1], c("Run", "N", "SR")], boxplot((1-(nUnmated/Nmales))~cut(meanFemaleAge, seq(1,2,.1)), xlab = "meanFemaleAge", ylab = "Proportion mated males"))

with(Checks, plot(log(nMalesSampled)~meanMaleTrait, 
                  col = ifelse(Rate==6, "red", "green")))

with(Checks, plot(log(nMalesSampled)~Nmales, 
                  col = ifelse(Rate==6, "red", "green")))

summary(glm(nMalesSampled~meanMaleTrait + Nmales,
            data = Checks, family = "poisson"))


with(Checks[,.SD[1], c("Run", "N", "SR")],
     plot(MorI~meanMaleTrait, main = "Random male encounter \n (N: 20-50, SR: .3-1)",
          xlab = "Mean male trait value", ylab = expression(I[delta])))
with(Checks[,.SD[1], c("Run", "N", "SR")],
     plot(MorI~meanFemaleAge, main = "Random male encounter \n (N: 20-50, SR: .3-1)",
          xlab = "Mean female age", ylab = expression(I[delta])))
# with(Checks[,.SD[1], c("Run", "N", "SR")],
#      plot(MorI~Nmales, main = "Random male encounter \n (N: 20-50, SR: .3-1)",
#           xlab = "SR", ylab = expression(I[delta])))


with(Checks, plot(MaleTraits~meanMaleTrait, 
                  col = ifelse(Rate==6, "red", "green")))
ggplot(data = Checks, 
       aes(y = Realised, x = MaleTraits, col = meanMaleTrait))+
  geom_point()+
  geom_smooth(method = glm, 
              method.args= list(family="binomial"),
              aes(fill=meanMaleTrait))
glm(Realised~MaleTraits,
    data = Checks, family =  binomial)
# 
# 
# plot(sapply(1:100, function(nMales) {mean(Pref(trait = rnorm(mean = 0, sd = 1, n = nMales), AssignMating = T,
#                                           Old = F, MidYoung = 0, RateYoung = 0))})~c(1:100),
#      add = T, col = "green", xlab = "Number of males", ylab = "Probability of mating with a given male", ylim = c(0,1),
#      main = "Mean male trait = 0")
# 
# points(sapply(1:100, function(x) {
#   mean(Pref(trait = rnorm(mean = 0, sd = 1, n = x), RateOld = 6, MidOld = 1, AssignMating = T, Old = T))
#   })~c(1:100), add = T, col = "red"
#        , xlab = "Number of males", ylab = "Probability of mating", main = "Preference curve")
# # points(sapply(1:100, function(x) {mean(Pref(trait = rnorm(mean = 0, sd = 1, n = x), 
# #                                             AssignMating = T, Old = F))})~c(1:100), add = T, col = "green")
legend("topright", bty = "n",
       pch = 1, cex = 1,
       legend = c("Young females", "Old females"),
       col =  c( "green","red")
)
# plot(sapply(1:100, function(x) {mean(Pref(trait = rnorm(mean = 1, sd = 1, n = x), AssignMating = T,
#                                           Old = F, MidYoung = 0, RateYoung = 0))})~c(1:100),
#      add = T, col = "green", xlab = "Number of males", ylab = "Probability of mating with a given male", ylim = c(0,1),
#      main = "Mean male trait = 1")
# 
# points(sapply(1:100, function(x) {
#   mean(Pref(trait = rnorm(mean = 1, sd = 1, n = x), RateOld = 6, MidOld = 1, AssignMating = T, Old = T))
# })~c(1:100), add = T, col = "red"
# , xlab = "Number of males", ylab = "Probability of mating", main = "Preference curve")
# # points(sapply(1:100, function(x) {mean(Pref(trait = rnorm(mean = 0, sd = 1, n = x), 
# #                                             AssignMating = T, Old = F))})~c(1:100), add = T, col = "green")
legend("topright", bty = "n",
       pch = 1, cex = 1,
       legend = c("Young females", "Old females"),
       col =  c( "green","red")
)
