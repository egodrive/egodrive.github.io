library(data.table)
nRuns = 25
tmax = 500

BetaRuns = lapply(c(.3,.75,.9), function(d){
  print(paste("Beta: ", d))
  nRuns = 25
  Runs = rbindlist(lapply(1:nRuns, function(run){
    print(run)
    set.seed(run)
    
    # Deviation from optimal
    #deviation =  sample(seq(-1,1,.5),1)
    deviation =  seq(-2,2,1)[ceiling(run/5)]
    
    set.seed(2)
    #Run time and population
    tmax = 1000
    nInd = 1
    
    ## Make grid
    maxD = 2500
    densP = 0.6 # tetthet av mat
    nP = floor(maxD*(densP)) # antall matlokasjoner
    Y = X = seq(0,maxD, by = 1)
    
    # Box size
    sizeP = 20
    Radii = rep(sizeP, nP)
    
    Grid = TRUE
    if(Grid==TRUE){
      uns = seq(0,maxD, Radii[1])+Radii[1]/2
      SampleGrid = expand.grid(uns,uns)
      plotsC = SampleGrid[sample(1:nrow(SampleGrid), nP),]
      #plotsC = plotsC+Radii[1]/2 # For å bruke midten av cellen
      #plot(1~1, type = "n", xlim = c(-20,maxD+20), ylim = c(-20, maxD+20))
      
      plots = lapply(1:nrow(plotsC), function(x){
        bbox = rbind(x = c(plotsC[x,1]-Radii[1]/2, plotsC[x,1]+Radii[1]/2),
                     y = c(plotsC[x,2]-Radii[1]/2, plotsC[x,2]+Radii[1]/2))
        #rect(plots[[x]][1,1], plots[[x]][2,1], plots[[x]][1,2], plots[[x]][2,2])
        
        bbox
      })
    }else{
      
      # Create initial kernels
      cP = 1
      plotsC = matrix(round(runif(1*2, min = 0, max = maxD)), ncol = 2)
      
      plots = list()
      plots[[cP]]<-rbind(x = c(plotsC[cP,1]-Radii[cP], plotsC[cP,1]+Radii[cP]),
                         y = c(plotsC[cP,2]-Radii[cP], plotsC[cP,2]+Radii[cP]))
      
      while(cP<nP){
        #print(cP)
        
        # Generer en ny en
        t1 = matrix(round(runif(1*2, min = 0, max = maxD)), ncol = 2)
        # Er den nye for nærme noen eksisterende?
        a = as.matrix(dist(rbind(t1,plotsC)))[-1,1]
        # Hvis nei, bygg boksen
        bnd = sqrt(Radii[1]^2 + Radii[1]^2)*2
        if(!any(a<bnd)){
          #print(min(a)<7.07)
          #Sys.sleep(.2)
          t2 = rbind(x = c(t1[1,1]-Radii[1], t1[1,1]+Radii[1]),
                     y = c(t1[1,2]-Radii[1], t1[1,2]+Radii[1]))
          
          cP <- cP + 1
          plots[[cP]]<-t2
          plotsC <-rbind(plotsC,t1)
        }
      }
    }
    
    plot(1~1, type = "n", xlim = c(-20,maxD+20), ylim = c(-20, maxD+20))
    
    lapply(seq_along(plots), function(x){
      rect(plots[[x]][1,1], plots[[x]][2,1], plots[[x]][1,2], plots[[x]][2,2])
    })
    
    ####
    # Exponent of the gain curve
    Betas = d#.75 #runif(nInd, min = 0.1, max = .90)
    
    Stay = rep(0, nInd)
    Arrival = rep(0, nInd)
    Residence = rep(0, nInd)
    
    Gains = rep(0, nInd)
    Gained = Gains
    
    # generate random start locations for the individuals
    set.seed(run)
    pos = t(sapply(1:nInd, function(x) round(runif(2, min = maxD*0.45, max = maxD*.55))))
    
    
    #pos = matrix(c(30,40), ncol = 2)
    DF = array(0,dim = c(tmax, 8, nInd),
               dimnames = list(c(1:tmax),
                               c("State", "Box", "ResidenceTime", "Traveltime","Gain","PosX", "PosY", "Dist"),
                               c(1:nInd)))
    DF[1,"Traveltime",]<-0# round(runif(nInd, min = 4, max = 10))
    class(DF[,"Traveltime",])<-"numeric"
    
    # Merk den gjennomsnittlige steg lengden er kortere enn lengden/bredden av en foraging box
    Moves = function(Mode = "Levy", alpha = 2, n = 1){
      t(sapply(1:n, function(x){if(Mode == "Levy"){
       dist = 0
         while(dist==0 | dist>500){# For å forhindre kjempelange steg
           #dist = runif(1)^(-1/alpha) 
           dist =   runif(1)^(1/(1 - alpha)) # Alternativ https://rdrr.io/cran/adehabitatLT/src/R/simm.levy.r
          }
        
         angle = runif(1)*2*pi
        
        dx = dist * sin(angle)
        dy = dist * cos(angle)
      }else{
        dx = rnorm(n = 1, mean = 0, sd = 5)
        dy = rnorm(n = 1, mean = 0, sd = 5)
      }
        c(dx,dy)}))
    }
    
    t = 1
    while(t < tmax){
      #print(t)
      #Sys.sleep(.1)
      inBox = apply(pos, 1, function(x){
        sapply(plots, function(P){ 
          x[1] %between% P[1,] & # x-coordinate
            x[2] %between% P[2,]# y-coordinate
        }) 
      }, simplify = T)
      
      # Evaluate whether individuals are in boxes or not
      for(cols in 1:ncol(inBox)){
        DF[t,c("PosX", "PosY"),cols]<-pos[cols,]
        #print("Forage evaluations")
        
        # in a box, and should stay there longer
        #if(any(inBox[,cols])){print("In box")}
        if(any(inBox[,cols]) & Arrival[cols]==0){
          #print("Box arrival")
          # First time landed in plot
          Arrival[cols] <- t
          # gets no food on arrival
          
          # How long should the individual remain in the box
          RT = round((as.numeric(DF[t,"Traveltime",cols])/(1 + (1/Betas[cols]))) +
                       deviation)
          DF[t,"ResidenceTime",cols] <- RT
          Stay[cols] <- t+as.numeric(DF[t,"ResidenceTime",cols])
          
          DF[t,"State",cols]<-"Arrival"
          # Which box is individual in?
          DF[t,"Box",cols]<-min(which(inBox[,cols]))
        }
        
        if(any(inBox[,cols]) & Arrival[cols]<t){ 
          # Is in plot from previous round, and we have passed the arrival time
          #print("Box stayed")
          
          # Which box is individual in?
          DF[t,"Box",cols]<-min(which(inBox[,cols]))
          
          # How long have the individual been in the plot
          ResidenceTime = as.numeric(DF[(t-1),"ResidenceTime",cols])
          # How much to gain?
          Gains[cols] <- Gains[cols]  + ((ResidenceTime^Betas[cols])-((ResidenceTime-1)^Betas[cols]))
          #print(Gains[cols])
          DF[t,"ResidenceTime",cols]<-as.numeric(DF[(t-1),"ResidenceTime",cols])-1
          DF[t,"State",cols]<-"Forage"
          
          # Reset travel time since you have stayed in a box
          DF[t,"Traveltime",cols]<-0
        }
        
        # Move individuals
        for(cols in 1:ncol(inBox)){
          #print(paste("Write gains: ",Gains[cols]))
          DF[t,"Gain", cols]<-(Gains[cols]) 
          #print(DF[t,"Gain", cols])
          #print("Gains should have been written")
          # plot positions
          points(pos[cols,2]~pos[cols,1],
                 cex = 2, 
                 pch = ifelse(DF[(t-1),"State",cols]=="Movement",15,17),
                 
                 col = ifelse(any(inBox[,cols]), "green", "red"))
          
          # Only move those who are outside of box or spent too much time in box
          if(!any(inBox[,cols]) | DF[t,"ResidenceTime",cols]==0){
            DF[t,"State",cols]<-"Movement"
            # Reset arrival as we have started to move again
            Arrival[cols] <- 0
            
            #print("Movement")
            
            # If in a box, move it to the edge to ensure movement away from box
            if(any(inBox[,cols])){
              
              tmp = plots[[min(which(inBox[,cols]))]]
              pos[cols,] <-round(c(runif(1,tmp[1,1],tmp[1,2]),
                                   runif(1,tmp[2,1],tmp[2,2])))
              
              # Include a while-loop, to run as long as the individual remains inside. this way to ensure movement away from the patch
              while(any(sapply(plots, function(P){ 
                pos[cols,][1] %between% P[1,] & # x-coordinate
                  pos[cols,][2] %between% P[2,]# y-coordinate
              }))){
                dist = Moves()
                pos[cols,] = round(pos[cols,] + dist)
              }
            }else{
              dist = Moves()
              pos[cols,] = round(pos[cols,] + dist)
            }
            DF[t,"Dist",cols]<-sqrt(sum(dist^2))
            
            # Increase travel time for those who are not in boxes
            DF[t+1,"Traveltime",cols]<-as.numeric(DF[t,"Traveltime",cols])+1
          }
        }
        # If moving outside the grid, relocate to new random place
        for(cols in 1:nInd){
          if(any(pos[cols,]>maxD) | any(pos[cols,]<0)){
            #print("Outside")
            pos[cols,]<-runif(2, min = maxD*0, max = maxD*1)
          }
        }
        #print(DF[t,"Gain", ])
        #Gained = rbind(Gained, Gains)
        t = t + 1
      }
    }
    DT = data.table(DF[,,1])
    DT[,Deviation:=deviation]
    
    DT[,Run:=run]
    DT[,Beta := Betas[1]]
    dd = dist(plotsC)
    DT[,DistanceBtwPlots_median:=median(dd)][,DistanceBtwPlots_var:=var(dd)]
    gc() # Clear memory
    DT
    #with(DT[State=="Arrival"], plot(ResidenceTime~Traveltime))
    #list(Gains = Gains, Deviation = deviation, Decisions = Decisions)
  }))
  gc() # Clear memory
  Runs
})

BetaRuns2 = rbindlist(BetaRuns)
BetaRuns2[,Time:=1:.N, c("Run", "Beta")]
BetaRuns2[,Traveltime:=as.numeric(Traveltime)]
BetaRuns2[,ResidenceTime:=as.numeric(ResidenceTime)]
BetaRuns2[,Gain := as.numeric(Gain)]
#BetaRuns2[, Rates:=as.numeric(Rates)]
BetaRuns2[State=="Arrival", Rates := Gain/Traveltime ]
library(ggplot2)
ggplot(data = BetaRuns2[State=="Arrival"], mapping = aes(col = factor(Beta)))+
  geom_point(mapping = aes(x = Traveltime, y = (ResidenceTime)))

ggplot(data = BetaRuns2[State=="Arrival"], mapping = aes(col = factor(Beta)))+
  geom_boxplot(mapping = aes(x = factor(Deviation), y = (Rates)))+
  facet_wrap(~Beta)

ggplot(data = BetaRuns2[,.SD[(.N-1)], c("Beta", "Run")], mapping = aes(col = factor(Beta)))+
  geom_boxplot(mapping = aes(x = factor(Deviation), y = (Gain)))+
  facet_wrap(~Beta)


Runs[,gDeviation:=cut(Deviation, seq(-5,5, .25))]
########FOR LOOP#########

ggplot(data = BetaRuns2[State=="Arrival"], mapping = aes(col = factor(Beta)))+
  geom_point(mapping = aes(x = Traveltime, y = ResidenceTime))+
  geom_line(mapping = aes(x = Traveltime, y = ResidenceTime))

ggplot(data = BetaRuns2[,.SD[which.max(Gain)], c("Beta", "Deviation", "Run")],
       mapping = aes(col = factor(Beta)))+
  geom_boxplot(mapping = aes(x = factor(Deviation), y = Gain))+
  facet_wrap(~Beta)

ggplot(data = BetaRuns2[,.SD[which.max(Gain)], c("Beta", "Deviation", "Run")],
       mapping = aes(col = factor(Beta)))+
  stat_smooth(mapping = aes(x = DistanceBtwPlots_var, 
                          y = Gain), method = "lm")

ggplot(data = BetaRuns2[,.SD[which.max(Gain)], c("Beta", "Deviation", "Run")],
       mapping = aes(col = factor(Beta)))+
  stat_smooth(mapping = aes(x = DistanceBtwPlots_clustering, 
                            y = Gain), method = "lm")

ggplot(data = BetaRuns2[,.SD[which.max(Gain)], c("Beta", "Deviation", "Run")],
       mapping = aes(col = factor(Beta)))+
  stat_smooth(mapping = aes(x = DistanceBtwPlots_clustering, 
                            y = ResidenceTime), method = "lm")


ggplot(data = BetaRuns2[,.SD[which.max(Gain)], c("Beta", "Deviation", "Run")],
       mapping = aes(col = factor(Beta)))+
  stat_smooth(mapping = aes(x = DistanceBtwPlots_median, 
                            y = Gain), method = "lm")

df = BetaRuns2[,.(Gain=max(Gain,na.rm = T), ResidenceTime = mean(ResidenceTime)),
               c("Beta", "Deviation", "Run", "DistanceBtwPlots_clustering","DistanceBtwPlots_var","DistanceBtwPlots_median")]
library(lme4);library(MuMIn)
df[,CV:=DistanceBtwPlots_var/DistanceBtwPlots_median]
summary(df)
dredge(lmer(log(Gain)~CV*Beta + 
              ResidenceTime*Beta + 
              DistanceBtwPlots_clustering*Beta + 
              log(DistanceBtwPlots_median)*Beta + 
              Deviation + 
              (1|Run), data = df, na.action = na.fail, REML = F)
       #, fixed = c("Deviation")
       )


