#### BROWNIAN MOTION
# https://www.stat.berkeley.edu/~aldous/Research/Ugrad/ZY1.pdf
library(data.table);library(ggplot2)
library(boot)
library(mipfp)
library(cowplot)
#6.Brownian Motion with Boundary Coding
DiseaseMotions = function(SimLength = 500, nBoxes = 4, nIndRange = 3:10, ArenaDimRange = 10:50,
                          DispRate = .05, Virulence = .25, InitialDisease = .2, VaccineCoverage = .1, VaccineStrategy = NULL, PrintState=TRUE, PrintTime = TRUE,TimePlot = TRUE){
  library(data.table);library(ggplot2)
  library(boot)
  library(mipfp)
  library(cowplot)
  # SimLength = 100; nBoxes = 4; nIndRange = 3:10; ArenaDimRange = 10:50;  DispRate = .05; Virulence = .25
  N = SimLength
  # Number of arenas
  nBoxes
  # Individuals per arenas
  nInd = sample(nIndRange,nBoxes, replace = T)
  
  
  # Create a matrix which contains info on which arena each player (column) is in at a given time (row)
  BoxStatus = matrix(nrow = N, ncol = sum(nInd))
  BoxStatus[1,]<-rep(1:nBoxes, nInd)
  
  # Dimensions of the arenas
  Dims = sample(ArenaDimRange, nBoxes)
  # set up arenas
  ub = Dims; lb = -Dims
  
  #Movements
  dx = array(rnorm(n = N*sum(nInd)), c(N,sum(nInd)))
  dy = array(rnorm(n = N*sum(nInd)), c(N,sum(nInd)))
  
  # Positions
  xpos = array(rep(1,N*sum(nInd)), c(N,sum(nInd))); ypos = array(rep(1,N*sum(nInd)), c(N,sum(nInd)))
  
  # Create a matrix which contains info on which infection status of each player (column) is in at a given time (row)
  IStatus = matrix(nrow = N, ncol = sum(nInd))
  # sample infected individuals, and how long they should be sick
  DiseaseDuration = 5
  Is = rep(0, sum(nInd))
  Is[sample(1:sum(nInd), round(sum(nInd)*InitialDisease), replace = F)]<-1
  Is[Is==1]<-rpois(length(which(Is==1)), DiseaseDuration)
  IStatus[1,]<-Is
  
  
  #print(unique(BoxStatus[1,]))
  
  for(Time in 1:(N-1)){
    if(PrintTime){print(paste("Time:",Time))}
    if(Time == 1){
      # For the first round, randomly place the players inside their respective arenas
      for(Box in 1:nBoxes){
        xpos[1,BoxStatus[1,]==Box]<-dx[1,BoxStatus[1,]==Box]<-runif(n = length(which(BoxStatus[1,]==Box)), min = lb[Box], max = ub[Box])
        ypos[1,BoxStatus[1,]==Box]<-dy[1,BoxStatus[1,]==Box]<-runif(n = length(which(BoxStatus[1,]==Box)), min = lb[Box], max = ub[Box])
      }
    }
    
    # Vaccination
    if(PrintState){print("State: Vaccination")}
    
    grouper <- function(balls= 13, individuals = 10) {
      # create a random number for each row
      base = floor(balls/individuals)
      rest <- balls- individuals * base
      dis = rep(base, individuals)
      idx = sample(1:individuals, rest, replace = F)
      dis[idx] <- dis[idx]+(rest/length(idx))
      dis = as.integer(dis)
      return(dis)  
    }
    nVaccines0 = round(sum(nInd)*VaccineCoverage)
    
    if(VaccineStrategy=="Democratic"){
      print(nVaccines0)
      nVaccinesBox = grouper(nVaccines0, nBoxes)
      
    }
    if(VaccineStrategy=="Density"){
      dis = floor(nVaccines0*prop.table(table(BoxStatus[Time,])))
      rest = nVaccines0-sum(dis)
      idx = sample(1:nBoxes, rest, replace = F)
      dis[idx] <- dis[idx]+(rest/length(idx))
      nVaccinesBox = as.integer(dis)
    }
    
    # if box is allocated more vaccines than redistribute these among the others
    t = 0
    while(any(nVaccinesBox-table(BoxStatus[Time,])>0)){
      # Which boxes have surplus vaccines?
      idx = which(nVaccinesBox-table(BoxStatus[Time,])>0)
      
      nVaccinesBoxRe = grouper(sum((nVaccinesBox-table(BoxStatus[Time,]))[idx]), nBoxes-(length(idx)))
      nu = nVaccinesBox
      nu[idx] = nu[idx]+(table(BoxStatus[Time,])-nVaccinesBox)[idx] # Remove from surplus
      nu[-idx] = nu[-idx] + nVaccinesBoxRe # add to receiving
      nVaccinesBox<-nu
      
    }
    for(Box in 1:nBoxes){
      print(Box)
      # Randomly select which individuals within each box to get vaccinated
      IStatus[Time,sample(which(BoxStatus[Time,]==Box), nVaccinesBox[Box], replace = F)]<--9999
    }
    
    # Should we move any individuals?
    if(nBoxes>1){
      if(PrintState){print("State: Dispersal")}
      for(Box in 1:nBoxes){
        # sample how many individuals that disperse from Box, ensuring that we don't disperse more than are actually present
        Dispersers = min(c(length(which(BoxStatus[Time,]==Box))-1,
                           rpois(1, lambda = DispRate*sqrt(ub[Box]^2))))
        Dispersers = max(c(Dispersers,0)) # Make sure to leave one individual behind
        
        # Sample target arena dependent on potential target areas, i.e. more likely to go to large patches than small patches
        Target = sample(1:(nBoxes-1),Dispersers,prob = (ub[-Box]^2)/sum(ub[-Box]^2), replace = T)
        BoxStatus[Time,sample(which(BoxStatus[Time,]==Box), Dispersers, replace = F)]<-Target
      }
      
    }
    if(PrintState){print("State: Movement")}
    # Start looping
    for(Ind in 1:sum(nInd)){
      BoxStatus[Time+1,Ind]<-BoxStatus[Time,Ind]
      if((xpos[Time,Ind] + dx[Time+1,Ind]) > ub[BoxStatus[Time,Ind]]){
        xpos[Time+1,Ind] = ub[BoxStatus[Time,Ind]] # x above, move to upper boundary
      }else{
        if((xpos[Time,Ind] + dx[Time+1,Ind]) < lb[BoxStatus[Time,Ind]]){
          xpos[Time+1,Ind] = lb[BoxStatus[Time,Ind]] # x below, move to lower boundary
        }else{
          xpos[Time+1,Ind] = xpos[Time,Ind] + dx[Time+1,Ind] # x within boundaries, move as planned
          if((ypos[Time,Ind] + dy[Time+1,Ind]) > ub[BoxStatus[Time,Ind]]){
            ypos[Time+1,Ind] = ub[BoxStatus[Time,Ind]] # y above, move to upper boundary
          }else{
            if((ypos[Time, Ind] + dy[Time + 1,Ind]) < lb[BoxStatus[Time,Ind]]){
              dy[Time+1,Ind] = lb[BoxStatus[Time,Ind]] # y below, move to lower boundary
            }else{
              ypos[Time+1,Ind] = ypos[Time,Ind] + dy[Time+1,Ind]# y within boundaries, move as planned
            }}}}
      #}
    }
    
    
    # # Infection status
    if(PrintState){print("State: Diseasespreading")}
    for(Box in 1:nBoxes){
      # First, ask if there are any sick individuals in the target arena
      # And there should be more than one individual present in the arena
      if(any(IStatus[Time,BoxStatus[Time,]==Box]>=Time) & length(which(BoxStatus[Time,]==Box))>1){
        # Calculate among individual distances
        tmp = as.matrix(dist(cbind(xpos[Time,BoxStatus[Time,]==Box],
                                   ypos[Time,BoxStatus[Time,]==Box]), diag = T, upper = T))
        
        
        # Who gets sick according to distance?
        mat = inv.logit(2 - Virulence*tmp)
        Infected = apply(as.matrix(mat), c(1,2), function(x){sample(1:0, 1, prob = c(x,1-x))})
        # Does any get infected?
        if(any(Infected[lower.tri(Infected)]>0)){
          # extract the ones who gets infected this round, i.e. remove those who are sick from before
          I2 = Infected[ which(IStatus[Time,BoxStatus[Time,]==Box]>=Time),]
          if(any(I2)>0){
            if(!is.null(dim(I2))){
              NuSick = which(sapply(which(apply(I2, 2, function(x) any(x>0))),
                                    function(x) !x %in% which(IStatus[Time,BoxStatus[Time,]==Box]>=Time)))
            }else{
              NuSick = which(sapply(which(I2>0),
                                    function(x) !x %in% which(IStatus[Time,BoxStatus[Time,]==Box]>=Time)))
            }
            # Locate the ones who gets infected, and insert the time until they are well again
            if(length(NuSick)>0){
              IStatus[Time,which(BoxStatus[Time,]==Box)[NuSick]]<-Time + rpois(n = length(NuSick), DiseaseDuration)
            }
          }
          
          
          
        }
      }
    }
    
    # Move the infection statuses forward
    IStatus[Time+1,]<-IStatus[Time,]
    
    # Should we plot per time step?
    if(TimePlot){
      par(mfrow = c(2,2), omi = c(0,0,0,0), mar = c(0,0,0,0))
      for(Box in 1:nBoxes){
        plot(x = xpos[Time,which(BoxStatus[Time,]==Box)],
             y = ypos[Time,which(BoxStatus[Time,]==Box)],
             xlim = c(-max(ArenaDimRange), max(ArenaDimRange)), ylim = c(-max(ArenaDimRange), max(ArenaDimRange)),
             xlab = "", ylab = "",
             col = ifelse(IStatus[Time,which(BoxStatus[Time,]==Box)]>=Time, "red",
                          ifelse(IStatus[Time,which(BoxStatus[Time,]==Box)]<0, "green", "black")),
             cex = 2,
             pch = 19)
        if(Box==1){
          legend("topleft",
                 c("Frisk", "Syk", "Vaksinert"),
                 pch = 19,
                 col = c("black", "red", "green"))
        }
        legend("bottomleft",
               c(paste("Vaksinert: ",length(which(IStatus[Time,]<0))),
                 paste("Syk: ",length(which(IStatus[Time,]<0))))
        )
        
        
      }
    }
  }
  
  #print(unique(BoxStatus[1,]))
  dt = rbindlist(lapply(1:sum(nInd), function(x){
    dt = data.table(Ind = x, posX = xpos[,x], posY = ypos[,x], Box = BoxStatus[,x])
    dt[,Time:=1:nrow(dt)]
    dt[,lb:=lb[BoxStatus[,x]]][,ub:=ub[BoxStatus[,x]]][,Sick:=IStatus[,x]]
    dt
  }))
  
  dt[,SickNow:=ifelse(Time<=Sick,"Yes", "No")]
  dt[,nInd:=uniqueN(Ind), c("Time", "Box")]
  dt[,Prevalence:=length(which(.SD$SickNow=="Yes"))/nInd, c("Time", "Box")]
  setorder(dt, Time, Ind)
  return(dt)
}

tmp1a = DiseaseMotions(SimLength = 30, nBoxes = 1, DispRate = 0.5, nIndRange = 10, 
                       VaccineCoverage = 0.1, VaccineStrategy = "Democratic")

tmp1b = DiseaseMotions(SimLength = 30, nBoxes = 4, DispRate = 0.5, nIndRange = 10, 
                       VaccineCoverage = 0.1, VaccineStrategy = "Density")




tmp1a = DiseaseMotions(SimLength = 200, nBoxes = 4, DispRate = 0.5, nIndRange = 10)
tmp1b = DiseaseMotions(SimLength = 200, nBoxes = 4, DispRate = 0.25, nIndRange = 10)

tmp2a = DiseaseMotions(SimLength = 200, nBoxes = 8, DispRate = 0.5, nIndRange = 10)
tmp2b = DiseaseMotions(SimLength = 200, nBoxes = 8, DispRate = 0.25, nIndRange = 10)

# ISSUES:
# we sometimes loose arenas!
# is it because everyone moves away immediately?
#
# and sometimes individuals tend to get lost in the south
library(ggplot2)
library(cowplot)
P1a = ggplot(data = tmp1a[,.SD[1], c("Time", "Box")], 
             mapping = aes(x = Time, y = Prevalence))+
  geom_line()+
  facet_grid(.~Box)+
  theme(axis.text.x = element_blank())

S1a = ggplot(data = tmp1a, mapping = aes(x = posX, y = posY, col = factor(Ind),xmin = lb, ymin = lb, xmax = ub, ymax = ub))+
  geom_point()+facet_grid(.~Box)+geom_rect(fill = NA, col = "black")+ theme(legend.position = "none")+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + xlab("") + ylab("")


P1b = ggplot(data = tmp1b[,.SD[1], c("Time", "Box")], 
             mapping = aes(x = Time, y = Prevalence))+
  geom_line()+
  facet_grid(.~Box)+
  theme(axis.text.x = element_blank())

S1b = ggplot(data = tmp1b, mapping = aes(x = posX, y = posY, col = factor(Ind),xmin = lb, ymin = lb, xmax = ub, ymax = ub))+
  geom_point()+facet_grid(.~Box)+geom_rect(fill = NA, col = "black")+ theme(legend.position = "none") + xlab("") + ylab("")

P2a = ggplot(data = tmp2a[,.SD[1], c("Time", "Box")], 
             mapping = aes(x = Time, y = Prevalence))+
  geom_line()+
  facet_grid(.~Box)+
  theme(axis.text.x = element_blank())

S2a = ggplot(data = tmp2a, mapping = aes(x = posX, y = posY, col = factor(Ind),xmin = lb, ymin = lb, xmax = ub, ymax = ub))+
  geom_point()+facet_grid(.~Box)+geom_rect(fill = NA, col = "black")+ theme(legend.position = "none")+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + xlab("") + ylab("")


P2b = ggplot(data = tmp2b[,.SD[1], c("Time", "Box")], 
             mapping = aes(x = Time, y = Prevalence))+
  geom_line()+
  facet_grid(.~Box)+
  theme(axis.text.x = element_blank())

S2b = ggplot(data = tmp2b, mapping = aes(x = posX, y = posY, col = factor(Ind),xmin = lb, ymin = lb, xmax = ub, ymax = ub))+
  geom_point()+facet_grid(.~Box)+geom_rect(fill = NA, col = "black")+ theme(legend.position = "none") + xlab("") + ylab("")


windows()
plot_grid(P1a,S1a, P1b, S1b, P2a, S2a, P2b, S2b, byrow = F, nrow = 2)






par(mfrow = c(1,nBoxes))
plot(1,1, type = "n", xlim = c(-20,20), ylim = c(-20,20))
rect(lb[1], lb[1], ub[1], ub[1])
lines(ypos[,1]~xpos[,1], type = "b", col = 1)


plot(1,1, type = "n", xlim = c(-20,20), ylim = c(-20,20))
rect(lb[2], lb[2], ub[2], ub[2])
lines(ypos[,2]~xpos[,2], type = "b", col = 2)
lines(ypos[,3]~xpos[,3], type = "b", col = 3)



par(mfrow = c(1,nBoxes))
for(i in 1:nBoxes){
  plot(1,1, type = "n", xlim = c(-20,20), ylim = c(-20,20))
  rect(lb[i], lb[i], ub[i], ub[i])
  
  
  for(N in which(BoxStatus[1,]==i)){
    lines(ypos[,BoxStatus[1,]==i]~xpos[,BoxStatus[1,]==i], type = "b", col = N)
  }
  
  
  
}




Sims = vector("list", length = 100)
for(x in 1:100){
  print(paste("Run: ",x,"a"))
  if(is.null(Sims[[x]]){
    a = Sys.time()
    tmp2a = DiseaseMotions(SimLength = 200, nBoxes = 6, DispRate = 0.5, nIndRange = 60)
    tmp2[,Run:=x][,DispRate := .5]
    
    tmp2b = DiseaseMotions(SimLength = 200, nBoxes = 6, DispRate = 0.25, nIndRange = 60)
    tmp2b[,Run:=x][,DispRate := .25]
    
    print(paste("Run: ",x,"b"))
    tmp2c = DiseaseMotions(SimLength = 200, nBoxes = 6, DispRate = 0.125, nIndRange = 60)
    tmp2c[,Run:=x][,DispRate := .125]
    tmp2d = DiseaseMotions(SimLength = 
                           200, nBoxes = 6, DispRate = 0.0625, nIndRange = 60)
    tmp2d[,Run:=x][,DispRate := .0625]
    
    dt = rbind(tmp2a, tmp2b, tmp2c, tmp2d, fill = T)
    Sims[[x]]<-dt
    print(Sys.time()-a)
  }
}


SimDT = rbindlist(Sims)

SimDT[,Density:=nInd/(ub^2)]
setorder(SimDT, Run, Box, Time)

#SimDT[,dtPrevalence:=NULL]
dt = SimDT[,.SD[1], c("Run", "Box","Time")][,c("Run", "Box","Time", "Prevalence"),with=F][,dtPrevalence:=c(NA,diff(Prevalence))]
SimDT = dt[SimDT, on = c("Run", "Box","Time")]
ggplot(data = SimDT[!is.na(DispRate)][,.SD[1], c("Run", "Box","Time")], aes(x = Density, y = dtPrevalence, col = factor(DispRate)))+
  geom_point()+
  stat_smooth(method = "lm")

ggplot(data = SimDT[!is.na(DispRate)][,.SD[1], c("Run", "Box","Time")], aes(x = Density, y = Prevalence, col = factor(DispRate)))+
  geom_point()+
  stat_smooth(mmethod="glm", family="binomial")+
  scale_x_log10()



# *******************************
# BROWNIAN MOTION SIMULATION
# December 2012 | Benjamin Tovar
# *******************************
#
#   REFERENCES
#   http://landshape.org/enm/r-code-for-brownian-motion/
#
#   According to Wikipedia the mathematical model for Brownian motion 
#   (also known as random walks) can also be used to describe many 
#   phenomena as well as the random movements of minute particles, 
#   such as stock market fluctuations and the evolution of physical 
#   characteristics in the fossil record. The simple form of the 
#   mathematical model for Brownian motion has the form:
#
#    S_t = eS_t-1
#
#    where e is drawn from a probability distribution.
#
#######################################################################
brownian <- function(n.times = 10, nboxes = 1){
  x <- y <- x.new <- y.new <- x.new.p <- y.new.p <- vector()
  # Patch sizes
  dims = sample(1:100, nboxes)
  # initial densities
  Dens = sample(2:3, nboxes)
  
  for(i in 1:n.times){
    
    
    for(p in 1:nboxes){
      # Initialize variables
      if(i ==1){
        x = runif(n = Dens[p], min = -dims[p], max = dims[p])
        y = runif(n = Dens[p], min = -dims[p], max = dims[p])
      }else{
        x <- rnorm(n = Dens[p])
        y <- rnorm(n = Dens[p])
      }
      # concatenate variables 
      # to increase the vector size
      x.new <- rbind(x.new,x)
      y.new <- rbind(y.new,y)
      # sum the vector numbers
      x.new.p <- colSums(x.new)
      y.new.p <- colSums(y.new)  
      
      # # plot the model
      # plot(x.new.p,y.new.p,type="b",
      #      main=paste("Brownian motion simulation in R\nTime =",i,sep=" "),
      #      xlab="x coordinates",ylab="y coordinates",
      #      col=c(rep("gray",i-1),"red"),
      #      pch=c(rep(20,i-1),1))    
    }
  }
  
}
melt(as.data.table(y.new)[,Time:=1:nrow(y.new)], id.vars = "Time", value.name = "posY", variable.name = "ID"
)[melt(as.data.table(x.new)[,Time:=1:nrow(x.new)], id.vars = "Time", value.name = "posX", variable.name = "ID")
  , on = c("ID", "Time")]

plot(1,1, type = "n", xlim = c(-20,20), ylim = c(-20,20))

# Test the function
brownian(10)
a
# ****************************************
# EXPORT BROWNIAN MOTION SIMULATION IMAGES
# ****************************************
export.brownian <- function(n.times){
  x <- y <- x.new <- y.new <- x.new.p <- y.new.p <- vector()
  for(i in 1:n.times){
    # Initialize variables
    x <- rnorm(1)
    y <- rnorm(1)
    # concatenate variables to increase the
    # vector size
    x.new <- c(x.new,x)
    y.new <- c(y.new,y)
    # sum the vector numbers
    x.new.p <- cumsum(x.new)
    y.new.p <- cumsum(y.new)  
    # plot the model
    #png(paste("image",i,"png",sep="."),width=600,height=600)
    plot(x.new.p,y.new.p,type="b",
         main=paste("Brownian motion simulation in R\nTime =",
                    i,sep=" "),
         xlab="x coordinates",ylab="y coordinates",
         col=c(rep("gray",i-1),"red"),
         pch=c(rep(20,i-1),1))
    #cat("image",i,"DONE",date(),"\n")
    #dev.off()
  }
}
# Test the function
a = export.brownian(10)
#######################
n = 1000
t = 100
No.Ex = 10
steps = seq(0,t,length=n+1)
A = replicate(No.Ex, {
  bm <- c(0, cumsum(rnorm(n,0,sqrt(t/n))))
}) 

cols = rainbow(No.Ex)
matplot(A, type = "l", col = cols, lty = 1)


Walk = function(n = 100, alpha = 2, fmin = .1, fmax = 100, type = "Levy"){
  alpha=2;
  x=rep(0,n)
  y=rep(0,n)
  
  for (i in 2:n){
    theta <- runif(n-1)*2*pi
    if(type=="Levy"){
      f <- runif(n-1, fmax^(-alpha), fmin^(-alpha))^(-1/alpha)
    }else if(type=="Brownian"){
      f <- rnorm(n-1)
    }
    
    x <- c(0, cumsum(f*cos(theta)))
    y <- c(0, cumsum(f*sin(theta)))
  }
  cbind(x, y)
}
df = Walk(n = 1000, type = "Levy")
df2 = Walk(n = 1000, type = "Brownian")
plot(df[,2]~df[,1], type ="b", xlim = range(c(df[,1], df2[,1])), ylim = range(c(df[,2], df2[,2])))
points(df2[,2]~df2[,1], type ="b", col = "red")
legend("topleft",
       c("Lévy walk", "Brownian motion"),
       col = c("black", "red"),
       pch = c(1,1))

########### Foraging motion #####
AClandscape = function(d, nLandscapes = 1, SILL = 0.025, RANGE = 5, ValueRange = NULL, Skew = 1, Model = "Gau", binary = F, binarythreshold = .4){
  # Binary-threshold: values below this gets set to zero
  #binary: divide into binary landscape
  require(gstat)
  require(sp)
  #cat(paste(r, x,";"))
  #set.seed(x) # to ensure the same seed for each scale
  xy <- expand.grid(1:d, 1:d)
  names(xy) <- c('x','y')
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                   model=vgm(psill=SILL, range=RANGE, model= Model), nmax=20)
  yy <- predict(g.dummy, newdata=xy, nsim=nLandscapes)
  yy[,3:(ncol(yy):nLandscapes)]<-yy[,3:(ncol(yy):nLandscapes)]^Skew
  gridded(yy) = ~x+y
  if(binary){yy@data[,1] = ifelse(yy@data[,1]<quantile(yy@data[,1],binarythreshold),0,1)}
  if(binarythreshold<1){yy@data[,1] = ifelse(yy@data[,1]<quantile(yy@data[,1],binarythreshold),0,yy@data[,1])}
  
  if(!is.null(ValueRange)){
    ys = (yy@data[,1]-min(yy@data[,1]))/(max(yy@data[,1])-min(yy@data[,1]))
    yy@data[,1] <- ValueRange[1] + diff(ValueRange)*ys
  }
  as.matrix(yy)
}


ForagingMotions = function(SimLength = 500, nIndRange = 3:10, AttackRate = 20, ArenaDimRange = 10:50, Torus = T,
                           DispRate = .05, nBoxes = 1, Virulence = .25, InitialDisease = .2, VaccineCoverage = .1, VaccineStrategy = NULL, PrintState=TRUE, PrintTime = TRUE,TimePlot = F){
  library(data.table)
  # Torus = F;SimLength = 100; nBoxes = 1; nIndRange = 3:10; ArenaDimRange = 10:50;  DispRate = .05; Virulence = .25
  N = SimLength
  # Number of arenas
  nBoxes = 1
  # Individuals per arenas
  nInd = sample(nIndRange,nBoxes, replace = T)
  
  # Create a matrix which contains info on which arena 
  # each player (column) is in at 
  # a given time (row)
  BoxStatus = matrix(nrow = SimLength, ncol = sum(nInd))
  BoxStatus[1,]<-rep(1:nBoxes, nInd)
  
  EnergyStatus = matrix(rep(0,nInd),nrow = SimLength, ncol = sum(nInd))
  # How long an individual should stay at a given point
  WaitTime = matrix(nrow = SimLength, ncol = sum(nInd))
  
  
  # Dimensions of the arenas
  Dims0 = sample(ArenaDimRange, nBoxes)
  if(!Torus){Dims <- Dims0*4}else{Dims <- Dims0}
  # set up arenas
  ub = round(Dims/2); lb = -round(Dims/2)
  # Create arena
  #Arena = matrix(rlnorm(n = Dims^2, meanlog = 2), ncol = Dims, nrow = Dims)
  #Arena = (Arena - min(Arena))/(max(Arena)-min(Arena))
  Arena = AClandscape(d = Dims, RANGE = 8, ValueRange = c(10,20))
  Arena0 = Arena
  AttackRate = AttackRate
  
  #Movements
  dx = array(rnorm(n = SimLength*sum(nInd)), c(SimLength,sum(nInd)))
  dy = array(rnorm(n = SimLength*sum(nInd)), c(SimLength,sum(nInd)))
  
  # Positions
  xpos = array(rep(1,SimLength*sum(nInd)), c(SimLength,sum(nInd)))
  ypos = array(rep(1,SimLength*sum(nInd)), c(SimLength,sum(nInd)))
  
  for(Time in 1:(SimLength-1)){
    if(PrintTime){print(paste("Time:",Time))}
    if(Time == 1){
      # For the first round, randomly place the players inside their respective arenas
      for(Box in 1:nBoxes){
        xpos[1,BoxStatus[1,]==Box]<-dx[1,BoxStatus[1,]==Box]<-runif(n = length(which(BoxStatus[1,]==Box)), min = -Dims0[Box], max = Dims0[Box])
        ypos[1,BoxStatus[1,]==Box]<-dy[1,BoxStatus[1,]==Box]<-runif(n = length(which(BoxStatus[1,]==Box)), min = -Dims0[Box], max = Dims0[Box])
      }
    }
    if(PrintState){print("State: Movement")}
    # Start looping
    for(Ind in 1:sum(nInd)){
      if(Torus){
        BoxStatus[Time+1,Ind]<-BoxStatus[Time,Ind]
        if((xpos[Time,Ind] + dx[Time+1,Ind]) > ub[BoxStatus[Time,Ind]]){
          xpos[Time+1,Ind] = ub[BoxStatus[Time,Ind]] # x above, move to upper boundary
        }else{
          if((xpos[Time,Ind] + dx[Time+1,Ind]) < lb[BoxStatus[Time,Ind]]){
            xpos[Time+1,Ind] = lb[BoxStatus[Time,Ind]] # x below, move to lower boundary
          }else{
            xpos[Time+1,Ind] = xpos[Time,Ind] + dx[Time+1,Ind] # x within boundaries, move as planned
            if((ypos[Time,Ind] + dy[Time+1,Ind]) > ub[BoxStatus[Time,Ind]]){
              ypos[Time+1,Ind] = ub[BoxStatus[Time,Ind]] # y above, move to upper boundary
            }else{
              if((ypos[Time, Ind] + dy[Time + 1,Ind]) < lb[BoxStatus[Time,Ind]]){
                dy[Time+1,Ind] = lb[BoxStatus[Time,Ind]] # y below, move to lower boundary
              }else{
                ypos[Time+1,Ind] = ypos[Time,Ind] + dy[Time+1,Ind]# y within boundaries, move as planned
              }}}}
        
        
      }else{
        xpos[Time+1,Ind] = xpos[Time,Ind] + dx[Time+1,Ind] 
        ypos[Time+1,Ind] = ypos[Time,Ind] + dy[Time+1,Ind]
        # Obtain energy
        
        
        if(dx[Time+1,Ind]>0 | dx[Time+1,Ind]>0){ # if you have moved, how long should you stay here?
          ny = round((Dims/2)-ypos[Time,Ind])
          nx =  round((Dims/2)-xpos[Time,Ind])
          # If inside the arena
          if(ny %between% c(0, Dims) & nx %between% c(0, Dims)){
            # Basert på Aastrøm 1990 eq. 1
            # G(t) = Si(1- (1/(1 + kt)))
            Si = Arena[ny,nx]
            
            ### MARGINAL VALUE THEOREM
            # Use the 5 latest travel distances (i.e. travel costs, 5 is arbitraily chosen)
            nT = min(c(Time,5))
            k = 1.5
            TravelCost = sum(sqrt(diff(dx[(Time-nT):Time, Ind])^2+diff(dy[(Time-nT):Time, Ind])^2))
            tStar = sqrt(TravelCost/k)
            
            ResidenceTime = tStar#1*round(rlnorm(n = 1, meanlog = 2)+0*Arena[ny,nx])
            ResidenceTime = min(c(Time + ResidenceTime, SimLength))
            
            dx[Time:ResidenceTime,Ind]<-0
            dy[Time:ResidenceTime,Ind]<-0
      
          }
        }
        ## EAT!!!
        ## FUNCTIONAL RESPONSE
        # Disc equation look-a-like, similar approach as Calcagno et al. 2014
        Foraged = Arena[ny,nx]^0.5
        Arena[ny,nx] = Arena[ny,nx] - Foraged
        EnergyStatus[Time + 1, Ind] <- EnergyStatus[Time, Ind] + Foraged
      }
    }
    
    #Metabolic cost
    EnergyStatus[Time+1,] <- EnergyStatus[Time+1,]-0.05
    
    # Arena regrowth
    Arena[Arena<0]<-0.0001
    Arena  = Arena + 1.1*((Arena0-Arena)/Arena0)
  }
  
  #print(unique(BoxStatus[1,]))
  dt = rbindlist(lapply(1:sum(nInd), function(x){
    dt = data.table(Ind = x, 
                    posX = xpos[,x], posY = ypos[,x], Energy = EnergyStatus[,x])
    dt[,Time:=1:nrow(dt)]
    dt[,lb:=lb[BoxStatus[,x]]][,ub:=ub[BoxStatus[,x]]]
    dt
  }))
  dt[,dist:=c(sqrt((posX-shift(posX))^2 + (posY-shift(posY))^2)), "Ind"]
  attr(dt, "Arena")<-Arena0
  setorder(dt, Time, Ind)
  return(dt)
}


tmp1a = ForagingMotions(SimLength = 3000, AttackRate = 100,
                        ArenaDimRange = 100,
                        nIndRange = 50:50, Torus = F)
summary(tmp1a)

ggplot(data = tmp1a, aes(x = Time, y = Energy, col = factor(Ind)))+
  geom_line()
tmp1a[,UniquecellsPerInd:=uniqueN(paste(round(posX), round(posY))), "Ind"]
ggplot(data = tmp1a[,.SD[.N], "Ind"], aes(x = Energy, y = log(UniquecellsPerInd)))+
  geom_point()+ stat_smooth(method = "lm")


tmp1a
Arena1a = attr(tmp1a, "Arena")
tmp1a[,Dist:=c(NA,sqrt(diff(posX)^2 + diff(posY)^2)), "Ind"]
tmp1a[, TravelCost:=frollsum(Dist, n = 5), "Ind"]

curve(sqrt(x/k), from = min(tmp1a$TravelCost, na.rm = T), to = max(tmp1a$TravelCost, na.rm = T), ylab = "ResidenceTime", xlab = "TravelCost")

library(ggplot2)
d2.df <- data.frame(x=rep(as.numeric(1:ncol(Arena1a)),each=ncol(Arena1a)),
                    y=rep(as.numeric(1:nrow(Arena1a)),times=nrow(Arena1a)),
                    z=as.vector(Arena1a))
scale = unique(na.omit(tmp1a$ub));scale
ggplot()+
  geom_tile(data=d2.df, mapping = aes(x=x,y=y,fill=z))+
  geom_point(data = tmp1a, mapping = aes(x = posX+scale, 
                                         y = posY+scale, 
                                         col = factor(Ind)))

tmp1a[,MSD:=sqrt((posX-posX[1])^2+(posY-posY[1])^2), "Ind"][,MSD0:=sqrt(Time), "Ind"]


# # Yes it is random
# quantile(tmp1a[,.SD[.N], "Ind"
# ][,.(RW = MSD-MSD0)]$RW, c(0.025,.5,.975))


ggplot(data = tmp1a, aes(x = Time, y = MSD, col = factor(Ind)))+
  geom_point()+stat_smooth(method = "lm")

ggplot()+geom_point(data = tmp1a, mapping = aes(x = posX, y = posY, 
                                      col = factor(Ind)))






# nTimes = 100
# DF = data.table(posX = 1, posY = 1, Time = 1:nTimes, Ind = letters[1])
# DF[,dx:=rnorm(n = nTimes)][,dy:=rnorm(n = nTimes)]
# DF[Time==1, posX:=rnorm(n = nrow(DF[Time==1]))]
# DF[Time==1, posY:=rnorm(n = nrow(DF[Time==1]))]
# DF[Time>1,posX := dx + shift(posX), "Ind"
#    ][Time>1,posY := dy + shift(posY), "Ind"]
# DF[,MSD:=sqrt((posX-posX[1])^2+(posY-posY[1])^2), "Ind"]
# with(DF, plot(MSD~Time))
