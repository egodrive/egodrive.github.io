fourPL = function(r = 1.01, K = 100, L = 50,tmax = 25){
  #K = asymptote
  # L  = lower asymptote
  times = seq(0,tmax, length.out = 100)
  
  P = K - L - M0
  # change in biomass
  dMdT = r*(M-L)*((K-M)/(K-L))
  #Biomass at time t
  biomass = L + (M0*(K-L))/(M0 + P*exp(-r*times))
  # AGR at time t
  AGR = (r*M0*(K-L)*P*exp(-r*times))/(M0 + P*exp(-r*times))^2
  par(mfrow = c(2,2))
  plot(biomass~times, type = "l")
  plot(AGR~biomass, type = "l")
  plot(dMdT~times, type = "l")
  plot(AGR~times, type = "l")
}
# Generalised logistic function
# https://en.wikipedia.org/wiki/Generalised_logistic_function
# http://www.metla.fi/silvafennica/full/sf33/sf334327.pdf
fourPL2 = function(r = 1.01,
                   InflectionPoint = 1,
                   K = 1, L = 0,
                   v = 1,
                   C = 1,
                   tmax = 10,
                   Mass = seq(0,1, length.out = 50),
                   out = "MDeltaM"){
  # K = asymptote
  # L  = lower asymptote
  times = seq(0,tmax, length.out = 1000)
  # K+(L-K)/(1+(x/InflectionPoint)^r)
  #Biomass at time t
  biomass =  L + ((K-L)/((C + InflectionPoint*exp(-r*times))^(1/v)))
  
  # absolute growth rate at time t
  # Found the derivative using 'deriv()'
  AGR = (K - L) * ((C + InflectionPoint * (exp(-r * times)))^((1/v) - 1) * 
                     ((1/v) * (InflectionPoint * ((exp(-r * times)) * r))))/((C + InflectionPoint * (exp(-r * times)))^(1/v))^2
  
  
  if(out=="MDeltaM"){
    # but what is AGR(biomass) = ?
    if(C ==1 & v==1){
      dMdT = function(M,...){r*(M-L)*((K-M)/(K-L))}
      MDeltaM = dMdT(M = Mass, r = r, K = K, L = L)
      MDeltaM
    }else{
      dMdT = function(M,...){
        #mass = seq(.1,25, length.out = 500)
        #M = mass;K = 26;L = 0.1;v = 1;C = 1; B = 1.01;Q = 8
        #M = biomass
        exp1 = (M - L)/(K - L)
        exp2 = (exp1^-v)-C
        exp3 = exp2/InflectionPoint
        time = log(exp3)/-r
        # time can now be inserted into the AGR function above
        
        EXP4 = InflectionPoint * exp(-r * log(((((M - L)/(K - L))^-v)-C)/InflectionPoint)/(-r))
        # EXP4 simplifies to, i.e. the inflection point has no influence
        EXP4 = ((((M-L)/(K-L))^-v)-C)
        
        AGR = (K - L) * ((C + EXP4)^((1/v) - 1) * 
                           ((1/v) * (EXP4 * r)))/((C + EXP4)^(1/v))^2
        # Further simplifications
        AGR = (EXP4*r*(K-L)*(EXP4 + C)^(-(1/v)-1))/v
        
      }
      MDeltaM = dMdT(M = Mass, r = r, K = K, L = L, C = C, v = v, InflectionPoint = InflectionPoint)
      
      
    }
  }
  
  par(mfrow = c(1,3))
  plot(biomass~times, type = "l")
  plot(AGR~biomass, type = "l")
  ls = list(times = times, AGR = AGR, Biomass = biomass)
  if(out == "MDeltaM"){
    plot(MDeltaM~Mass, type = "l")
    attr(ls, "MDeltaM")<-MDeltaM
  }
  ls
}
library(data.table)

data.table(v = seq(0,1, length.out = 50))

datas = rbindlist(lapply(seq(.25,2, length.out = 10),function(x){
  aa = fourPL2(v = x,out = "MDeltaM")
  data.table(v = x,
             M = seq(0,1, length.out = 50), 
             dM = attr(aa, "MDeltaM"))
}))

datas2 = rbindlist(lapply(seq(.25,2, length.out = 10),function(x){
  aa = fourPL2(v = x)
  data.table(v = x,
             Time = aa$times, 
             Mass = aa$Biomass)
}))
datasR = rbindlist(lapply(seq(1.01,2, length.out = 10),function(x){
  aa = fourPL2(r = x,out = "MDeltaM")
  data.table(r = x,
             M = seq(0,1, length.out = 50), 
             dM = attr(aa, "MDeltaM"))
}))

library(ggplot2);library(cowplot)
qplot(data = datas, x = M, y = dM, col = v, group = v, geom = "line")
qplot(data = datasR, x = M, y = dM, col = r, group = r, geom = "line")
qplot(data = datas2, x = Time, y = Mass, col = v, group = v, geom = "line")


### PREAMBLE #####
getwd()
.libPaths(c("C:/Bibliotek_jobb/SpatialScaleInfluence/RScripts/RPakker",.libPaths()))
#install.packages(c("gstat","raster","MuMIn","minqa", "nloptr","KernSmooth", "spam", "fields", "MASS", "mgcv", "lme4", "MuMIn", "lattice", "gamm4", "gee", "sp", "spatstat", "adehabitatHR", "adehabitatLT", "data.table", "cowplot", "stringr","raster", "rgdal", "sp", "maptools", "circular","gstat"))
#install.packages(c("cowplot", "data.table"))

# Libraries
library(KernSmooth)
library(spam)
library(fields)
library(MASS)
library(mgcv)
library(lme4)
library(lattice)
library(gamm4)
library(gee)
library(sp)
library(spatstat)
library(stringr)
library(circular)
library(gstat)

#####        GRASS-HERBIVORE SIMULATION         #####
MoveSim = function(xy,
                   CenterRelease = 3, # grid cells away from center animasl are releasd
                   CentreOffset = c(0,0), # should we move the centre that we relase the animals from? x, y position that move the centre with
                   tmax = 100,
                   tburning = 1,
                   N0 = 2,
                   ss = 5,
                   concentration = 0,
                   stepDistribution = c(1.0001,3),
                   d = 20, 
                   P0 = 50, 
                   TakeShelter = TRUE,
                   ForageRate = 50,# 
                   BMR = 10, #Energy loss per time step
                   telp = 1, telq = 1,failp=0, 
                   MoistureLayer = NULL,
                   AcidLayer = NULL,
                   Boundary = "torus",
                   VarGrowthRate = NULL, # If "Spatial", make the growth rate vary spatially btw 1.001-1.5
                   ...){
  #tmax = 500; N0 = 1; SL = c(1.1,3.5); d = 50; Boundary = "torus"; ForageRate = 50
  #tmax: Time units
  #tburnin: Initialisation time units
  #N0: total population size of roe deer
  #P0: initial grass density
  # stepDistribution = c(1.0001,3)# range of values which describes the variation in step length among animals. Only matters when stepCost is non-linear
  # telp = 1;telq = 1;failp=0;PerceptionLength = 2;concentration = 0;MemoryDecay = 0.00001; d = 20; tmax = 100; N0 = 1; ss = 5; SL = c(1.0001,3); P0 = 50; ForageRate = 50;Boundary="torus";MoistureLayer = NULL;AcidLayer = NULL;TakeShelter = TRUE
  library(KernSmooth)
  library(spam)
  library(fields)
  library(MASS)
  library(mgcv)
  library(lme4)
  library(lattice)
  library(gamm4)
  library(gee)
  library(sp)
  library(spatstat)
  library(stringr)
  library(circular)
  library(gstat)
  
  # Helper function for torusoid
  TorusOider = function(Mat, coords){
    Cntr = matrix(coords, ncol = 2)
    dd = dim(Mat)
    cc = round(median(c(1,dd[1])))
    Shftr = cc-Cntr[1,1]
    Rows = 1:dd[1]
    Rows = binhf::shift(Rows,abs(Shftr), dir = ifelse(Shftr<0,"left", "right"))
    
    ShftC = cc-Cntr[1,2]
    Cols = 1:dd[2]
    Cols = binhf::shift(Cols,abs(ShftC), dir = ifelse(ShftC<0,"left", "right"))
    
    Mat[Rows,Cols]
  }
  
  
  # Function to generate autocorrealted landscape
  rng<-(d+1):(2*d)
  AClandscape = function(d, nLandscapes = 1, SILL = 0.025, RANGE = 5, ValueRange = NULL, Skew = 1, Model = "Gau", binary = F, binarythreshold = .9){
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
    if(!is.null(ValueRange)){
      ys = (yy@data[,1]-min(yy@data[,1]))/(max(yy@data[,1])-min(yy@data[,1]))
      yy@data[,1] <- ValueRange[1] + diff(ValueRange)*ys
    }
    as.matrix(yy)
  }
  
  if(is.null(MoistureLayer)){ # Allow to use predesignated env. layers in case we want repeated runs in similar environments
    M = AClandscape(d, RANGE = 8, ValueRange = c(10,50))
  }else{
    M = MoistureLayer
    M = matrix(M[which(!is.na(M), arr.ind = T)], ncol = d, nrow = d)
  }
  if(is.null(AcidLayer)){
    pH = AClandscape(d, RANGE = 8, ValueRange = c(4,8))
    pH = matrix(pH[which(!is.na(pH), arr.ind = T)], ncol = d, nrow = d)
  }else{pH = AcidLayer}
  
  if(is.null(VarGrowthRate)){
    GrowthRate = 1.0001
  }else if(VarGrowthRate=="Spatial"){
    GrowthRate = AClandscape(d, RANGE = 8, ValueRange = c(1.0001,2))
    GrowthRate = matrix(GrowthRate[which(!is.na(GrowthRate), arr.ind = T)], ncol = d, nrow = d)
  }
  
  
  
  # FUNCTION: Updates grass density according to a local growth model
  grass<-function(P, X1 = M, X2 = pH, waterlog=F, r = GrowthRate)
  {
    # r: Intrinsic growth rate for grass
    a0<- -1
    a1<- 0.05
    a2<- 1.8
    a3<- -0.16
    d<-dim(P)[1]
    K<-exp(a0+a1*X1+a2*X2+a3*X2^2) # Carrying capacity as a function of moisture and pH
    
    PN<-matrix(rpois(length(P), P*exp(r*(1-P/K))), d, d)
    #sum(exp(r*(1-P/K)))
    if(waterlog==T) PN<-PN*matrix(exp(-18+0.5*X1)/(1+exp(-18+0.5*X1))<0.85, d, d)
    return(PN)
  }
  
  P<-matrix(P0, d, d) # Map of standing grass density
  P<-grass(P,M,pH, waterlog=F) # Grass growth for given time unit
  
  MemT = MemE = array(rep(NA, d*d), dim = c(d,d,N0))
  N<-matrix(0, d, d) # Map of animal usage
  
  if(Boundary=="no"){
    Buffer = PerceptionLength*20 # adds 10 to each side
    Add1 = function(Mat, Buffer){
      Arena2 = cbind(matrix(NA, ncol = Buffer, nrow = nrow(Mat)), Mat,matrix(NA, ncol = Buffer, nrow = nrow(Mat)))
      Arena2 = rbind(matrix(NA, ncol = ncol(Arena2), nrow = Buffer),Arena2,matrix(NA, ncol = ncol(Arena2), nrow = Buffer))
      Arena2
    }
    pH = Add1(pH, Buffer = Buffer)
    M = Add1(M, Buffer = Buffer)
    P = Add1(P, Buffer = Buffer)
    MemT = MemE = array(rep(NA, (Buffer*2+d)*(2*Buffer + d)), dim = c((d + 2*Buffer),(d + Buffer*2),N0))
    N = Add1(N, Buffer = Buffer)
  }
  P_hist = array(P0, dim = c(ncol(P),ncol(P), tmax))
  
  Co = ifelse(TakeShelter,.5,1)
  
  bx = min(c(0,min(which(!is.na(P), arr.ind = T)[,1])))
  bx = ifelse(bx==1,0,bx)
  by = min(c(0,min(which(!is.na(P), arr.ind = T)[,2])))
  by = ifelse(by==1,0,by)
  
  Center = rep(round(d/2),2) + CentreOffset
  
  if(!is.null(CenterRelease)){
    xs = Center[1] + c(CenterRelease,-CenterRelease)
    xs = round(runif(N0,min = min(xs),max = max(xs))) + bx
    ys = Center[2] + c(CenterRelease,-CenterRelease)
    ys = round(runif(N0,min = min(ys),max = max(ys))) + by
    
  }else{
    xs = round(runif(N0, 1,d)) + bx # The last term to account for potential buffer, if we do not have any boundary
    ys = round(runif(N0, 1,d)) + by # The last term to account for potential buffer, if we do not have any boundary
  }
  
  Ns=Ns2<-cbind("x" = xs,
                "y" = ys,
                "en"=runif(N0, BMR,100), 
                "be"=rbinom(N0, 1, Co), 
                "stepLength" = sample(seq(stepDistribution[1],
                                          stepDistribution[length(stepDistribution)],length.out  = 100), N0),
                "relAngle" = rvonmises(n = N0, mu = circular(pi/2), kappa = 0), # a random direction from previous location to start with
                "concentration" = sample(seq(0,concentration,length.out  = 100), N0))
  
  Coords <- expand.grid(x = 1:ncol(P), y = 1:ncol(P))
  
  
  # Specifications of telemetry study design 
  #telp: Proportion of population being tagged for telemetry observation
  #telq: Max proportion of total time over which animals will be observed
  #failp: Daily failure probability for each telemetry tag
  
  telP<-min(N0, ceiling(N0*telp)) # Number of tagged animals
  telQ<-min(tmax, ceiling(tmax*telq)) # Number of observation time units
  telQI<-tmax-telQ+1 # First observation to be recorded for tagged animals
  telD<-as.data.frame(matrix(NA, nrow=telQ*telP, ncol=2+ncol(Ns))) # Data structure for positions of tagged individuals
  #telD<-as.data.frame(matrix(NA, nrow=tmax*N0, ncol=2+ncol(Ns))) # Data structure for positions of tagged individuals
  names(telD)<-c("time","id", colnames(Ns))
  telL<-rep(0, telP)
  teln<-1 # Row number for entry of telemetry data
  
  ### The movement function
  move<-function(x,y,
                 X1,
                 MemoryTime,
                 MemoryEnv,
                 #MemoryDecay = 0.001,
                 TimeStep, 
                 PerceptionLength = 2,
                 stepDistribution = 2, 
                 stepCost = "linear",
                 concentration = 0, 
                 relAngle = NULL, 
                 Boundary = c("reflect", "torus", "no"),
                 ss,
                 d = NA, 
                 Coords = NULL
  ){
    #x/y = coordinates
    # X1 = env. layer
    # MemoryTime = a matrix containing values on WHEN the cell was last visited
    # MemoryEnv = a matrix containging values on WHAT the cell contained when visited
    # Coords: a matrix containing all possible coordinates. if not present it will be calculated by can also be given to avoid extra computaion time
    # TimeStep: an integer marking time. necessary for memory movement
    # ss = noise in perception
    # d = arena dimension
    # stepDistribution = mean of a log-normal
    # concentration = to model any correlated random walk using Von Mises distribution
    # ReflectiveBoundary: are animals on the arena boundary reflected or sent across to the other side?
    
    # k = 0.01;X1 = P; x = 100;y = 25; PerceptionLength = 2; stepDistribution = 2;concentration = 0;Boundary = "torus";ss = 5;MemoryTime=MemT[,,1]; MemoryEnv = MemT[,,1]
    # x = Ns[i,"x"];y = Ns[i,"y"];X1 = P;ss=ss;d=d;Coords = Coords;TimeStep = time;MemoryEnv = MemE[,,i];MemoryTime = MemT[,,i];stepDistribution = Ns[i,"stepLength"];relAngle = Ns[i, "relAngle"];concentration = Ns[i, "concentration"]
    d = ncol(X1)
    if(is.null(Coords)){Coords <- expand.grid(x = 1:d, y = 1:d)}
    if(Boundary=="reflect"){
      # Neighborhood definition and boundary checks.
      ys = (y-PerceptionLength):(y+PerceptionLength)
      ys<-ifelse(ys>d-1, d-1,ys)
      ys<-ifelse(ys<2, 2, ys)
      xs = (x-PerceptionLength):(x+PerceptionLength)
      xs<-ifelse(xs>d-1, d-1,xs)
      xs<-ifelse(xs<2, 2, xs)
      xs = unique(xs);ys = unique(ys)
    }else if(Boundary=="torus"){
      ys = (y-PerceptionLength):(y+PerceptionLength)
      ys = ifelse(ys>=d, 1 + (ys-d),ys)
      ys = ifelse(ys<1,d+ys, ys)
      
      xs = (x-PerceptionLength):(x+PerceptionLength)
      xs = ifelse(xs>=d, 1 + (xs-d),xs)
      xs = ifelse(xs<1,d+xs, xs)
      
    }else if(Boundary=="no"){# Allow the animal to move outside the edge. This area has no resources so the animal will through memory navigate its way bay
      ys = (y-PerceptionLength):(y+PerceptionLength)
      xs = (x-PerceptionLength):(x+PerceptionLength)
    }
    
    xs = unique(xs);ys = unique(ys)
    
    # Coordinates of where the animal can see
    # matrix rows are x's and columns are y's
    idx = as.matrix(expand.grid(row = xs,col = ys))
    
    #MemT = MemE = array(rep(NA, d*d), dim = c(d,d,N0));MemoryEnv = MemE[,,i];MemoryTime = MemT[,,i]
    MemEnvNu = MemoryEnv
    MemEnvNu[idx]<-X1[idx] # Update the cells you perceive with the perceived values
    
    MemTimeNu = MemoryTime
    MemTimeNu[idx]<-TimeStep
    
    # Memory weight function
    # Based on that from Merkle et al. 2014, which is similar to that of McNamara and Houston 1985
    Weights = function(k,visited = MemoryTime, time = TimeStep){
      mat = 1/(1 + k*(time-visited))
      mat[visited==0 | is.na(visited)]<-0
      mat # Unvisited patches set to 0
    }
    
    # Use the summer values foundin merkle 2014
    kPrevVis = 0.00000001 #low values suggest that individuals are equally likely to choose recently visited patches as those visisted long time ago
    PrevVis = Weights(k = kPrevVis, visited = MemTimeNu, time = TimeStep)
    
    kRelRefPoint = 0.008 # weight of the most recent foraged patches
    rrWeights = Weights(k = kRelRefPoint, visited = MemoryTime, time = TimeStep)
    RelRefPoint = (sum(rrWeights*MemoryEnv, na.rm = T)/sum(rrWeights, na.rm = T))- # The mean of previously visited patches compared to the source (the one currently inhabited)
      MemEnvNu[x,y] # Positive values indicate that previously visited patches are better than the current patch
    RelRefPoint = ifelse(TimeStep==1,mean(MemEnvNu[idx], na.rm = T),RelRefPoint) # Have no experience to start with
    
    kRefPoint = 0.000001
    rrWeights = Weights(k = kRefPoint, visited = MemTimeNu)
    RefPoint = sum(rrWeights*MemEnvNu, na.rm = T)/sum(rrWeights, na.rm = T) # The weighted mean of visited patches
    RefPoint = ifelse(TimeStep==1,mean(MemEnvNu[idx], na.rm = T),RefPoint) # Have no experience to start with
    
    kProfitPrevVis = 0.004 # memory decay of the PROFITABILITY of previously visisted patches, i.e. how important it is to remember the quality of the patch
    MemEnvNu[is.na(MemEnvNu)]<-0
    # Fill the non-perceptect patches with expected values
    ExpProfit = RefPoint + (MemEnvNu - RefPoint)*Weights(k = kProfitPrevVis, visited = MemTimeNu) #unvisisted patches are assigned with the ref point value
    #ExpectedProfit = MemoryEnv*(1/(1 + k*(TimeStep-MemoryTime)))
    #ExpectedProfit[idx]<-X1[idx] # Update the cells you perceive with the perceived values
    
    
    if(Boundary=="torus"){
      Distances = sqrt((round(median(range(Coords[,1])))-Coords[,1])^2 +
                         (round(median(range(Coords[,2]))) - Coords[,2])^2)
    }else{Distances = sqrt((x-Coords[,1])^2 + (y - Coords[,2])^2)
    }
    MinDistance = 0.49
    Distances[Distances==0]<-MinDistance
    Distances = matrix(Distances, 
                       ncol = ncol(MemoryEnv), nrow = nrow(MemoryEnv))
    
    #library(poweRlaw); MovementWeights = matrix(dplcon(Distances, MinDistance, stepDistribution), ncol = d, nrow = d) # Based on a power-law distribution
    # Cost of movement
    if(stepCost=="normal"){
      MovementWeights = matrix(dnorm(Distances, sd = stepDistribution), 
                               ncol = ncol(MemoryEnv), nrow = nrow(MemoryEnv))
    }else if(stepCost=="linear"){
      MovementWeights = matrix(1/Distances, 
                               ncol = ncol(MemoryEnv), nrow = nrow(MemoryEnv))
    }
    
    
    # Correalted walk is not UPDATED
    if(concentration>0){ #if there is a correlated random walk 
      # relAngle = Ns[i, "relAngle"]
      if(is.null(relAngle)){print("Needs input on angle from previous location")}else{
        
        
        # weighting according to direction
        mat = sapply(1:length(ys), function(Y){
          sapply(1:length(xs), function(X){
            atan2(ys[Y]-y,xs[X]-x)
          })})  
        wsR = apply(mat, c(1,2), function(x) circular::dvonmises(x, mu=circular(relAngle), kappa=concentration))
        wsR = wsR/sum(wsR)
        MovementWeights = MovementWeights*wsR
        
      }}
    
    # The sum of them all
    # by multiplying we assume that when a cell is void of food, the animal will not be "forced" to move there due to step length preferences
    # For now we use the beta coefficients for winter found in Merkle 2014 
    
    Profit =  (1.948*ExpProfit +# More likely to visit above-average profitable patches
                 1.284*PrevVis + # Likely to go back to known location (e.g. risk-sensitive?)
                 0.183*RelRefPoint*PrevVis) # More likely to visit previous patches if past experience is bettern than current location
    
    if(Boundary=="torus"){
      MovementKernel = (TorusOider(Profit, coords = c(x,y)) + -0.153*Distances) # The last term is already "torusoided" if needed
    }else{
      MovementKernel = (Profit + -0.153*Distances)}
    MovementKernel = MovementKernel/sum(MovementKernel, na.rm = T)
    
    # Insert 0's for NA-values
    MovementKernel[is.na(MovementKernel)]<-0
    
    pos = which(MovementKernel == max(MovementKernel, na.rm = T), arr.ind = TRUE) # this will collect position of the desired cell
    
    if(Boundary == "torus"){
      cc = round(apply(apply(Coords, 2, range),2, median))
      pos = c(x,y) + (pos-cc)
      pos = ifelse(pos<1,d,ifelse(pos>d,1,pos)) # Move to other edge if neccessary
    }
    
    #pos = c(as.numeric(rownames(MovementKernel)[jump[1]]),
    #as.numeric(colnames(MovementKernel)[jump[2]]))
    
    attr(pos, "StepLength")<-sqrt((x-pos[1])^2 + (y-pos[2])^2)
    attr(pos, "relAngle")<-atan2(pos[2]-y,pos[1]-x)
    MemoryEnv[idx] = rowMeans(cbind(ifelse(MemoryEnv[idx]==0,
                                           NA,
                                           MemoryEnv[idx]),X1[idx]), na.rm = T) # Update experience so that previously visisted cells will get the mean of old and newly observed. Since unvisited cells are NA, this will mean that unvisisted cells will get the value perceived. 
    MemoryTime[idx] = TimeStep # Update the memory on when the cell was last visited
    attr(pos, "MemoryEnv") <- MemoryEnv
    attr(pos, "MemoryTime") <- MemoryTime
    
    return(pos)
  }
  
  
  
  print("Beginning loop")
  for (time in 1:(tmax)) # Time loop
  {
    P<-grass(P,M,pH, waterlog=F) # Grass growth for given time unit
    P_hist[,,time]<-P
    if(time/tmax*10==round(time/tmax*10)) print(paste("Completed :", time/tmax*100, "%"))
    for(i in 1:N0) # Animal loop
    {
      #print(paste("Timestep",time," - Animal",i))
      if(Ns[i,"be"]==1) # If the animal is foraging
      {
        pos<-move(x = Ns[i,"x"],
                  y = Ns[i,"y"],
                  X1 = P,# Move animal according to grass density
                  ss=ss, # noise in perception
                  d=d,
                  Coords = Coords,
                  TimeStep = time,
                  MemoryEnv = MemE[,,i],
                  MemoryTime = MemT[,,i],
                  stepDistribution = Ns[i,"stepLength"],
                  relAngle = Ns[i, "relAngle"],
                  concentration = Ns[i, "concentration"])
        
        Ns[i,"x"]<-pos[1]
        Ns[i, "y"]<-pos[2]
        pos=ceiling(pos) # Position to use for matrix extraction
        Ns[i, "relAngle"]<-attr(pos, "relAngle") # update the relative angle from previous position
        
        MemE[,,i] <-attr(pos, "MemoryEnv") # Update memory components
        MemT[,,i] <-attr(pos, "MemoryTime")
        
        
        P[pos[1],pos[2]]<-P[pos[1],pos[2]]-ForageRate#round(ForageRate*P[pos[1],pos[2]]) # Grass depletion caused by single animal
        Ns[i,"en"]<-Ns[i,"en"]+ForageRate#0.01*P[pos[1],pos[2]] # Gain in energetic state through grazing
        #Ns[i,"en"]<-Ns[i,"en"]-attr(pos, "StepLength")*.75 # Gain in energetic state through grazing - cost of moving
        
        if(TakeShelter & Ns[i,"en"]>80) Ns[i,"be"]<-0 # Sets the animal to sheltering mode if its energy is high
      }
      
      if(Ns[i,"be"]==0) # If the animal is sheltering
      {
        pos<-move(x = Ns[i,"x"],
                  y = Ns[i,"y"],
                  X1 = pH, # Move animal according to "conifer" density
                  ss=ss, # noise in perception
                  d=d,
                  Coords = Coords,
                  TimeStep = time,
                  MemoryEnv = MemE[,,i],
                  MemoryTime = MemT[,,i],
                  stepDistribution = Ns[i,"stepLength"],
                  relAngle = Ns[i, "relAngle"],
                  concentration = Ns[i, "concentration"])
        Ns[i,"x"]<-pos[1]
        Ns[i, "y"]<-pos[2]
        pos=ceiling(pos) # Position to use for matrix extraction
        Ns[i, "relAngle"]<-attr(pos, "relAngle") # update the relative angle from previous position
        MemE[,,i] <-attr(pos, "MemoryEnv") # Update memory components
        MemT[,,i] <-attr(pos, "MemoryTime")
        
        if(Ns[i,"en"]<20) Ns[i,"be"]<-1 # Sets the animal to foraging mode if its energy is low
      }
      
      Ns[i,"en"]<-Ns[i,"en"]-BMR # Basal metabolic cost
      
      if(time>0.2*tmax) N[pos[1],pos[2]]<-N[pos[1],pos[2]]+1 # If initial settling is done, increment usage by one unit for this animal
      
      #Recording of telemetry data
      if(time>=telQI && i<=telP && telL[i]!=1)
      {
        telD[teln,]<-c(time, i,
                       Ns[i,]#+runif(2,-0.5,0.5) # Observation error in telemetry points
        )
        teln<-teln+1 # Keeps track of number of telemetry points
        telL[i]<-rbinom(1,1,failp) # Malfunction in tags as a Bernoulli process
      }
    }
  }
  
  telD$N = N0
  #with(telD, plot(y~x))
  ls = list(Data = telD, N = as.matrix(N), P = P_hist,
            M = as.matrix(M), pH  =as.matrix(pH))
  if(!is.null(VarGrowthRate)){
    ls = c(ls, GrowthRate = as.matrix(GrowthRate))
  }
  ls
}

# Find K, and vary the N0 around this
aa = MoveSim(tmax = 1000, 
             CenterRelease = 2,
             CentreOffset = c(5,5),
             N0 = 1, d = 15, Boundary = "torus",
              stepCost = "linear",
              ForageRate = 80, BMR = 10)

#aa = aa2
library(data.table);library(cowplot)
TotData = data.table(aa$Data)
TotData
#TotData[,Bound := diff(x)>40 | diff(y)>40,"id"]

qplot(data = TotData, x = x, y = y
      , col = factor(id), size = sqrt(time))


# Graphical output
par(mfrow=c(1,3), bty = "L")
plot(log(aa$N)~aa$P[,,100])
#with(datas, plot(en~time))
image(aa$P[,,100], main="Grass density", col=terrain.colors(100))
image(aa$N, main="Animal usage", col=terrain.colors(100))
par(mfrow=c(1,1))


AllData = CJ(x = unique(TotData$x), y = unique(TotData$y), time = unique(TotData$time))
AllData = TotData[,.(nIndXY=uniqueN(id)), c("x","y", "time")][AllData, on = c("x", "y", "time")]
AllData[,nIndXY := ifelse(is.na(nIndXY),0,nIndXY)]
AllData[,VarXY:=var(nIndXY)]
AllData[,SumXY:=sum(nIndXY), c("x", "y")]
summary(AllData)

# Variation in energy per patch
MeanP = apply(aa$P, c(1,2), mean, na.rm = T)
colnames(MeanP)<-1:ncol(MeanP)
MeanP = melt(data.table(x = 1:nrow(MeanP), MeanP), id.var = "x")
names(MeanP) = c("x", "y", "MeanP")
MeanP[,x:=paste(x)][,y :=paste(y)]

VarP = apply(aa$P, c(1,2), var, na.rm = T)
colnames(VarP)<-1:ncol(VarP)
VarP = melt(data.table(x = 1:nrow(VarP), VarP), id.var = "x")
names(VarP) = c("x", "y", "VarP")
VarP[,x:=paste(x)][,y :=paste(y)]

AllData[,x:=paste(x)][,y :=paste(y)]
AllData = VarP[AllData, on = c("x","y")]
AllData = MeanP[AllData, on = c("x","y")]

qplot(data = AllData[,.SD[1], c("x","y")], x = MeanP, y = VarP)
qplot(data = AllData[,.SD[1], c("x","y")], x = VarP, y = VarXY)

TotData = CJ(x = unique(TotData$x), y = unique(TotData$y))[,cellid:=1:.N][TotData, on = c("x", "y")]
TotData[,ID:=id]

BurnIn = 1
TotData[,TAC:=sapply(BurnIn:max(time), function(aa){uniqueN(.SD[time<=aa & time >=BurnIn]$cellid)}),c("ID")]

TotData

library(cowplot)
ggplot(data = TotData[time>=BurnIn & TAC>1], aes(x = time, y = TAC)) + stat_smooth()


## Influence of density on variance in energy per patch####
nRuns = 5
TIME = 500
aa = vector("list", nRuns)
aa2 = vector("list", nRuns)
aa3 = vector("list", nRuns)
aa4 = vector("list", nRuns)

a1 = lapply(which(sapply(aa, "class")=="NULL"), function(x){
  print(x)
  if(x==1){aa[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5), d = 50, Boundary = "torus", 
                               ForageRate = .5)
  }else{
    aa[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5), d = 50, Boundary = "torus", 
                        ForageRate = .5, 
                        AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M)$Data
  }
})

# Make this loops so all simualtions are done on the same landscape.
a2 = lapply(which(sapply(aa2, "class")=="NULL"), function(x){ print(x)
  aa2[[x]] <<- MoveSim(tmax = TIME, N0 = 5, tburning = 0, SL = c(1.1,3.5),
                       AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M,
                       d = 50, Boundary = "torus", ForageRate = .5)$Data
})

a3 = lapply(which(sapply(aa3, "class")=="NULL"), function(x){ print(x)
  aa3[[x]] <<- MoveSim(tmax = TIME, N0 = 10, tburning = 0, SL = c(1.1,3.5),
                       AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M, 
                       d = 50, Boundary = "torus", ForageRate = .5)$Data
})

a4 = lapply(which(sapply(aa4, "class")=="NULL"), function(x){ print(x)
  aa4[[x]] <<- MoveSim(tmax = TIME, N0 = 20, tburning = 0, SL = c(1.1,3.5),
                       AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M, 
                       d = 50, Boundary = "torus", ForageRate = .5)$Data
})



library(data.table)
TotData1 = rbind(rbindlist(lapply(which(sapply(aa, "class")=="data.frame"), function(ii) {
  df = data.table(aa[[ii]])
  df$Run=ii
  df$N0 = 1
  df
})),
rbindlist(lapply(which(sapply(aa2, "class")=="data.frame"), function(ii) {
  df = data.table(aa2[[ii]])
  df$Run=ii
  df$N0 = 5
  df
})),
rbindlist(lapply(which(sapply(aa3, "class")=="data.frame"), function(ii) {
  df = data.table(aa3[[ii]])
  df$Run=ii
  df$N0 = 10
  df
})),
rbindlist(lapply(which(sapply(aa4, "class")=="data.frame"), function(ii) {
  df = data.table(aa4[[ii]])
  df$Run=ii
  df$N0 = 20
  df
}))
)

TotData1 = CJ(x = unique(TotData1$x), y = unique(TotData1$y))[,cellid:=1:.N][TotData1, on = c("x", "y")]
TotData1[,ID:=paste(Run, N,N0, id, sep = "_")]

BurnIn = 1000
TotData1[,TAC:=sapply(BurnIn:max(time), function(aa){
  uniqueN(.SD[time<=aa & time >=BurnIn]$cellid)
}),c("ID")]
library(cowplot)
ggplot(data = TotData1[time>=BurnIn], aes(x = time, y = TAC, col = (N0), group = factor(N0))) + stat_smooth()



# TotData[,TimeIntervals:=cut(time, seq(0,max(time),500), include.lowest = T)
#         ][,c("cx", "cy"):=list(mean(x), mean(y)), c("ID","TimeIntervals")
#           ][,cShift:=sqrt(diff(unique(cx, shift(cx)))^2 + diff(unique(cy, shift(cy)))^2),"ID"
#             ]
# TotData[id==1 & time %in% c(500,501)]
# ggplot(data = TotData[,.SD[.N], c("ID", "TimeIntervals")],
#        aes(x = time, y = cShift, col = factor(Decay))) + stat_smooth()




### DECAY RATE AT DENS == 1 ####
nRuns = 20
TIME = 5000
aa = vector("list", nRuns)
aa2 = vector("list", nRuns)
aa3 = vector("list", nRuns)
aa4 = vector("list", nRuns)

a1 = lapply(which(sapply(aa, "class")=="NULL"), function(x){
  print(x)
  if(x==1){aa[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5), d = 50, Boundary = "no", 
                               ForageRate = .5, MemoryDecay = 0.0001)
  }else{
    aa[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5), d = 50, Boundary = "no", 
                        ForageRate = .5, MemoryDecay = 0.0001,
                        AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M)$Data
  }
})

a2 = lapply(which(sapply(aa2, "class")=="NULL"), function(x){ print(x)
  aa2[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5),
                       AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M,
                       d = 50, Boundary = "no", ForageRate = .5, MemoryDecay = 0.00001)$Data
})

a3 = lapply(which(sapply(aa3, "class")=="NULL"), function(x){ print(x)
  aa3[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5),
                       AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M, 
                       d = 50, Boundary = "no", ForageRate = .5, MemoryDecay = 0.000001)$Data
})

a4 = lapply(which(sapply(aa4, "class")=="NULL"), function(x){ print(x)
  aa4[[x]] <<- MoveSim(tmax = TIME, N0 = 1, tburning = 0, SL = c(1.1,3.5),
                       AcidLayer = aa[[1]]$pH, MoistureLayer = aa[[1]]$M, 
                       d = 50, Boundary = "no", ForageRate = .5, MemoryDecay = 0.0000001)$Data
})

TotData = rbind(rbindlist(lapply(which(sapply(aa, "class")=="data.frame"), function(ii) {
  df = data.table(aa[[ii]])
  df$Run=ii
  df$Decay = 0.0001
  df
})),
rbindlist(lapply(which(sapply(aa2, "class")=="data.frame"), function(ii) {
  df = data.table(aa2[[ii]])
  df$Run=ii
  df$Decay = 0.00001
  df
})),
rbindlist(lapply(which(sapply(aa3, "class")=="data.frame"), function(ii) {
  df = data.table(aa3[[ii]])
  df$Run=ii
  df$Decay = 0.000001
  df
})),
rbindlist(lapply(which(sapply(aa4, "class")=="data.frame"), function(ii) {
  df = data.table(aa4[[ii]])
  df$Run=ii
  df$Decay = 0.0000001
  df
}))
)

TotData = CJ(x = unique(TotData$x), y = unique(TotData$y))[,cellid:=1:.N][TotData, on = c("x", "y")]
TotData[,ID:=paste(Run, N, Decay, id, sep = "_")]

BurnIn = 1000
TotData[,TAC:=sapply(BurnIn:max(time), function(aa){
  uniqueN(.SD[time<=aa & time >=BurnIn]$cellid)
}),c("ID")]
ggplot(data = TotData[time>=BurnIn], aes(x = time, y = TAC, col = (Decay), group = factor(Decay))) + stat_smooth()
TotData[,TimeIntervals:=cut(time, seq(0,max(time),500), include.lowest = T)
        ][,c("cx", "cy"):=list(mean(x), mean(y)), c("ID","TimeIntervals")
          ][,cShift:=sqrt(diff(unique(cx, shift(cx)))^2 + diff(unique(cy, shift(cy)))^2),"ID"
            ]
TotData[id==1 & time %in% c(500,501)]
ggplot(data = TotData[,.SD[.N], c("ID", "TimeIntervals")],
       aes(x = time, y = cShift, col = factor(Decay))) + stat_smooth()

#test1 = MoveSim(tmax = 2000, N0 = 10, tburning = 500, SL = c(1,3), d = 100)
#test2 = MoveSim(tmax = 2000, N0 = 25, tburning = 500, SL = c(1,3), d = 100)
#test3 = MoveSim(tmax = 2000, N0 = 50, tburning = 500, SL = c(1,3), d = 100)
#test4 = MoveSim(tmax = 2000, N0 = 100, tburning = 500, SL = c(1,3), d = 100)
#test5 = MoveSim(tmax = 2000, N0 = 150, tburning = 500, SL = c(1,3), d = 100)
MoveSim(tmax = 100, N0 = 1, tburning = 0, SL = c(1.1,3.5), d = 50, Boundary = "no", ForageRate = .5, MemoryDecay = 0.0001)

n1 = lapply(1:10,function(x){
  test2 = MoveSim(tmax = 2000, N0 = 1, tburning = 0, SL = c(1.1,3.5), d = 50, Boundary = "no", ForageRate = .5, MemoryDecay = 0.0001)
})
n5 = lapply(1:10,function(x){
  test2 = MoveSim(tmax = 2000, N0 = 5, tburning = 0, SL = c(1.1,3.5), d = 50, ForageRate = .5)
})
n10 = lapply(1:10,function(x){
  test2 = MoveSim(tmax = 2000, N0 = 10, tburning = 0, SL = c(1.1,3.5), d = 50, ForageRate = .5)
})

n20 = lapply(1:10,function(x){
  test2 = MoveSim(tmax = 2000, N0 = 20, tburning = 0, SL = c(1.1,3.5), d = 50, ForageRate = .5)
})

TotData = rbind(rbindlist(lapply(seq_along(n1), function(ii) data.table(n1[[ii]]$Data)[,Run:=paste0("n1_",ii)][,N:=1])),
                rbindlist(lapply(seq_along(n5), function(ii) data.table(n5[[ii]]$Data)[,Run:=paste0("n5_",ii)][, N:=5])),
                rbindlist(lapply(seq_along(n10), function(ii) data.table(n10[[ii]]$Data)[,Run:=paste0("n10_",ii)][,N:=10])),
                rbindlist(lapply(seq_along(n20), function(ii) data.table(n20[[ii]]$Data)[,Run:=paste0("n20_",ii)][,N:=20])))

#MemoryDecay
library(ggplot2);library(data.table)
#datas = data.table(test2$Data)
TotData = CJ(x = unique(TotData$x), y = unique(TotData$y))[,cellid:=1:.N][TotData, on = c("x", "y")]
TotData[,uniqueN(time), c("Run","N","id")]

TotData[,TAC:=sapply(1:max(time), function(aa){
  uniqueN(.SD[time<=aa]$cellid)
}),c("id", "N", "Run")]
setorder(TotData, Run, id, time)
ggplot(data = TotData[N==1], aes(x = time, y = TAC, col = as.factor(id), linetype = factor(N))) + geom_line()

P = test2$P[,,1]
P_var = apply(test2$P, c(1,2), function(x) var(x))
P = apply(test2$P, c(1,2), function(x) mean(x))
N = test2$N
pH = test2$pH


datas[,EnergyDiff:=c(NA,diff(en)), "id"]
datas[,NextDist:=c(sqrt(diff(x)^2 + diff(y)^2), NA), "id"]
datas[,Displacement:=sqrt((x-x[1])^2 + (y-y[1])^2), c("id")]
with(datas[id==2], plot(Displacement~time))

datas[,c("MSD11","MSD21", "R"):={
  df1 = .SD[time %between% c(max(time)*.8,max(time)*.9)]
  df2 = .SD[time %between% c(max(time)*.9,max(time)*1)]
  mx = mean(df1$x); my = mean(df1$y)
  mx2 = mean(df2$x); my2 = mean(df2$y)
  MSD11 = mean((sqrt((df1$x-mx)^2 + (df1$y-my)^2))^2)
  MSD21 = mean((sqrt((df2$x-mx)^2 + (df2$y-my)^2))^2)
  MSD22 = mean((sqrt((df2$x-mx2)^2 + (df2$y-my2)^2))^2)
  list(MSD11,
       MSD21,
       R = abs(MSD21- MSD11)/(MSD11 + MSD22)
  )
},"id"]

with(datas[,.SD[1], "id"], plot(R~stepLength))
with(datas[,.SD[1], "id"], summary(R))


with(datas[,(mean(NextDist, na.rm = T)), c("id", "stepLength")], cor.test(V1,stepLength))
# Graphical output
x<-1:ncol(P)
y<-1:nrow(P)
par(mfrow=c(2,3), bty = "L")
plot(log(N)~P)
with(datas[,(mean(NextDist, na.rm = T)), c("id", "stepLength")], plot(V1~stepLength))
#with(datas, plot(en~time))
image(x, y, P, main="Grass density", col=terrain.colors(100))
image(x, y, pH, main="Cover density", col=terrain.colors(100))
image(x, y, N, main="Animal usage", col=terrain.colors(100))
par(mfrow=c(1,1))

with(datas[id==1 & time >950],
     points(y = y, x = x, col = ifelse(time==1, "blue", ifelse(time == max(datas$time), "red", "black")), lty = datas$id, type = "l"))

datas[,id:=paste(id)]
library(adehabitatHR)
datas = data.frame(datas)
HRs = lapply(unique(datas$id), function(x){
  #x = unique(datas$id)[10]
  df = datas[datas$time>1000 & datas$id==x,]
  coordinates(df)<- ~x + y
  ud = kernelUD(df)
  ud
})
HRa = lapply(HRs, function(x) {
  udHR = kernel.area(x, percent = 95)
  udHR
})
library(raster)
P1 = raster(P, xmn = 1, xmx = max(datas$x), ymn = 1, ymx = max(datas$y), crs = CRS("+proj=utm +zone=32 +datum=WGS84"))
pH1 = raster(pH, xmn = 1, xmx = max(datas$x), ymn = 1, ymx = max(datas$y), crs = CRS("+proj=utm +zone=32 +datum=WGS84"))
P2 = as(P1, "SpatialPixelsDataFrame")

resources = (t(sapply(HRs, function(x){
  #print(x)
  xx = getverticeshr(x, percent = 50)
  proj4string(xx) <- proj4string(P2)
  c(HR=xx$area,
    food = sum(extract(P1, xx)[[1]]),
    cover = sum(extract(pH1, xx)[[1]]))
})))

datas2 = data.table(ID=unique(datas$id), resources)
par(mfrow = c(2,1))
with(datas2, plot(HR~food))
with(datas2, plot(y = HR, x = cover, col = "red"))


