###### EXAMPLE #####
library(ggplot2)
th = theme(
  legend.position = "bottom",
  legend.text = element_text(color="black", size=18),
  plot.title = element_text(color="black", size=20, face="bold"),
  axis.title = element_text(color="black", size=18,face = "bold"),
  axis.text = element_text(color="black", size=18),
  strip.text = element_text(size = 18, colour = "black")
)

####
# Genereal form of the gain function
# [1] Charnov EL, Parker GA. 1995 Dimensionless Invariants from Foraging Theory Marginal Value Theorem. Proc. Natl. Acad. Sci. USA 92, 1446-1450

Gm = .5 # Proportionate to patch size
g = .55 # Constant reflecting initial capture rate of prey
c = g/Gm # rate of prey capture
GeneralGainForm = function(Gmax, c, t){
  Gmax*(1-exp(-c*t))
}
GeneralGainForm(Gmax = .5, c = 1.1, t =1)

GeneralGainFormDeriv = function(mu, t, Gmax, c){
  GeneralGainForm(Gmax = Gmax, c = c, t = t)/(mu+t)
}
GeneralGainFormDeriv(mu = 2, Gmax = .5, c = 1.1, t =1)

curve(Gm*(1-exp(-5*x)), from = 0, to = 10)
curve(Gm*(1-exp(-2*x)), from = 0, to = 10,add = T)
curve(Gm*(1-exp(-1*x)), from = 0, to = 10,add = T)
curve(Gm*(1-exp(-.5*x)), from = 0, to = 10,add = T)
curve(Gm*(1-exp(-.4*x)), from = 0, to = 10,add = T)
curve(Gm*(1-exp(-.3*x)), from = 0, to = 10,add = T)
curve(Gm*(1-exp(-.1*x)), from = 0, to = 10,add = T)
curve(Gm*(1-exp(-.0000000000001*x)), from = 0, to = 100, add = T)
curve(GeneralGainFormDeriv(mu = 2, Gmax = .5, c = 1.1, t =x), from=0,to = 10,add=T, col="green")

# For a given travel/search time mu
mu = 2
Toptim = (1/c)*(1+(1/(c*mu + 1)))*log(c*mu+1)
GeneralGainForm(Gmax = .5, c = 1.1, t =Toptim)

####################
asymptote = 1
LogRate = 0.5
yInt = 0
MVT_ = function(aymptote = 1, Lograte = 0.5, yInt = 0, HandlingTime = seq(0,10,.25)){
  # aymptote = 1; Lograte = 0.5; yInt = 0; HandlingTime = seq(0,10,.25)
  df = data.frame(Asymptote = asymptote,
                  Rate = LogRate,
                  Intercept = yInt, HandlingTime = HandlingTime)
  
  df$Gain = asymptote +(yInt - asymptote)*exp(-exp(LogRate)*HandlingTime)
  
  # To find the tangents
  x = seq(0,10,.1)
  # Computes the derivatives
  fGras = deriv(~asymptote +(yInt - asymptote)*exp(-exp(LogRate)*x),"x")
  fGras1 = deriv(~2*asymptote +(yInt - 2*asymptote)*exp(-exp(LogRate)*x),"x")
  fGras2 = deriv(~asymptote +(yInt - asymptote)*exp(-exp(2*LogRate)*x),"x")
  fGras3 = deriv(~2*asymptote +(yInt - 2*asymptote)*exp(-exp(2*LogRate)*x),"x")
  bsG = eval(fGras)
  bsG1 = eval(fGras1)
  bsG2 = eval(fGras2)
  bsG3 = eval(fGras3)
  
  # The y - values along the x-values
  yG = bsG[seq_along(x)]
  yG1 = bsG1[seq_along(x)]
  yG2 = bsG2[seq_along(x)]
  yG3 = bsG3[seq_along(x)]
  
  # The derivats
  TangentG = attr(bsG,"grad")
  TangentG1 = attr(bsG1,"grad")
  TangentG2 = attr(bsG2,"grad")
  TangentG3 = attr(bsG3,"grad")
  
  xInterceptG = yG -TangentG*x
  xInterceptG1 = yG1 -TangentG1*x
  xInterceptG2 = yG2 -TangentG2*x
  xInterceptG3 = yG3 -TangentG3*x
  
  TravelTimeG= abs((0-xInterceptG)/TangentG)
  TravelTimeG1= abs((0-xInterceptG1)/TangentG1)
  TravelTimeG2= abs((0-xInterceptG2)/TangentG2)
  TravelTimeG3= abs((0-xInterceptG3)/TangentG3)
  
  par(mfrow =c(2,1))
  plot(yG~x, xlim = c(0,5),ylim = c(0,2), main = "Gain curve", type="l", xlab ="Time in patch", ylab = "Gain")
  lines(yG1~x, xlim = c(0,5), col = "red")
  lines(yG2~x, xlim = c(0,5), col = "green")
  lines(yG3~x, xlim = c(0,5), col = "blue")
  
  
  plot(log(x)~log(TravelTimeG), xlim = c(0,5), main = "", type="l",
       xlab = "Travel distance (ln-dist)", ylab = "Opt. patch residence time (ln-time)")
  lines(log(x)~log(TravelTimeG1), xlim = c(0,5), col = "red",lty="dotted")
  
  lines(log(x)~log(TravelTimeG2), xlim = c(0,5), col = "green", main = "", type="l")
  lines(log(x)~log(TravelTimeG3), xlim = c(0,5), col = "blue",lty="dotted")
  # Doubling of rate shifts the elevation, whereas, shifting the asymptote does nothing
  
  }
#curve(asymptote +(yInt - asymptote)*exp(-exp(LogRate)*x), from = 1, to = 10, "red")
plot(MVT_()~seq(0,10,.25), type = "l")



# To find the tangents
x = seq(0,10,.1)
# Computes the derivatives
fGras = deriv(~asymptote +(yInt - asymptote)*exp(-exp(LogRate)*x),"x")
fGras1 = deriv(~2*asymptote +(yInt - 2*asymptote)*exp(-exp(LogRate)*x),"x")
fGras2 = deriv(~asymptote +(yInt - asymptote)*exp(-exp(2*LogRate)*x),"x")
fGras3 = deriv(~2*asymptote +(yInt - 2*asymptote)*exp(-exp(2*LogRate)*x),"x")
bsG = eval(fGras)
bsG1 = eval(fGras1)
bsG2 = eval(fGras2)
bsG3 = eval(fGras3)

# The y - values along the x-values
yG = bsG[seq_along(x)]
yG1 = bsG1[seq_along(x)]
yG2 = bsG2[seq_along(x)]
yG3 = bsG3[seq_along(x)]

# The derivats
TangentG = attr(bsG,"grad")
TangentG1 = attr(bsG1,"grad")
TangentG2 = attr(bsG2,"grad")
TangentG3 = attr(bsG3,"grad")

xInterceptG = yG -TangentG*x
xInterceptG1 = yG1 -TangentG1*x
xInterceptG2 = yG2 -TangentG2*x
xInterceptG3 = yG3 -TangentG3*x

TravelTimeG= abs((0-xInterceptG)/TangentG)
TravelTimeG1= abs((0-xInterceptG1)/TangentG1)
TravelTimeG2= abs((0-xInterceptG2)/TangentG2)
TravelTimeG3= abs((0-xInterceptG3)/TangentG3)

par(mfrow =c(2,1))
plot(yG~x, xlim = c(0,5),ylim = c(0,2), main = "Gain curve", type="l", xlab ="Time in patch", ylab = "Gain")
lines(yG1~x, xlim = c(0,5), col = "red")
lines(yG2~x, xlim = c(0,5), col = "green")
lines(yG3~x, xlim = c(0,5), col = "blue")


plot(log(x)~log(TravelTimeG), xlim = c(0,5), main = "", type="l",
     xlab = "Travel distance (ln-dist)", ylab = "Opt. patch residence time (ln-time)")
lines(log(x)~log(TravelTimeG1), xlim = c(0,5), col = "red",lty="dotted")

lines(log(x)~log(TravelTimeG2), xlim = c(0,5), col = "green", main = "", type="l")
lines(log(x)~log(TravelTimeG3), xlim = c(0,5), col = "blue",lty="dotted")
# Doubling of rate shifts the elevation, whereas, shifting the asymptote does nothing


Gm = 1
CaptureRateGras = 0.20
fGras = deriv(~Gm*(1-exp(-CaptureRateGras*x)),"x")
bsG = eval(fGras)
yG = bsG[seq_along(x)]
TangentG = attr(bsG,"grad")
xInterceptG = yG -TangentG*x
TravelTimeG= abs((0-xInterceptG)/TangentG)

Herbm = 3*Gm
CaptureRateHerb = 0.35
fHerb = deriv(~Herbm*(1-exp(-CaptureRateHerb*x)),"x")
Herbm = 1
bsHerb = eval(fHerb)
yHerb = bsHerb[seq_along(x)]
TangentHerb = attr(bsHerb,"grad")
xInterceptHerb = yHerb -TangentHerb*x
TravelTimeH= abs((0-xInterceptHerb)/TangentHerb)
plot(bsHerb~TravelTimeH, type = "l")

Zero = 3*Gm
CaptureRateZero = 1
fZero = deriv(~Zero*(x),"x")
Zerom = 1
bsZero = eval(fZero)
yZero = bsZero[seq_along(x)]
TangentZero = attr(bsZero,"grad")
xInterceptZero = yZero -TangentZero*x
TravelTimeZ= abs((0-xInterceptZero)/TangentZero)
plot(bsZero~TravelTimeZ, type="l")

par(mfrow=c(1,2), mai=c(1.02,1,0.82,0.42))
curve(Herbm*(1-exp(-CaptureRateHerb*x)),from =0, to =100,
      #plot(eval(fGras),  type="l",ylim=c(0,1),
      ylab=expression(Gain^- time), xlab="Time",cex.lab = 1.5,cex.axis = 1.5)
curve(Gm*(1-exp(-CaptureRateGras*x)),from =0, to =100,
      add= T, col = "green")
curve(Zero*(1-exp(-CaptureRateZero*x)),from =0, to =100,
      add= T, col = "red")

plot(log(TravelTimeG),log(x), type="l", ylab="Residence time (ln)",
     xlab="Travel time (ln)",cex.lab = 1.5,cex.axis = 1.5)
lines(log(TravelTimeH),log(x), col ="green", type="l")
lines(log(TravelTimeH),log(x), col ="green", type="l")



# plot(eval(fGras),  type="l")
# for(i in 1:length(TravelTimeG)){
#      abline(b = TangentG[i], a = xInterceptG[i])
# }

# for(i in 1:length(TravelTimeH)){
#      abline(b = TangentHerb[i], a = xInterceptHerb[i],col="green")
# }

getwd()
png("MVT_illustration.png",
    width = 480*2, height = 480,)

par(mfrow=c(1,2), mai=c(1.02,1,0.82,0.42))
plot(eval(fGras),  type="l",ylim=c(0,1),
     ylab=expression(Gain^- time), xlab="Time",cex.lab = 1.5,cex.axis = 1.5)
lines(eval(fHerb), col = "green")

plot(log(TravelTimeG),log(x), type="l", ylab="Residence time (ln)",
     xlab="Travel time (ln)",xlim=c(0,3),cex.lab = 1.5,cex.axis = 1.5)
lines(log(TravelTimeH),log(x), col ="green", type="l")
dev.off()

# ### Following stephens 1993
# TravelTime = travel time
# m = rate parameter of the gain function
# LogArea = 10.228
# LogPrev = 4.355
# LogNext=4.316
# pframe[(LogNext<4.39 & LogNext>4.3) &
#             (LogPrev >4.33 & LogPrev < 4.4) &
#             (LogArea>10.2 & LogArea<10.25) ]
tc = 10# Time required to extract a certain proportion of enery
m = 1/tc
m2 = 2*m
TravelTime = 2
k = 1

toptim = (1/m)*(log(m*TravelTime+1) + ((log(m*TravelTime)+1)/(m*TravelTime + 1)))
sgain = k*(1-exp(-m*time))
curve(k*(1 - exp(-m*x)), from = , to = 100, ylab = "Sgain")


par(mfrow=c(1,2))
curve(k*(1-exp(-m*x)), from = 0, to = 30,
      xlab = "Residence time", ylab = "Gain")
curve(k*(1-exp(-m2*x)), from = 0, to = 30,
      xlab = "Residence time", ylab = "Gain", add = T,
      col = "green")
curve(.03*(x), from = 0, to = 30,
      xlab = "Residence time", ylab = "Gain", add = T,
      col = "red")

# For traveltime Tt
Tt = 10
Rt = (1/m)*(log(m*Tt + 1) +
              (log(m*Tt + 1)/(m*Tt + 1)))
curve((1/m)*(log(m*x + 1) +
               (log(m*x + 1)/(m*x + 1))),
      from = 0, to = 100,
      xlab="Travel time",
      ylab = "Optimal residence time")
curve((1/m2)*(log(m2*x + 1) +
                (log(m2*x + 1)/(m2*x + 1))),
      from = 0, to = 100, add = T, col = "green")

####
# For a given travelling time and two gain functions
# ÃstrÃ¸m et al. 1990
Thi = .5 #Handlint time/time in patch
k1 = .5# acceleration towards asympote
si = 2# patch size
Gain = si*(1-(1/(1+k1*Thi)))
Gain = si*(1-(1/(1+(k1*Thi/si))))
# Ts = travel time
Toptim = sqrt(Ts/k1)
k1a = .5
k1b = .2
Ts = 10

MVT_Aastroem = function(Thi = seq(.5,10,.25), rate = .5, si = 2, Ts = NULL, Tangent = TRUE, Toptim=NULL){
  # Thi = seq(.5,10,.25); rate = .5; si = 3; Ts = seq(2,12,2); Tangent = TRUE; Toptim=NULL
  
  # Thi = handling time / time in patch
  # k1 = rate = acceleration towards asymptote
  # si = patch size
  # Ts = average inter-patch travel time
  
  if(!is.null(Toptim) & !is.null(Ts)){
    # If we don't give...?
    k1 = Ts/(Toptim^2)
    k1 = unique(round(k1,3))
    if(length(k1)==1){
      Times = Thi
      Gain = (1-(1/(1+k1*Times)))
      attr(Gain, "k1")<-k1
      attr(Gain, "tstar")<-Times
    }else{print("No unique solution")}
    
  }else{
    lapply(rate, function(k1){
      Gain = data.frame(HandlingTime = Thi,
                        PatchSize = si,
                        Rate = k1)
      
      if(Method == "Aastrom"){
        Gain$Gain = with(Gain, PatchSize*(1-(1/(1+(Rate*HandlingTime))))) # Eq. 1 in the paper)
        if(!is.null(Ts)){
          # For a gain function with a given travel time
          # The optimal residence time for a given patch size is
          
          Gain = do.call("rbind", lapply(Ts, function(x) {
            df = Gain
            df$TravelTime = x
            df}))
          
          Gain$tstar = with(Gain, sqrt(TravelTime/Rate))# Eq. 3 in paper,
          
          if(Tangent){
            tmp = t( sapply(1:nrow(Gain), function(i){
              si = Gain$PatchSize[i]
              sa = mean(Gain$PatchSize)
              k1 = Gain$Rate[i]
              tOptim = Gain$tstar[i]
              
              derivateda = deriv(~si*(1-(1/(1+k1*tOptim))),"tOptim")
              
              sqrt((si*Ts)/(sa*k1)) + (1/k1)(sqrt(si/sa)-1) # eq. 6
              
              # Slope in this point is
              slope = attr(eval(derivateda),"gradient")
              gain = eval(derivateda)[1]
              # Where the tangent crosses the y-axis, given it should
              c(Slope = slope, Intercept = Gain$TravelTime[i]*slope, Gain = gain)
            }))
            Gain = cbind(Gain, tmp)
          }}
        
      }else if(Method == "Aastrom2"){
        Gain$Gain2 = with(Gain, PatchSize*(1-(1/(1+(Rate*HandlingTime*PatchSize))))) # Eq. 2 in the paper)
        
        if(!is.null(Ts)){
          # For a gain function with a given travel time
          # The optimal residence time for a given patch size is
          
          Gain = do.call("rbind", lapply(Ts, function(x) {
            df = Gain
            df$TravelTime = x
            df}))
          
          Gain$tstar = with(Gain, sqrt(TravelTime/Rate*mean(Gain$PatchSize)))# Eq. 4 in paper,
          
          if(Tangent){
            tmp = t( sapply(1:nrow(Gain), function(i){
              si = Gain$PatchSize[i]
              sa = mean(Gain$PatchSize)
              k1 = Gain$Rate[i]
              tOptim = Gain$tstar[i]
              derivateda = deriv(~si*(1-(1/(1+k1*tOptim*si))),"tOptim") # Eq. 7
              
              # Slope in this point is
              slope = attr(eval(derivateda),"gradient")
              gain = eval(derivateda)[1]
              # Where the tangent crosses the y-axis, given it should
              # intersect at x = -10, for y = 0.
              
              #Gain$slope<-slope
              #Gain$intercept<-Ts*slope
              #Gain$gain<-gain
              c(Slope = slope, Intercept = Gain$TravelTime[i]*slope, Gain = gain)
            }))
            Gain = cbind(Gain, tmp)
          }}
        
      }
      
      

      Gain
    }
    )
    
  }
  Gain
}

Times = seq(.5,10,.5)
ys = MVT_Aastroem(Thi = Times, si = 3, Ts = seq(2,12,2))

par(mfrow = c(2,2))
with(ys, plot(Gain~HandlingTime ,type = "l", xlim = c(-10,10)))

abline(v = attr(ys, "tstar"))
abline(a = attr(ys, "intercept"), b = attr(ys, "slope"))
plot(attr(ys, "TravelTimes")~attr(ys, "tstar"),type = "l")

a2 = MVT_Aastroem(Thi = Times, k1 = c(.25,.5))

#### ALTNERNATIVE
library(data.table);library(RColorBrewer);library(stringr)

MVT_loglog = function(BetaTimeSpent=c(.15,.30),x.values = Times){
  x = x.values
  lists =  lapply(BetaTimeSpent, function(i){
    fx = deriv(expr = y~x^i,"x")
    bs = eval(fx)
    y = bs[seq_along(x)]
    tangents = attr(bs,"grad")
    ints = y -tangents*x
    return(data.table(x, y, Slope = tangents, Intercept = ints))})
  datas =    rbindlist(lists)
  names(datas)[3]="Tangent"
  names(datas)[4]="xIntercept"
  betas = rep(BetaTimeSpent, each=length(x))
  datas[,TimeBetas:=betas]
  
  datas[,TravelTime := abs((0-xIntercept)/Tangent)]
  datas = datas[is.finite(Tangent) & is.finite(TravelTime) ]
  
  library(RColorBrewer)
  set =  row.names(subset(brewer.pal.info,colorblind==TRUE & category=="qual" & maxcolors>length(BetaTimeSpent)))[1]
  cols =  brewer.pal(length(BetaTimeSpent), set)
  cols = if(length(BetaTimeSpent)<3){cols[c(1,length(cols))]}else{cols}
  datas[,colr:=cols[as.numeric(as.factor(TimeBetas))]]
  # So, what are we doing here?
  d = table(cut(datas$TravelTime,seq(min(datas$TravelTime),max(datas$TravelTime),.5)),
            datas$TimeBetas)
  # For which intervals do we have estimates for all betas?
  d =(apply(d, c(1,2),sign))
  names(which(rowSums(d)>=length(BetaTimeSpent))[1])
  c = as.numeric(unlist(str_split(gsub("\\(|\\]","",
                                       names(which(rowSums(d)>=length(BetaTimeSpent))[1])),",")))
  c = datas[TravelTime>=c[1] & TravelTime<=c[2]]
  
  # Relationship btw cummulative gain, travel time and time spent in patch
  #par(mfrow=c(1,2))
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  # gain function
  par(mar=c(0, 4.1, 1, 2.1))
  with(datas,
       plot(c(1:5),c(1:5),xaxt="n",
            xlim = c(-(c$TravelTime[1] +20),max(x)),
            #xlim=c(-max(TravelTime),max(x)),
            ylim=c(0,max(y)),type="n",
            ylab = "", xlab ="TravelTime/TimeSpentInPatch"))
  mtext(cex=2, side = 2, text = "Gain", line =2)
  legend("topleft", legend = BetaTimeSpent, lty = c(1,1), col = cols,cex = 1.5, lwd = 2)
  text(-42, 2.47, expression(paste("Gain = ", Time^beta, sep = '')), pos = 2,cex=2)
  abline(v = 0, lty="dashed")
  for(i in 1:nrow(datas)){
    try(
      with(datas,abline(b = Tangent[i], a =xIntercept[i],
                        col= colr[i],lwd =2))
      ,silent=F)
  }
  
  # Illustrating how residence time with increasing quality for a given travel time
  par(mar=c(5.1, 4.25, 1, 2.1))
  with(c,
       plot(c(1:5),c(1:5), xlim=c(-max(TravelTime),max(x)+2),ylim=c(0,max(y)),type="n",
            ylab = "",
            xlab ="",
            xaxt="n"))
  mtext(cex=2, side = 2, text = "Gain", line =2)
  mtext(cex=1.5, side = 1, text = "TravelTime/TimeSpentInPatch", line =2.5)
  axis(side = 1, at = round(seq(-max(c$TravelTime),(max(x)+2), by = 10)), labels = abs(round(seq(-max(c$TravelTime),(max(x)+2), by = 10))))
  
  abline(v = 0, lty="dashed")
  for(i in seq_along(BetaTimeSpent)){
    print(i)
    curve(x^BetaTimeSpent[i], from = 0, to = max(x),add = T, col = cols[i],lwd =2)
    with(c[TimeBetas==BetaTimeSpent[i]],
         abline(a = xIntercept, b = Tangent, col = colr,lwd =2))
    with(c[TimeBetas==BetaTimeSpent[i]],
         abline(v = x, col = colr,lwd =2))
  }
  
  #abline(a = 0, b = mean(datas$Tangent),lty="dashed")
  # Relationship between travel time and patch time
  par(mar=c(5.1, 4.4, 1, 2.1))
  with(datas,plot(log(x)~log(TravelTime), type="n",
                  ylab = "", xlab="",
                  ylim=c(min(log(x)), max(log(x))), xlim=range(na.omit(log(TravelTime)))))
  for(i in seq_along(BetaTimeSpent)){
    with(datas[TimeBetas==BetaTimeSpent[i]],
         lines(y = log(x), x = log(TravelTime), col = colr,lwd =2))     }
  mtext(cex=2, side = 2, text = "Time spent \n in patch (ln)", line =2)
  mtext(cex=1.5, side = 1, text = "Travel time (ln)", line =2.5)
  
  par(mfrow=c(1,1))
  return(datas)
}
Times = seq(.5,10,.5)


png("foo.png", width = 800, height = 800)
a = MVT_loglog(BetaTimeSpent = c(.15,.3), x.values = Times)
dev.off()

with(a,
     plot(c(1:5),c(1:5), xlim=c(-max(TravelTime),max(x)),ylim=c(0,max(y)),type="n",
          ylab = "Gain", xlab ="TravelTime/TimeSpentInPatch"))
abline(v = 0)
a[TravelTime >38 & TravelTime <40]
for(i in 1:nrow(a)){
  try(
    with(a,abline(b = Tangent[i], a =xIntercept[i],
                  col= colr[i]))
    ,silent=F)
}

# The different shapes will result in different linear relationships
# btw patch residency and travel time

with(a, plot(log(x)~log(TravelTime), col = colr, type = "p", pch = 19))
with(a, plot(log(x)~log(TravelTime), col = colr, type = "p", pch = 19))
