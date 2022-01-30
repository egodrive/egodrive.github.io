### DATA #####
library(data.table);library(ggplot2)
MVTm = data.table(read.csv2("C:/Users/endre/Downloads/MVT.csv", header=F))
names(MVTm)<-c("Distance", "Residence_time")
setorder(MVTm, Distance, Residence_time)
MVTm[,Species:=ifelse(Residence_time<140, "Myrmica_ruginodis", "Lasius_niger")]
MVTm[,Species:=ifelse(Residence_time>100 & Distance<1.5, "Lasius_niger", Species)]
MVTm$Distance = round(MVTm$Distance, 2)
MVTm$Residence_time = round(MVTm$Residence_time, 2)

setorder(MVTm, Species, Distance, Residence_time)
MVTm[,Measure:=rep(c("M", "H"), .N/2), Species]

#Manual corrections
MVTm[Distance<1.3 & Residence_time>155]$Measure<-c("H","M")
MVTm[Distance<2.3 & Residence_time>240]$Measure<-c("M","H")

setorder(MVTm, Species, Distance, Measure)
MVTm$idx = rep(1:20, each = 2)

MVTm[Residence_time %between% c(100,150) & Distance %between% c(.9,1.4), idx := 3]
MVTm[Residence_time %between% c(150,200) & Distance %between% c(.9,1.4), idx := 2]

MVTm[Residence_time %between% c(150,240) & Distance %between% c(1.6,2.4)]$idx = 5
MVTm[Residence_time %between% c(240,280) & Distance %between% c(1.6,2.4)]$idx = 6


MVTm[Residence_time %between% c(150,200) & Distance %between% c(1.6,2.4), idx := 2]

MVTm[,StdError:=max(Residence_time)-min(Residence_time), "idx"]
MVTm[,Distance_2:=Distance^2]

MVTm = MVTm[Measure=="M"]

write.csv2(MVTm, file = "F:/egodrive.github.io/Bonser_et_al_1998_MVT.csv", row.names = F)

ggplot(data = MVTm, aes(y = Residence_time, x = Distance, col = Measure, shape = Species))+geom_point(size = 3)+
  geom_errorbar(data = MVTm[Measure=="M"],aes(ymax = Residence_time+StdError, ymin = Residence_time-StdError), width = 0.2)


library(MuMIn)
dredge(lm(Residence_time~Distance*Species + Distance_2*Species, data = MVTm, weights = 1/StdError, na.action = na.fail))
dredge(lm(Residence_time~log(Distance)*Species + log(Distance_2)*Species, data = MVTm, weights = 1/StdError, na.action = na.fail))


summary((lm(log(Residence_time)~log(Distance) + Species, data = MVTm, weights = 1/StdError)))

summary((lm(log(Residence_time)~log(Distance) * Species, data = MVTm, weights = 1/StdError)))

step((lm(log(Residence_time)~log(Distance*2) * Species, data = MVTm, weights = 1/StdError)))
step((lm((Residence_time)~(Distance) * Species, data = MVTm, weights = 1/StdError)))


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
# General form of the gain function
# [1] Charnov EL, Parker GA. 1995 Dimensionless Invariants from Foraging Theory Marginal Value Theorem. Proc. Natl. Acad. Sci. USA 92, 1446-1450

MVT_1 = function(Method = "Standard", Asymptote = c(.5,1), Rate = c(0.5,1), HalfTime = 4, yInt = c(1,2), HandlingTime = seq(0,10,.25)){
  # Asymptote = c(.5,1); Rate = c(0.5,1.2); yInt = 0; HandlingTime = seq(0,10,.25);HalfTime = 4;Method = "Sigmoid"
  library(gridExtra);library(data.table)
  
  df = CJ(Asymptote = Asymptote,
                   Rate = Rate,
                   Intercept = yInt, HandlingTime = HandlingTime)
  
  # To find the tangents
  # Computes the derivatives
  x = df$HandlingTime
  df$x = x
  if(Method == "Standard"){
    GainFunction = deriv(~Asymptote +(Intercept - Asymptote)*exp(-exp(Rate)*x),"x")
  }
  if(Method == "CharnovParker1995"){
    #Asymptote: Gm = .5 # Proportionate to patch size
    #
    #Rate: 
    #     g = .55 # Constant reflecting initial capture rate of prey
    #     c = g/Gm # rate of prey capture
    # Handling Time = t
    # GeneralGainForm = function(Gmax, c, t){  Gmax*(1-exp(-c*t))}
    
    GainFunction = deriv(~Asymptote*(1-exp(-Rate*x)),"x")
  }
  if(Method=="Aastrom_1"){
    # Asymptote = si = Patch size
    # Rate = k1
    # Handling time = Thi
    # Gain = si*(1-(1/(1+(k1*Thi)))) # Eq. 1 in the paper)
    
    GainFunction = deriv(~Asymptote*(1-(1/(1+Rate*x))),"x") # Original
    GainFunction = deriv(~Asymptote -(Asymptote/(1+Rate*x)),"x") # Omskrevet
    
  }
  if(Method=="Aastrom_2"){
    # Asymptote = si = Patch size
    # Rate = k1
    # Handling time = Thi
    # Gain = si*(1-(1/(1+((k1*Thi)/si)))) # Eq. 2 in the paper)
    
    GainFunction = deriv(~Asymptote*(1-(1/(1+(Rate*x/Asymptote)))),"x") # original fra paper
    GainFunction = deriv(~Asymptote - (Asymptote^2)/(Asymptote+Rate*x),"x") # omskrevet
    
  }
  if(Method=="LogLog"){
    GainFunction = deriv(expr = ~exp(Intercept)*x^Rate,"x")
  }
  if(Method=="Sigmoid"){
    # Gompertz
    # Asymptote = 1 # asympmtote, opprinnelig a
    # b: half-way, x displacement
    # Rate: c
    #HalfWay = mean(HandlingTime)
    # print("Set halftime of the function to be equal to mean of handling times under consideration")
    # Logistic function
    # Asymptote/(1 + exp(-Rate*(x-HalfTime)))
    
    GainFunction = deriv(~Asymptote/(1 + exp(-Rate*(x-HalfTime))),"x")
    
  }
  
  # The y - values along the x-values
  bsG = eval(GainFunction, envir = df)
  df$Gain = bsG[seq_along(x)]
  
  # The derivatives
  df$TangentG = c(attr(bsG,"grad"))
  # Second derivative
  df[,SecDer:= c(NA,diff(TangentG)), c("Rate", "Asymptote")]
  
  # And where the tangent line crosses the x-axis
  df$xInterceptG = with(df, Gain -TangentG*HandlingTime)
  
  df$TravelTimeG = with(df, abs((0-xInterceptG)/TangentG))
  
  
  ### PLOTS
  df = data.table(df)
  df[,Rate := as.factor(Rate)]
  df[,Asymptote := as.factor(Asymptote)]
  library(ggpubr)
  p1 = ggplot(data = df, aes(x = HandlingTime , y = Gain, col = Rate, linetype = Asymptote, group = interaction(Rate, Asymptote, Intercept)))+
    geom_line()
  
  p2 = ggplot(data = df, aes(x = HandlingTime , y = TangentG, col = Rate, linetype = Asymptote, group = interaction(Rate, Asymptote, Intercept)))+
    geom_line() + ylab("Functional response")
  
  
  p3 = ggplot(data = df, aes(x = TravelTimeG , y = HandlingTime, col = Rate, linetype = Asymptote, group = interaction(Rate, Asymptote, Intercept)))+
    geom_line()
  
  p4 = ggplot(data = df, aes(x = log(TravelTimeG) , y = log(HandlingTime), col = Rate, linetype = Asymptote, group = interaction(Rate, Asymptote, Intercept)))+
    geom_line()
  
  print(ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom"))
  return(as.data.table(df))
}
MVT_1(HandlingTime = seq(0,3,.05), Rate = c(1.2,1.2), yInt = c(1,2), HalfTime = 1.5, Method = "LogLog")



MVT_1(HandlingTime = seq(0,3,.05), Rate = c(1.2,.75), Asymptote = c(1,2), HalfTime = 1.5, Method = "Aastrom_1")
MVT_1(HandlingTime = seq(0,3,.05), Rate = c(1.2,.75), Asymptote = c(1,2), HalfTime = 1.5, Method = "Aastrom_2")


head(a)


#### including variability: https://www.sciencedirect.com/science/article/abs/pii/S000334720580226X
