library(data.table)
nSimulations = 10
tmax = 200
ParameterSpace = expand.grid(nIterations = 1:nSimulations,
                             #RangeØyer = seq(40,60,length.out = 1),
                             RangeStartProp = seq(.35,.6,length.out = 1),
                             RangeLandskap = seq(3.4,3.7,length.out = 1),
                             RangeKonstantC = seq(.24,.75,length.out =  1),# .23->.24
                             RangeKonstantE = seq(60,100,length.out = 1),# 325->320
                             eNoise = seq(0,3,length.out = 1))
setorder(ParameterSpace, nIterations)

ParameterSpace$SnittPrevalens = NA
ParameterSpace$SdPrevalens = NA
ParameterSpace$MinPrevalens = NA
ParameterSpace$AndelUtvikling = NA
iteration = 0

RangeØyer = 8
nØyer = RangeØyer # ParameterSpace$RangeØyer[i] #i

TidsArray = array(rep(0, tmax*nØyer), 
                  dim = c(tmax, nØyer,  nrow(ParameterSpace)))

TidsserieKoloniseringsrater = array(rep(0, tmax*nØyer), 
                                    dim = c(tmax, nØyer, nrow(ParameterSpace)))

GeografiArray = array(rep(NA, 3*nØyer),
                      dim = c(3, nØyer, nrow(ParameterSpace)))


for(i in 1:nrow(ParameterSpace)){
  #print(paste(i , "av", nrow(ParameterSpace)))
  # Antall øyer/habitat i systemet
  
  # lengde på tidsserie
  tmax = 100
  
  iteration <- ParameterSpace$nIterations[i]
  
  AntallBebodd = ParameterSpace$RangeStartProp[i]*nØyer #i
  set.seed(iteration) # slik at ulike miljøstøy-nivå vil testes på de samme øy-konfigurasjonene
  Bebodd = sample(1:nØyer,AntallBebodd)
  
  
  Tilstede = rep(0, nØyer)
  Tilstede[Bebodd]<-1
  
  TidsserieTilstede = matrix(rep(rep(NA, nØyer), tmax), ncol = nØyer)
  TidsserieTilstede[1,]<-Tilstede
  
  TidsArray[1,,iteration]<-Tilstede
  TidsserieKoloniseringsrater[1,,iteration]<-Tilstede
  
  # tilfeldig plassering av øyene
  ArenaStørrelse = ParameterSpace$RangeLandskap[i]
  set.seed(iteration) # slik at ulike miljøstøy-nivå vil testes på de samme øy-konfigurasjonene
  xPos = runif(n = nØyer, min = 0, max = ArenaStørrelse)
  set.seed(-iteration) # slik at ulike miljøstøy-nivå vil testes på de samme øy-konfigurasjonene
  yPos = runif(n = nØyer, min = 0, max = ArenaStørrelse)
  
  # Tilfeldig størrelse
  set.seed(iteration) # slik at ulike miljøstøy-nivå vil testes på de samme øy-konfigurasjonene
  ØyStr = rlnorm(n = nØyer, meanlog = .1, sdlog = .2) #Ai
  
  GeografiArray[,,i]<-rbind(ØyStr, xPos, yPos)
  # Regn ut mellom-øy avstand
  MellomØyAvstand = dist(cbind(xPos, yPos), diag = T, upper = T) # dij
  Alpha = -1#1/mean(MellomØyAvstand)
  # Konstanter
  c_konstant = ParameterSpace$RangeKonstantC[i]
  e_konstant = ParameterSpace$RangeKonstantE[i]
  
  eNoise = ParameterSpace$eNoise[i]
  # Utryddelsesrater er konstant over tid
  UtryddelsesRateØy_i = e_konstant/ØyStr
  
  # fra Hanski, Ilkka, et al. "The quantitative incidence function model and persistence of an endangered butterfly metapopulation." Conservation Biology 10.2 (1996): 578-590.
  # ØystrEffekt = 0.952
  # muTick = 0.158
  
  # AjAi = sapply(ØyStr, function(x) x*ØyStr)
  # mij = exp(-Alpha*as.matrix(MellomØyAvstand))*AjAi # number of immigrants
  
  # Tidssteg 
  t = 1
  
  NoiseE = rnorm(n = tmax, sd = eNoise)
  while(t<tmax){
    
    KoloniseringsRateØy_i = sapply(1:nØyer, function(i){
      Aj = ØyStr[-i]
      Aj = Aj + NoiseE
      Pj = Tilstede[-i]
      dij = as.matrix(MellomØyAvstand)[-i,i]
      c_konstant*sum(Pj*Aj*
                       exp(-Alpha*dij))
      }) # Si
    
    
    TidsserieKoloniseringsrater[t,,i]<-KoloniseringsRateØy_i
    
    
    # Sannsynligheten for bebodd øy
    Pi = KoloniseringsRateØy_i/(KoloniseringsRateØy_i + UtryddelsesRateØy_i) # Også kjent som Ji
    Pi = ifelse(Pi>1,1,ifelse(Pi<0,0,Pi))
    
    Kolonisert = sapply(Pi, function(x) sample(0:1, 1, prob = c(1-x,x)))
    
    # Oppdater bebodd-status
    Tilstede<-Kolonisert
    TidsserieTilstede[t,]<-Tilstede
    TidsArray[t,,i]<-Tilstede
    
    # Oppdater tidssteg
    t = t+1
  }
  
  ParameterSpace$SnittPrevalens[i] = median(rowMeans(TidsserieTilstede), na.rm = T)
  ParameterSpace$SdPrevalens[i] = sd(rowMeans(TidsserieTilstede), na.rm = T)
  ParameterSpace$MinPrevalens[i] = min(rowMeans(TidsserieTilstede), na.rm = T)
  ParameterSpace$AndelUtvikling[i] = coef(lm(rowMeans(TidsserieTilstede)~I(1:length(rowMeans(TidsserieTilstede)))))[2]
  # if(i==nrow(ParameterSpace)){
  #   saveRDS(list(ParameterSpace = ParameterSpace,
  #                Geografi = GeografiArray,
  #                TidsArray = TidsArray,
  #                TidsserieKoloniseringsrater = TidsserieKoloniseringsrater),
  #           "ParameterSpace_MetaPopStoch.rds")
  # }
}



library(ggplot2);library(cowplot);library(scales);library(gridExtra)
# p1 = ggplot(data = ParameterSpace, aes(y = MinPrevalens, x = eNoise, group = factor(eNoise)))+geom_boxplot()+theme_cowplot()+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+ylab("Andel (minimum)")
# p2 = ggplot(data = ParameterSpace, aes(y = SnittPrevalens, x = eNoise, group = factor(eNoise)))+geom_boxplot()+theme_cowplot()+
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank())+ylab("Andel (gj.snitt)")
# 
# p3 = ggplot(data = ParameterSpace, aes(y = AndelUtvikling*200, x = eNoise, group = factor(eNoise)))+geom_hline(yintercept = 0, linetype="dashed")+theme_cowplot()+
#   geom_boxplot()+ylab("Andelsvekst per år x 200") + xlab("Miljøstokastisitet")#+
#scale_y_continuous(breaks = seq(-.2,.2,.1))+coord_cartesian(ylim = c(-0.1,0.1))

Geografi = GeografiArray#aa$Geografi
TidsArray = TidsArray#aa$TidsArray
TidsserieKoloniseringsrater = TidsserieKoloniseringsrater#aa$TidsserieKoloniseringsrater

#### Centrality ####
# Eigenvector centrality er et mål på hvor mange noder som peker til deg. 
# Per parametersetting
Centrality = lapply(1:nrow(Geografi), function(x){
  print(x)
  # Per landskap
  dij = as.matrix(dist(cbind(Geografi[2,,x], Geografi[3,,x]), diag = T, upper = T))
  
  Centrality = sapply(1:nrow(dij), function(i){
    js = sapply(1:ncol(dij), function(j){
      exp(dij[i,j])*Geografi[1,i,x]*Geografi[1,j,x]
    })
    js
  })
  diag(Centrality)<-0 # set i = j <- 0
  Eigens = eigen(Centrality)
  EVec = Eigens$vectors[,1]^2
  EVec = EVec/sum(EVec)
  
  list(Eigenvector = EVec,
       EigenValue = Eigens$value[1])
}) 

ParameterSpace$EigenvalueCentrality<-lapply(Centrality, function(x) x$EigenValue)

## Connectivity ####
Alpha = 1/1 # Alpha beskriver overlevelse for spredningen, dvs. 1/gjennomsnittlig spredningsdistanse
# Per parametersetting
Connectivity = lapply(1:dim(Geografi)[3], function(x){
  # Per landskap
  dij = as.matrix(dist(cbind(Geografi[2,,x], Geografi[3,,x]), diag = T, upper = T))
  # Geographic aspect
  Geo = exp(-Alpha*dij)*Geografi[1,,x]*dij 
  # Per tidssteg
  Si_per_t = t(sapply(1:nrow(TidsArray), function(T){
    # per øy
    Si = sapply(1:ncol(TidsArray), function(Ø){
      Si = sum(Alpha*TidsArray[T,,x]*Geo[Ø,])
    })
    Si
  }))
}) 

length(Connectivity)

Connectivities = sapply(Connectivity, function(x){
  colMeans(x)
})


### Eksempel
df = as.data.frame(t(Geografi[,,1]))
names(df)<- c("Areal", "xPos", "yPos")

#df$Connectivity = colMeans(t(Connectivities))
df$Connectivity = t(Connectivities)[1,]
df$EigenCentrality = Centrality[[1]]$Eigenvector^2*Centrality[[1]]$EigenValue

sd(t(Connectivities)[1,])/mean(t(Connectivities)[1,])

ggplot()+
  geom_point(data = df, aes(x = xPos, y = yPos, size = Connectivity, col = Connectivity))
ggplot()+
  geom_point(data = df, aes(x = xPos, y = yPos, size = log(EigenCentrality), col = log(EigenCentrality)))

ParameterSpace = as.data.table(ParameterSpace)
ParameterSpace[,AverageConnectivity:=colMeans(Connectivities)]

#plot(ParameterSpace[,6:12])

# Hettemåke foraging flight 4.6-11.8 km flight range # https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.6291
# sildemåke 30.9 #https://www.researchgate.net/publication/306119852_Terrestrial_and_Marine_Foraging_Strategies_of_an_Opportunistic_Seabird_Species_Breeding_in_the_Wadden_Sea/link/57b2ef8508aeaf239baefb27/download
ggplot(mapping = aes(x = log(df$EigenCentrality), y = colMeans(TidsserieKoloniseringsrater[,,1])))+
  geom_point()+
  stat_smooth(method="lm")

ggplot(mapping = aes(x = colMeans(Connectivity[[1]])*df$Areal, y = colMeans(TidsserieKoloniseringsrater[,,1])))+
  geom_point()+
  stat_smooth(method="lm")

with(df, sum(Connectivity*Areal)/(sum(Areal)))


### Konkretiser ####
library(data.table)
VO = fread("F:/egodrive.github.io/Verneområder_centroider.csv")
#VO = VO[-1,]
#Mosvatnet
531039,61; N: 6539879 Ø: 310963
# Klostervågen
847287,21; N: 6556056 Ø: 304732

Geografi = VOgeo = as.matrix(VO[,c("Areal", "PosX", "PosY")])
VOgeo1 = rbind(VOgeo, c(847287.21, 304732, 6556056)) # Klostervågen
VOgeo2 = rbind(VOgeo, c(531039.61, 310963, 6539879)) # Mosvatnet
VOgeo3 = rbind(VOgeo2, c(847287.21, 304732, 6556056)) # Begge to

Centrality = function(Geografi, DispersalDistance = 30.9, Alpha = NULL, unit = "m"){
  # Geografi = dataramme med areal, x og y koordinate
  Areal = Geografi[,1]
  if(unit=="m"){
    Areal = Areal/1000000
    }
  dij = as.matrix(dist(cbind(Geografi[,2], Geografi[,3]), diag = T, upper = T))
  if(unit == "m"){
    dij = dij/1000# transform from meter to kilometer
    
  }
  
  Alpha = 1/DispersalDistance
  
  Si = sapply(1:nrow(Geografi), function(i){
    Geo = sum(exp(-Alpha*dij[i,-i])*Geografi[-i,1]) 
  })
  
  M = sapply(seq_along(Areal), function(i){
    sapply(seq_along(Areal), function(j){
      exp(-Alpha*dij[i,j])*Areal[i]*Areal[j]
    })
  })
  
  diag(M)<-0
  ev = eigen(M)$vectors[,1]^2
  ev = ev/(sum(ev))
  we = log(ev*eigen(M)$value[1])
  
  ConnectedHabitat = (sum(Si*Areal)/(sum(Areal)))*(sum(Areal^2)/sum(Areal)) # Hanski 1999
  
  ls = list(Centrality = Si,
            Eigenvectors = we,
            ConnectedHabitat = ConnectedHabitat,
            EigenValue = eigen(M)$value[1]
            )
  
  attr(ls, "Areal")<-Areal
  attr(ls, "xPos")<-Geografi[,2]
  attr(ls, "yPos")<-Geografi[,3]
  ls
}

Measures = Centrality(Geografi = VOgeo)

cor(Measures$Centrality, Measures$Eigenvectors)

p0_c = ggplot()+
  geom_point(mapping = aes(x = VO$PosX, y = VO$PosY,
                           size = Centrality(Geografi = VOgeo)$Centrality,
                           col =Centrality(Geografi = VOgeo)$Centrality))+
  theme(legend.position = "none")

p1_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo1[,2], y = VOgeo1[,3],
                           size = scale(Centrality(Geografi = VOgeo1)$Centrality),
                           col = scale(Centrality(Geografi = VOgeo1)$Centrality)))+
  theme(legend.position = "none")+
  geom_point(data = as.data.frame(tail(VOgeo1,1)), 
             mapping = aes(x = PosX, y = PosY),
             col = "red")

p2_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo2[,2], y = VOgeo2[,3],
                           size = scale(Centrality(Geografi = VOgeo2)$Centrality),
                           col = scale(Centrality(Geografi = VOgeo2)$Centrality)))+
  theme(legend.position = "none")+
  geom_point(data = as.data.frame(tail(VOgeo2,1)), mapping = aes(x = PosX, y = PosY), col = "red")

p3_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo3[,2], y = VOgeo3[,3],
                           size = scale(Centrality(Geografi = VOgeo3)$Centrality),
                           col = scale(Centrality(Geografi = VOgeo3)$Centrality)))+
  theme(legend.position = "none")+
  geom_point(data = as.data.frame(tail(VOgeo3,2)), mapping = aes(x = PosX, y = PosY), col = "red")

library(gridExtra)
grid.arrange(p0_c, p1_c, p2_c, p3_c)




#### EIGENVECTORS ####
p0_c = ggplot()+
  geom_point(mapping = aes(x = VO$PosX, y = VO$PosY,
                           size = log(Centrality(Geografi = VOgeo)$Eigenvectors),
                           col =log(Centrality(Geografi = VOgeo)$Eigenvectors)))+
  geom_label(mapping = aes(x = min(VOgeo3[,2]), y = max(VOgeo3[,3]),
                          label = paste("Metapopulasjonskapasitet = ",
                                        round(Centrality(Geografi = VOgeo)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

p1_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo1[,2], y = VOgeo1[,3],
                           size = log(Centrality(Geografi = VOgeo1)$Eigenvectors),
                           col = log(Centrality(Geografi = VOgeo1)$Eigenvectors)))+
  geom_point(data = as.data.frame(tail(VOgeo1,1)), mapping = aes(x = PosX, y = PosY), col = "red")+
  geom_label(mapping = aes(x = min(VOgeo1[,2]), y = max(VOgeo1[,3]), 
                          label = paste("Metapopulasjonskapasitet = ",
                                        round(Centrality(Geografi = VOgeo1)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

p2_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo2[,2], y = VOgeo2[,3],
                           size = log(Centrality(Geografi = VOgeo2)$Eigenvectors),
                           col = log(Centrality(Geografi = VOgeo2)$Eigenvectors)))+
  geom_point(data = as.data.frame(tail(VOgeo2,1)), mapping = aes(x = PosX, y = PosY), col = "red")+
  geom_label(mapping = aes(x = min(VOgeo2[,2]), y = max(VOgeo2[,3]),
                          label = paste("Metapopulasjonskapasitet = ",
                                        round(Centrality(Geografi = VOgeo2)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

p3_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo3[,2], y = VOgeo3[,3],
                           size = log(Centrality(Geografi = VOgeo3)$Eigenvectors),
                           col = log(Centrality(Geografi = VOgeo3)$Eigenvectors)))+
  geom_point(data = as.data.frame(tail(VOgeo3,2)), mapping = aes(x = PosX, y = PosY), col = "red")+
  geom_label(mapping = aes(x = min(VOgeo3[,2]), y = max(VOgeo3[,3]),
                          label = paste("Metapopulasjonskapasitet = ",
                                        round(Centrality(Geografi = VOgeo3)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

th = theme_cowplot()+theme(axis.text=element_blank(),axis.title=element_blank(),legend.position = "none")
library(gridExtra)
grid.arrange(p0_c+th, p1_c+th, p2_c+th, p3_c+th)

library(ggmap)
library(mapproj)
library(cowplot)

Plotting = function(Geografi = VOgeo3, AntallNye = 0){
  InputData = as.data.frame(Geografi)
  
  bb = bbox(SpatialPoints(coords = InputData[,2:3], 
                          proj4string = CRS("+init=epsg:32632")))
  
  VOs = data.frame(coordinates(spTransform(SpatialPoints(coords = InputData[,2:3], 
                                                         proj4string = CRS("+init=epsg:32632")),
                                           CRS("+init=epsg:4326"))))
  VOs$EV = log(Centrality(Geografi = InputData)$Eigenvectors)
  
  
  ## Plot kartet
  aaa = get_map(c(left = min(VOs[,1])-.2, 
                  right = max(VOs[,1])+.2,
                  bottom = min(VOs[,2])-.2, 
                  top = max(VOs[,2])+.2))
  
  Labs = data.frame(Y = max(VOs$PosY),#as.numeric(attr(aaa, "bb")[3]),
                    X = min(VOs$PosX),#as.numeric(attr(aaa, "bb")[2]), 
                    MC = paste("Metapopulasjonskapasitet = ",
                               round(Centrality(Geografi = 
                                                  InputData)$EigenValue,2)))
  
  ggmap(aaa)+
    xlab("Lengdegrad")+ ylab("Breddegrad")+theme_cowplot()+
    geom_point(data= VOs, 
               mapping = aes(x = PosX, y = PosY,
                             size = EV, col = EV))+
    geom_point(data = as.data.frame(tail(VOs,AntallNye)),
               mapping = aes(x = PosX, y = PosY), col = "red")+
    geom_label(data = Labs, 
               mapping = aes(y = Y, x = X,
                             label = MC),
               fontface = "bold", fill = alpha("white",0.5), hjust = 0.05, vjust = -.7)+
    th
}

Plotting(VOgeo1, AntallNye = 2)


grid.arrange(Plotting(VOgeo)+coord_fixed(),
             Plotting(VOgeo1, AntallNye = 1)+coord_fixed(),
             Plotting(VOgeo2, AntallNye = 1)+coord_fixed(),
             Plotting(VOgeo3, AntallNye = 2)+coord_fixed())


### Konkretiser ####
library(data.table)
a = read.csv("Verneområder_centroider.csv", encoding = "UTF-8")
VO = as.data.table(a)

#VO = VO[-1,]
#Mosvatnet
#531039,61; N: 6539879 Ø: 310963
# Klostervågen
#847287,21; N: 6556056 Ø: 304732

Geografi = VOgeo = as.matrix(VO[,c("Areal", "PosX", "PosY")])
VOgeo1 = rbind(VOgeo, c(847287.21, 304732, 6556056)) # Klostervågen
VOgeo2 = rbind(VOgeo, c(531039.61, 310963, 6539879)) # Mosvatnet
VOgeo3 = rbind(VOgeo2, c(847287.21, 304732, 6556056)) # Begge to

VOGeo4 = rbind(VOgeo3, c(847287.21, 304532, 6556256)) # Fiktiv
VOGeo4 = cbind(VOGeo4, Ny = c(rep(0,nrow(VOgeo)), rep(1,3)))



Measures = Centrality(Geografi = VOgeo)

cor(Measures$Centrality, log(Measures$Eigenvectors))
library(ggplot2)
# p0_c = ggplot()+
#   geom_point(mapping = aes(x = VO$PosX, y = VO$PosY,
#                            size = Centrality(Geografi = VOgeo)$Centrality,
#                            col =Centrality(Geografi = VOgeo)$Centrality))+
#   theme(legend.title = element_blank(),legend.position = "none")
# 
# p1_c = ggplot()+
#   geom_point(mapping = aes(x = VOgeo1[,2], y = VOgeo1[,3],
#                            size = scale(Centrality(Geografi = VOgeo1)$Centrality),
#                            col = scale(Centrality(Geografi = VOgeo1)$Centrality)))+
#   theme(legend.title = element_blank(),legend.position = "none")+
#   geom_point(data = as.data.frame(tail(VOgeo1,1)), mapping = aes(x = PosX, y = PosY), col = "red")
# 
# p2_c = ggplot()+
#   geom_point(mapping = aes(x = VOgeo2[,2], y = VOgeo2[,3],
#                            size = scale(Centrality(Geografi = VOgeo2)$Centrality),
#                            col = scale(Centrality(Geografi = VOgeo2)$Centrality)))+
#   theme(legend.title = element_blank(),legend.position = "none")+
#   geom_point(data = as.data.frame(tail(VOgeo2,1)), mapping = aes(x = PosX, y = PosY), col = "red")
# 
# p3_c = ggplot()+
#   geom_point(mapping = aes(x = VOgeo3[,2], y = VOgeo3[,3],
#                            size = scale(Centrality(Geografi = VOgeo3)$Centrality),
#                            col = scale(Centrality(Geografi = VOgeo3)$Centrality)))+
#   theme(legend.title = element_blank(),legend.position = "none")+
#   geom_point(data = as.data.frame(tail(VOgeo3,2)), mapping = aes(x = PosX, y = PosY), col = "red")
# 
# library(gridExtra);library(cowplot)
# grid.arrange(p0_c, p1_c, p2_c, p3_c)
# 
# EIGENVECTORS
p0_c = ggplot()+
  geom_point(mapping = aes(x = VO$PosX, y = VO$PosY,
                           size = log(Centrality(Geografi = VOgeo)$Eigenvectors),
                           col =log(Centrality(Geografi = VOgeo)$Eigenvectors)))+
  geom_label(mapping = aes(x = min(VOgeo3[,2]), y = max(VOgeo3[,3]),
                           label = paste("Metapopulasjonskapasitet = ",
                                         round(Centrality(Geografi = VOgeo)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

p1_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo1[,2], y = VOgeo1[,3],
                           size = log(Centrality(Geografi = VOgeo1)$Eigenvectors),
                           col = log(Centrality(Geografi = VOgeo1)$Eigenvectors)))+
  geom_point(data = as.data.frame(tail(VOgeo1,1)), mapping = aes(x = PosX, y = PosY), col = "red")+
  geom_label(mapping = aes(x = min(VOgeo1[,2]), y = max(VOgeo1[,3]), 
                           label = paste("Metapopulasjonskapasitet = ",
                                         round(Centrality(Geografi = VOgeo1)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

p2_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo2[,2], y = VOgeo2[,3],
                           size = log(Centrality(Geografi = VOgeo2)$Eigenvectors),
                           col = log(Centrality(Geografi = VOgeo2)$Eigenvectors)))+
  geom_point(data = as.data.frame(tail(VOgeo2,1)), mapping = aes(x = PosX, y = PosY), col = "red")+
  geom_label(mapping = aes(x = min(VOgeo2[,2]), y = max(VOgeo2[,3]),
                           label = paste("Metapopulasjonskapasitet = ",
                                         round(Centrality(Geografi = VOgeo2)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

p3_c = ggplot()+
  geom_point(mapping = aes(x = VOgeo3[,2], y = VOgeo3[,3],
                           size = log(Centrality(Geografi = VOgeo3)$Eigenvectors),
                           col = log(Centrality(Geografi = VOgeo3)$Eigenvectors)))+
  geom_point(data = as.data.frame(tail(VOgeo3,2)), mapping = aes(x = PosX, y = PosY), col = "red")+
  geom_label(mapping = aes(x = min(VOgeo3[,2]), y = max(VOgeo3[,3]),
                           label = paste("Metapopulasjonskapasitet = ",
                                         round(Centrality(Geografi = VOgeo3)$EigenValue,2))), fontface = "bold", fill = alpha("white",0.5), hjust = -.1)

th = theme_cowplot()+theme(axis.text=element_blank(),axis.title=element_blank(),legend.position = "none")
library(gridExtra)
grid.arrange(p0_c+th, p1_c+th, p2_c+th, p3_c+th)

library(ggmap)
library(mapproj)
library(cowplot)
library(splancs)

Plotting = function(Geografi = VOGeo4, AntallNye = 1:2, Iterer = T){
  
  InputData = as.data.frame(Geografi)
  
  bb = bbox(SpatialPoints(coords = InputData[,2:3], 
                          proj4string = CRS("+init=epsg:32632")))
  
  Results = rbindlist(lapply(AntallNye, function(i){
    
    aa = combn(1:nrow(InputData[InputData$Ny==1,]),m = i)
    Results = rbindlist(lapply(1:ncol(aa), function(ii) {
      InputData2 = rbind(InputData[InputData$Ny==0,],
                         InputData[InputData$Ny==1,][aa[,ii],])
      VOs = data.frame(coordinates(spTransform(SpatialPoints(coords = InputData2[,2:3], 
                                                             proj4string = CRS("+init=epsg:32632")),
                                               CRS("+init=epsg:4326"))))
      VOs$EV = log(Centrality(Geografi = InputData2)$Eigenvectors)
      
      Results = data.frame(t(aa[,ii]),Centrality = round(Centrality(Geografi = 
                                                                      InputData2)$EigenValue,2))
      names(Results)[1:length(unique(aa[,ii]))]<-unique(aa[,ii])
      Results[1:length(unique(aa[,ii]))]<-"X"
      data.table(Results)
    }), fill = T, use.names = T)
    
  }), use.names = T)
  Results[,nNye:=apply(as.matrix(Results),1, function(x){length(unique(which(x=="X")))})]
  Results = cbind(Results[,-"Centrality"]
                  ,Results[,"Centrality"])
  setorder(Results, -Centrality)
  Results
  
  ## Plot kartet
  InputData[]
  
  GEOS0 = rbind(InputData[InputData$Ny==0,],
                InputData[InputData$Ny==1,][as.numeric(names(Results)[which(Results[1,]=="X")]),])
  
  GEOS = data.frame(coordinates(spTransform(SpatialPoints(coords = GEOS0[,2:3], 
                                                          proj4string = CRS("+init=epsg:32632")),
                                            CRS("+init=epsg:4326"))))
  
  GEOS$EV = log(Centrality(Geografi = GEOS0)$Eigenvectors)
  
  
  aaa = get_map(c(left = min(GEOS[,1])-.2, 
                  right = max(GEOS[,1])+.2,
                  bottom = min(GEOS[,2])-.2, 
                  top = max(GEOS[,2])+.2))
  
  Labs = data.frame(Y = max(VOs$PosY),#as.numeric(attr(aaa, "bb")[3]),
                    X = min(VOs$PosX),#as.numeric(attr(aaa, "bb")[2]), 
                    MC = paste("Metapopulasjonskapasitet = ",
                               round(Centrality(Geografi = 
                                                  InputData)$EigenValue,2)))
  
  ggmap(aaa)+
    xlab("Lengdegrad")+ ylab("Breddegrad")+theme_cowplot()+
    
    geom_point(data= GEOS, 
               mapping = aes(x = PosX, y = PosY,
                             size = EV, col = EV))+
    geom_point(data = as.data.frame(tail(GEOS,length(which(Results[1,]=="X")))),
               mapping = aes(x = PosX, y = PosY), col = "red")+
    geom_label(data = Labs, 
               mapping = aes(y = Y, x = X,
                             label = MC),
               fontface = "bold", 
               fill = alpha("white",0.5), 
               hjust = 0.05, vjust = -.7)+
    th
  list(Results,GEOS)
}

aa = Plotting(VOGeo4, AntallNye = 3)


### DRAGEHODE #####
#install.packages("rgbif")
library(rgbif)
#library(scrubr)
library(maps)
library(data.table)
library(sp)
library(raster)
library(ggplot2)
library(cowplot)

myspecies <- c("Meligethes norvegicus","Dracocephalum ruyschiana")
gbif_data <- occ_data(scientificName = myspecies, 
                      hasCoordinate = TRUE,
                      country = "NO", 
                      limit = 20000)

Host = data.table( gbif_data$`Dracocephalum ruyschiana`$data)
Beetle = data.table(gbif_data$`Meligethes norvegicus`$data)

Joint = rbind(Beetle, Host, fill = T)

Joint = cbind(data.table(coordinates(spTransform(
  SpatialPoints(coords = Joint[,c("decimalLongitude", "decimalLatitude")],
                proj4string = CRS("+init=epsg:4326")),
  CRS("+init=epsg:32632"))))
  , Joint)
names(Joint)[1:2]<-c("UTMx", "UTMy")
DispDistanceM = 1500

# Create clusters, aka. populations
ssp = Joint[,c(1:2)]
chc <- hclust(dist(ssp), method="complete")

# Distance with a 1500 threshold  
chc.d <- cutree(chc, h=DispDistanceM) 
Joint[,Populasjon:=chc.d]

# Lage et grid (kan brukes til å definere populasjon)
Joint[,c("gridY", "gridX"):=list(cut(UTMy, seq(from = min(Joint$UTMy, na.rm = T),
                                               to  = max(Joint$UTMy, na.rm = T),
                                               by = DispDistanceM)),
                                 cut(UTMx, seq(from = min(Joint$UTMx, na.rm = T),
                                               to  = min(Joint$UTMx, na.rm = T),
                                               by = DispDistanceM)))]
# Finne sentrum hvor hver populasjon
Joint[,c("centrY", "centrX"):=list(mean(UTMy, na.rm = T),
                                   mean(UTMx, na.rm = T)),
      c("Populasjon")]

Joint[,nHost := (max(.SD[species=="Dracocephalum ruyschiana"]$individualCount, na.rm = T)),
      c("centrY", "centrX")]

Joint = Joint[decimalLatitude>0 & decimalLongitude>10 & 
                decimalLongitude<11]
setorder(Joint, species)

ggplot(data = Joint[nHost>0],
       aes(y = centrY, x = centrX, size = nHost,  col = species))+
  geom_point()+
  labs(xlab="UTM (km)")

Joint[,BeetlePresent:=uniqueN(species)-1, c("centrX", "centrY")]

tmp = (Joint[nHost>0,.SD[which.max(nHost)], 
             c("species", "centrX", "centrY", "nHost")
][,c("centrX", "centrY", "nHost","BeetlePresent")]
)

# Erstatter rene tilstedeværelse 
#(dvs nHost = 1, med det laveste observerte ikke-1-verdi)

tmp[,nHost:=ifelse(nHost==1,min(.SD[nHost>1]$nHost), nHost)]
tmp[,3]<-log(tmp[,3])
table(tmp$BeetlePresent)

Centrality = function(Geografi, DispersalDistance = 1.5, Alpha = NULL, unitAreal = "km2", unitDistance = "m"){
  # Geografi = dataramme med areal, x og y koordinate
  if(ncol(Geografi)>2){
    Areal = (Geografi[,3])
  }else{
    Geografi = cbind(Geografi,1)
    Areal = (Geografi[,3])
  }
  
  if(unitAreal=="m2"){
    Areal = Areal/1000000
  }
  
  dij = as.matrix(dist(cbind(Geografi[,1], Geografi[,2]), diag = T, upper = T))
  
  if(unitDistance == "m"){
    dij = dij/1000# transform from meter to kilometer
    
  }
  Alpha = 1/DispersalDistance
  
  Si = sapply(1:nrow(Geografi), function(i){
    Geo = sum(exp(-Alpha*dij[i,-i])*Geografi[-i,3]) 
  })
  
  M = sapply(seq_along(Areal), function(i){
    sapply(seq_along(Areal), function(j){
      exp(-Alpha*dij[i,j])*Areal[i]*Areal[j]
    })
  })
  
  diag(M)<-0
  
  ev = eigen(M)$vectors[,1]^2
  ev = ev/(sum(ev))
  we = ev*eigen(M)$value[1]
  
  
  list(Centrality = Si,
       EigenValue = eigen(M)$value[1],
       Eigenvectors = we)
  
}

library(ggmap)
library(mapproj)
library(cowplot)
library(splancs)

Plotting = function(Geografi = tmp[,1:3], Nye = tmp[,4], AntallNye = 1, Iterer = T){
  
  InputData = as.data.frame(cbind(Geografi, Ny = Nye))
  InputData$Id = 1:nrow(InputData)
  bb = bbox(SpatialPoints(coords = InputData[,2:3], 
                          proj4string = CRS("+init=epsg:32632")))
  
  Results = rbindlist(lapply(AntallNye, function(i){
    aa = combn(1:nrow(InputData[InputData$Ny==1,]),m = i)
    
    Results = rbindlist(lapply(1:ncol(aa), function(ii) {
      print(paste(ii,ncol(aa), sep = " av "))
      
      InputData2 = rbind(InputData[InputData$Ny==0,],
                         InputData[InputData$Ny==1,][aa[,ii],])
      
      #VOs = data.frame(coordinates(
      #  spTransform(SpatialPoints(coords = InputData2[,2:3],
      #                            proj4string = CRS("+init=epsg:32632")),
      #              CRS("+init=epsg:4326"))))
      #VOs$EV = log(Centrality(Geografi = InputData2)$Eigenvectors)
      
      Results = data.frame(t(aa[,ii]),
                           Centrality = round(Centrality(Geografi = 
                                                           InputData2)$EigenValue,2))
      names(Results)[1:length(unique(aa[,ii]))]<-unique(aa[,ii])
      Results[1:length(unique(aa[,ii]))]<-"X"
      data.table(Results)
    }), fill = T, use.names = T)
  }), use.names = T)
  
  Results[,nNye:=apply(as.matrix(Results),1, function(x){length(unique(which(x=="X")))})]
  
  
  NullModell = round(Centrality(Geografi = InputData[InputData[,4]==0,1:3])$EigenValue,2)
  
  Results = cbind(Results[,-"Centrality"]
                  ,Results[,"Centrality"])
  
  setorder(Results, -Centrality)
  Results$Centrality<- Results$Centrality-NullModell
  Results
  
  # hvilke nye områder er inkludert?
  as.numeric(names(Results)[which(Results[1,]=="X")])
  
  ## Plot kartet
  GEOS0 = rbind(InputData[InputData$Ny==0,],
                InputData[InputData$Ny==1,
                ][as.numeric(names(Results
                )[which(Results[1,]=="X")]),])
  
  GEOS = data.frame(coordinates(
    spTransform(SpatialPoints(coords = GEOS0[,c("centrX", "centrY")],
                              proj4string = CRS("+init=epsg:32632")),
                CRS("+init=epsg:4326"))))
  
  GEOS$EV = log(Centrality(Geografi = GEOS0)$Eigenvectors)
  names(GEOS)[1:2]<-c("PosX", "PosY")
  
  GEOS$Bille<-GEOS0[,4]
  # Alle lokaliteter
  aaa = get_map(c(left = min(GEOS[,1])-.2, 
                  right = max(GEOS[,1])+.2,
                  bottom = min(GEOS[,2])-.2, 
                  top = max(GEOS[,2])+.2))
  
  # Området akkurat rundt den nye lokaliteten
  aaaSub = get_map(c(left = tail(GEOS,1)[,1]-.3, 
                     right = tail(GEOS,1)[,1]+.3,
                     bottom = tail(GEOS,1)[,2]-.3, 
                     top = tail(GEOS,1)[,2]+.3))
  
  #ggmap(aaaSub)
  
  Labs = data.frame(Y = max(GEOS$PosY),#as.numeric(attr(aaa, "bb")[3]),
                    X = min(GEOS$PosX),#as.numeric(attr(aaa, "bb")[2]), 
                    MC = paste("Metapopulasjonskapasitet = ",
                               round(Centrality(Geografi = 
                                                  InputData)$EigenValue,2)))
  FinPlot = ggmap(aaaSub)+
    xlab("Lengdegrad")+ ylab("Breddegrad")+theme_cowplot()+
    
    geom_point(data= GEOS, 
               mapping = aes(x = PosX, y = PosY,
                             size = EV, col = EV))+
    geom_point(data = as.data.frame(tail(GEOS,length(which(Results[1,]=="X")))),
               mapping = aes(x = PosX, y = PosY), col = "red")+
    geom_label(data = Labs, 
               mapping = aes(y = Y, x = X,
                             label = MC),
               fontface = "bold", 
               fill = alpha("white",0.5), 
               hjust = 0.05, vjust = -.7)+
    th
  
  list(Inputs = GEOS, Ranking = Results, Plot = FinPlot)
}

Results = Plotting(AntallNye = 1)
Results[[2]]

# Man har funnet at centrality-målet korresponderer bra 
# med sannsynligheten for forekomst. 
names(Results)[which(Results[,1]=="X")]

Results$Plot

ggmap(aaa)+geom_path(data = GEOS[c((combn(1:42,2))),], mapping = aes(x = PosX, y = PosY))


