library(data.table)
library(ggfortify)
library(readxl)
polls <- data.table(read_excel("C:/Users/endre/OneDrive/Dokumenter/Pollinator_popularitet_planter.xlsx",
                               col_types ="guess"))

polls = polls[,1:11]
names(polls)[2]<-"Native"
names(polls)[3]<-"Life_cycle"
Polls_num = apply(polls[,c(4:11)], c(1,2), as.numeric)
Polls_num = apply(Polls_num, 2, function(x){(x-min(x))/(max(x)-min(x))})
navn = sapply(as.list(strsplit(polls$`Plant cultivar`,'“')), function(x) x[[1]])
Gr2 = sapply(navn, function(x){
  if(x %in% RL_2021$`Vitenskapelig navn`){
    RL_2021$`Kategori 2021`[which(x==RL_2021$`Vitenskapelig navn`)]
  }else{
    if(x %in% FA$`Vitenskapelig navn`){
      FA$`Risikonivå 2018`[which(x==FA$`Vitenskapelig navn`)]
    }else{
      "-"
    }
  }
})
PopNavn = sapply(navn, function(x){
  if(x %in% RL_2021$`Vitenskapelig navn`){
    RL_2021$Populærnavn[which(x==RL_2021$`Vitenskapelig navn`)]
  }else{
    if(x %in% FA$`Vitenskapelig navn`){
      FA$Populærnavn[which(x==FA$`Vitenskapelig navn`)]
    }else{
      "-"
    }
  }
})


row.names(Polls_num)<-PopNavn

RL_2021 <- data.table(read_excel("C:/Users/endre/OneDrive/ArterLR/rødliste-2021.xlsx",
                                 col_types ="guess"))
RL_2021 = RL_2021[Vurderingsområde=="Norge"]
FA = fread("C:/Users/endre/OneDrive/ArterLR/Fremmedartslista_2018_karplanter_2.csv")


Polls_2 = cbind(Polls_num[,c(1:3,5:6)],
                Gruppe =paste(unlist( kmeans(Polls_num[,c(1:3,5:6)], 3)[1])),
                Gr2)

windows()
autoplot(prcomp(Polls_num[,c(1:3,5:6)]),
         data = Polls_2,
         colour = 'Gr2',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

windows()
autoplot(kmeans(Polls_num[,c(1:3,5:6)], 3),
         data = Polls_num[,c(1:3,5:6)], 
         #size = `Bumblebees/m2`,
         label = TRUE, label.size = 3)


