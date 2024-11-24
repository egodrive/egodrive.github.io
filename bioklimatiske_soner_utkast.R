---
title: "Bioklimatiske Soner"
description: |
  Bioklimatiske soner er tøffe!
author:
  - name: Endre Grüner Ofstad
url: {}
date: 2024-03-29
output:
  distill::distill_article:
  self_contained: false

---
  
  Her kommer det en tekst om bioklimatiske soner. 
```{r Grunnlagskode for å hente og sammenstille data, eval = F, echo = F, cache = T}
########################
#devtools::install_github("mdsumner/distancetocoast")
#devtools::install_github("ropensci/rnaturalearthhires")
#### OTHER ####
# https://chelsa-climate.org/bioclim/

### Forests ####
# Skogkart
# https://forobs.jrc.ec.europa.eu/GFC
# tiles = list.files("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/", "JRC_GFC",
#                    full.names = T)
# 
# ext <- extent(4, 34, 57, 72)
# forests = lapply(tiles, function(x){
#   rast(x)
# })

#forests_m = do.call("sprc", forests)
#forests_m  = merge(forests_m)
#forests_m = do.call("merge", forests)
#image(forests_m[[1]])
#### Klima variabler ####
#"https://geodata.ucdavis.edu/climate/worldclim/2_1/tiles/tile/"
#ff = list.files("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/", "tif", full.names = F)
#ff = sapply(ff, function(x){
#  paste(c(unlist(strsplit(x,"_"))[1:4], ""), collapse = "_")
#})

tiles = list.files("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/", "tile_",
                   full.names = F)
imgs = data.table(file = list.files("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/", "tile_",
                                    full.names = T),
                  var = sapply(tiles, function(x){
                    gsub(".tif","",unlist(strsplit(x,"_"))[5], fixed = T)
                  }),
                  tiles = sapply(tiles, function(x){
                    gsub(".tif","",unlist(strsplit(x,"_"))[2], fixed = T)
                  }))
imgs = imgs[var!="srad"]
imgs$file

res = lapply(c(7,8,19,18), function(x){
  print(x)
  idx = which(imgs$tiles==x)
  
  # https://archive.org/details/methodsinclimato032726mbp/page/n217/mode/2up?q=humid
  t2 = diff(range(rast(imgs$file[which(imgs$tiles==x & imgs$var == "tavg")])))
  Theta = rast(matrix(crds(t2, na.rm = F)[,2], byrow = T, ncol = ncol(t2), nrow = nrow(t2)))
  
  Theta = sin((Theta*(pi/180))+(10*(pi/180)))
  ext(Theta) <- ext(t2);crs(Theta) <- crs(t2)
  t2b = (t2*1.7)
  t3c = (t2b/Theta)#-14
  names(t3c)<-"Conrad"
  
  zones = lapply(idx, function(x){
    rast(imgs$file[x])
  })
  names(zones)<-imgs$var[idx]
  
  zones[length(zones)+1]<-t3c
  #rm(t3c, t2b, t2, Theta)
  
  # merge all layers
  zone = do.call("c", zones)
  names(zone) = gsub("tile_7_wc2.1_30s_","",names(zone))
  names(zone)
  image(zone[[1]])
  #zone = crop(zone, extent(4, 34, 57, 72))
  
  zone_sc = do.call("c", zones)
  rm(zones, r2, r2b, r2c)
  gc()
  zone_sc
}
)
imgsMerged = do.call("merge", res)
writeRaster(imgsMerged, filename = "C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/Fennoscandia.tif")

```

https://hastie.su.domains/Papers/gap.pdf

```{r setup, cache = T, echo = F}
library(sp);library(raster);library(data.table);library(rgdal);library(distancetocoast);library(raster);library(rnaturalearth);library(rnaturalearthhires);library(terra);library(data.table);library(sp);library(factoextra);library(moments);library(geodata)
library(factoextra);library(cluster);library(rgbif);library(maps);library(gridExtra);library(tidyterra)
knitr::opts_chunk$set(fig.width=12, fig.height=8) 

```




```{r Sampler data, cache = T, echo = F, results = "hide"}

imgsMerged = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/Fennoscandia.tif")

#### SAMPLE POINTS####
countries <- world(resolution = 5, path = "maps")  # you may choose a smaller (more detailed) resolution for the polygon borders, and a different folder path to save the imported map
cntry_codes <- country_codes()
# add this table to the countries map attributes:
countries <- merge(countries, cntry_codes, by.x = "GID_0", by.y = "ISO3", all.x = TRUE)

#dd1 = countries[countries$UNREGION1 %in% c("Northern Europe","Western Europe") & !(countries$NAME_0 %in% c("Iceland","Svalbard and Jan Mayen")),]
dd1 = countries[countries$NAME_0 %in% c("Norway", "Sweden", "Denmark"),]
tmp = terra::crop(imgsMerged[[1]],dd1, mask = T)

set.seed(10)
dd = na.omit(spatSample(tmp,20000000, xy = T, method = "regular"))[,1:2]
dd$ID = 1:nrow(dd)
sampleData = extract(imgsMerged,dd[,c("x","y")])

Temp_lavest = apply(sampleData[,grepl("tavg", colnames(sampleData))], 1, min)
sampleData = cbind(sampleData, Temp_lavest)
rm(Temp_lavest)

## OAK ####
# 
# fileT = "C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/BroadLeavedTrees_Norway.tif"
# 
# if(!file.exists(fileT)){
#   myspecies <- c("Quercus robur","Quercus petraea","Fagus sylvatica","Fraxinus excelsior","Ulmus glabra")
#   # download GBIF occurrence data for this species; this takes time if there are many data points!
#   gbif_data <- occ_data(scientificName = myspecies, country = "NO",
#                         hasCoordinate = TRUE, limit = 20000)
#   ocdata = rbindlist(lapply(gbif_data, function(x) {data.table(x$data)}), fill = T)
#   tmp = ocdata[,c("genusKey","scientificName", "decimalLatitude", "decimalLongitude"),with = F]
#   # Round of the coordinates so they to a large extent reflect different individuals/populations
#   tmp$decimalLatitude = round(tmp$decimalLatitude,4)
#   tmp$decimalLongitude = round(tmp$decimalLongitude,4)
#   tmp = unique(tmp)
#     
#   r = r2 = (imgsMerged[[1]])
#   rFact = 2
#   r2 = disagg(r2, fact = rFact)# 2-double the resolution
#   r2 = raster(r2)
#   r2[]<-0
#   # get the cell index for each point and make a table:
#   counts = table(cellFromXY(r,tmp[,c("decimalLongitude","decimalLatitude"),with =F]))
#   # fill in the raster with the counts from the cell index:
#   r2[as.numeric(names(counts))] = ifelse(counts==0, 0, log(counts))
#   r2 = aggregate(r2, fact = rFact)# 2-double the resolution
#   # Write a copy so we don't have to run this part every time
#   writeRaster(r2, filename = fileT)
#   
# }
# r2 = rast(fileT)
# oaks = extract(r2,dd[,c("x","y")])
# 
# names(oaks)[2]="mOaks"
# sampleData = cbind(sampleData, Broadleaved = (oaks$mOaks))

# COAST #####
# Distance to coast, calculate distance to neareas NA-cell i.e. the coast
D2C = extract(rast(distance_to_coastline_10),dd[,c("x","y")])
sampleData = cbind(sampleData, D2C = log(D2C$layer))
rm(D2C)
### Evapotranspiration ####
# # https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/5

PET = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/et0_v3_yr.tif")
PET = extract(PET,dd[,c("x","y")])
sampleData = cbind(sampleData, PET = PET[,2])

lf = list.files("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/Global-ET0_v3_monthly", ".tif", full.names = T)
library(stringr)
lf = lf[str_count(lf, "\\.")==1]
pets = sapply(lf, function(x){
  pt = rast(x)
  pet = extract(pt,dd[,c("x","y")])
  pet[,2]
})
colnames(pets)<-paste0("PET_",str_sub(colnames(pets),-6,-5))

PETsd = apply(pets, 1, function(x) sd(x, na.rm = T))*100

sampleData = cbind(sampleData, pets, PETsd)
rm(PET);rm(pets);rm(PETsd)
### Forests ####
#forests = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/GEOCARBON_AGB_Map.tif")
#forests = extract(forests,dd[,c("x","y")])
#sampleData = cbind(sampleData, forestBiomass = forests[,2])

#### Growing season and snow cover ####
GS = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/CHELSA_gsl_1981-2010_V.2.1.tif")
GS = extract(GS,dd[,c("x","y")])

GSP = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/CHELSA_gsp_1981-2010_V.2.1.tif")
GSP = extract(GSP,dd[,c("x","y")])
GSP$`CHELSA_gsp_1981-2010_V.2.1` = ifelse(GSP$`CHELSA_gsp_1981-2010_V.2.1`>10000,0,GSP$`CHELSA_gsp_1981-2010_V.2.1`)

SC = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/CHELSA_scd_1981-2010_V.2.1.tif")
SC = extract(SC,dd[,c("x","y")])

GDD5 = rast("C:/Users/endre/OneDrive/DiverseProsjekt/Bioklimatiske soner/CHELSA_gdd5_1981-2010_V.2.1.tif")
GDD5 = extract(GDD5,dd[,c("x","y")])

sampleData = cbind(sampleData, GS = GS[,2], GSP = GSP[,2], SC = SC[,2], GDD5 = GDD5[,2])
rm(GSP);rm(GDD5);rm(SC);rm(GS);gc()
colnames(sampleData) = gsub("tile_7_wc2.1_30s_","",
                            colnames(sampleData))
# remove what skew we can
sampleData[,"elev"] = ifelse(is.na(sampleData[,"elev"])|sampleData[,"elev"]<0.1,1,sampleData[,"elev"])
sampleData[,"elev"] = log(sampleData[,"elev"])

## Normalize to (0,1) #####
sampleData = apply(sampleData, 2, function(x){
  vals = x
  vals = scale(vals)
  y2 = (vals-min(vals,na.rm = T))/(max(vals, na.rm = T)-min(vals,na.rm = T))
  y2
}
)

# add id column
sampleData[,"ID"] <- as.integer(dd$ID)

#sampleData[,"Broadleaved"] = ifelse(is.nan(sampleData[,"Broadleaved"]),0,sampleData[,"Broadleaved"])

sampleData = sampleData[,colnames(sampleData)!="GSP"] # minst variasjon
sampleData = sampleData[,colnames(sampleData)!="GS"] # minst variasjon


#sampleData2 = sampleData[,-1]

# Fjerne de punktene som ligger helt på kystlinja
#ii = which(colnames(sampleData)=="D2C")
#sampleData = sampleData[sampleData[,ii]>0,]
summary(sampleData)


sampleData = na.omit(sampleData)
ids2 = sampleData[,"ID"]
sampleData2 = sampleData[,-1]
# rm(sampleData)
# Test
# dd2 = data.table(ID = ids2, 
#            Conrad = sampleData2[,1])[data.table(dd),
#                                      on  = "ID"]
# dd_sp = vect(x = dd2, geom = c("x","y"))
# 
# r0 = imgsMerged[[1]]
# r0 = crop(r0, dd_sp)
# r0[]<-0
# TestRaster = terra::rasterize(x = dd_sp, r0, field = "Conrad")
# image(TestRaster)

```
# Prepper figurer

```{r antall klynger}



```



```{r Figurer output, cache = T, echo = F}
library(ggplot2);library(terra)
###  PCA ####
pca.data = prcomp(sampleData2)

load = as.matrix(pca.data$rotation[,1:2])

##fviz_pca_var(pca.data, addEllipses = T)

### rotation test ####

# Try rotate axis 1 according to Elev
# https://rstudio-pubs-static.s3.amazonaws.com/742270_9e031b38243a43b49b316a9a9bd65346.html
#xx = load[row.names(load)=="elev",]
#xx = atan(xx[2]/xx[1])/pi*180

psi <- -50*pi/180
PSI_x <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), byrow = T, nrow = 2)

# Try rotate axis to according to Tavg 12,1,2
#xx = load[row.names(load)=="elev",]
#xx = atan(xx[2]/xx[1])/pi*180

psi <- -52*pi/180
PSI_y <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), byrow = T, nrow = 2)

(rotDF <- cbind(load %*% PSI_x[,1],
                load %*% PSI_y[,2]))
load2 = data.frame(load)
load2$var = row.names(load)

rotDF = data.frame(rotDF)
rotDF$var = row.names(rotDF)



NewData = as.matrix(sampleData2)%*%as.matrix(rotDF[,1:2])

```
```{r env space, cache = T}
tmp = abs(rotDF[,1:2])
tmp = tmp/apply(tmp, 2, max)
tmp = cbind(tmp,idx = apply(tmp, 1, function(x) which.max(x)), nidx = apply(tmp, 1, function(x) length(which(x[1:2]>.5))))

# Fjerner de som loader i begge dimensjoner ("skaper støy")
tmp = tmp[!tmp$nidx>1,]

tmp = tmp[apply(tmp[,1:2], 1, function(x) any(x>.5)),]
tmp$x2Ord = order(tmp$X2)
tmp$x1Ord = order(tmp$X1)
rbind(tmp[tmp$X1<1/3 & tmp$X2>1/3,],
      tmp[tmp$X2<1/3 & tmp$X1>1/3,])

Soner = row.names(tmp[tmp$X1<1/3 & tmp$X2>1/3,])
Seksjoner = row.names(tmp[tmp$X2<1/3 & tmp$X1>1/3,])

Orig = ggplot(data = load2, aes(x = PC1, y = PC2, label = var))+
  geom_point(col = "grey30")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_text()

rotDF$Gruppering = with(rotDF, ifelse(var %in% Soner, "Sone",
                                      ifelse(var %in% Seksjoner, "Seksjon", "Annet")))

New = ggplot(data = rotDF, aes(x = X1, y = X2, label = var,col =  Gruppering))+
  geom_point(col = "grey30")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_text()+ 
  theme(legend.position="none")

grid.arrange(Orig, New)


fviz_screeplot(pca.data, addlabels = TRUE, ylim = c(0, 50))

```

# Results

```{r clustering, echo = T, cache = T}


#### clustering #####

#calculate gap statistic for each number of clusters (up to 25 clusters)
#gapstats = "C:/Users/endre/OneDrive/Tekster/egodrive.github.io/_posts/2024-03-29-bioklimatiske-soner/hclust_gap_stats.rds"
gapstats = "C:/Users/endre/OneDrive/Tekster/egodrive.github.io/_posts/2024-03-29-bioklimatiske-soner/kmeans_gap_stats.rds"
set.seed(10);idx = sample(1:nrow(sampleData2),30000)

if(!file.exists(gapstats)){
  
  gap_stat_k_soner <- clusGap(sampleData2[idx,Soner], FUN = kmeans, K.max = 20, B = 50)
  gap_stat_k_seksjoner <- clusGap(sampleData2[idx,Seksjoner], FUN = kmeans, K.max = 20, B = 50)
  
  saveRDS(list(gap_stat_k_soner = gap_stat_k_soner,
               gap_stat_k_seksjoner = gap_stat_k_seksjoner)
          ,"C:/Users/endre/OneDrive/Tekster/egodrive.github.io/_posts/2024-03-29-bioklimatiske-soner/kmeans_gap_stats.rds") 
  #gap_stat_hc_soner <- clusGap(sampleData2[idx,Soner], FUN = hcut, K.max = 23, B = 40)
  #gap_stat_hc_seksjoner <- clusGap(sampleData2[idx,Seksjoner], FUN = hcut, K.max = 23, B = 40)
  
  
  # saveRDS(list(gap_stat_hc_soner = gap_stat_hc_soner,
  #              gap_stat_hc_seksjoner = gap_stat_hc_seksjoner,
  #              gap_stat_k_soner = gap_stat_k_soner,
  #              gap_stat_k_seksjoner = gap_stat_k_seksjoner)
  #         ,gapstats) 
}
lf = readRDS(gapstats)

gap_stat_k_soner = lf$gap_stat_k_soner
gap_stat_k_seksjoner = lf$gap_stat_k_seksjoner
#tmp = Seksjoner; Seksjoner = Soner; Soner = tmp

library(factoextra)
k1 = fviz_gap_stat(gap_stat_k_soner, maxSE = list(method = "Tibs2001SEmax", SE.factor = 1));k1 # 10 når man bruker soner >.5
k2 = fviz_gap_stat(gap_stat_k_seksjoner, maxSE = list(method = "Tibs2001SEmax", SE.factor = 1));k2 # 14 når man bruker seksjoner >0.5
nKsoner_0 = (which(with(k1$data,(gap+SE.sim)-(c(gap[-1],NA)-SE.sim))>0))
nKsoner_0 = nKsoner_0[nKsoner_0>1]
nKsoner = min(nKsoner_0[nKsoner_0>1])
nKseksjoner_0 = which(with(k2$data, (gap+SE.sim)-(c(gap[-1],NA)-SE.sim))>0)
nKseksjoner_0 = nKseksjoner_0[nKseksjoner_0>1]
nKseksjoner = min(nKseksjoner_0[nKseksjoner_0>1])
# h1 = fviz_gap_stat(gap_stat_hc_soner, maxSE = list(method = "Tibs2001SEmax", SE.factor = 1));h1
# h2 = fviz_gap_stat(gap_stat_hc_seksjoner, maxSE = list(method = "Tibs2001SEmax", SE.factor = 1));k2
# nHsoner = min(which(with(h1$data,(gap+SE.sim)-(c(gap[-1],NA)-SE.sim))>0))
# nHseksjoner = min(which(with(h2$data, (gap+SE.sim)-(c(gap[-1],NA)-SE.sim))>0))



clusters = "C:/Users/endre/OneDrive/Tekster/egodrive.github.io/_posts/2024-03-29-bioklimatiske-soner/clusters.rds"
if(!file.exists(clusters)){
  
  HC_soner = hcut(sampleData2[idx,Soner], nKsoner)
  K_soner = kmeans(sampleData2[,Soner], iter.max = 100, nstart = 5)
  Kseksjoner = kmeans(sampleData2[,Seksjoner], nKseksjoner)
  HCseksjoner = hcut(sampleData2[idx,Seksjoner], nKseksjoner)
  
  saveRDS(list(HC_soner = HC_soner,
               K_soner = K_soner,
               Kseksjoner = Kseksjoner,
               HCseksjoner = HCseksjoner
               )
          ,clusters) 
}

# Soner
ls = readRDS(clusters)
HC_soner = ls$HC_soner
K_soner = ls$K_soner
Kseksjoner = ls$Kseksjoner
HCseksjoner = ls$HCseksjoner


# Merge the clusters back on the spatial data
tt = data.table(cbind(ID = ids2,
                      Sone = K_soner$cluster),
                x = NewData[,1], y = NewData[,2])
tmp = data.table(ID = ids2[idx]
                 , hClust = cutree(HC_soner, k = nKsoner)
                 , hClust6 = cutree(HC_soner, k = 6)
)
tt = tmp[tt, on = "ID"]

# hClust ble bare gjort på et subset pga kapasitetsbegrensningser
# gjør det derfor slik at per sone, estimert fra kmeans, så blir andre punkt i den sone
# lagt til den klyngen som forekommer mest i den kmeans-sonen
tt[,hClust:=names(which.max(table(.SD$hClust))),"Sone"]
tt[,hClust6:=names(which.max(table(.SD$hClust6))),"Sone"]
names(tt)<-c("ID", "hClust", "hClust6", "Sone", "x", "y")
BK_soner_data = tt

tmp = BK_soner_data[,median(y), "Sone"]
tmp = arrange(tmp, -V1);tmp[,Grp:=1:.N]
BK_soner_data = tmp[BK_soner_data, on = "Sone"]
setnames(BK_soner_data, c("Sone", "Grp"), c("Grp", "Sone"))
BK_soner_data$Grp<-NULL
setnames(BK_soner_data, c("x", "y"), c("PCAx", "PCAy"))

tmp = BK_soner_data[,median(PCAy), "hClust"]
tmp = arrange(tmp, -V1);tmp[,Grp:=1:.N]
BK_soner_data = tmp[BK_soner_data, on = "hClust"]
setnames(BK_soner_data, c("hClust", "Grp"), c("Grp", "hClust"))
BK_soner_data$V1<-NULL;BK_soner_data$i.V1<-NULL
tmp = BK_soner_data[,median(PCAy), "hClust6"]
tmp = arrange(tmp, -V1);tmp[,Grp:=1:.N]
BK_soner_data = tmp[BK_soner_data, on = "hClust6"]
setnames(BK_soner_data, c("hClust6", "Grp"), c("Grp", "hClust6"))
BK_soner_data$V1<-NULL;BK_soner_data$i.V1<-NULL

dd_1 = data.frame(dd)
dd2 = BK_soner_data[dd_1, on = "ID"]

NewData = data.table(NewData)
NewData[,ID:=1:.N]
gc()

# Seksjoner

# Merge the clusters back on the spatial data
tt = data.table(ID = ids2,Seksjon = Kseksjoner$cluster, 
                x = NewData[,1], y = NewData[,2])
names(tt)<-c("ID", "Seksjon", "x", "y")
tmp = data.table(ID = ids2[idx]
                 , Seksjon_hClust = cutree(HCseksjoner, k = nKseksjoner)
                 , Seksjon_hClust6 = cutree(HCseksjoner, k = 6)
)
tt = tmp[tt, on = "ID"]

# hClust ble bare gjort på et subset pga kapasitetsbegrensningser
# gjør det derfor slik at per seksjon, estimert fra kmeans, så blir andre punkt i den seksjonen
# lagt til den klyngen som forekommer mest i den kmeans-seksjonen
tt[,Seksjon_hClust:=names(which.max(table(.SD$Seksjon_hClust))),"Seksjon"]
tt[,Seksjon_hClust6:=names(which.max(table(.SD$Seksjon_hClust6))),"Seksjon"]
BK_seksjoner_data = tt

tmp = BK_seksjoner_data[,median(x), "Seksjon"]
tmp = arrange(tmp, -V1);tmp[,Grp:=1:.N]
BK_seksjoner_data = tmp[BK_seksjoner_data, on = "Seksjon"]
setnames(BK_seksjoner_data, c("Seksjon", "Grp"), c("Grp", "Seksjon"))
BK_seksjoner_data$Grp<-NULL;BK_seksjoner_data$V1<-NULL
setnames(BK_seksjoner_data, c("x", "y"), c("PCAx", "PCAy"))

tmp = BK_seksjoner_data[,median(PCAx), "Seksjon_hClust"]
tmp = arrange(tmp, -V1);tmp[,Grp:=1:.N]
BK_seksjoner_data = tmp[BK_seksjoner_data, on = "Seksjon_hClust"]
setnames(BK_seksjoner_data, c("Seksjon_hClust", "Grp"), c("Grp", "Seksjon_hClust"))
BK_seksjoner_data$Grp<-NULL;BK_seksjoner_data$V1<-NULL
tmp = BK_seksjoner_data[,median(PCAx), "Seksjon_hClust6"]
tmp = arrange(tmp, -V1);tmp[,Grp:=1:.N]
BK_seksjoner_data = tmp[BK_seksjoner_data, on = "Seksjon_hClust6"]
setnames(BK_seksjoner_data, c("Seksjon_hClust6", "Grp"), c("Grp", "Seksjon_hClust6"))
BK_seksjoner_data$Grp<-NULL;BK_seksjoner_data$V1<-NULL

dd2 = BK_seksjoner_data[,c("ID", "Seksjon","Seksjon_hClust","Seksjon_hClust6")][dd2, on = "ID"]


dd_sp = vect(x = dd2, geom = c("x","y"))
r0 = imgsMerged[[1]]
r0 = crop(r0, dd_sp)
r0[]<-0
SoneRaster = terra::rasterize(x = dd_sp, r0, field = "Sone")
#SoneRaster_hc = terra::rasterize(x = dd_sp, r0, field = "hClust")
#SoneRaster_hc6 = terra::rasterize(x = dd_sp, r0, field = "hClust6")

SeksjonRaster = terra::rasterize(x = dd_sp, r0, field = "Seksjon")
#SeksjonRaster_hc = terra::rasterize(x = dd_sp, r0, field = "Seksjon_hClust")
#SeksjonRaster_hc6 = terra::rasterize(x = dd_sp, r0, field = "Seksjon_hClust6")

```

```{r flere figurer, echo = T, cache = T}
idx = sample(1:nrow(BK_seksjoner_data), 30000)

# Fra: #538FCB
# Via: #81C343
# Via: #F2EC46
# Til: #F064A4
SoneFarger = colorRampPalette(colors = c("#538FCB", "#81C343", "#F2EC46", "#F064A4"))

BK_soner_data[,Sone:=as.factor(Sone)]
BK_soner = qplot(x = PCAx, y = PCAy, col = Sone, 
                 data = BK_soner_data[idx,])+
  discrete_scale("col", scale_name = "personal", palette = SoneFarger)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)

# Fra:  #3853A3
# Via: #98D7DC
# Til: #E7F6F8
SeksjonsFarger = colorRampPalette(colors = c("#E7F6F8","#98D7DC","#3853A3"))
BK_seksjoner_data[,Seksjon:=as.factor(Seksjon)]
BK_seksjoner = qplot(x = PCAx, y = PCAy, col = Seksjon,
                     data = BK_seksjoner_data[idx,])+
  discrete_scale("col", scale_name = "personal", palette = SeksjonsFarger)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)

grid.arrange(BK_soner, BK_seksjoner)

sc <- scale_fill_gradientn(colours = SoneFarger(100),
                           limits=c(1, nKsoner))
Sone = ggplot()+
  geom_spatraster(data = SoneRaster, aes(fill = last))+
  sc+ theme(legend.position="none")

sc <- scale_fill_gradientn(colours = SeksjonsFarger(100), limits=c(1, nKseksjoner))
Seksjon = ggplot()+
  geom_spatraster(data = SeksjonRaster, aes(fill = last))+
  sc+ theme(legend.position="none")

grid.arrange(Sone+
               theme(plot.margin= unit(c(0, 0, 0, 0), "lines"),
                     axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)),
             Seksjon+theme(axis.text.y=element_blank()),
             ncol = 2)

```

stats = BK_soner_data[,c("ID", "Sone")][data.table(ID = ids2, sampleData2[,Soner]), on = "ID"]
stats = melt(stats, id.vars = c("ID", "Sone"))
stats = stats[variable %in% c("SC", "GDD5", "elev", "Conrad", "bio17", "PET")]
stats$Sone = as.integer(stats$Sone)
stats[,variable:=as.factor(variable)]
setorder(stats,-variable)

a = ggplot(data = stats[,.SD[sample(1:.N,10)], c("Sone", "variable")], 
       aes(x = value, y = Sone,  group = interaction(Sone, variable)
           , fill = Sone))+
  geom_boxplot(position = position_dodge())+
  scale_fill_gradientn(colours = SoneFarger(100), limits=c(1, nKsoner))+
  facet_wrap(variable~.)

library(ggpattern)
# How do the zones cluster
tmp0 = sapply(nKsoner:1, function(x){cutree(HC_soner, x)})
tmp0 = unique(tmp0)
row.names(tmp0)<-tmp0[,1]
tmp = hcut(tmp0)
tmp$height<- tail(sort(HC_soner$height), nKsoner-1)
#plot(tmp)
library(ggdendro);library(plotly)
library(ggraph)
#<-letters[1:4]
tmp2 = as.dendrogram(tmp)
tmp2 = as.phylo(tmp)

tmp2$order<-c(6,5,4,3,2,1)

b = (ggdendrogram(tmp2, rotate = T, size = 2))
grid.arrange(a, b, ncol = 2)



```{r leaflet kart, echo = T, cache = T}
library(raster);library(leaflet)
so = raster(SoneRaster)
SonePal <- colorNumeric(c("#538FCB", "#81C343", "#F2EC46", "#F064A4"), values(so),
                    na.color = "transparent")

leaflet() %>% addTiles() %>%
  addRasterImage(so, colors = SonePal, opacity = 0.8, group = "Soner") %>%
  addRasterImage(se, colors = SeksjonsPal, opacity = 0.8, group = "Seksjoner")%>%
addLayersControl(
  overlayGroups = c("Soner", "Seksjoner"),
  options = layersControlOptions(collapsed = FALSE)
)


%>%
  addLegend(pal = SonePal, values = values(so),
            title = "Bioklimatisk sone")

se = raster(SeksjonRaster)
SeksjonsPal = colorNumeric(c("#E7F6F8","#98D7DC","#3853A3"), values(se),
                           na.color = "transparent")

a <- leaflet() %>%
  # Base groups
  addTiles(group = "OSM (default)")

a %>%
  addProviderTiles(providers$Stadia.StamenToner, group = "Toner") %>%
  addProviderTiles(providers$Stadia.StamenTonerLite, group = "Toner Lite") %>%
  # Overlay groups
  addRasterImage(so, colors = SonePal, opacity = 0.8, group = "Soner") %>%
  addRasterImage(se, colors = SeksjonsPal, opacity = 0.8, group = "Seksjoner") %>%
  # Layers control
  addLayersControl(
    baseGroups = c("OSM (default)", "Toner", "Toner Lite"),
    overlayGroups = c("Soner", "Seksjoner"),
    options = layersControlOptions(collapsed = FALSE)
  )


```
Sys.time()
kk_Sone <- kmH(x = sampleData2[idx,Soner],maxclus = 15,
               desired.ncores = detectCores()-1,verbose = TRUE)
Sys.time()
saveRDS(kk_Sone, file="KMH_Soner.rds")
Sys.time()
kk_seksjon <- kmH(x = sampleData2[idx,Seksjoner],maxclus = 15,
                  desired.ncores = detectCores()-1,verbose = TRUE)
Sys.time()
saveRDS(kk_seksjon, file="KMH_Seksjoner.rds")
