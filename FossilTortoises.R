### packages ####
library(dplyr)
library(ggplot2)
library(paleoTS)
library(stats)
library(moments)
library(pgirmess)
library(vegan)
#library(plotly)

#### Data basis - import data set from Excel ####
setwd("C:/Users/Jule/Documents/Uni/MA")
tidyCL<-read.csv("tortoises_tidy.csv", sep=";", header=TRUE)

## change to:
# setwd("C:/Users/Jule/Documents/Uni/MA/PublicationTortoises")
# tidyCL<-read.csv("tortoises_tidy_pub.csv", sep=";", header=TRUE)
# C:\Users\Jule\Documents\Uni\MA\PublicationTortoises



#####prepare for analysis (fix column names in .csv-file after converting table from excel####
colnames(tidyCL)[6] <- "MAmin"
colnames(tidyCL)[7] <- "Mamax"
colnames(tidyCL)[17] <- "CL"
colnames(tidyCL)[18] <- "PL"
colnames(tidyCL)[21] <- "estimated"

tidyCL <-  tidyCL %>%
  mutate(Age= ((MAmin)+(as.numeric(Mamax)))/2)  #%>%
# filter(estimated=="m" | estimated=="mo"| estimated=="mf") #%>%
# filter(!is.na(CL))

####### import extant data ####
extant <- read.csv("MFN_testudinidae.csv", sep=";", header=TRUE, dec=".", na.strings = "NA", stringsAsFactors=FALSE)  # file: MFN_testudinidae.csv

colnames(extant)[5] <- "SCL"
colnames(extant)[11] <- "PL"
colnames(extant)[12] <- "PLmid"


Extant <- extant %>%
  mutate(CL = SCL * 10, PL=PL*10, PLmid=PLmid*10) 

#### estimating CL from PL ####

#table shows CL/PL ratios for all fossil taxa that have both measurements available
CLPLtidy <- tidyCL %>%
  filter(!is.na(CL) & !is.na(PL)) %>%
  dplyr::select(Taxon, CL, PL, size, Age, Island, Continent, Genus) %>%
  mutate(ratio=CL/PL) %>%
  group_by(Taxon) %>% #to show ratios per Taxon, leave out to get a total ratio
  dplyr::summarise(meanRatio=round(mean(ratio),2), sdRatio=round(sd(ratio),2), n=n(), min=min(ratio), max=max(ratio))


#table shows CL/PL ratios for all extant taxa that have both measurements available
CLPLextant <- Extant %>%
  dplyr::filter(!is.na(CL) & !is.na(PL)) %>%
  dplyr::select(Taxon=Species, CL, PL, PLmid, Island, Continent, Genus) %>%
  mutate(ratio=CL/PL, ratioMid=CL/PLmid) #%>%
#  group_by(Taxon) %>% #to show ratios per Taxon, leave out to get a total ratio
#  summarise(meanRatio=round(mean(ratio),2), sdRatio=round(sd(ratio),2), n=n(), min=min(ratio), max=max(ratio))


#kruskal.test(CLPLextant$meanRatio, CLPLextant$Taxon)
#KruskalMC <- data.frame(kruskalmc(CLPLextant$ratio, CLPLextant$Taxon))

Ratio <- CLPLextant %>%
  dplyr::summarise(meanRatio=round(mean(ratio),2), sdRatio=round(sd(ratio),2), n=n(), min=min(ratio), max=max(ratio))

RatioSpecies <- CLPLextant %>%
  group_by(Taxon) %>% #to show ratios per Taxon, leave out to get a total ratio
  dplyr::summarise(meanRatio=round(mean(ratio),2), sdRatio=round(sd(ratio),2), n=n(), min=min(ratio), max=max(ratio))


fossil <- tidyCL %>%
  mutate(Taxon=as.character(Taxon), Island=as.character(Island), Continent=as.character(Continent), estimated=as.character(estimated), Genus=as.character(Genus)) %>%
  dplyr::select(Taxon, CL, PL, size, estimated, Age, Island, Continent, Genus)

modern <- Extant %>%
  dplyr::select(Taxon=Species, CL, PL, estimated, Age, Island, Continent, Genus)

All <- bind_rows(modern, fossil)

testRatio <- All %>% #tidyCL
  dplyr::select(Taxon, CL, PL, size, estimated, Age, Island, Continent, Genus) %>%
  mutate(extraCL = PL*Ratio$meanRatio) %>%
  dplyr::select(Taxon, CL, extraCL, PL, size, estimated, Age, Island, Continent, Genus)

testRatio$CL[is.na(testRatio$CL)] <- testRatio$extraCL[is.na(testRatio$CL)]



##### Time bins (stratigraphic stages) ######
# Bin data, smaller bins

Miocene <- testRatio %>% #testRatio or tidyCL 
  filter(Age < 23.000)

# PleiPlio$bin <- cut(PleiPlio$Age, c(0, 0.000001, 0.0117, 0.126, 0.781, 2.588, 3.6, 5.332, 11.608))
# EpochBins <- as.vector(c("Modern", "Modern", "Upper Pleistocene", "Middle Pleistocene", "Lower Pleistocene", "Upper Pliocene", "Lower Pliocene", "Upper Miocene"))
# MeanBins <- as.vector(c((0+0.000001)/2, (0+0.0117)/2, (0.0117+0.126)/2, (0.126+0.781)/2, (0.781+2.588)/2,(2.588+3.6)/2, (3.6+5.332)/2, (5.332+11.608)/2))
# new cuts: 1.806, 4.466, 7.424, 9.516, 13.789, 18

Miocene$bin <- cut(Miocene$Age, c(0, 0.0117, 0.126, 0.781, 1.806, 2.588, 3.6,
                                  5.332, 7.246, 11.608, 13.82, 15.97, 23.03))


EpochBins <- as.vector(c("Modern", "Upper Pleistocene", "Middle Pleistocene", "Lower Pleistocene", "Gelasian", "Piacencian", "Zanclean", "Messinian","Tortonian", "Serravallian","Langhian",
                         "Burdigalian/Aquitanian"))

Stages <- as.vector(c("Modern", "Upper Pleistocene", "Middle Pleistocene", "Lower Pleistocene", "Lower Pleistocene", "Upper Pliocene", "Lower Pliocene", "Upper Miocene","Upper Miocene",
                      "Middle Miocene","Middle Miocene", "Lower Miocene"))


MeanBins <- as.vector(c((0+0.0117)/2, (0.0117+0.126)/2, (0.126+0.781)/2, (0.781+1.806)/2, (1.806+2.588)/2,(2.588+3.6)/2, (3.6+5.332)/2, (5.332+7.246)/2, (7.246+11.608)/2, (11.608+13.82)/2, (13.82+15.97)/2, (15.97+23.03)/2))

#1.806, 4.466, 7.424, 9.516, 13.789, 18

BinsMio <- Miocene %>%
  select(bin) %>%
  group_by(bin) %>%
  dplyr::summarise(nIndividuals=n())

BinsSpeciesMio <- Miocene %>%
  select(bin, Taxon) %>%
  group_by(bin, Taxon) %>%
  dplyr::summarise(nSpecies=n()) %>%
  dplyr::summarise(nSpecies=n())


BinsGeneraMio <- Miocene %>%
  select(bin, Genus) %>%
  group_by(bin, Genus) %>%
  dplyr::summarise(nGenera=n()) %>%
  dplyr::summarise(nGenera=n())

bin <- as.vector(unique(BinsMio$bin))



BINSMio <- data.frame(bin, EpochBins, Stages, MeanBins) %>%
  merge(BinsMio) %>%
  merge(BinsSpeciesMio) %>%
  merge(BinsGeneraMio) %>%
  arrange(MeanBins)

BINSMio$EpochBins = factor(BINSMio$EpochBins, 
                           levels= c("Modern", "Upper Pleistocene", "Middle Pleistocene", "Lower Pleistocene", "Gelasian", "Piacencian", "Zanclean", "Messinian","Tortonian", "Serravallian","Langhian",
                                     "Burdigalian/Aquitanian" ))

#BINS <- read.table("timebins.txt", sep="\t", header=TRUE)

#kable(Bins, caption="Time bins with corresponding sample sizes (individuals)")

#kable(BINSMio, caption="Smaller time bins with age range, epoch name, mean age and corresponding sample sizes (on individual, species and genus level)")

PleiPlioCL <- Miocene %>%
  merge(BINSMio) %>%
  filter(!is.na(CL))

# SUM <- PleiPlioCL %>%
#   group_by(MeanBins) %>%
#   summarise(n=n())

#na <- PleiPlioCL$EpochBins[!complete.cases(PleiPlioCL$EpochBins)]



#### Data overview - Scatterplot CL~Time ####

Continent <- as.vector(unique(PleiPlioCL$Continent))
Con <- as.vector(c("Asia", "Africa", "Europe", "America", "America", "America", "Europe"))#
Continents <- data.frame(Continent, Con)

scatter <- testRatio %>%
  filter(!is.na(CL)) %>%
  filter(Age < 23.000) %>%
  merge(Continents) %>%
  mutate(Insularity=ifelse(Island=="y","insular", "continental")) %>%
  select(CL, Age, Insularity, Con, estimated) %>%
  ggplot(aes(-Age, CL))+ # for time to start at the left: -Age
  geom_vline(xintercept = c(0, -0.0117, -0.126, -0.781, -1.806, -2.588, -3.6,
                            -5.332, -7.246, -11.608, -13.82, -15.97, -23.03))+ 
  geom_jitter() + #aes(shape=Insularity, colour= Con)
  theme_classic() +
  # geom_vline(xintercept= c(0, 0.0117, 0.126, 0.781, 2.588, 3.6, 5.332, 11.608, 15.97, 23.03)) +
  geom_vline(xintercept= -20.44, linetype="dashed")  +
  scale_color_manual(values=c("#000000", "#E69F00",  "#FF00FF", "#009E73", "#56B4E9","#F0E442", "#0072B2", "#D55E00", "#CC79A7") , name="Continents")+
  # scale_shape_manual(values= c(16, 17), name="Insularity", breaks=c("continental", "insular")) +
  xlab("Time [mya]") + ylab("Carapace length [mm]")#+ #+ coord_trans(x="log2") +
#scale_x_log10()
# labs(shape="Insularity", colour="Continents")
# 


scatter
# fig.cap="Scatterplot of carapace length over time, indicating insular (triangle) and continental (circles)
# nd colour indicating continents. Lines indicate stratigraphic stages which were used as time bins, the dashed
# line is the border between the two stages of the Lower Miocene, which were consideres as one time bin.


#ggsave("OvervieData.png", height=5, width=8, units='in', dpi=800)

#length(scatter$CL)


#### Maps - fossil occurences of testudinidae #####

CL <- tidyCL %>%
  #  filter(Age < 11.000) %>%
  dplyr::select(Locality, Genus, Taxon, Latitude, Longitude, Country, Age) %>%
  group_by(Locality) %>%
  mutate(count= n())

length(unique(CL$Locality))

FossilOccurrences<-read.csv("tortoises_occurrences.csv", sep=";", header=TRUE)

FossOcc <- FossilOccurrences %>%
  mutate(Age= (as.numeric(as.character(MA.min)) + as.numeric(as.character(Ma.max)))/2) %>%
  select(Locality, Country, Latitude, Longitude, Age, Genus, Taxon, Clavailability) %>%
  #  merge(tidyCL) %>%
  group_by(Locality) %>%
  mutate(count= n(), Longitude=as.numeric(as.character(Longitude)))


Occurrences <- FossilOccurrences %>%
  mutate(Age= (as.numeric(as.character(MA.min)) + as.numeric(as.character(Ma.max)))/2) %>%
  select(Locality, Country, Latitude, Longitude, Age, Genus, Taxon, Author) %>%#, comment, Clavailability
  arrange(Age)

#kable(Occurrences, row.names = TRUE, caption="fossil occurrences")
#na <- FossOcc[!complete.cases(FossOcc),]

# OccMap <-  merge(Map, FossOcc, by="Locality", incomparables=TRUE) %>%
#  distinct()

# unique(CL$Locality) #wie viele localities
#length(unique(FossOcc$Locality[which(FossOcc$Clavailability == "yes")]))


mapWorld <- borders("world", colour="azure3", fill="azure3") # create a layer of borders, run line at the beginning (before loading plotly)

cbbPalette <- c("#000000", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

OwnPalette <- c("#000000", "#000000", "#56B4E9") #alternativ: c("#000000","#000000", "#E69F00", "#56B4E9")#

mapOc <- FossOcc %>%
  ggplot(aes(Longitude, Latitude)) +# borders("world", ylim = c(-60, 90)) +
  mapWorld +
  theme_classic() +
  geom_point(data=FossOcc, alpha=1/4,(aes(as.numeric(Longitude), Latitude, colour=FossOcc$Clavailability))) +  # alpha= Transparency
  geom_point(data=CL, alpha=1/2,(aes(Longitude, Latitude, colour="yes"))) +
 # geom_point(data=pdbd, alpha=1/6,(aes(Longitude, Latitude, colour="pdbd"))) +
  theme(legend.position="none") +
  scale_colour_manual(values=OwnPalette) +#cbbPalette +
  scale_x_continuous(breaks=c(-180, -100, 0, 100, 180)) +
  scale_y_continuous(breaks=c(-90, -50, 0, 50, 90))

#unique(FossOcc$CL)

mapOc
# fig.cap="Map displaying all fossil occurrences of testudinids, with color indicating whether relevant 
# literature was available (black if not) and if it was, whether body size data was available or not 
# (yes and no, respectively).


#ggsave("MapOccurrences.png", height=5, width=8, units='in', dpi=800)

# require(plotly)
# ggplotly(mapOc)

# length(CL$Locality)
# length(unique(CL$Locality)) #wie viele localities
# length(unique(FossOcc$Locality[which(FossOcc$Clavailability == "yes")]))

##### compare with PDBD data! #####
#pbdb_testudinidae.csv
setwd("C:/Users/Jule/Documents/Uni/MA/PublicationTortoises")
PDBD<-read.csv("pbdb_testudinidae.csv", sep=",", header=TRUE, skip=15)
names(PDBD)
pdbd <- PDBD %>%
  select(collection_no, identified_name, accepted_name, early_interval, late_interval, max_ma, min_ma,
         reference_no, Longitude=lng, Latitude=lat, cc, state, county, Locality=geogcomments)  %>%
  mutate(Age= (as.numeric(as.character(min_ma)) + as.numeric(as.character(max_ma)))/2) 

length(unique(pdbd$Locality)) # geogcomments = Locality in my table, I guess


mapWorld <- borders("world", colour="azure3", fill="azure3") # create a layer of borders, run line at the beginning (before loading plotly)

#cbbPalette <- c("#000000", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#OwnPalette <- c("#000000", "#000000", "#56B4E9") #alternativ: c("#000000", "#E69F00", "#56B4E9")

mapPDBD <- pdbd %>%
  ggplot(aes(Longitude, Latitude)) +# borders("world", ylim = c(-60, 90)) +
  mapWorld +
  theme_classic() +
  geom_point(data=pdbd, alpha=1/4,(aes(as.numeric(Longitude), Latitude, colour="#000000"))) +
 # geom_point(data=CL, alpha=1/2,(aes(Longitude, Latitude, colour="yes"))) +
  theme(legend.position="none") +
 # scale_colour_manual(values=OwnPalette) +#cbbPalette +
  scale_x_continuous(breaks=c(-180, -100, 0, 100, 180)) +
  scale_y_continuous(breaks=c(-90, -50, 0, 50, 90))

#unique(FossOcc$CL)

mapPDBD

##### check sample size (?) for Miocene!#####

CLpdbdCom <- tidyCL %>%
  filter(Age < 23.000 & Age > 5.330) %>%
  dplyr::select(Locality, Genus, Taxon, Latitude, Longitude, Country, Age) %>%
  group_by(Locality) %>%
  mutate(count= n())

FossOccPdbdCom <- FossilOccurrences %>%
  mutate(Age= (as.numeric(as.character(MA.min)) + as.numeric(as.character(Ma.max)))/2) %>%
  filter(Age < 23.000 & Age > 5.330) %>%
  select(Locality, Country, Latitude, Longitude, Age, Genus, Taxon, Clavailability) %>%
  #  merge(tidyCL) %>%
  group_by(Locality) %>%
  mutate(count= n(), Longitude=as.numeric(as.character(Longitude)))

Pdbd <- pdbd %>%
  filter(Age < 23.000 & Age > 5.330) %>%
  group_by(Locality) %>%
  mutate(count= n())

mapWorld <- borders("world", colour="azure3", fill="azure3") # create a layer of borders, run line at the beginning (before loading plotly)

cbbPalette <- c("#000000", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

OwnPalette <- c("#000000","#000000", "#E69F00", "#56B4E9")#c("#000000", "#000000", "#56B4E9") #alternativ: 

mapOcPDBDMiocene <- FossOccPdbdCom %>%
  ggplot(aes(Longitude, Latitude)) +# borders("world", ylim = c(-60, 90)) +
  mapWorld +
  theme_classic() +
  geom_point(data=FossOccPdbdCom, alpha=1/4,(aes(as.numeric(Longitude), Latitude, colour="FossOcc"))) +  # alpha= Transparency
  geom_point(data=CLpdbdCom, alpha=1/2,(aes(Longitude, Latitude, colour="CL"))) +
  geom_point(data=Pdbd, alpha=1/6,(aes(Longitude, Latitude, colour="PDBD"))) +
#  theme(legend.position="none") +
  scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9")
                      #,labels=c("FossOcc", "CL", "PDBD")
                      ) +#cbbPalette +
  scale_x_continuous(breaks=c(-180, -100, 0, 100, 180)) +
  scale_y_continuous(breaks=c(-90, -50, 0, 50, 90))

#unique(FossOcc$CL)

mapOcPDBDMiocene


#ggsave("MapMioceneComparison.png", height=5, width=8, units='in', dpi=800)
# require(plotly)
# ggplotly(mapOcPDBDMiocene)

##### PDBD data #####
pdbd
unique(pdbd$accepted_name)
length(unique(pdbd$accepted_name))

unique(pdbd$identified_name)
length(unique(pdbd$identified_name))

PDBDMiocene <- pdbd %>%
  filter(Age < 23.000 & Age > 5.330)

speciesPDBDMio <- as.vector(unique(PDBDMiocene$identified_name))
length(speciesPDBDMio)


CLMiocene <- tidyCL %>%
  filter(Age < 23.000 & Age > 5.330)

speciesCLMiocene <- as.vector(unique(CLMiocene$Taxon))
length(speciesCLMiocene)
  

compareSpeciesMio <- setdiff(speciesPDBDMio, speciesCLMiocene)



write.table(compareSpeciesMio, file="PDBDSpeciesMioNeeded.txt", sep="\t", row.names = FALSE)

# do the same for period younger than miocene, just to be sure

PDBDQ <- pdbd %>%
  filter(Age < 5.330)

speciesPDBDQ <- as.vector(unique(PDBDQ$identified_name))
length(speciesPDBDQ)


CLQ <- tidyCL %>%
  filter(Age < 5.330)

speciesCLQ <- as.vector(unique(CLQ$Taxon))
length(speciesCLQ)


compareSpeciesQ <- setdiff(speciesPDBDQ, speciesCLQ)

write.table(compareSpeciesQ, file="PDBDSpeciesHolPleiPlio_justToCheck.txt", sep="\t", row.names = FALSE)



#### Map - Body size of testudinidae ####

tidyCL$bin <- cut(tidyCL$Age, c(0, 0.0117, 0.126, 0.781, 1.806, 2.588, 3.6,
                                5.332, 7.246, 11.608, 13.82, 15.97, 20.44, 23.03)
)

#                c(0, 0.000001, 0.0117, 0.126, 0.781, 2.588, 3.6, 5.332, 11.608, 15.97, 23.03,50.000))


MapCL <- tidyCL %>%
  merge(BINSMio) %>%
  #  filter(Age < 23.000) %>%
  #  filter(Latitude != "-") %>%
  dplyr::select(Genus, Taxon, Latitude, Longitude, Country, CL, PL, MeanBins, size) %>%
  group_by(Latitude) %>%
  mutate(count= n()) 


mapWorld <- borders("world", colour="azure3", fill="azure3") # create a layer of borders, run line at the beginning (before loading plotly)


mapCL <- MapCL %>%
  ggplot(aes(Longitude, Latitude)) +# borders("world", ylim = c(-60, 90)) +
  mapWorld +
  theme_classic() +
  geom_point(aes(Longitude, Latitude,size=count, colour=MeanBins)) +  #colour=size, 
  scale_colour_gradientn(name="Age [mya]", colors=c("orange", "red", "purple", "blue", "green", "yellow"))+
  scale_x_continuous(breaks=c(-180, -100, 0, 100, 180)) +
  scale_y_continuous(breaks=c(-90, -50, 0, 50, 90)) +
  scale_size_continuous(name="n")

#values=c("#000000", "#E69F00",  "#FF00FF", "#009E73", "#56B4E9","#F0E442", "#0072B2", "#D55E00", "#CC79A7")


mapCL
# fig.cap="Map displaying all localities for which body size data for testudinids was available in the 
# literature. Size of points denotes sample size, color denotes approximate age.


#ggsave("MapBodysize.png", height=5, width=8, units='in', dpi=800)


#### Get overview over sample sizes #####
OverviewSpecies <- PleiPlioCL %>%
  filter(EpochBins != "Modern") %>%
  group_by(EpochBins, Taxon) %>%
  dplyr::summarise(n=n(), meanCL = mean(CL))

#kable(OverviewSpecies, caption="Overview over fossil species per time bin, with sample size and mean CL.")

OverviewSpeciesFossil <- PleiPlioCL %>%
  filter(EpochBins != "Modern") %>%
  group_by(Taxon) %>%
  dplyr::summarise(n=n(), meanCL = mean(CL))

#write.table(OverviewSpecies,file="OverviewSpeciesSampleSize.txt", sep="\t", row.names = FALSE) 
#kable(OverviewSpeciesFossil, caption="General overview over fossil species, with sample size and mean CL")

OverviewGeneraTime <- PleiPlioCL %>%
  group_by(EpochBins, Genus) %>%
  dplyr::summarise(n=n(), meanCL = mean(CL))

#kable(OverviewGeneraTime, caption="Overview over genera (modern and fossil) per time bin, with sample sizes and mean CL.")


OverviewGenera <- PleiPlioCL %>%
  group_by(Genus) %>%
  dplyr::summarise(n=n(), meanCL = mean(CL))

#kable(OverviewGenera, caption="General overview over genera, with sample sizes and mean CL.")

# plot(tidyCL$CL~tidyCL$Latitude)
# plot(tidyCL$CL~tidyCL$Longitude)


##### Sampling Accumulation Curves #####
# SACSpecies

#Species Accumulation Curve
veganSL <- tidyCL %>%
  dplyr::select(Locality, Taxon) %>%
  group_by(Locality, Taxon) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Taxon, n, fill=0)


#par(mfrow = c(2, 1))

veganSL=veganSL[,-1]
veganspSL=specaccum(veganSL,method="rarefaction", permutations=1000)
plot(veganspSL,xlab="No. of Localities",ylab="Species Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2 )


veganSR <- tidyCL %>%
  dplyr::select(Reference, Taxon) %>%
  group_by(Reference, Taxon) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Taxon, n, fill=0)

veganSR=veganSR[,-1]
veganspSR=specaccum(veganSR,method="rarefaction", permutations=1000)
#jpeg('SACSpecies.jpg')
plot(veganspSR,xlab="No. of References",ylab="Species Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#dev.off()
#fig.cap="Species Accumulation Curve of fossil species per Locality and reference



#SACGenera

veganGL <- tidyCL %>%
  dplyr::select(Locality, Genus) %>%
  group_by(Locality, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


#par(mfrow = c(2, 1))


veganGL=veganGL[,-1]
veganspGL=specaccum(veganGL,method="rarefaction", permutations=1000)
plot(veganspGL,xlab="No. of Localities",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2, main="Fossil genera, CL, per Locality" )
# --> appendix!

veganGR <- tidyCL %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)

#library(vegan)

veganGR=veganGR[,-1]
veganspGR=specaccum(veganGR,method="rarefaction", permutations=1000)
#jpeg('SACGenera.jpg')
plot(veganspGR,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#dev.off()
#fig.cap="Sampling Accumulation Curve of fossil genera per reference


#SACGEurasia

#Species Accumulation Curve with Genera, Eurasia

overviewVegan <- PleiPlioCL %>%
  merge(Continents) %>%
  filter(EpochBins != "Modern") %>%
  group_by( Con) %>%  #Continent,
  # filter(Con=="Europe" | Con=="Asia" )%>%
  filter(!is.na(CL)) %>%
  dplyr::summarise(meanCL=mean(CL), sdCL= sd(CL), n=n(), meanAge=mean(Age))


veganEA <- tidyCL %>%
  merge(Continents) %>%
  filter(Con=="Europe" | Con=="Asia") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


veganEA=veganEA[,-1]
veganspEA=specaccum(veganEA,method="rarefaction", permutations=1000)
plot(veganspEA,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, Eurasia


#SACGEurope
#Species Accumulation Curve with Genera, Europe

veganEu <- tidyCL %>%
  merge(Continents) %>%
  filter(Con=="Europe") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


veganEu=veganEu[,-1]
veganspEu=specaccum(veganEu,method="rarefaction", permutations=1000)
plot(veganspEu,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, Europe


#SACGAfrica
#Species Accumulation Curve with Genera, Africa

veganAf <- tidyCL %>%
  merge(Continents) %>%
  filter(Con=="Africa") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)

veganAf=veganAf[,-1]
veganspAf=specaccum(veganAf,method="rarefaction", permutations=1000)
plot(veganspAf,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, Africa


#SACGAmerica

#Species Accumulation Curve with Genera, America

veganAm <- tidyCL %>%
  merge(Continents) %>%
  filter(Con=="America") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


veganAm=veganAm[,-1]
veganspAm=specaccum(veganAm,method="rarefaction", permutations=1000)
plot(veganspAm,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, America


#SACGNAmerica
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, N-America"}
#Species Accumulation Curve with Genera, N-America


veganNA <- tidyCL %>%
  merge(Continents) %>%
  filter(Continent=="N-America") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


veganNA=veganNA[,-1]
veganspNA=specaccum(veganNA,method="rarefaction", permutations=1000)
plot(veganspNA,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, N-America


#SACGSAmerica
#Species Accumulation Curve with Genera, S-America

veganSA <- tidyCL %>%
  merge(Continents) %>%
  filter(Continent=="S-America") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


veganSA=veganSA[,-1]
veganspSA=specaccum(veganSA,method="rarefaction", permutations=1000)
plot(veganspSA,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, S-America


#SACGAsia
#Species Accumulation Curve with Genera, Asia

veganAs <- tidyCL %>%
  merge(Continents) %>%
  filter(Con=="Asia") %>%
  dplyr::select(Reference, Genus) %>%
  group_by(Reference, Genus) %>%
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)


veganAs=veganAs[,-1]
veganspAs=specaccum(veganAs,method="rarefaction", permutations=1000)
plot(veganspAs,xlab="No. of References",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2)
#fig.cap="Sampling Accumulation Curve of fossil genera per reference, Asia


#SAC fossil occurences
SACFossOc <- FossilOccurrences %>%
  dplyr::select(Locality, Genus) %>%  #, Taxon
  group_by(Locality, Genus) %>%  #, Taxon
  dplyr::summarise(n=n()) %>%
  tidyr::spread(Genus, n, fill=0)  #Taxon

#par(mfrow = c(2, 1))

SACFossOc=SACFossOc[,-1]
SACFossOcSL=specaccum(SACFossOc,method="rarefaction", permutations=1000)
plot(SACFossOcSL,xlab="No. of Localities",ylab="Genera Richness", xvar="individuals", ci.type="line", ci.lty=2, ci.col="grey", col="deepskyblue4", lwd=2 )


#####Histograms#####


#HistAll
#Histograms of body size data, all

HistCL <- PleiPlioCL %>%
  filter( !is.na(CL)) %>%
  mutate(CL=log(CL)) %>%
  #  filter(EpochBins != "Modern") %>%
  ggplot(aes(CL)) + 
  #geom_histogram( col="black", fill="gray") + 
  geom_histogram( aes(y=..density..),col="black", fill="gray")+
  #geom_density() +
  geom_density(alpha=.2, fill="darkslategray1") +
  #facet_wrap(~EpochBins, scales="free_x")+   #binwidth=10,
  theme_classic()+ 
  #theme_gray() +
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  xlab("Log Carapace length [mm]")

#HistCLModern
HistCL  
##fig.cap="Distribution of body size data, logtransformed, all data.
#ggsave("HistAll.png", height=5, width=8, units='in', dpi=800)



StatsAll <- PleiPlioCL %>%
  filter(!is.na(CL)) %>%
  #  group_by(Island) %>%
  dplyr::summarise(nCL=length(CL), #range=range(CL),
                   min=min(CL), max=max(CL), 
                   var=var(CL), mean=round(mean(CL),1), logm=round(mean(log10(CL)),1),
                   med=round(median(CL),1), logmed=round(median(log10(CL)),1), #CLmode=mode(CL), CLlogmode=mode(log10(CL)),
                   skew=round(skewness(CL),2), logsk=round(skewness(log10(CL)),2),
                   kurt=round(kurtosis(CL),2), logku=round(kurtosis(log10(CL)),2) ) %>%
  unique()

#normalDistribution
qqnorm(PleiPlioCL$CL); qqline(PleiPlioCL$CL, col=2)
qqnorm(log10(PleiPlioCL$CL)); qqline(log10(PleiPlioCL$CL), col=2)



## per time bin
#HistBins
#Histograms of body size data, per time bin


HistCLBins <- PleiPlioCL %>%
  # merge(BINS) %>%
  filter( !is.na(CL)) %>%
  mutate(CL=log(CL)) %>%
  arrange(MeanBins) %>%
  #  filter(EpochBins != "Modern") %>%
  ggplot(aes(CL)) + 
  geom_histogram( aes(y=..density..),col="black", fill="gray")+
  #geom_density() +
  geom_density(alpha=.2, fill="darkslategray1") +
  facet_wrap(~EpochBins, scales="free_x")+   #binwidth=10,
  theme_classic()+ 
  #theme_gray() +
  theme(panel.border = element_rect(colour = "black", fill=NA))+ 
  xlab("Log Carapace length [mm]")


HistCLBins
#fig.cap="Distribution of body size data per time bin, logtransformed.


StatsBins <- PleiPlioCL %>%
  filter(!is.na(CL)) %>%
  group_by(EpochBins) %>%
  dplyr::summarise(nCL=length(CL), min=min(CL), max=max(CL), var=var(CL), mean=round(mean(CL),1), logm=round(mean(log10(CL)),1),
                   med=round(median(CL),1), logmed=round(median(log10(CL)),1), #CLmode=mode(CL), CLlogmode=mode(log10(CL)),
                   skew=round(skewness(CL),2), logsk=round(skewness(log10(CL)),2),
                   kurt=round(kurtosis(CL),2), logku=round(kurtosis(log10(CL)),2) ) %>%
  unique()



## modern vs. fossil
#HistFosMo
#Histograms of body size data, modern vs. fossil


Epoch <- as.vector(c("Modern", "Pleistocene", "Pleistocene", "Pleistocene", "Pleistocene", 
                     "Pliocene", "Pliocene",
                     "Miocene", "Miocene", "Miocene",  "Miocene", "Miocene"))
EPOCH <- as.vector(c("Modern", "Fossil", "Fossil", "Fossil"
                     , "Fossil", "Fossil", "Fossil", "Fossil", "Fossil", "Fossil"
                     , "Fossil", "Fossil"))
EpochBins <- as.vector(unique(PleiPlioCL$EpochBins))
Epochs <- data.frame(Epoch,EpochBins, EPOCH)


Epochs$EPOCH = factor(Epochs$EPOCH,levels=c("Modern", "Fossil"))

HistCLFossil <- PleiPlioCL %>%
  merge(Epochs) %>%
  filter( !is.na(CL)) %>%
  mutate(CL=log(CL)) %>%
  #  filter(EpochBins != "Modern") %>%
  ggplot(aes(CL)) + 
  geom_histogram( aes(y=..density..),col="black", fill="gray")+
  #geom_density() +
  geom_density(alpha=.2, fill="darkslategray1") +
  facet_wrap(~EPOCH)+   #binwidth=10,, scales="free_x"
  theme_classic()+ 
  #theme_gray() +
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  xlab("Log Carapace length [mm]")


HistCLFossil
#ggsave("HistModernFossil.png", height=5, width=8, units='in', dpi=800)


StatsFossil <- PleiPlioCL %>%
  merge(Epochs) %>%
  filter(!is.na(CL)) %>%
  group_by(EPOCH) %>%
  dplyr::summarise(nCL=length(CL), min=min(CL), max=max(CL), var=var(CL), mean=round(mean(CL),1), logm=round(mean(log10(CL)),1),
                   med=round(median(CL),1), logmed=round(median(log10(CL)),1), #CLmode=mode(CL), CLlogmode=mode(log10(CL)),
                   skew=round(skewness(CL),2), logsk=round(skewness(log10(CL)),2),
                   kurt=round(kurtosis(CL),2), logku=round(kurtosis(log10(CL)),2) ) %>%
  unique()


## modern vs. fossil, continental vs. insular
#Histograms of body size data, modern vs. fossil, continental vs. insular

HistCLFossilIsland <- PleiPlioCL %>%
  merge(Epochs) %>%
  filter( !is.na(CL)) %>%
  mutate(CL=log(CL)) %>%
  mutate(Insularity=ifelse(Island=="y","insular", "continental")) %>%
  #  filter(EpochBins != "Modern") %>%
  ggplot(aes(CL)) + 
  geom_histogram( aes(y=..density..),col="black", fill="gray")+
  #geom_density() +
  geom_density(alpha=.2, fill="darkslategray1") + 
  facet_grid(EPOCH~Insularity)+   #binwidth=10,, scales="free_x"
  theme_classic()+ 
  #theme_gray() +
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  xlab("Log Carapace length [mm]")


HistCLFossilIsland

ggsave("HistModernFossilInsularContinental.png", height=5, width=8, units='in', dpi=800)
#fig.cap="Distribution of body size data modern vs. fossil, continental vs. insular logtransformed.

StatsFossilIsland <- PleiPlioCL %>%
  merge(Epochs) %>%
  filter(!is.na(CL)) %>%
  group_by(EPOCH, Island) %>%
  dplyr::summarise(nCL=length(CL), min=min(CL), max=max(CL), var=var(CL), mean=round(mean(CL),1), logm=round(mean(log10(CL)),1),
                   med=round(median(CL),1), logmed=round(median(log10(CL)),1), #CLmode=mode(CL), CLlogmode=mode(log10(CL)),
                   skew=round(skewness(CL),2), logsk=round(skewness(log10(CL)),2),
                   kurt=round(kurtosis(CL),2), logku=round(kurtosis(log10(CL)),2) ) %>%
  unique()



## continental vs. insular

##HistCI", echo=FALSE, 
#Histograms of body size data, continental vs. insular

HistCLIsland <- PleiPlioCL %>%
  filter( !is.na(CL)) %>%
  mutate(CL=log(CL)) %>%
  mutate(Insularity=ifelse(Island=="y","insular", "continental")) %>%
  #  filter(EpochBins != "Modern") %>%
  ggplot(aes(CL)) + 
  geom_histogram( aes(y=..density..),col="black", fill="gray")+
  #geom_density() +
  geom_density(alpha=.2, fill="darkslategray1") + 
  facet_wrap(~Insularity)+   #binwidth=10,, scales="free_x"
  theme_classic()+ 
  #theme_gray() +
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  xlab("Log Carapace length [mm]")

#HistCLModern
HistCLIsland

#fig.cap="Distribution of body site data of continental (n) and insular(y) species, logtransformed.


#ggsave("HistInsularContinental.png", height=5, width=8, units='in', dpi=800)


StatsIsland <- PleiPlioCL %>%
  filter(!is.na(CL)) %>%
  group_by(Island) %>%
  dplyr::summarise(nCL=length(CL), min=min(CL), max=max(CL), var=var(CL), mean=round(mean(CL),1), logm=round(mean(log10(CL)),1),
                   med=round(median(CL),1), logmed=round(median(log10(CL)),1), #CLmode=mode(CL), CLlogmode=mode(log10(CL)),
                   skew=round(skewness(CL),2), logsk=round(skewness(log10(CL)),2),
                   kurt=round(kurtosis(CL),2), logku=round(kurtosis(log10(CL)),2) ) %>%
  unique()

## continents
##"HistCon", echo=FALSE,
#Histograms of body size data, split by continents

Continent <- as.vector(unique(PleiPlioCL$Continent))
Con <- as.vector(c("Asia", "Africa", "Europe", "America", "America", "America", "Europe"))#
Continents <- data.frame(Continent, Con)

HistCLContinents <- PleiPlioCL %>%
  filter( !is.na(CL)) %>%
  mutate(CL=log(CL)) %>%
  merge(Continents) %>%
  #  filter(EpochBins != "Modern") %>%
  ggplot(aes(CL)) + 
  geom_histogram( aes(y=..density..),col="black", fill="gray")+
  #geom_density() +
  geom_density(alpha=.2, fill="darkslategray1") +
  facet_wrap(~Con)+   #binwidth=10,, scales="free_x"
  theme_classic()+ 
  #theme_gray() +
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  xlab("Log Carapace length [mm]")


HistCLContinents
#fig.cap="Distribution of body site data per continent, logtransformed.



StatsContinents <- PleiPlioCL %>%
  merge(Continents) %>%
  filter(!is.na(CL)) %>%
  group_by(Con) %>%
  dplyr::summarise(nCL=length(CL), min=min(CL), max=max(CL), var=var(CL), mean=round(mean(CL),1), logm=round(mean(log10(CL)),1),
                   med=round(median(CL),1), logmed=round(median(log10(CL)),1), #CLmode=mode(CL), CLlogmode=mode(log10(CL)),
                   skew=round(skewness(CL),2), logsk=round(skewness(log10(CL)),2),
                   kurt=round(kurtosis(CL),2), logku=round(kurtosis(log10(CL)),2) ) %>%
  unique()



#### Descriptive statistics ####
#Tabel Stats Distribution CL", echo=FALSE}
Stats <- bind_rows(StatsAll, StatsBins, StatsFossil, StatsIsland, 
                   StatsFossilIsland, StatsContinents) %>%
  select(-EpochBins, -Island, -Con, -EPOCH)



Stats$Variable <- c("all", "Modern", "Upper Pleistocene", "Middle Pleistocene", "Lower Pleistocene", "Gelasian", "Piacencian", "Zanclean", "Messinian","Tortonian", "Serravallian","Langhian",
                    "Burdigalian/Aquitanian", "Modern","Fossil", "continental",
                    "insular","modern-con", "modern-ins", "fossil-con", "fossil-ins",
                    "Africa", "America", "Asia", "Europe")


#kable(Stats, caption="General statistics of body size data: all, per time bin, insular and continental, per continent (all referring to CL: min, max, variance, mean, logmean, median, logmedian, skewness, logskewness, kurosis, logkurtosis")


####### Boxplots #####
## genera per time bins

#BPGBins", echo=FALSE, warning=FALSE, 


IndGenera <- PleiPlioCL %>%
  group_by(EpochBins, Genus) %>%
  filter(!is.na(CL)) %>%
  #  filter(Island == "y") %>%
  dplyr::summarise(GenusMean=mean(CL), GenusSD=sd(CL), n=n()) 

IndGenera[is.na(IndGenera)] <-0

# givesample sizes
give.n <- function(x){
  return(c(y = 2500, label = length(x)))
}



BinGenera <- PleiPlioCL %>%
  group_by(EpochBins, Genus) %>%
  filter(!is.na(CL)) %>%
  #  filter(Island == "y") %>%
  dplyr::summarise(GenusMean=mean(CL), GenusSD=sd(CL), n=n()) %>%
  #  summarise(TimeBinMean=mean(GenusMean), TimeBinSD=sd(GenusMean), n=n()) %>%
  ungroup() %>%#, CL=CL
  ggplot(aes(EpochBins, GenusMean)) + #, colour=Genus
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot() +
  geom_pointrange(data=IndGenera, position="jitter", aes(x=EpochBins, y=GenusMean, colour=EpochBins, ymin = GenusMean-GenusSD, ymax = GenusMean+GenusSD)) +
  stat_summary(fun.data = give.n, geom = "text") +
  theme(legend.background = element_rect(colour = 'black'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  scale_colour_manual(values=c("#999999", "#969696", "#999999", "#999999", "#999999", "#969696",
                               "#999999", "#999999", "#999999", "#969696", "#999999",
                               "#999999", "#969696", "#999999", "#999999", "#999999",
                               "#999999", "#999999")) +
  theme(legend.position="none") + ylab("Carapace length [mm]") + xlab("Stratigraphic Stages")+
  annotate("text", x=0.7, y=2500,label= "n=")

BinGenera
#fig.cap="Boxplots of mean CL per time bin, including mean and sd CL for each genus (as pointrange).


#ggsave("BoxplotBins.png", height=5, width=8, units='in', dpi=800)

####generic level
# Multiple comparison after Kruskal-Wallis Test



names(IndGenera)
library(pgirmess)
kruskalmc(IndGenera$GenusMean,IndGenera$EpochBins)
MCK <- data.frame(kruskalmc(IndGenera$GenusMean,IndGenera$EpochBins))


# resp<-c(0.44,0.44,0.54,0.32,0.21,0.28,0.7,0.77,0.48,0.64,0.71,0.75,0.8,0.76,0.34,0.80,0.73,0.8)
# categ<-as.factor(rep(c("A","B","C"),times=1,each=6))
# kruskalmc(resp, categ)


#### 
names(PleiPlioCL)

kruskalmc(PleiPlioCL$CL,PleiPlioCL$EpochBins)
MCKW <- data.frame(kruskalmc(PleiPlioCL$CL,PleiPlioCL$EpochBins))




########## paleoTS analysis #########

## all (continental and insular)

### genera (all)
##paleoTSAll", echo=FALSE, include=FALSE, 
#paleoTS plot with genus mean, including island species

GenusMean <- PleiPlioCL %>%
  filter(EpochBins != "Modern") %>%
  filter(!is.na(CL)) %>%
  group_by(EpochBins, Genus) %>%
  #  filter(CL < 999 & Island =="y") %>%
  dplyr::summarise(meanCL = mean(CL), n=n()) %>%
  #  filter(EpochBins != "Aquitanian") %>%
  #  ungroup() %>%
  #  group_by(EpochBins) %>%
  dplyr::summarise(mm = mean(meanCL), nn=n(), vv=var(meanCL)) %>%
  merge(BINSMio) %>%
  select(mm, nn, vv, tt=MeanBins)

#na <- PleiPlioCL[!complete.cases(PleiPlioCL),]

GenusMean[is.na(GenusMean)] <- 0


GenusModern <- PleiPlioCL %>%
  filter(EpochBins == "Modern") %>%
  filter(!is.na(CL)) %>%
  #  filter(CL < 999& Island =="y") %>%
  group_by(Genus) %>%
  dplyr::summarise(meanCL=mean(CL), sdCL=sd(CL), n=n(), Age=mean(Age), EpochBins=unique(EpochBins), MeanBins=unique(MeanBins)) #%>%
#  dplyr::summarise(mm=mean(meanCL), nn=n(), vv=var(meanCL), tt=mean(MeanBins))



sumTort <- read.csv("tortoises_summary.csv", sep=";", header=TRUE)
#colnames(sumTort)[1] <- "Taxon"
colnames(sumTort)[7] <- "meanCL"
colnames(sumTort)[8] <- "sdCL"
#sumTort$EpochBins <- "Modern"
#sumTort$MeanBins <- 0.0000005


sumTortoises <- sumTort %>%
  mutate(Genus= as.character(Genus)) %>%
  mutate(MeanBins=(Mamin+Mamax)/2) %>%
  select(Genus, MeanBins,  meanCL, sdCL, n, Island) %>% #n,
  bind_rows(GenusModern)


SumTort <- sumTortoises %>%
  group_by(Genus) %>%
  #  filter(meanCL < 999 & Island =="y") %>%
  mutate(tt=MeanBins, vv=sdCL^2, nn=n, mm=meanCL) %>%
  dplyr::select(mm, nn, vv, tt, Genus) %>%
  mutate(nx = nn*mm) %>%
  mutate(mmall=sum(nx)/sum(nn)) %>%
  mutate(SD=sqrt(nx), d=mm-mmall) %>%
  mutate(nsd=((nx^2+d^2)*nn)) %>%
  mutate(varall=sum(nsd)/sum(nn), n=sum(nn)) %>%
  dplyr::select(mm=mmall, vv=varall, nn=n, tt, Genus) %>%
  unique() %>%
  dplyr::select(CL=mm, n=nn, var=vv, tt, Genus)

#write.table(SumTort,file="SumTortModernGenus.txt", sep="\t", row.names = FALSE)

#kable(SumTort, caption="Overview over body size means per time bin on genus level.")

#boxplot(SumTort$CL, caption="Body size distribution in time bin 'modern'.")


modernMeanGenus <- SumTort %>%
  ungroup() %>%
  select(CL, n, tt) %>%
  filter( !is.na(CL)) %>%
  #group_by(tt) %>%
  dplyr::summarise(mm = mean(CL), vv=var(CL), nn=n(), tt=5.85e-03)

#bis hier alle modernen taxa zusammengefasst (summarised and MFN_testudinidae)


GenusPaleo <- modernMeanGenus %>%
  bind_rows(GenusMean)%>%
  arrange(tt) %>%
  select(tt,nn,mm,vv) %>%
  filter(nn!=0)

GenusPaleo$vv[is.na(GenusPaleo$vv)] <- 0


#kable(GenusPaleo,caption="paleoTS object, all data")


paleoGen <-as.paleoTS(GenusPaleo$mm, GenusPaleo$vv, GenusPaleo$nn, GenusPaleo$tt, MM = NULL, genpars = NULL, reset.time=TRUE)

paleoGen$tt = -paleoGen$tt


#jpeg('paleoTSAll.jpg', width=800, height=500)
plot(paleoGen)
#fig.cap="paleoTS plot with genus mean, all

# plot(GenusPaleo$tt, GenusPaleo$mm, type="b", xlab="Time", ylab="Trait Mean")
# arrows(GenusPaleo$tt,GenusPaleo$mm-GenusPaleo$vv,GenusPaleo$tt,GenusPaleo$mm+GenusPaleo$vv, code=3, length=0.02, angle = 90)

abline(h=mean(GenusPaleo$mm), lty=5)

points(x=c(-2.59, -5.33), y=c(245, 245), pch=17)
#dev.off()


####Model-fitting, genera ######

PaleoGenFit <- (fit3models(paleoGen, silent=FALSE, method="AD", pool=FALSE))





#________________________________________________________________________________________________
#######try Juans paleoTS-Code ####


paleo_data <- paleoGen  #from Juan
plot(paleo_data)   
ou <- opt.GRW(paleo_data, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE)  #from Juan
bm <- opt.URW(paleo_data, pool = TRUE, cl = list(fnscale=-1), meth = "L-BFGS-B", hess = FALSE)    #from Juan
st <- opt.Stasis(paleo_data, pool = TRUE, cl = list(fnscale=-1), meth = "L-BFGS-B", hess = FALSE) #from Juan
#pun <- fitGpunc(paleo_data,oshare=F,method= "AD" , pool=F)   ## not working!! #from Juan

#compareModels(bm, ou,st, silent = FALSE)#pun,   #from Juan

cat(ou$value, bm$value, st$value, "\n")

#documentation:
## generate data for a directional sequence 
## y <- sim.GRW(ns=30, ms=1, vs=1) 
# plot(y) 
# m.rw<- opt.GRW(y) 
# m.rwu<- opt.URW(y) 
# m.sta<- opt.Stasis(y)
## print log-likelihoods; easier to use function fit3models() 
## cat(m.rw$value, m.rwu$value, m.sta$value, "\n")




####NOT WORKING YET

# library(MuMIn)
# aic.w <- Weights(c(st$AICc,pun$AICc))
# cbind(c("stasis", "shift"), aic.w)
# aic.w[2]/aic.w[1]
# shift.time <- paleo_data$tt[pun$parameters[5]]
# 
# 
# 
# par(mfrow=c(2,2))
# plot(paleo_data, modelFit= bm, pch=19, lwd=1.5, ylim=NULL, las=1)
# plot(paleo_data, modelFit= ou, pch=19, lwd=1.5, ylim=NULL, las=1)
# plot(paleo_data, modelFit= st, pch=19, lwd=1.5, ylim=NULL, las=1)
# plot(paleo_data, modelFit=pun, pch=19, lwd=1.5, ylim=NULL, las=1)

