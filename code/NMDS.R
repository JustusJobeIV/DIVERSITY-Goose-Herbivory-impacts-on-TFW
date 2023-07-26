###Packages
library(vegan)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

NACE_Herb<-read.csv("~/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/raw_data/NACE herbivory study 2009_2021 veg data 21jul21.csv")

##Clear out N/A
NACE_Herb<-na.omit(NACE_Herb)
## remove comment column
NACE_Herb<-NACE_Herb[,-9]

# Make June 2009 its own dataframe
NACE_June2009=NACE_Herb =filter(NACE_Herb, (month %in% c("JUN")))
NACE_June2009= filter(NACE_June2009,(year %in% "2009"))

#### Remove data from 2016-2021####

NACE_Herb =filter(NACE_Herb, (year %in% c("2009", "2010","2011","2012","2013","2014","2015")))
NACE_Herb =filter(NACE_Herb, (month %in% c("AUG")))

## Make June 2009 its own dataframe
NACE_June2009=filter(NACE_Herb, (month %in% c("JUN")))
NACE_June2009= filter(NACE_June2009,(year %in% "2009"))

##Combine both dataframes
NACE_HerbCombined<- full_join(NACE_Herb, NACE_June2009)
write.csv(NACE_HerbCombined,"NACE_Final_Data.csv", row.names = TRUE)

# casting data into a site by species matrix
NACE_Herb.wide_sum <- dcast(NACE_HerbCombined,
                            year + month + module + trt + habitat ~ spp,
                            value.var = "percov",
                            fun.aggregate = sum, fill = 0)

write_csv(NACE_Herb.wide_sum, file='/Users/justusjobe/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/processed data/NACE_HERB.SppMAX.csv')

NACE_Herb_spp.mat<-NACE_Herb.wide_sum[6:36]

### Removing Trash
NACE_Herb_spp.mat<-NACE_Herb_spp.mat[-28]



### transforming to relative abundance

NACE.spp.rel <-decostand(NACE_Herb_spp.mat, method = "total")

NACE.spp_NMS <- metaMDS(NACE.spp.rel, 
                        distance = "bray",
                        k=2, 
                        maxit=999,
                        trymax=1000, 
                        wascores =TRUE)
#### No Convergence on k=2
NACE.spp_NMS <- metaMDS(NACE.spp.rel, 
                        distance = "bray",
                        k=3, 
                        maxit=999,
                        trymax=1000, 
                        wascores =TRUE)

### Convergence on k=3 #### 

# Shepards test/goodness of fit
goodness(NACE.spp_NMS) # Produces a results of test statistics for goodness of fit for each point
stressplot(NACE.spp_NMS) # Produces a Shepards diagram


plot(NACE.spp_NMS$points)

NACE_xy<-data.frame(NACE.spp_NMS$points)

NACE_xy$year<-NACE_Herb.wide_sum$year
NACE_xy$month<-NACE_Herb.wide_sum$month
NACE_xy$module<-NACE_Herb.wide_sum$module
NACE_xy$trt<-NACE_Herb.wide_sum$trt
NACE_xy$habitat<-NACE_Herb.wide_sum$habitat

### Habitat X Treatment column
NACE_xy1<-NACE_xy

NACE_xy1$Habxtrt <- NA ## Blank column

## If then statement 

NACE_xy1<-NACE_xy1%>%mutate(Habxtrt= case_when(
  trt=="C" & habitat=="veg" ~ "Veg_Control",
  trt=="C" & habitat=="unveg" ~ "Unveg_Control",
  trt=="X" & habitat=="veg" ~ "Veg_Fenced",
  trt=="X" & habitat=="unveg" ~ "Unveg_Fence",
  TRUE ~ as.character(Habxtrt)
))


##### Preliminary Graphing ##### 
ggplot(NACE_xy1, aes(MDS1, MDS2, color = year)) + geom_point() + theme_bw()

nmds.NACE<-ggplot(NACE_xy1, aes(x=MDS1, y=MDS2))+
  geom_point(aes(MDS1, MDS2, colour = factor(NACE_xy1$Habxtrt), shape = factor(NACE_xy1$year)), size = 2)

nmds.NACE

### Pull the June 2009 data out from the data frame now that we have the coordinates for it. 
NACEJune_mds<- filter(NACE_xy1, (month %in% c("JUN")))
NACEJune_mds$HTY <- NA



###If then to summarize groups within June 

NACEJune_mds<-unite(NACEJune_mds,HTY,sep = "_",c(month,year,trt,habitat))


#### The code is broken here for some reason, the group_by function is not working. #################
### Create Centroids for June data. 
June_centroid<-NACEJune_mds %>% 
  group_by(HTY) %>%
  summarise(MDS1=mean(MDS1),
            MDS2=mean(MDS2),
            MDS3=mean(MDS3), .groups="drop")

June_centroid<-separate(June_centroid,col="HTY", sep = "_", into=c("month","year","trt","habitat"))

write_csv(June_centroid, file='/Users/justusjobe/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/processed data/Junecentroid.csv')




#### August Centroids #### 
NACE_xy2<-filter(NACE_xy1, (month %in% c("AUG")))
NACE_xy2$HTY <- NA

NACE_xy2<-unite(NACE_xy2,HTY,sep = "_",c(month,year,trt,habitat))

## Creation of centroids 
centroid<-NACE_xy2 %>% 
  group_by(HTY) %>%
  summarise(MDS1=mean(MDS1),
            MDS2=mean(MDS2),
            MDS3=mean(MDS3), .groups="drop")

centroid<-separate(centroid,col="HTY", sep = "_", into=c("month","year","trt","habitat"))

write_csv(centroid, file='/Users/justusjobe/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/processed data/centroid.csv')


#### Joining them back together ### 
June_centroid$Month<-"June"
centroid$Month<-"August"
Centroids_Final<-full_join(centroid, June_centroid)



### Habitat X Treatment column
Centroids_Final$Habxtrt<-NA

## If then statement 

Centroids_Final<-Centroids_Final%>%mutate(Habxtrt= case_when(
  trt=="C" & habitat=="veg" ~ "Vegetated Unfenced",
  trt=="C" & habitat=="unveg" ~ "Unvegetated Unfenced",
  trt=="X" & habitat=="veg" ~ "Vegetated Fenced",
  trt=="X" & habitat=="unveg" ~ "Unvegetated Fenced",
  TRUE ~ as.character(Habxtrt)
))

Centroids_Final$Habxtrt<-factor(Centroids_Final$Habxtrt, levels = c("Unvegetated Unfenced","Unvegetated Fenced","Vegetated Unfenced","Vegetated Fenced"))
df$var2 <- factor(df$var2, levels=c("levelB", "levelC", "levelA"))


Centroids_Final$year<-as.numeric(Centroids_Final$year)


#### I used both theme_minimal and theme_dark
ppi=600

png(filename ="figures/NMDS.png", width = 8*ppi, height= 8*ppi, res = ppi)
q<-ggplot(NACE_xy1, aes(x=MDS1, y=MDS2))+
  geom_point(shape=22, color="grey")+
  geom_point(data =Centroids_Final, mapping = aes(x=MDS1,y=MDS2, colour=year, shape=Habxtrt, size=Month)) +
  scale_colour_gradientn(colours = brewer.pal(7,"Blues"),name="Year")+
  scale_shape_manual(values=seq(15,18),name="Treatment")+scale_size_manual(values=c(4,2))+theme_dark()+theme(aspect.ratio =1/1 )

q+guides(shape = guide_legend(override.aes = list(size = 4)))

dev.off()

## Adding Significant Species Vectors 

NACE.spp.fit<-envfit(NACE.spp_NMS,NACE.spp.rel, permutations = 999)
spp.scrs<- as.data.frame(scores(NACE.spp.fit, display = "vectors"))
spp.scrs<- cbind(spp.scrs, Species= rownames(spp.scrs))
spp.scrs<- cbind(spp.scrs, pval=NACE.spp.fit$vectors$pvals)
sig.spp.scrs<- subset(spp.scrs, pval<=0.05)
sig.spp.scrs<-sig.spp.scrs[-c(14),]
sig.spp.scrs<-sig.spp.scrs[-c(6),]




#### This code works better than the one above ### 


png(filename = "figures/NMDS_sppvectors.png", width = 8*ppi, height= 8*ppi, res = ppi)
q+geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+guides(shape = guide_legend(override.aes = list(size = 4)))
dev.off()

#PERMANOVA ### 

NACE_ENV<-read.csv("~/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/raw_data/NACE_ENV.csv")
NACE_A<-NACE.spp.rel
NACE_A.dist<-vegdist(NACE_A)
attach(NACE_ENV)
NACE_ENV$year=as.numeric(NACE_ENV$year)

adonis(NACE_A.dist~trt*habitat*year, strata = NACE_ENV$module, data=NACE_ENV)






