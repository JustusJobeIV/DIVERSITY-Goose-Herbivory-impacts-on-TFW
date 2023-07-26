library(vegan)
library(tidyverse)
library(reshape2)
library(tibble)
library(plyr)
library(ggpubr)
library(dplyr)
library(emmeans)


NACE_Herb<-read.csv("~/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/raw_data/NACE herbivory study 2009_2021 veg data 21jul21.csv")

##Clear out N/A
NACE_Herb<-na.omit(NACE_Herb)
## remove comment column
NACE_Herb<-NACE_Herb[,-9]

#### Pull June out ### 
NACE_June2009=filter(NACE_Herb, (month %in% c("JUN")))
NACE_June2009= filter(NACE_June2009,(year %in% "2009"))

#### Remove data from 2016-2021####
NACE_Herb_Div = filter(NACE_Herb, (year %in% c("2009", "2010","2011","2012","2013","2014","2015")))
NACE_Herb_Div =filter(NACE_Herb_Div, (month %in% c("AUG")))

##Combine both dataframes
NACE_HerbCombined_Div<- full_join(NACE_Herb_Div, NACE_June2009)

# casting data into a site by species matrix
NACE_Herb.sppmax <- dcast(NACE_HerbCombined_Div,
                          year + month + module + trt + habitat ~ spp,
                          value.var = "percov",
                          fun.aggregate = sum, fill = 0)

NACE_Herb.wide_sum<-NACE_Herb.sppmax

write_csv(NACE_Herb.wide_sum, file='/Users/justusjobe/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/processed data/NACE_HERB.SppMAX.csv')


NACE_Herb.sppmax<-NACE_Herb.sppmax[6:36]

#transforming to relative abundance
NACE.spp.rel.Div <-decostand(NACE_Herb.sppmax, method = "total")

NACE_Div2<-NACE.spp.rel.Div


## Adding back columns

NACE_Div2$year<-NACE_Herb.wide_sum$year
NACE_Div2$month<-NACE_Herb.wide_sum$month
NACE_Div2$module<-NACE_Herb.wide_sum$module
NACE_Div2$trt<-NACE_Herb.wide_sum$trt
NACE_Div2$habitat<-NACE_Herb.wide_sum$habitat

NACE_Div2 <- NACE_Div2 %>% relocate(year, .before = AMACAN)
NACE_Div2 <- NACE_Div2 %>% relocate(month, .before = AMACAN)
NACE_Div2 <- NACE_Div2 %>% relocate(module, .before = AMACAN)
NACE_Div2 <- NACE_Div2 %>% relocate(trt, .before = AMACAN)
NACE_Div2 <- NACE_Div2 %>% relocate(habitat, .before = AMACAN)


##### I can only have one column that is not relative abundance.

NACE_Div2$HTY <- NA
NACE_Div2 <- NACE_Div2 %>% relocate(HTY, .before = year)

NACE_Div2<-unite(NACE_Div2,HTY,sep = "_",c(month,year,module,trt,habitat))
NACE_Div2 <- NACE_Div2 %>% relocate(HTY, .before = AMACAN)

### Remove Trash, Bare, Deadwood, Detritus ###
NACE_Div2<-NACE_Div2[,-c(3,9,29)]


###Richness
RichnessF_NACE<-ddply(NACE_Div2,~HTY,function(x) {
  data.frame(RICHNESSF_NACE=sum(x[-1]>0))
})


##Shannon-Wiener Index (H')
Shannon_Wiener_NACEF<- ddply(NACE_Div2,~HTY,function(x) {
  +         data.frame(SHANNON=diversity(x[-1], index="shannon"))
})



###Merging everything together
NACE_DiversityF_1<-merge(RichnessF_NACE,Shannon_Wiener_NACEF, by.x ="HTY", by.y = "HTY" )

NACE_DiversityF_1<-separate(NACE_DiversityF_1,col="HTY", sep = "_", into=c("month","year","module","trt","habitat"))

write_csv(NACE_DiversityF_1, file='/Users/justusjobe/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/processed data/NACE_Div2.csv')


NACE_Div3<-NACE_DiversityF_1
### Habitat X Treatment column
NACE_Div3$Habxtrt<-NA

## If then statement 

NACE_Div3<-NACE_Div3%>%mutate(Habxtrt= case_when(
  trt=="C" & habitat=="veg" ~ "Vegetated Unfenced",
  trt=="C" & habitat=="unveg" ~ "Unvegetated Unfenced",
  trt=="X" & habitat=="veg" ~ "Vegetated Fenced",
  trt=="X" & habitat=="unveg" ~ "Unvegetated Fenced",
  TRUE ~ as.character(Habxtrt)
))

NACE_Div3$Year<-NA

### I do this in order to separate June 2009 and August 2009 data so that it is easier to graph using ggpubr and I can ignore the month. I will come back and fix the labeling issue in Illustrator when I am doing my final touches to the figure.
NACE_Div3<-NACE_Div3%>%mutate(Year= case_when(
  year=="2009" & month=="JUN" ~ "2008",
  year=="2009" & month=="AUG" ~ "2009",
  year=="2010" & month=="AUG" ~ "2010",
  year=="2011" & month=="AUG" ~ "2011",
  year=="2012" & month=="AUG" ~ "2012",
  year=="2013" & month=="AUG" ~ "2013",
  year=="2014" & month=="AUG" ~ "2014",
  year=="2015" & month=="AUG" ~ "2015",
  TRUE ~ as.character(Year)
))

NACE_Div3$Year<-as.numeric(NACE_Div3$Year)

### Graphing and pulling them off at a higher resolution### 
### Richness ###
ppi<-600
png(filename ="figures/RichnessWJ.png", width = 8*ppi, height= 8*ppi, res = ppi)
p=ggline(NACE_Div3, x = "Year", y = "RICHNESSF_NACE",color = "Habxtrt",shape="Habxtrt", add="mean_se", palette = "Paired", point.size = 1.1)
s=ggpar(p,ylab="Richness",xlab ="Year",legend ="bottom",legend.title = "Treatments")
s+scale_shape_manual(values = c(16,15,18,17))
dev.off()


## Shannon_Wiener 
png(filename ="figures/ShannonWJ.png", width = 8*ppi, height= 8*ppi, res = ppi)
p=ggline(NACE_Div3, x = "Year", y = "SHANNON",color = "Habxtrt",shape= "Habxtrt",add = "mean_se", palette = "Paired", point.size = 1.1)
x=ggpar(p, ylab="Shannon-Wiener Index",xlab ="Year",legend ="bottom",legend.title = "Treatments", yticks.by = 0.25,xticks.by = 1)
x+scale_shape_manual(values = c(16,15,18,17))
dev.off()



### Calculating Differences ### Arrange the rows so they line up with year and module 

##Richness
NACE_Div3$logRich= log(NACE_Div3$RICHNESSF_NACE+1)
NACE_DiffRich<-NACE_Div3 %>%
  arrange(desc(Year), desc(module), desc(trt)) %>% 
  group_by(Year, module) %>%
  mutate(DiffRichness=lag(logRich)-logRich)

## Removing rows that don't contain the difference value
NACE_DiffRich<-filter(NACE_DiffRich, DiffRichness !="NA")

### Shannon-Wiener 
NACE_Div3$logShannon= log(NACE_Div3$SHANNON+1)
NACE_DiffShannon<-NACE_Div3 %>%
  arrange(desc(Year), desc(module),desc(trt)) %>% 
  group_by(Year, module) %>%
  mutate(DiffShannon=lag(logShannon)-logShannon)

## Removing rows that don't contain the difference value
NACE_DiffShannon<-filter(NACE_DiffShannon, DiffShannon !="NA")

library(rstatix)
library(nlme)

###Richness
NACE_DiffRich %>% 
  group_by(habitat,Year) %>% 
  identify_outliers(DiffRichness)


NACE_DiffRich %>% 
  group_by(Year) %>%
  shapiro_test(DiffRichness)

ggqqplot(NACE_DiffRich,"DiffRichness", facet.by = "Year")


NACE_DiffRich$Year<-as.factor(NACE_DiffRich$Year)
NACE_DiffRich$time<-NA

NACE_DiffRich<-NACE_DiffRich%>%mutate(time= case_when(
  Year=="2008"  ~ "1",
  Year=="2009" ~ "2",
  Year=="2010" ~ "3",
  Year=="2011" ~ "4",
  Year=="2012" ~ "5",
  Year=="2013" ~ "6",
  Year=="2014" ~ "7",
  Year=="2015" ~ "8",
  TRUE ~ as.character(time)
))


NACE_DiffRich$time<-as.numeric(NACE_DiffRich$time)
NACE_DiffRich$Year<-as.factor(NACE_DiffRich$Year)
NACE_DiffRich$module<-as.numeric(NACE_DiffRich$module)





### Unstructured
USRich<-gls(DiffRichness~habitat*Year,correlation = corSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_DiffRich)
USRich
anova(USRich)

##Autoregressive
AGRich<-gls(DiffRichness~habitat*Year,correlation = corAR1(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_DiffRich)
AGRich
anova(AGRich)

##Compound Symmetry
CSRich<-gls(DiffRichness~habitat*Year,correlation=corCompSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time), data=NACE_DiffRich )
CSRich
anova(CSRich)

AIC(USRich,AGRich,CSRich)
BIC(USRich,AGRich,CSRich)




##Shannon_Wiener

NACE_DiffShannon %>% 
  group_by(habitat,year) %>% 
  identify_outliers(DiffShannon)

NACE_DiffShannon %>% 
  group_by(year) %>%
  shapiro_test(DiffShannon)

ggqqplot(NACE_DiffShannon,"DiffShannon", facet.by = "Year")



NACE_DiffShannon$year<-as.factor(NACE_DiffShannon$year)

NACE_DiffShannon$time<-NA

NACE_DiffShannon<-NACE_DiffShannon%>%mutate(time= case_when(
  Year=="2008"  ~ "1",
  Year=="2009" ~ "2",
  Year=="2010" ~ "3",
  Year=="2011" ~ "4",
  Year=="2012" ~ "5",
  Year=="2013" ~ "6",
  Year=="2014" ~ "7",
  Year=="2015" ~ "8",
  TRUE ~ as.character(time)
))

NACE_DiffShannon$time<-as.numeric(NACE_DiffShannon$time)
NACE_DiffShannon$Year<-as.factor(NACE_DiffShannon$Year)
NACE_DiffShannon$module<-as.numeric(NACE_DiffShannon$module)


### Unstructured
USShan<-gls(DiffShannon~habitat*Year,correlation = corSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_DiffShannon)
USShan
anova(USShan)

##Autoregressive
AGShannon<-gls(DiffShannon~habitat*Year,correlation = corAR1(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_DiffShannon)
AGShannon
anova(AGShannon)

##Compound Symmetry
CSShannon<-gls(DiffShannon~habitat*Year,correlation=corCompSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time), data=NACE_DiffShannon )
CSShannon
anova(CSShannon)

### Comparing the model. 
AIC(USShan,AGShannon,CSShannon)
BIC(USShan,AGShannon,CSShannon)





###### T-test for the effect of caging 
library(tidyverse)
library(rstatix)

Veg_Div<-subset(NACE_Div3,habitat=="veg")
Unveg_Div<-subset(NACE_Div3,habitat==("unveg"))


## RICHNESS

stat.test.Richness <- Veg_Div %>%
  group_by(Year) %>%
  pairwise_t_test(
    RICHNESSF_NACE ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.Richness


stat.test.RICHNESS <- Unveg_Div %>%
  group_by(Year) %>%
  pairwise_t_test(
    RICHNESSF_NACE ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.RICHNESS



##Shannon_Wiener

stat.test.Shannon <- Veg_Div %>%
  group_by(Year) %>%
  pairwise_t_test(
    SHANNON ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.Shannon


stat.test.SHANNON <- Unveg_Div %>%
  group_by(Year) %>%
  pairwise_t_test(
    SHANNON ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.SHANNON


###Getting the Averages###

NACE_RICHNESS_Mean<-NACE_Div3 %>% 
  group_by(year,Habxtrt) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),RICHNESSF_NACE )

NACE_RICHNESS_Mean2<-NACE_Div3 %>% 
  group_by(Habxtrt) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),RICHNESSF_NACE )


NACE_SW_Mean<-NACE_Div3 %>% 
  group_by(year,Habxtrt) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),SHANNON)

NACE_SW_Mean2<-NACE_Div3 %>% 
  group_by(Habxtrt) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),SHANNON)



Rich_DIFF_Mean<-NACE_DiffRich %>% 
  group_by(year,) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),DiffRichness )

SW_DIFF_Mean<-NACE_DiffShannon %>% 
  group_by(year) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),DiffShannon )









