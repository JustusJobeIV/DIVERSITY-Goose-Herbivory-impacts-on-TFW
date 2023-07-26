## Packages
library(tidyverse)
library(reshape2)
library(tibble)
library(plyr)
library(ggpubr)
library(nlme)
library(rstatix)


NACE_Herb<-read.csv("~/Desktop/R Projects/NACE PROJECT/Diversity_Manuscript_RCode/raw_data/NACE herbivory study 2009_2021 veg data 21jul21.csv")

##Clear out N/A
NACE_Herb<-na.omit(NACE_Herb)
## remove comment column
NACE_Herb<-NACE_Herb[,-9]

#### Remove data from 2016-2021####
#### Pull June out ### 
NACE_June2009_Per=filter(NACE_Herb, (month %in% c("JUN")))
NACE_June2009_Per= filter(NACE_June2009_Per,(year %in% "2009"))

#### Remove data from 2016-2021####
NACE_Herb_per = filter(NACE_Herb, (year %in% c("2009", "2010","2011","2012","2013","2014","2015")))
NACE_Herb_per =filter(NACE_Herb_per, (month %in% c("AUG")))

##Combine both dataframes
NACE_HerbCombined_per<- full_join(NACE_Herb_per, NACE_June2009_Per)


# casting data into a site by species matrix
NACE_Herb.wide_sum_percov <- dcast(NACE_HerbCombined_per,
                                   year + month + module + trt + habitat ~ spp,
                                   value.var = "percov",
                                   fun.aggregate = sum, fill = 0)

NACE_Herb_per_cov<-NACE_Herb.wide_sum_percov
## Removing undesired columns 
## Getting rid of bare
NACE_Herb_per_cov<-NACE_Herb_per_cov[,-7]
## Getting rid of detritus
NACE_Herb_per_cov<-NACE_Herb_per_cov[,-12]
### Getting rid of trash
NACE_Herb_per_cov<-NACE_Herb_per_cov[,-31]


### Summing all of the columns by row other than the columns describing the plot/treatment
NACE_Herb_per_cov$percov= rowSums(NACE_Herb_per_cov[,c(-1,-2,-3,-4,-5)])


### Habitat X Treatment column
NACE_Herb_per_cov$Habxtrt<-NA

## If then statement 

NACE_Herb_per_cov<-NACE_Herb_per_cov%>%mutate(Habxtrt= case_when(
  trt=="C" & habitat=="veg" ~ "Vegetated Unfenced",
  trt=="C" & habitat=="unveg" ~ "Unvegetated Unfenced",
  trt=="X" & habitat=="veg" ~ "Vegetated Fenced",
  trt=="X" & habitat=="unveg" ~ "Unvegetated Fenced",
  TRUE ~ as.character(Habxtrt)
))

NACE_Herb_per_cov$Year<-NA

### I do this in order to separate June 2009 and August 2009 data so that it is easier to graph using ggpubr and I can ignore the month. I will come back and fix the labeling issue in Illustrator when I am doing my final touches to the figure.
NACE_Herb_per_cov<-NACE_Herb_per_cov%>%mutate(Year= case_when(
  year=="2009" & month=="JUN" ~ "2008",
  year=="2009" & month=="AUG" ~ "2009",
  year=="2010"  ~ "2010",
  year=="2011"  ~ "2011",
  year=="2012"  ~ "2012",
  year=="2013"  ~ "2013",
  year=="2014"  ~ "2014",
  year=="2015"  ~ "2015",
  TRUE ~ as.character(Habxtrt)
))

NACE_Herb_per_cov$Year<-as.numeric(NACE_Herb_per_cov$Year)

ppi=600

png(filename = "figures/PerCov.png", width = 8*ppi, height= 8*ppi, res = ppi)
r=ggline(NACE_Herb_per_cov, x = "Year", y = "percov",color = "Habxtrt",shape="Habxtrt",add = "mean_se", palette = "Paired",point.size =1.1)
t=ggpar(r,ylab="Total Percent Cover",xlab ="Year",legend ="bottom",legend.title = "Treatments")
t+scale_y_continuous(breaks = get_breaks(by =20, from = 0),limits = c(0, 200))+scale_shape_manual(values = c(16,15,18,17))
dev.off()

NACE_Herb_per_cov$logpercov= log(NACE_Herb_per_cov$percov+1)
NACE_Herb_per_cov2<-NACE_Herb_per_cov%>%
  arrange(desc(Year), desc(module)) %>% 
  group_by(Year, module)%>%
  mutate(Diffpercov=logpercov-lag(logpercov))

NACE_Herb_per_cov2<-filter(NACE_Herb_per_cov2, Diffpercov !="NA")


NACE_Herb_per_cov2 %>% 
  group_by(habitat,year) %>% 
  identify_outliers(Diffpercov)

NACE_Herb_per_cov2 %>% 
  group_by(Year) %>%
  shapiro_test(Diffpercov)

ggqqplot(NACE_Herb_per_cov2,"Diffpercov", facet.by = "Year")

NACE_Herb_per_cov2<-ungroup(NACE_Herb_per_cov2)

NACE_Herb_per_cov2$time<-NA

NACE_Herb_per_cov2<-NACE_Herb_per_cov2%>%mutate(time= case_when(
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

NACE_Herb_per_cov2$time<-as.numeric(NACE_Herb_per_cov2$time)
NACE_Herb_per_cov2$Year<-as.factor(NACE_Herb_per_cov2$Year)
NACE_Herb_per_cov2$module<-as.factor(NACE_Herb_per_cov2$module)
NACE_Herb_per_cov2$Diffpercov<-as.numeric(NACE_Herb_per_cov2$Diffpercov)
################################
## In these models time needs to be numeric since it is the covariate. 
### Unstructured 
USpercov<-gls(Diffpercov~habitat*Year,correlation = corSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_Herb_per_cov2)
USpercov
anova(USpercov)



##Autoregressive
AGpercov<-gls(Diffpercov~habitat*Year,correlation = corAR1(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_Herb_per_cov2)
AGpercov
anova(AGpercov)

###Compound Symmetry
CSpercov<-gls(Diffpercov~habitat*Year,correlation=corCompSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time), data=NACE_Herb_per_cov2)
CSpercov
anova(CSpercov)

### Comparing the model. 
AIC(USpercov,AGpercov,CSpercov)
BIC(USpercov,AGpercov,CSpercov)



### Multiple Comparisons
library(emmeans)
Posthoc_Percov <- emmeans(AGpercov, ~ habitat | Year)
Posthoc_Percov
pairs(Posthoc_Percov)


###### T-test for the effect of caging 

library(tidyverse)
library(rstatix)

Veg_percov<-subset(NACE_Herb_per_cov,habitat=="veg")
Unveg_percov<-subset(NACE_Herb_per_cov,habitat==("unveg"))



stat.test.veg <- Veg_percov %>%
  group_by(Year) %>%
  pairwise_t_test(
    percov ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.veg




stat.test.unveg <- Unveg_percov %>%
  group_by(Year) %>%
  pairwise_t_test(
    percov ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.unveg




###Getting the Averages###

PERCOV_Mean<-NACE_Herb_per_cov %>% 
  group_by(Year,Habxtrt) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),percov )

PERCOV_DIFF_Mean<-NACE_Herb_per_cov2 %>% 
  group_by(Year,Habxtrt) %>%
  summarise_each(funs(mean,sd,se=sd(.)/sqrt(n())),Diffpercov )



### Percent Coverage of Annual Species
NACE_Herb_annual<-NACE_Herb_per_cov

NACE_Herb_annual<-NACE_Herb_annual %>%
  rowwise() %>%
  mutate(annual=sum(c(BIDCON,BIDFRO,BIDLAE,BIDSPP,ECHSPP,ECHMUR,ECLPRO,IMPCAP,PERHYD,PERPUN,PERSPP,ZIZAQU)))


NACE_Herb_annual$AP= ((NACE_Herb_annual$annual/NACE_Herb_annual$percov))

NACE_Herb_annual[is.na(NACE_Herb_annual)]<-0


#### Correct graph for ANNUALS####

png(filename="figures/Annual.png", width = 8*ppi, height= 8*ppi, res = ppi)
x=ggline(NACE_Herb_annual, x = "Year", y = "annual",color = "Habxtrt",shape  = "Habxtrt",add = "mean_se", palette = "Paired", point.size = 1.1)
v=ggpar(x,ylab="Proportion of Percent Cover Comprised of Annual Plant Species",xlab ="Year",legend ="bottom",legend.title = "Treatments")
v+scale_shape_manual(values = c(16,15,18,17))
dev.off()



##### Creating a separate data frame for t-tests #### 
NACE_Herb_annual2<-NACE_Herb_annual

### Run the stats
NACE_Herb_annual$logannual= log(NACE_Herb_annual$annual+1)
NACE_Herb_annual<-NACE_Herb_annual%>%
  arrange(desc(Year), desc(module)) %>% 
  group_by(Year, module)%>%
  mutate(Diffannual=logannual-lag(logannual))



NACE_Herb_annual<-filter(NACE_Herb_annual, Diffannual !="NA")

NACE_Herb_annual<-NACE_Herb_annual%>%mutate(Year= case_when(
  year=="2009" & month=="JUN" ~ "2008",
  year=="2009" & month=="AUG" ~ "2009",
  year=="2010"  ~ "2010",
  year=="2011"  ~ "2011",
  year=="2012"  ~ "2012",
  year=="2013"  ~ "2013",
  year=="2014"  ~ "2014",
  year=="2015"  ~ "2015",
  TRUE ~ as.character(Habxtrt)
))





NACE_Herb_annual$time<-NA

NACE_Herb_annual<-NACE_Herb_annual%>%mutate(time= case_when(
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





NACE_Herb_annual %>% 
  group_by(habitat,Year) %>% 
  identify_outliers(Diffannual)

NACE_Herb_annual %>% 
  group_by(Year) %>%
  shapiro_test(Diffannual)

ggqqplot(NACE_Herb_annual,"Diffannual", facet.by = "Year")






NACE_Herb_annual$time<-as.numeric(NACE_Herb_annual$time)
NACE_Herb_annual$Year<-as.factor(NACE_Herb_annual$Year)
NACE_Herb_annual$module<-as.factor(NACE_Herb_annual$module)


### Unstructured ### 
USannual<-gls(Diffannual~habitat*Year,correlation = corSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_Herb_annual)
USannual
anova(USannual)

##Autoregressive
AGannual<-gls(Diffannual~habitat*Year,correlation = corAR1(form=~time|module),weights = nlme::varIdent(form = ~1 | time),data=NACE_Herb_annual)
AGannual
anova(AGannual)

###Compound Symmetry
CSannual<-gls(Diffannual~habitat*Year,correlation=corCompSymm(form=~time|module),weights = nlme::varIdent(form = ~1 | time), data=NACE_Herb_annual)
CSannual
anova(CSannual)

### Comparing the model. 
AIC(USannual,AGannual,CSannual)
BIC(USannual,AGannual,CSannual)

Posthoc_annual <- emmeans(AGannual, ~ Year*habitat)
Posthoc_annual
pairs(Posthoc_annual)



Veg_annual<-subset(NACE_Herb_annual2,habitat=="veg")
Unveg_annual<-subset(NACE_Herb_annual2,habitat==("unveg"))


###### T-test for the effect of caging ###### 
stat.test.annual <- Veg_annual %>%
  group_by(Year) %>%
  pairwise_t_test(
    annual ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.annual





stat.test.annual2 <- Unveg_annual %>%
  group_by(Year) %>%
  pairwise_t_test(
    annual ~ trt, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) 
stat.test.annual2










