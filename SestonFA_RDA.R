library(xlsx)
library(tidyr)
library(dplyr)
library("stringr")
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(Hmisc)
library(ggpubr)

#-----In here we try to build the relationship between 
                   # 1).environmental parameters and Seston FA and the 
                   # 2). seston FA and phytoplankton biomass. 
#(Seston FA determines by both factors; then we will do variation partitioning to identify the effect of environmental factors and phytoplankton biomass on seston FA)

#1)--------------running RDA for Seston FA composition and environmental variables--------------------

Sys.setlocale(category = "LC_ALL", locale = "greek") # To identify the Greek Symbols

#importing seston FA concentrations.

# import the file Seston_coastal and page 5

C<-read.xlsx(file.choose(),header=TRUE,5,check.names=FALSE,encoding = 'UTF-8')
C$Month = plyr::revalue(C$Month,c("may"="May","july"="July",'june'="June","august"="August",'september'="September"))
C<-filter(C, !Site=="ORE" & !Site=="OFF"  & !Site=="HORNE"& !Site=="DE")
C<-filter(C, !Station=="B")
C<-filter(C,!filtrate_sum_mg=="NA")
head(C)

#Importing physio-chemical parameters

# import the file COAST2018 Data matrix 20210629 and sheet 1

A<-read.xlsx(file.choose(),header=TRUE,1,check.names=FALSE,encoding = 'UTF-8') 
names(A)
attach(A)
sest_ENV<-filter(A,!Site=="degerfjarden" & !Site=="degerfjaerden" & !Site=="ORE" & !Site=="OFF"  & !Site=="HORNE")
sest_ENV<-na.omit(sest_ENV)
colnames(sest_ENV)
sest_ENV<-dplyr::select(sest_ENV,Site,Month,Station,Salinity,`TDP (mg.m-3)`,`TDN (mg.m-3)`,`Humic substances (mg.m-3)`,`Integrated primary production (mgC.m-2.d-1)`,`Integrated bacteria production (mgC.m-2.d-1)`,`DOC (g.m-3)`,`NH4-N (mg.m-3)`,`NO3-N (mg.m-3)`,`NO2-N (mg.m-3)`,"In situ Temp (°C)",`PO4-P (mg.m-3)`,"pH")

# change the character values in columns into numeric ones

sest_ENV[sest_ENV== "<0.6"] <- 0.6
sest_ENV[sest_ENV== "<0.9"] <- 0.9
sest_ENV[sest_ENV== "<0.4"] <- 0.4
sest_ENV[sest_ENV== "<0.3"] <- 0.3

str(sest_ENV)
sest_ENV$`NH4-N (mg.m-3)`<-as.numeric(sest_ENV$`NH4-N (mg.m-3)`)
sest_ENV$`NO3-N (mg.m-3)`<-as.numeric(sest_ENV$`NO3-N (mg.m-3)`)
sest_ENV$`PO4-P (mg.m-3)`<-as.numeric(sest_ENV$`PO4-P (mg.m-3)`)
sest_ENV$`NO2-N (mg.m-3)`<-as.numeric(sest_ENV$`NO2-N (mg.m-3)`)

str(sest_ENV)

sest_ENV<-sest_ENV%>%
  rowwise()%>%
  mutate(DIN=`NH4-N (mg.m-3)`+`NO3-N (mg.m-3)`,`NO2-N (mg.m-3)`)

sest_ENV # environmental variables
C        # seston FA

#adding above to merged dataframe

Seston_FA_ENV<-merge(sest_ENV,(C%>%select(Month,Site,Station,EPA:`??6 PUFA`)),by=c("Month","Site","Station")) # merge seston FA and env.
Seston_FA_ENV

Seston_FA_ENV[Seston_FA_ENV== 0] <- 0.0000001 # replace 0 values to small value
Seston_FA_ENV<-select(Seston_FA_ENV,Site,Station,Month,Salinity,`Humic substances (mg.m-3)`,DIN, `DOC (g.m-3)`,`In situ Temp (°C)`,`PO4-P (mg.m-3)`,DIN: `??6 PUFA`)
log_Seston_FA_ENV<-apply(Seston_FA_ENV[,4:16],2, function(x) (log10(x)))

#binding the site, station and Month

log_Seston_FA_ENV<-as.data.frame(log_Seston_FA_ENV)
log_Seston_FA_ENV$Site<-Seston_FA_ENV$Site
log_Seston_FA_ENV$Station<-Seston_FA_ENV$Station
log_Seston_FA_ENV$Month<-Seston_FA_ENV$Month
head(log_Seston_FA_ENV)

##################################################################################################################################################################################################################

# Running the RDA#
# before RDA correlation check of explanatory variables; variables with high correlations should be removed. 

#--------------correlation check to decide which variables to retain--------------

cor_data_seston_FA_ENV = cor((as.matrix(log_Seston_FA_ENV%>%select(Salinity:`PO4-P (mg.m-3)`))))
cor_data_seston_FA_ENV # display correlation coefficient matrix
corrplot::corrplot(cor_data_seston_FA_ENV , method="number",is.corr = TRUE) #visualize correlation matrix
corrplot::corrplot(cor_data_seston_FA_ENV *cor_data_seston_FA_ENV,method='number') # display R2 value between variables

# as per the results, DOC salinity and the Humic substances correlate each other. Therefore remove salinity and humic substances.

# Environmental factors matrix

Seston_FA_EV_dataframe<-(log_Seston_FA_ENV%>%filter(!Month=="May" & !Month=="June")%>%select(DIN:`PO4-P (mg.m-3)`)) #after first run I removed Salinity, DOC and Humic substances that cause multicolinearily by VIF and corr plot
names(Seston_FA_EV_dataframe)<-c("DIN","DOC","Temperature(°C)","Phospate")

# FA matrix

Seston_FA_FA_dataframe<-(log_Seston_FA_ENV%>%filter(!Month=="May" & !Month=="June")%>%select(EPA:`??6 PUFA`)) #if want selected month: (logmerged%>%filter(!Month=="May")%>%select("EPA":"PUFA"))


###running global model with all variables for seston FA and environmental parametes#######

tb_rda.EV_seston_FA<- rda (Seston_FA_FA_dataframe~.,Seston_FA_EV_dataframe,scale = TRUE)

# check their significance; three methods exists; 

#i) by global model
anova (tb_rda.EV_seston_FA)

#ii) by term
anova (tb_rda.EV_seston_FA,by="term")

#iii) by axis
anova (tb_rda.EV_seston_FA,by="axis")

# check summary
summary(tb_rda.EV_seston_FA)
summaryEV<-summary(tb_rda.EV_seston_FA)

tb_rda.EV_seston_FA$call


sitescorestb_rda.EV_seston_FA<-as.data.frame(tb_rda.EV_seston_FA$CCA$u) # This is location data; u is site data
speciescorestb_rda.EV_seston_FA<-as.data.frame(tb_rda.EV_seston_FA$CCA$v) # This is species data; v is species data; https://rdrr.io/rforge/vegan/man/cca.object.html
envscoresEVtb_rda.EV_seston_FA<-tb_rda.EV_seston_FA$CCA$biplot              # This is environment scores

#attaching location and month data for both location and species scores
sitescorestb_rda.EV_seston_FA$site<-(log_Seston_FA_ENV%>%filter(!Month=="May" & !Month=="June"))$Site #if want selected site: (logmerged%>%filter(!Month=="May"))$Site
sitescorestb_rda.EV_seston_FA$Month<-(log_Seston_FA_ENV%>%filter(!Month=="May" & !Month=="June"))$Month # if want selected month: (logmerged%>%filter(!Month=="May"))$Month



# plotting the ggplot RDA Plot

SesontFA_ENV_RDA_GRAPH<-ggplot(sitescorestb_rda.EV_seston_FA,aes(x=RDA1,y=RDA2)) + 
  geom_point(aes(fill=Month,
                 shape=site),
             size=4) +
  scale_fill_brewer(palette="Paired")+
  scale_shape_manual(values = c(21,22,25,24))+ # use fill able shapes 21-25
  guides(shape = guide_legend(title = "Site"),
         fill= guide_legend(title = "Month",override.aes=list(shape=21)))+
  geom_segment(data = as.data.frame(envscoresEVtb_rda.EV_seston_FA), 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
               alpha = 0.75, color = 'blue')+
  #adding phytoplankton data as vectors
  geom_text(data = as.data.frame(envscoresEVtb_rda.EV_seston_FA), 
            aes(x = RDA1, y = RDA2, 
                label = rownames(envscoresEVtb_rda.EV_seston_FA)), 
            col = 'red',
            size=4,
            vjust="inward",
            hjust='inward')+
  #geom_point(data = speciescoresEV, 
  #aes(x = RDA1, y = RDA2), 
  #col = 'red')+
  geom_segment(data = as.data.frame(speciescorestb_rda.EV_seston_FA), 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.2,"cm"),type = "closed"),
               alpha = 0.75, color = 'red',size=0.1)+
  ggrepel::geom_text_repel(data = as.data.frame(speciescorestb_rda.EV_seston_FA), 
                           aes(x = RDA1, y = RDA2, 
                               label = rownames(speciescorestb_rda.EV_seston_FA)), 
                           col = 'black',
                           size=4) +
  labs(shape="Month", colour="Site")+
  labs(x=paste("RDA 1 (", format(100 *summaryEV$cont[[1]][2,1], digits=4), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *summaryEV$cont[[1]][2,2], digits=4), "%)", sep=""))+
  theme_minimal()+
  theme(axis.text.y   = element_text(size=12),
        axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = 'top',
        panel.border = element_rect(colour = "black", fill=NA, size=1))

SesontFA_ENV_RDA_GRAPH


###########################################################################################################
#2) Understand the relationship between Seston FA and PB biomass
###########################################################################################################

#importing phytoplankton dataframe

# import COAST2018 Phytoplankton taxonomy data TBH file

PB<-read.xlsx(file.choose(),header=TRUE,1,check.names=FALSE,encoding = 'UTF-8')
attach(PB)
names(PB)
head(PB)
levels(factor(PB$ClassName))
PB<-filter(PB, !PB$Site=="ORE" & !PB$Site=="OFF"  & !PB$Site=="HORNE")
head(PB)

# calculating sumbiomass

Biomass<-PB%>%
  dplyr::group_by(Month,Site,Station,ClassName)%>% 
  summarise(sumbiomass=sum(`ug C/m3`)) # this backticks occur since there is space between words
attach(Biomass)
names(Biomass)
levels(factor(Biomass$ClassName))
Biomass<-filter(Biomass,!ClassName=="No Class")
levels(factor(Biomass$ClassName))
Biomass<-Biomass%>%
  pivot_wider(names_from=ClassName,values_from = sumbiomass,values_fill=0)
Biomass

#Few Classes not found abundatntly; so we have to remove them
Biomass<-select(Biomass,-c("EUGLENOPHYCEAE","CHAROPHYCEAE","EBRIIDEA","LITOSTOMATEA"))#; Litostomatea is cilliate; Ebriddea is type of phytoplankton; Charophyceae not present almost all stations; Euglenophyceace present in all stations, in later months.


#importing seston FA dataframe, which is C.

C
Biomass

#adding above to merged dataframe
Seston_FA_phyto<-merge(Biomass,(C%>%select(Month,Site,Station,EPA:`??6 PUFA`)),by=c("Month","Site","Station")) # merge seston FA and PB
Seston_FA_phyto

Seston_FA_phyto[Seston_FA_phyto== 0] <- 0.0000001
Seston_FA_phyto<-select(Seston_FA_phyto,Site,Station,Month,CHLOROPHYCEAE:`PRASINOPHYCEAE (MICROMONADOPHYCEAE)`,EPA: `??6 PUFA`)
log_Seston_FA_phyto<-apply(Seston_FA_phyto[,4:18],2, function(x) (log10(x)))

#binding the site, station and Month

log_Seston_FA_phyto<-as.data.frame(log_Seston_FA_phyto)
log_Seston_FA_phyto$Site<-Seston_FA_phyto$Site
log_Seston_FA_phyto$Station<-Seston_FA_phyto$Station
log_Seston_FA_phyto$Month<-Seston_FA_phyto$Month
head(log_Seston_FA_phyto)

Seston_FA_PB_dataframe<-(log_Seston_FA_phyto%>%filter(!Month=="May" & !Month=="June")%>%
                           select(CHLOROPHYCEAE:`PRASINOPHYCEAE (MICROMONADOPHYCEAE)`))
names(Seston_FA_PB_dataframe)<-c("Chlorophytes","Chrysophytes","Cryptophytes","Diatoms","Dinophytes","Cyanophytes","Prymnesiophytes","Prasinophytes")
Seston_FA_FA_phyto_dataframe<-(log_Seston_FA_phyto%>%filter(!Month=="May" & !Month=="June")%>%select(EPA:`??6 PUFA`))


# run global model with all variables
tb_rda.seston_FA_phyto<- rda (Seston_FA_FA_phyto_dataframe~.,Seston_FA_PB_dataframe)

# check their significance; three methods exists; 

#i) by global model
anova (tb_rda.seston_FA_phyto)

#ii) by term
anova (tb_rda.seston_FA_phyto,by="term")

#iii) by axis
anova (tb_rda.seston_FA_phyto,by="axis")

# check summary
summary(tb_rda.seston_FA_phyto)

summary_phyto<-summary(tb_rda.seston_FA_phyto)

tb_rda.seston_FA_phyto$call

#mulcolinearity of the variables can be checked by using vif factors. if sqrt(vit)>2 high multicolinearity. We can remove those variables and re-run the analysis again. 
sqrt(vif.cca(tb_rda.seston_FA_phyto))

#explained varinace can be searched by adjusted R2; since we have more variables it should be adjusted. 
adjR2.tb_rda.seston_FA_phyto <- RsquareAdj (tb_rda.seston_FA_phyto)$adj.r.squared
adjR2.tb_rda.seston_FA_phyto

#Preparation of data for plotting. 

sitescorestb_rda.seston_FA_phyto<-as.data.frame(tb_rda.seston_FA_phyto$CCA$u) # This is location data; u is site data
speciescorestb_rda.seston_FA_phyto<-as.data.frame(tb_rda.seston_FA_phyto$CCA$v) # This is species data; v is species data; https://rdrr.io/rforge/vegan/man/cca.object.html
envscoresEVtb_rda.seston_FA_phyto<-tb_rda.seston_FA_phyto$CCA$biplot              # This is environment scores

#attaching location and month data for both location and species scores
sitescorestb_rda.seston_FA_phyto$site<-(log_Seston_FA_phyto%>%filter(!Month=="May" & !Month=="June"))$Site #if want selected site: (logmerged%>%filter(!Month=="May"))$Site
sitescorestb_rda.seston_FA_phyto$Month<-(log_Seston_FA_phyto%>%filter(!Month=="May" & !Month=="June"))$Month # if want selected month: (logmerged%>%filter(!Month=="May"))$Month



# plotting the ggplot

SesontFA_Phyto_RDA_GRAPH<-ggplot(sitescorestb_rda.seston_FA_phyto,aes(x=RDA1,y=RDA2)) + 
  geom_point(aes(fill=Month,
                 shape=site),
                 size=4) +
  scale_fill_brewer(palette="Paired")+
  scale_shape_manual(values = c(21,22,25,24))+ # use fill able shapes 21-25
  guides(shape = guide_legend(title = "Site"),
         fill= guide_legend(title = "Month",
                            override.aes=list(shape=21)))+
  geom_segment(data = as.data.frame(envscoresEVtb_rda.seston_FA_phyto), 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
               color = 'red')+
  #adding phytoplankton class names as vectors
  geom_text(data = as.data.frame(envscoresEVtb_rda.seston_FA_phyto), 
            aes(x = RDA1, y = RDA2, 
                label = rownames(envscoresEVtb_rda.seston_FA_phyto)), 
            col = 'red',
            size=4,
            vjust="inward",
            hjust='inward')+
  #geom_point(data = speciescoresEV, 
  #aes(x = RDA1, y = RDA2), 
  #col = 'red')+
  #adding fatty acid vectors
  geom_segment(data = as.data.frame(speciescorestb_rda.seston_FA_phyto), 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
               color = 'black')+
  #adding fatty acid names
  ggrepel::geom_text_repel(data = as.data.frame(speciescorestb_rda.seston_FA_phyto), 
                           aes(x = RDA1, y = RDA2, 
                               label = rownames(speciescorestb_rda.seston_FA_phyto)), 
                           col = 'black',
                           size=4) +
  labs(shape="Month", colour="Site")+
  labs(x=paste("RDA 1 (", format(100 *summary_phyto$cont[[1]][2,1], digits=3), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *summary_phyto$cont[[1]][2,2], digits=3), "%)", sep=""))+
  theme_classic()+
  theme(axis.text.y   = element_text(size=12),
        axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=12),
        axis.title.x  = element_text(size=12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        #legend.position = 'top',
        panel.border = element_rect(colour = "black", fill=NA, size=1))
SesontFA_Phyto_RDA_GRAPH


# -------------------------Doing the variation partitioning of the data------------------------------------------------

# doing variation partitioning for checking Environmental parameters or phyto-composition determines the variability.


#We have different number of rows for seston FA, phyto-biomass and environmental variables. We will merge these together and Run varpart

Seston_FA_ENV # 52 observations and 16 variables.
Seston_FA_phyto # 53 observations and 18 variables.
C #  53 observations and 61 variables
Merged_partioning<-merge(Seston_FA_phyto%>%select(Site,Month,Station,"CHLOROPHYCEAE":`PRASINOPHYCEAE (MICROMONADOPHYCEAE)`), (C%>%select(Site,Month,Station,EPA:'??6 PUFA')),by=c("Site","Month","Station")) %>%
                   merge((Seston_FA_ENV%>%select(Site,Station,Month,DIN:`PO4-P (mg.m-3)`)),by=c("Site","Month","Station"))



VP<- vegan::varpart ((log10(Merged_partioning[,c(12:18) ])), (log10(Merged_partioning[,c(4:11) ])) , (log10(Merged_partioning[,c(19:22) ])) ,scale = TRUE)
VP$part

#log10(Merged_partitioning[,c(12:18) ]=phytoplankton
#(log10(Merged_partitioning[,c(19:22) ]))= Environmental

plot(VP,
     Xnames = c("Phytoplankton", "Environmental"), # name the partitions
     bg = c("seagreen3", "mediumpurple"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)

#plot cant combine with ggarrange; therefore we have to customize it as below and add to ggarrange. 

# ggforce to draw ven diagrams
library(ggforce)# ggforce to draw ven diagrams
df.venn <- data.frame(x = c(3, 1),y = c(1, 1),labels = c('2.8%', '20%'))
p <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, 
                        fill = df.venn$labels)) +
                        geom_circle(alpha = .5, size = 1, colour = 'black',show.legend = FALSE ) +
                        coord_fixed()+
                        annotate("text", x = df.venn$x , y = df.venn$y,label=df.venn$labels ,size = 5)+
                        annotate("text", x = 2 , y =1,label="15.8%" ,size = 5)+
                        annotate("text", x = 4 , y =-0.5,label="Residual=61%" ,size = 5)+
                        annotate("text", x = 0 , y =2.5,label="Phytoplankton" ,size = 5)+
                        annotate("text", x = 4 , y =2.5,label="Environmental variables" ,size = 5)+
                        theme_classic()+
                        theme(axis.text.y   = element_blank(),
                              axis.text.x   = element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.ticks.y=element_blank(),
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              panel.border = element_rect(colour = "black", fill=NA, size=1))
p

# arranging three plots in one file

Final_RDA_plots<-ggpubr::ggarrange(SesontFA_Phyto_RDA_GRAPH, SesontFA_ENV_RDA_GRAPH,
                      common.legend = TRUE,
                      labels = c("A", "B"),
                      nrow = 1,ncol=2)
Final_RDA_plots
ggsave("Final_RDA_plots.eps",width = 140, height = 70,units = "mm",scale = 1.9,device=cairo_ps)

# saving plot in elsevier format; not done here, it cause absurd format. We will do in photoshop.



###########################################################################################################

#3). The other relationship we can try is running the RDA for phytoplankton composition and environmental parameters. 
##########################################################################################################

#--------------------#running RDA for phytoplankton composition and environmental variables----------

#Importing physio-chemical parameters

# import the file COAST2018 Data matrix 20210629

A<-read.xlsx(file.choose(),header=TRUE,1,check.names=FALSE,encoding = 'UTF-8') 
names(A)
attach(A)
D<-filter(A,!Site=="degerfjarden" & !Site=="degerfjaerden" & !Site=="ORE" & !Site=="OFF"  & !Site=="HORNE")
D<-na.omit(D)
EV<-dplyr::select(D,Site,Month,Station,Salinity,`TDP (mg.m-3)`,`TDN (mg.m-3)`,`Humic substances (mg.m-3)`,`Integrated primary production (mgC.m-2.d-1)`,`Integrated bacteria production (mgC.m-2.d-1)`,`DOC (g.m-3)`,`NH4-N (mg.m-3)`,`NO3-N (mg.m-3)`,"In situ Temp (°C)",`PO4-P (mg.m-3)`,"pH") #Sonia told remove TDP and TDN since they are not much useful in the phyto, here included we can remove in later selection.
EV
colnames(EV)
D2<-dplyr::select(D,Site,Month,Station,Salinity,`TDP (mg.m-3)`,`TDN (mg.m-3)`,`Humic substances (mg.m-3)`,`Integrated primary production (mgC.m-2.d-1)`,`Integrated bacteria production (mgC.m-2.d-1)`,`DOC (g.m-3)`,`NH4-N (mg.m-3)`,`NO3-N (mg.m-3)`,`NO2-N (mg.m-3)`,"In situ Temp (°C)",`PO4-P (mg.m-3)`,"pH")

#some of the environmental parameters were below limits; we have to change those. 

D2[D2== "<0.6"] <- 0.6
D2[D2== "<0.9"] <- 0.9
D2[D2== "<0.4"] <- 0.4
D2[D2== "<0.3"] <- 0.3

#check whether those are numeric
str(D2)
# change chr columns to numeric

D2$`NH4-N (mg.m-3)`<-as.numeric(D2$`NH4-N (mg.m-3)`)
D2$`NO3-N (mg.m-3)`<-as.numeric(D2$`NO3-N (mg.m-3)`)
D2$`PO4-P (mg.m-3)`<-as.numeric(D2$`PO4-P (mg.m-3)`)
D2$`NO2-N (mg.m-3)`<-as.numeric(D2$`NO2-N (mg.m-3)`)

str(D2)

D2<-D2%>%
  rowwise()%>%
  mutate(DIN=`NH4-N (mg.m-3)`+`NO3-N (mg.m-3)`,`NO2-N (mg.m-3)`)

D2


#Importing phytoplankton data

# import the file COAST2018 Phytoplankton taxonomy data TBH

PB<-read.xlsx(file.choose(),header=TRUE,1,check.names=FALSE,encoding = 'UTF-8')
attach(PB)
names(PB)
head(PB)
levels(factor(PB$ClassName))
PB<-filter(PB, !PB$Site=="ORE" & !PB$Site=="OFF"  & !PB$Site=="HORNE")
head(PB)

# calculating sumbiomass
Biomass<-PB%>%
  dplyr::group_by(Month,Site,Station,ClassName)%>% 
  summarise(sumbiomass=sum(`ug C/m3`)) # this backticks occur since there is space between words
attach(Biomass)
names(Biomass)
levels(factor(Biomass$ClassName))
Biomass<-filter(Biomass,!ClassName=="No Class")
levels(factor(Biomass$ClassName))
Biomass<-Biomass%>%
  pivot_wider(names_from=ClassName,values_from = sumbiomass,values_fill=0)
Biomass
#Few Classes not found abundatntly; so we have to remove them
Biomass<-select(Biomass,-c("EUGLENOPHYCEAE","CHAROPHYCEAE","EBRIIDEA","LITOSTOMATEA"))#; Litostomatea is cilliate; Ebriddea is type of phytoplankton; Charophyceae not present almost all stations; Euglenophyceace present in all stations, in later months.
names(Biomass)<-c("Month","Site","Station","Chlorophytes","Chrysophytes","Cryptophytes","Diatoms","Dinophytes","Cyanophytes","Prymnesiophytes","Prasinophytes")

MergedEV_phyto<-merge(Biomass,D2,by=c("Month","Site","Station"))
MergedEV_phyto
attach(MergedEV_phyto)


#log conversion of all data

MergedEV_phyto[MergedEV_phyto== 0] <- 0.0000001
logmergedEV_phyto<-apply(MergedEV_phyto[,4:25],2, function(x) (log10(x)))
logmergedEV_phyto<-as.data.frame(logmergedEV_phyto)
logmergedEV_phyto$Site<-MergedEV_phyto$Site
logmergedEV_phyto$Station<-MergedEV_phyto$Station    
logmergedEV_phyto$Month<-MergedEV_phyto$Month
head(logmergedEV_phyto)



#--------------correlation check to decide which variables to retain--------------

cor_data = cor((as.matrix(logmergedEV_phyto%>%select(Salinity:DIN))))
cor_data # display correlation coefficient matrix
corrplot::corrplot(cor_data, method="number",is.corr = TRUE) #visualize correlation matrix
corrplot::corrplot(cor_data*cor_data,method='number') # display R2 value between variables
p_values<-rcorr (as.matrix(logmergedEV_phyto%>%select(Salinity:`Integrated bacteria production (mgC.m-2.d-1)`)))
print(p_values)

# we can do easily in excel by using correlation in data analysis

#it seems that there is high correlation between salinity; humic substances and the DOC; therefore I will remove the salinity and humic substances and keep DOC in the model

EV_phyto_dataframe<-((logmergedEV_phyto%>%filter(!Month=="May" & !Month=="June"))%>%select(`DOC (g.m-3)`,`In situ Temp (°C)`,`PO4-P (mg.m-3)`:DIN)) #after first run I removed Salinity, DOC and Humic substances that cause multicolinearily by VIF and corr plot

phytoplankton_dataframe<-((logmergedEV_phyto%>%filter(!Month=="May" & !Month=="June"))%>%select("Chlorophytes":"Prasinophytes")) #if want selected month: (logmerged%>%filter(!Month=="May")%>%select("EPA":"PUFA"))

# run global model with all variables
tb_rda.EV_phyto<- rda (phytoplankton_dataframe~.,EV_phyto_dataframe)
# check their significance; three methods exists; i) by global model
anova (tb_rda.EV_phyto)
#ii) by term
anova (tb_rda.EV_phyto,by="term")
#iii) by axis
anova (tb_rda.EV_phyto,by="axis")
# check summary
summary(tb_rda.EV_phyto)

tb_rda.EV_phyto$call

#mulcolinearity of the variables can be checked by using vif factors. if sqrt(vit)>2 high multicolinearity. We can remove those variables and re-run the analysis again. 
sqrt(vif.cca(tb_rda.EV_phyto))

#explained varinace can be searched by adjusted R2; since we have more variables it should be adjusted. 
adjR2.tbrdaEV_phyto <- RsquareAdj (tb_rda.EV_phyto)$adj.r.squared
adjR2.tbrdaEV_phyto

# since the model consist of high multicolinearity variables we have removed it DOC; Salinity, and Humic substances and re-ran the model

summaryEV_phyto<-summary(tb_rda.EV_phyto)

sitescoresEV_phyto<-as.data.frame(tb_rda.EV_phyto$CCA$u) # This is location data; u is site data
speciescoresEV_phyto<-as.data.frame(tb_rda.EV_phyto$CCA$v) # This is species data; v is species data; https://rdrr.io/rforge/vegan/man/cca.object.html
envscoresEV_phyto<-tb_rda.EV_phyto$CCA$biplot              # This is environment scores

#attaching location and month data for both location and species scores
sitescoresEV_phyto$site<-(logmergedEV_phyto%>%filter(!Month=="May" & !Month=="June"))$Site #if want selected site: (logmerged%>%filter(!Month=="May"))$Site
sitescoresEV_phyto$Month<-(logmergedEV_phyto%>%filter(!Month=="May" & !Month=="June"))$Month # if want selected month: (logmerged%>%filter(!Month=="May"))$Month

# plotting the ggplot

Phyto_env_RDA_graph<-ggplot(sitescoresEV_phyto,aes(x=RDA1,y=RDA2)) + 
                     geom_point(aes(fill=Month,
                                    shape=site),
                                    size=4) +
  scale_fill_brewer(palette="Paired")+
  scale_shape_manual(values = c(21,22,25,24))+ # use fill able shapes 21-25
  guides(shape = guide_legend(title = "Site"),
         fill= guide_legend(title = "Month",
                            override.aes=list(shape=21)))+
  #adding Environmental variables text 
  ggrepel::geom_text_repel(data = as.data.frame(envscoresEV_phyto), 
            aes(x = RDA1, y = RDA2, 
                label = rownames(envscoresEV_phyto)), 
                col = 'black')+
  #adding environmental variables vectors
  geom_segment(data = as.data.frame(envscoresEV_phyto), 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
               color = 'black')+
  #adding phytoplankton class variables text
  ggrepel::geom_text_repel(data = as.data.frame(speciescoresEV_phyto), 
                           aes(x = RDA1, y = RDA2, 
                               label = rownames(speciescoresEV_phyto)), 
                           col = 'red') +
  #geom_point(data = speciescoresEV, 
             #aes(x = RDA1, y = RDA2), 
             #col = 'red')+
  #adding phytoplankton classes as vectors
  geom_segment(data = as.data.frame(speciescoresEV_phyto), 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
               color = 'red')+
  labs(shape="Month", colour="Site")+
  labs(x=paste("RDA 1 (", format(100 *summaryEV_phyto$cont[[1]][2,1], digits=3), "%)", sep=""),
       y=paste("RDA 2 (", format(100 *summaryEV_phyto$cont[[1]][2,2], digits=3), "%)", sep=""))+
  theme_classic()+   
  theme(axis.text.y   = element_text(size=12),
                          axis.text.x   = element_text(size=12),
                          axis.title.y  = element_text(size=12),
                          axis.title.x  = element_text(size=12),
                          panel.background = element_blank(),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          #legend.position = 'top',
                          panel.border = element_rect(colour = "black", fill=NA, size=1))

arrange<-ggarrange(Phyto_env_RDA_graph,SesontFA_Phyto_RDA_GRAPH,common.legend = TRUE,legend = 'bottom',labels = c('a','b'))
ggsave("arrangedplot.eps", arrange,width = 10,height = 6,units = 'in',device = cairo_ps)




#####################################################################################################################################
#           This section is deprecated; it includes code chunks to reproduce forward selection methods. IF want we can select from this.CTRL+SHIFT+C TO BULK CODE INACTIVAE ANDDEACTIVATE#####################################################################################################################

#--------------RDA with forward selection method-----------------------------------------------------------

# #Since our case global model was significant; we could run the forward selection model
# # To run stepwise regression-the global model should be significant. If not we couldn't run it. We have to stick to the global model with all variables. 
# #full model
# # FA_dataframe2<-(logmergedEV%>%select("EPA":`??6 PUFA`)) #if want selected month: (logmerged%>%filter(!Month=="May")%>%select("EPA":"PUFA"))
# # full_model_RDA_EV <- rda (FA_dataframe2 ~ 1, data = EV_dataframe,scale = TRUE) # No constraints; only intercept; this is equal to pca
# #full_model_RDA
# Model_all_variable_EV<- rda (FA_dataframe2 ~ ., data = EV_dataframe,scale = TRUE) # with all variables/ constraints
# sel.osR2_EV <- ordiR2step (full_model_RDA_EV, scope = formula (Model_all_variable_EV), R2scope = adjR2.tbrdaEV, direction = 'forward', permutations = 999)
# anova(sel.osR2_EV ) # To check whether model is significant; out from the result of forward selection; go ahead with those variables and run RDA again. 
# 
# 
# #odel_all_variable12_EV<- rda (FA_dataframe2 ~`Bacteria production (mgC.m-3.d-1)`, data = EV_dataframe, scale = TRUE)
# #odel_all_variable12_EV
# 
# anova (sel.osR2_EV )
# #ii) by term
# anova (sel.osR2_EV ,by="term")
# #iii) by axis
# anova (sel.osR2_EV ,by="axis")
# 
# 
# RsquareAdj (sel.osR2_EV)$adj.r.squared
# sqrt(vif.cca(sel.osR2_EV))
# anova (odel_all_variable12_EV)
# 
# 
# #run summary of the model
# summary_forward_EV<-summary(sel.osR2_EV)
# summary_forward_EV
# 
# #Extracting components for drawing RDA plots
# 
# sitescores_forward_EV<-as.data.frame(sel.osR2_EV$CCA$u) # This is location data; u is site data
# speciescores_forward_EV<-as.data.frame(sel.osR2_EV$CCA$v) # This is species data; v is species data; https://rdrr.io/rforge/vegan/man/cca.object.html
# envscores_forward_EV<-sel.osR2_EV$CCA$biplot              # This is environment scores
# 
# #attaching location and month data for both location and species scores
# sitescores_forward_EV$site<-logmergedEV$Site #if want selected site: (logmerged%>%filter(!Month=="May"))$Site
# sitescores_forward_EV$Month<-logmergedEV$Month # if want selected month: (logmerged%>%filter(!Month=="May"))$Month
# 
# attach(sitescores_forward_EV)
# attach(speciescores_forward_EV)
# 
# 
# # plotting the ggplot
# ggplot(sitescores_forward_EV) + 
#   geom_point(aes(x = RDA1, y = RDA2, 
#                  col = factor(site),
#                  shape = factor(Month)),
#              size=5,
#              alpha=0.5) +
#   #adding phytoplankton data as vectors
#   geom_text(data = as.data.frame(envscores_forward_EV), 
#             aes(x = RDA1, y = RDA2, 
#                 label = rownames(envscores_forward_EV)), 
#             col = 'red') +
#   geom_segment(data = as.data.frame(envscores_forward_EV), 
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
#                alpha = 0.75, color = 'blue')+
#   geom_point(data = speciescores_forward_EV, 
#              aes(x = RDA1, y = RDA2), 
#              col = 'red')+
#   geom_segment(data = as.data.frame(speciescores_forward_EV), 
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow=arrow(length=unit(0.2,"cm"),type = "closed"),
#                alpha = 0.75, color = 'red',size=0.1)+
#   ggrepel::geom_text_repel(data = as.data.frame(speciescores_forward_EV), 
#                            aes(x = RDA1, y = RDA2, 
#                                label = rownames(speciescores_forward_EV)), 
#                            col = 'black') +
#   labs(shape="Month", colour="Site")+
#   labs(x=paste("RDA 1 (", format(100 *summary_forward_EV$cont[[1]][2,1], digits=4), "%)", sep=""),
#        y=paste("RDA 2 (", format(100 *summary_forward_EV$cont[[1]][2,2], digits=4), "%)", sep=""))+
#   theme_bw() 
# 
# #deprecated...
# 
# # Importing physio-chemical parameters
# 
# A<-read.xlsx(file.choose(),header=TRUE,1,check.names=FALSE,encoding = 'UTF-8') #physio-chemical parameters
# names(A)
# attach(A)
# B<-filter(A,!Site=="degerfjarden" & !Site=="degerfjaerden" & !Site=="ORE" & !Site=="OFF"  & !Site=="HORNE")
# B<-na.omit(B)
# B<-select(B,Month,Site,Station,`DOC (g.m-3)`,`Humic substances (mg.m-3)`)
# head(B)
# 
# #adding above to merged dataframe
# Merged<-merge(Merged,B,by=c("Month","Site","Station"))
# Merged
# 
# 
# #mulcolinearity of the variables can be checked by using vif factors. if sqrt(vit)>2 high multicolinearity. We can remove those variables and re-run the analysis again. 
# sqrt(vif.cca(tb_rda.EV_seston_FA))
# 
# #explained varinace can be searched by adjusted R2; since we have more variables it should be adjusted. 
# adjR2.tb_rda.EV_seston_FA <- RsquareAdj (tb_rda.EV_seston_FA)$adj.r.squared
# adjR2.tb_rda.EV_seston_FA
# 
# # since the model consist of high multicolinearity variables we have removed it DOC; Salinity, and Humic substances and re-ran the model
# 
# summaryEV<-summary(tb_rda.EV_seston_FA)
# 
# 
# #1)....Try to identify the relationship between Seston FA and phytoplankton biomass------------------------
# 
# 
# 
# #Importing FA values of seston
# 
# C<-read.xlsx(file.choose(),header=TRUE,5,check.names=FALSE,encoding = 'UTF-8')
# C$Month = plyr::revalue(C$Month,c("may"="May","july"="July",'june'="June","august"="August",'september'="September"))
# C<-filter(C, !Site=="ORE" & !Site=="OFF"  & !Site=="HORNE"& !Site=="DE")
# C<-filter(C, !Station=="B")
# C<-filter(C,!filtrate_sum_mg=="NA")
# head(C)
# 
# 
# #Importing phytoplankton data
# 
# PB<-read.xlsx(file.choose(),header=TRUE,1,check.names=FALSE,encoding = 'UTF-8')
# attach(PB)
# names(PB)
# head(PB)
# levels(factor(PB$ClassName))
# PB<-filter(PB, !PB$Site=="ORE" & !PB$Site=="OFF"  & !PB$Site=="HORNE")
# head(PB)
# 
# # calculating sumbiomass
# Biomass<-PB%>%
#   dplyr::group_by(Month,Site,Station,ClassName)%>% 
#   summarise(sumbiomass=sum(`ug C/m3`)) # this backticks occur since there is space between words
# attach(Biomass)
# names(Biomass)
# levels(factor(Biomass$ClassName))
# Biomass<-filter(Biomass,!ClassName=="No Class")
# levels(factor(Biomass$ClassName))
# Biomass<-Biomass%>%
#   pivot_wider(names_from=ClassName,values_from = sumbiomass,values_fill=0)
# Biomass
# #Few Classes not found abundatntly; so we have to remove them
# Biomass<-select(Biomass,-c("EUGLENOPHYCEAE","CHAROPHYCEAE","EBRIIDEA","LITOSTOMATEA"))#; Litostomatea is cilliate; Ebriddea is type of phytoplankton; Charophyceae not present almost all stations; Euglenophyceace present in all stations, in later months.
# 
# 
# #Merging two dataframes; seston FA and phyto biomass
# 
# Merged<-merge(C,Biomass,by=c("Month","Site","Station"))
# Merged
# attach(Merged)
# Merged<-select(Merged,-c("filter_ID":"Type"))
# 
# #log conversion of all data
# 
# Merged[Merged== 0] <- 0.0000001
# logmerged<-apply(Merged[,4:49],2, function(x) (log10(x)))
# 
# 
# #z score standaridazation; we can go ahead with this if we opt not to use scale function later. 
# #mvals<-apply(logmerged, 2, mean)
# #sdvals<-apply(logmerged, 2, sd)
# #logmerged<-(logmerged - mvals) %*% diag(1 /sdvals)
# #logmerged
# #logmerged<-as.data.frame(logmerged)
# #colnames(logmerged)<-names(Merged[,4:50])
# 
# 
# logmerged<-as.data.frame(logmerged)
# logmerged$Site<-Merged$Site
# logmerged$Station<-Merged$Station
# logmerged$Month<-Merged$Month
# head(logmerged)
# 
# 
# #--------------correlation check to decide which variables to retain--------------
# cor_data_phyto = cor((as.matrix(logmerged%>%select(CHLOROPHYCEAE:`PRASINOPHYCEAE (MICROMONADOPHYCEAE)`))))
# cor_data_phyto# display correlation coefficient matrix
# corrplot::corrplot(cor_data_phyto, method="number",is.corr = TRUE) #visualize correlation matrix
# corrplot::corrplot(cor_data_phyto*cor_data_phyto,method='number') # display R2 value between variables
# p_values<-rcorr (as.matrix(logmerged%>%select(CHLOROPHYCEAE:`PRASINOPHYCEAE (MICROMONADOPHYCEAE)`)))
# print(p_values)
# 
# # we can do easily in excel by using correlation in data analysis
# 
# 
# 
# FA_dataframe<-(logmerged%>%select("EPA":`??6 PUFA`)) #if want selected month: (logmerged%>%filter(!Month=="May")%>%select("EPA":"PUFA"))
# PB<-dataframe<-(logmerged%>%select("CHLOROPHYCEAE":"DOC (g.m-3)"))
# 
# # run global model with all variables
# tb_rda.all <- rda (FA_dataframe~.,PB,scale = TRUE) #scale true since we have different measurement for DOC and phyto..
# # check their significance; three methods exists; i) by global model
# anova (tb_rda.all)
# #ii) by term
# anova (tb_rda.all,by="term")
# #iii) by axis
# anova (tb_rda.all,by="axis")
# # check summary
# summary(tb_rda.all)
# 
# #mulcolinearity of the variables can be checked by using vif factors. if sqrt(vit)>2 high multicolinearity. We can remove those variables and re-run the analysis again. 
# sqrt(vif.cca(tb_rda.all))
# 
# #explained varinace can be searched by adjusted R2; since we have more variables it should be adjusted. 
# adjR2.tbrda <- RsquareAdj (tb_rda.all)$adj.r.squared
# adjR2.tbrda
# 
# summary<-summary(tb_rda.all)
# 
# sitescores<-as.data.frame(tb_rda.all$CCA$u) # This is location data; u is site data
# speciescores<-as.data.frame(tb_rda.all$CCA$v) # This is species data; v is species data; https://rdrr.io/rforge/vegan/man/cca.object.html
# envscores<-tb_rda.all$CCA$biplot              # This is environment scores
# 
# #attaching location and month data for both location and species scores
# sitescores$site<-logmerged$Site #if want selected site: (logmerged%>%filter(!Month=="May"))$Site
# sitescores$Month<-logmerged$Month # if want selected month: (logmerged%>%filter(!Month=="May"))$Month
# 
# 
# # plotting the ggplot
# ggplot(sitescores,aes(RDA1,y=RDA2)) + 
#   geom_point(aes(shape=site,
#                  fill=Month),
#              size=5,
#   ) +
#   scale_fill_brewer(palette="Paired")+
#   scale_shape_manual(values = c(21,22,25,24))+ # use fill able shapes 21-25
#   guides(shape = guide_legend(title = "Site"),
#          fill= guide_legend(title = "Month",override.aes=list(shape=21)))+
#   #adding phytoplankton data as vectors
#   geom_text(data = as.data.frame(envscores), 
#             aes(x = RDA1, y = RDA2, 
#                 label = rownames(envscores)), 
#             col = 'red') +
#   geom_segment(data = as.data.frame(envscores), 
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
#                alpha = 0.75, color = 'blue')+
#   #geom_point(data = speciescores, 
#   #aes(x = RDA1, y = RDA2), 
#   #col = 'red')+
#   geom_segment(data = as.data.frame(speciescores), 
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow=arrow(length=unit(0.2,"cm"),type = "closed"),
#                alpha = 0.75, color = 'red',size=0.1)+
#   ggrepel::geom_text_repel(data = as.data.frame(speciescores), 
#                            aes(x = RDA1, y = RDA2, 
#                                label = rownames(speciescores)), 
#                            col = 'black') +
#   labs(shape="Month", colour="Site")+
#   labs(x=paste("RDA 1 (", format(100 *summary$cont[[1]][2,1], digits=4), "%)", sep=""),
#        y=paste("RDA 2 (", format(100 *summary$cont[[1]][2,2], digits=4), "%)", sep=""))+
#   theme_bw() 
# 
# 
# #----------if the global model is significant we could run the forward selection model----------------------
# 
# #Since our case global model was significant; we could run the forward selection model
# # To run stepwise regression-the global model should be significant. If not we couldn't run it. We have to stick to the global model with all variables. 
# #full model
# full_model_RDA <- rda (FA_dataframe ~ 1, data = PB,scale = TRUE) # No constraints; only intercept; this is equal to pca
# #full_model_RDA
# Model_all_variable<- rda (FA_dataframe ~ ., data = PB,scale = TRUE) # with all variables/ constraints
# sel.osR2 <- ordiR2step (full_model_RDA, scope = formula (Model_all_variable), R2scope = adjR2.tbrda, direction = 'forward', permutations = 999)
# anova(sel.osR2 ) # To check whether model is significant; out from the result of forward selection; go ahead with those variables and run RDA again. 
# odel_all_variable12<- rda (FA_dataframe ~ CHLOROPHYCEAE + `NOSTOCOPHYCEAE (CYANOPHYCEAE)`
#                            + `DIATOMOPHYCEAE (BACILLARIOPHYCEAE)` + CRYPTOPHYCEAE, data = PB, scale = TRUE)
# odel_all_variable12
# RsquareAdj (odel_all_variable12)$adj.r.squared
# sqrt(vif.cca(odel_all_variable12))
# anova (odel_all_variable12)
# 
# 
# #run summary of the model
# summary_forward<-summary(sel.osR2)
# summary_forward
# 
# #Extracting components for drawing RDA plots
# 
# sitescores_forward<-as.data.frame(sel.osR2$CCA$u) # This is location data; u is site data
# speciescores_forward<-as.data.frame(sel.osR2$CCA$v) # This is species data; v is species data; https://rdrr.io/rforge/vegan/man/cca.object.html
# envscores_forward<-sel.osR2$CCA$biplot              # This is environment scores
# 
# #attaching location and month data for both location and species scores
# sitescores_forward$site<-logmerged$Site #if want selected site: (logmerged%>%filter(!Month=="May"))$Site
# sitescores_forward$Month<-logmerged$Month # if want selected month: (logmerged%>%filter(!Month=="May"))$Month
# 
# attach(sitescores_forward)
# attach(speciescores_forward)
# 
# 
# # plotting the ggplot
# ggplot(sitescores_forward) + 
#   geom_point(aes(x = RDA1, y = RDA2, 
#                  col = factor(site),
#                  shape = factor(Month)),
#              size=5,
#              alpha=0.5) +
#   #adding phytoplankton data as vectors
#   geom_text(data = as.data.frame(envscores_forward), 
#             aes(x = RDA1, y = RDA2, 
#                 label = rownames(envscores_forward)), 
#             col = 'red') +
#   geom_segment(data = as.data.frame(envscores_forward), 
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow=arrow(length=unit(0.3,"cm"),type = "closed"),
#                alpha = 0.75, color = 'blue')+
#   geom_point(data = speciescores_forward, 
#              aes(x = RDA1, y = RDA2), 
#              col = 'red')+
#   geom_segment(data = as.data.frame(speciescores_forward), 
#                aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow=arrow(length=unit(0.2,"cm"),type = "closed"),
#                alpha = 0.75, color = 'red',size=0.1)+
#   ggrepel::geom_text_repel(data = as.data.frame(speciescores_forward), 
#                            aes(x = RDA1, y = RDA2, 
#                                label = rownames(speciescores_forward)), 
#                            col = 'black') +
#   labs(shape="Month", colour="Site")+
#   labs(x=paste("RDA 1 (", format(100 *summary_forward$cont[[1]][2,1], digits=4), "%)", sep=""),
#        y=paste("RDA 2 (", format(100 *summary_forward$cont[[1]][2,2], digits=4), "%)", sep=""))+
#   theme_bw() 






