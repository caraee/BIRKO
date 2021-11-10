library("here")
library("tidyverse")
library("extrafont")
library("scales")
library("grid")
library("gtable")
library("data.table")
library("Amelia")
library("brms")
library("future")
library("gt")

#### Setup ####
#loadfonts()
GTT<-read_delim("GTT.txt",
                delim="\t",col_names=T)%>%pivot_longer(7:12,
                                                       names_to = "Time",values_to = "BG.Start",
                                                       values_drop_na = T)
GTT$BG.End<-case_when(GTT$BG.Start==33.3~147,
                      TRUE~GTT$BG.Start)
GTT$Group<-paste(GTT$Geno,GTT$Sex,GTT$Diet,sep=".")
GTT$Group<-factor(GTT$Group)
GTT$Geno<-factor(GTT$Geno)
GTT$Sex<-factor(GTT$Sex)
GTT$Diet<-factor(GTT$Diet)
GTT$Time<-as.numeric(GTT$Time)
GTT$BG<-ifelse(GTT$BG.End==GTT$BG.Start,GTT$BG.End,NA)

GTT<-data.frame(GTT)

GenoShapes<-c(22,3,24,21)
names(GenoShapes)<-levels(GTT$Geno)

c("#F0E442","#009E73","#0072B2","#D55E00","#CC79A7","#000000","#E69F00","#56B4E9")
GenoPal<-c("#D81B60","#004D40","#1E88E5","#FFC107")
names(GenoPal)<-levels(GTT$Geno)

GenoLines<-c(2,1,3,4)
names(GenoLines)<-levels(GTT$Geno)

#### Imputations ####
pr <- matrix(c(0,14,33.3,147,0.99999), byrow=T, ncol = 5) #Amelia priors for blood glucose above LOD

save.image("./preAmelia.RData")
set.seed(12345)
GTT.imp <- amelia(GTT, m = 10, p2s=1, idvars = c("Genotype","BG.Start","BG.End","Cohort", #m=10 only 1% of data is missing
                                                 "Group"),
                  ts="Weeks",cs="AnimalID",intercs = F,
                  polytime=1, #linear effect of time
                  ords=c("Time"),
                  noms = c("Geno","Sex","Diet","Dose"),
                  priors=pr,
                  multicore = 6)

plot(GTT.imp,which.vars = 14)
GTT.imp <- transform(GTT.imp, WeeksNom = factor(Weeks))
GTT.imp <- transform(GTT.imp, TimeNom = factor(Time))

#### Bayesian modeling ####
summary(GTT.imp)
imps_GTT<-lapply(GTT.imp$imputations, function(im){
  im
}) #change to list of dataframes

imps_GTT_LFD<-lapply(imps_GTT, function(im){
  im<-filter(im,Diet=="LFD")
  im
})

imps_GTT_HFD<-lapply(imps_GTT, function(im){
  im<-filter(im,Diet=="HFD")
  im
})

save.image("./prebrms.RData")

#run on compute cluster
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)

load("./prebrms.RData")

fit_GTT_LFD <- brm_multiple(bf(BG~Group*(WeeksNom/TimeNom)+(1|AnimalID)),
                            family="gaussian",
                            iter=8000,
                            warmup=4000,
                            thin=1,
                            prior=c(set_prior("normal(0,20)", class = "b")),
                            data = imps_GTT_LFD,
                            future=T,
                            silent=0)

save(fit_GTT_LFD,file="fit_GTT_LFD_orig.RData")

fit_GTT_HFD <- brm_multiple(bf(BG~Group*(WeeksNom/TimeNom)+(TimeNom|AnimalID)),
                            family="gaussian",
                            iter=10000,
                            warmup=5000,
                            thin=1,
                            prior=c(set_prior("normal(0,50)", class = "b")),
                            data = imps_GTT_HFD,
                            future=T,
                            silent=0)

save(fit_GTT_HFD,file="fit_GTT_HFD_orig.RData")

load("./fit_GTT_HFD_orig.RData")
load("./fit_GTT_LFD_orig.RData")

pp_check(fit_GTT_HFD)
pp_check(fit_GTT_LFD)

bayestestR::check_prior(fit_GTT_HFD)
bayestestR::check_prior(fit_GTT_LFD)

summary(fit_GTT_HFD)
max(fit_GTT_HFD$rhats) #no issues

summary(fit_GTT_LFD)
max(fit_GTT_LFD$rhats) #no issues

#### Figures ####
newdata = data.frame(Group = rep(levels(GTT$Group),6*6)|>str_sort(),
                     WeeksNom = rep(c(4,9,21,25,39,54),#unique(GTT$Weeks),
                                    each=6)|>rep(16),
                     TimeNom = rep(unique(GTT$Time),16*6))
newdata$Diet<-str_split_fixed(newdata$Group,"[.]",3)[,3]
newdata$Sex<-str_split_fixed(newdata$Group,"[.]",3)[,2]
newdata$Geno<-str_split_fixed(newdata$Group,"[.]",3)[,1]
newdata$Time<-newdata$TimeNom
newdata$Weeks<-newdata$WeeksNom

newdataLFD<-filter(newdata,Diet=="LFD")|>filter(Weeks!=54)|>filter(Weeks!=25)
newdataHFD<-filter(newdata,Diet=="HFD")|>filter(Weeks!=4)

fitLFD <- fitted(fit_GTT_LFD,
                 newdata = newdataLFD, 
                 re_formula = NA) # extract the full MCMC
fitHFD <- fitted(fit_GTT_HFD,
                 newdata = newdataHFD, 
                 re_formula = NA) # extract the full MCMC

ffLFD <- fitLFD |>
  as_tibble() |>
  bind_cols(newdataLFD)

ffHFD <- fitHFD |>
  as_tibble() |>
  bind_cols(newdataHFD)

ffHFD$Group<-factor(ffHFD$Group,levels=levels(GTT$Group))
ffHFD$Sex<-factor(ffHFD$Sex,levels=levels(GTT$Sex))
ffHFD$Geno<-str_split_fixed(ffHFD$Group,"[.]",3)[,1] |> factor(levels=levels(GTT$Geno))
ffHFD$Geno<-factor(ffHFD$Geno,levels=levels(GTT$Geno))
ffHFD<-filter(ffHFD,Geno!="InsR")|>filter(!(Sex=="Female"&Diet=="HFD"&Weeks==25))|>
  filter(!(Sex=="Male"&Diet=="HFD"&Weeks==21))
ffHFD$Weeks[ffHFD$Sex=="Male"&ffHFD$Diet=="HFD"&ffHFD$Weeks==25]<-21

#not interested in "InsR" genotype
#Week 21 in male mice fed HFD was mostly above LOD and not performed in all mice, so repeated at 25 weeks
#Use as 3rd GTT
#Week 25 in female mice fed HFD was 1 g/kg (too low)
GTT<-filter(GTT,Geno!="InsR")|>filter(!(Sex=="Female"&Diet=="HFD"&Weeks==25))|>
  filter(!(Sex=="Male"&Diet=="HFD"&Weeks==21))
GTT$Weeks[GTT$Sex=="Male"&GTT$Diet=="HFD"&GTT$Weeks==25]<-21

ffLFD$Group<-factor(ffLFD$Group,levels=levels(GTT$Group))
ffLFD$Sex<-factor(ffLFD$Sex,levels=levels(GTT$Sex))
ffLFD<-filter(ffLFD,Geno!="InsR")
ffLFD$Geno<-str_split_fixed(ffLFD$Group,"[.]",3)[,1]
ffLFD$Geno<-factor(ffLFD$Geno,levels=levels(GTT$Geno))

write_csv(ffHFD,file="./GTT_HFD_fitted_orig.csv")
write_csv(ffLFD,file="./GTT_LFD_fitted_orig.csv")

labsGeno<-c("HetCre"=expression(paste(italic("Insr"^"f/wt"),italic("Ins1"^"cre/wt"),"nTnG"^"+/wt",sep="")),
            "InsRCre"=expression(paste(italic("Insr"^"f/f"),italic("Ins1"^"cre/wt"),"nTnG"^"+/wt",sep="")),
            "InsR"=expression(paste(italic("Insr"^"f/f"),italic("Ins1"^"wt/wt"),"nTnG"^"+/wt",sep="")),
            "WTCre"=expression(paste(italic("Insr"^"wt/wt"),italic("Ins1"^"cre/wt"),"nTnG"^"+/wt",sep="")))

limits<-c(-5,125)
breaks=c(0,15,30,60,90,120)
pd <- position_dodge(width=2)

dat_fem_HFD <- data.frame(
  label = c("♀"),
  Weeks   = c(9),
  Sex = c("Female"),
  x     = c(15),
  y     = c(3)
)

dat_male_HFD <- data.frame(
  label = c("♂"),
  Weeks   = c(9),
  Sex = c("Male"),
  x     = c(15),
  y     = c(3)
)

dat_weeks_HFD<-data.frame(
  label = c("9","21","39","54"),
  Weeks   = c(9,21,39,54),
  Sex = c("Male"),
  x     = c(60),
  y     = c(1)
)

pHFD<-ggplot(filter(GTT,Diet=="HFD"),
          aes(x=Time,y=BG.Start))+
  facet_grid(Sex~Weeks,scales = "free",drop=T)+
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,group=AnimalID,colour=Geno),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.9)+
  scale_colour_manual(breaks=names(GenoPal)[-2],
                      name="Geno", values=GenoPal,labels=labsGeno)+
  scale_shape_manual(breaks=names(GenoShapes)[-2],name="Geno", values=GenoShapes,labels=labsGeno)+ #,labels=labsGen
  scale_fill_manual(breaks=names(GenoPal)[-2],
                    name="Geno", values=GenoPal,labels=labsGeno)+
  scale_linetype_manual(breaks=names(GenoLines)[-2],
                        name="Geno", values=GenoLines,labels=labsGeno)+
  geom_text(data    = dat_fem_HFD,
            mapping = aes(x = x, y = y, label = label),
            colour="#D81B60",
            size=10)+
  geom_text(data    = dat_male_HFD,
            mapping = aes(x = x, y = y, label = label),
            colour="#1014D4", 
            size=10)+
  geom_text(data    = dat_weeks_HFD,
            mapping = aes(x = x, y = y, label = label),
            colour="#666666", 
            size=4)+
  ggtitle("HFD")+
  theme(plot.title = element_text(family = "Arial",color="#666666", size=12,vjust = -5,hjust=0.01))+
  geom_line(aes(colour=Geno,group=AnimalID,linetype=Geno),position=pd,alpha=0.9)+
  geom_hline(data=data.frame(yint = 33.3),aes(yintercept=yint),colour="#666666",linetype=3)+
  labs(x="",y="Blood Glucose (mM)")+
  coord_cartesian(ylim=c(0, 35))+
  scale_x_continuous(breaks = breaks)+
  theme(axis.title = element_text(family = "Arial", color="black", size=10),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=10))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=10))+
  theme(legend.text = element_text(family = "Arial",color="black",size=10), 
        legend.title = element_blank())+
  theme(strip.text = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position ="bottom")+
  geom_smooth(data = ffHFD,
              aes(y = Estimate, ymin = Q2.5, 
                  ymax = Q97.5,
                  fill = Geno,
                  colour=Geno,group=Geno),
              stat = "identity",
              alpha = 0.4, size = 1/2) +
  geom_point(data=ffHFD,aes(y=Estimate,shape=Geno,
                            colour=Geno,group=Geno,fill=Geno),position=pd,size=1.75,
             alpha=1)
pHFD

ggsave("GTT_BG_HFD_Bayes_orig.png",
       path="./Figures",width = 8, height = 5, units = "in")
ggsave("GTT_BG_HFD_Bayes_orig.svg",
       path="./Figures",width = 8, height = 5, units = "in")

dat_fem_LFD <- data.frame(
  label = c("♀"),
  Weeks   = c(4),
  Sex = c("Female"),
  x     = c(15),
  y     = c(3)
)

dat_male_LFD <- data.frame(
  label = c("♂"),
  Weeks   = c(4),
  Sex = c("Male"),
  x     = c(15),
  y     = c(3)
)

dat_weeks_LFD<-data.frame(
  label = c("4","9","21","39"),
  Weeks   = c(4,9,21,39),
  Sex = c("Male"),
  x     = c(60),
  y     = c(1)
)

pLFD<-ggplot(GTT[GTT$Diet=="LFD",],
          aes(x=Time,y=BG.Start))+
  facet_grid(Sex~Weeks,scales = "free")+
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,
                     group=AnimalID,colour=Geno),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.9)+
  scale_colour_manual(breaks=names(GenoPal)[-2],
    name="Geno", values=GenoPal,labels=labsGeno)+
  scale_shape_manual(breaks=names(GenoShapes)[-2],name="Geno", values=GenoShapes,labels=labsGeno)+ #,labels=labsGen
  scale_fill_manual(breaks=names(GenoPal)[-2],
                    name="Geno", values=GenoPal,labels=labsGeno)+
  scale_linetype_manual(breaks=names(GenoLines)[-2],
                        name="Geno", values=GenoLines,labels=labsGeno)+
  geom_text(data    = dat_fem_LFD,
            mapping = aes(x = x, y = y, label = label),
            colour="#D81B60",
            size=10)+
  geom_text(data    = dat_male_LFD,
            mapping = aes(x = x, y = y, label = label),
            colour="#1014D4", 
            size=10)+
  geom_text(data    = dat_weeks_LFD,
            mapping = aes(x = x, y = y, label = label),
            colour="#666666", 
            size=4)+
  ggtitle("LFD")+
  theme(plot.title = element_text(family = "Arial",color="#666666", size=12,vjust = -5,hjust=0.01))+
  geom_line(aes(colour=Geno,group=AnimalID,linetype=Geno),position=pd,alpha=0.9)+
  geom_hline(data=data.frame(yint = 33.3),aes(yintercept=yint),colour="#666666",linetype=3)+
  labs(x="",y="Blood Glucose (mM)")+
  coord_cartesian(ylim=c(0, 35))+
  scale_x_continuous(breaks = breaks)+
  theme(axis.title = element_text(family = "Arial", color="black", size=10),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=10))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=10))+
  theme(legend.text = element_text(family = "Arial",color="black",size=10), 
        legend.title = element_blank())+
  theme(strip.text = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position ="bottom")+
  geom_smooth(data = ffLFD,
              aes(y = Estimate, ymin = Q2.5, 
                  ymax = Q97.5,
                  colour=Geno,
                  fill = Geno,
                  group=Geno),
              stat = "identity",
              alpha = 0.4, size = 1/2) +
  geom_point(data=ffLFD,aes(y=Estimate,shape=Geno,
                            colour=Geno,fill=Geno,
                            group=Geno),position=pd,size=1.75,
             alpha=1)
pLFD

ggsave("GTT_BG_LFD_Bayes_orig.png", path="./Figures",width = 8, height = 5, units = "in")
ggsave("GTT_BG_LFD_Bayes_orig.svg", path="./Figures",width = 8, height = 5, units = "in")

#### Comparisons ####
######LFD#####
#Female
#Reference is GroupHetCre.Female.LFD, WeeksNom4, TimeNom0
#no between sex comparisons
LFDtab<-rbind(
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD)-(Intercept+GroupWTCre.Female.LFD)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD)-(Intercept)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept)-(Intercept+GroupWTCre.Female.LFD)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom15+GroupWTCre.Female.LFD:WeeksNom4:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom15)-
           (Intercept+WeeksNom4:TimeNom15)>0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom4:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom15+GroupWTCre.Female.LFD:WeeksNom4:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom30+GroupWTCre.Female.LFD:WeeksNom4:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom30)-
           (Intercept+WeeksNom4:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom4:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom30+GroupWTCre.Female.LFD:WeeksNom4:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom60+GroupWTCre.Female.LFD:WeeksNom4:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom60)-
           (Intercept+WeeksNom4:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom4:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom60+GroupWTCre.Female.LFD:WeeksNom4:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom90+GroupWTCre.Female.LFD:WeeksNom4:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom90)-
           (Intercept+WeeksNom4:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom4:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom90+GroupWTCre.Female.LFD:WeeksNom4:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom120+GroupWTCre.Female.LFD:WeeksNom4:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom4:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom4:TimeNom120)-
           (Intercept+WeeksNom4:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom4:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom4:TimeNom120+GroupWTCre.Female.LFD:WeeksNom4:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9)-(Intercept+WeeksNom9)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom9)-(Intercept+GroupWTCre.Female.LFD+WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom15)-
           (Intercept+WeeksNom9+WeeksNom9:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom9+WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom30)-
           (Intercept+WeeksNom9+WeeksNom9:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom9+WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom60)-
           (Intercept+WeeksNom9+WeeksNom9:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom9+WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom90)-
           (Intercept+WeeksNom9+WeeksNom9:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom9+WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom9+GroupInsRCre.Female.LFD:WeeksNom9:TimeNom120)-
           (Intercept+WeeksNom9+WeeksNom9:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom9+WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupWTCre.Female.LFD:WeeksNom9+GroupWTCre.Female.LFD:WeeksNom9:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21)-(Intercept+WeeksNom21)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom21)-(Intercept+GroupWTCre.Female.LFD+WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom15)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom30)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom60)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom90)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom21+GroupInsRCre.Female.LFD:WeeksNom21:TimeNom120)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Female.LFD:WeeksNom21+GroupWTCre.Female.LFD:WeeksNom21:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39)-(Intercept+WeeksNom39)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom39)-(Intercept+GroupWTCre.Female.LFD+WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom15)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom30)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom60)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom90)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Female.LFD:WeeksNom39+GroupInsRCre.Female.LFD:WeeksNom39:TimeNom120)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Female.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Female.LFD:WeeksNom39+GroupWTCre.Female.LFD:WeeksNom39:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD)-(Intercept+GroupWTCre.Male.LFD)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom15+GroupWTCre.Male.LFD:WeeksNom4:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom30+GroupWTCre.Male.LFD:WeeksNom4:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom60+GroupWTCre.Male.LFD:WeeksNom4:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom90+GroupWTCre.Male.LFD:WeeksNom4:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom120+GroupWTCre.Male.LFD:WeeksNom4:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD)-(Intercept+GroupHetCre.Male.LFD)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom15)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom15+GroupHetCre.Male.LFD:WeeksNom4:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom30)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom30+GroupHetCre.Male.LFD:WeeksNom4:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom60)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom60+GroupHetCre.Male.LFD:WeeksNom4:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom90)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom90+GroupHetCre.Male.LFD:WeeksNom4:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom4:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom4:TimeNom120)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom120+GroupHetCre.Male.LFD:WeeksNom4:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD)-(Intercept+GroupWTCre.Male.LFD)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom15+GroupHetCre.Male.LFD:WeeksNom4:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom15+GroupWTCre.Male.LFD:WeeksNom4:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom30+GroupHetCre.Male.LFD:WeeksNom4:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom30+GroupWTCre.Male.LFD:WeeksNom4:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom60+GroupHetCre.Male.LFD:WeeksNom4:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom60+GroupWTCre.Male.LFD:WeeksNom4:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom90+GroupHetCre.Male.LFD:WeeksNom4:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom90+GroupWTCre.Male.LFD:WeeksNom4:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom4:TimeNom120+GroupHetCre.Male.LFD:WeeksNom4:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom4:TimeNom120+GroupWTCre.Male.LFD:WeeksNom4:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom9+GroupInsRCre.Male.LFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom15+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom30+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom60+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom90+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupHetCre.Male.LFD:WeeksNom9+GroupHetCre.Male.LFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom9+WeeksNom9:TimeNom120+GroupWTCre.Male.LFD:WeeksNom9+GroupWTCre.Male.LFD:WeeksNom9:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom21+GroupInsRCre.Male.LFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupHetCre.Male.LFD:WeeksNom21+GroupHetCre.Male.LFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Male.LFD:WeeksNom21+GroupWTCre.Male.LFD:WeeksNom21:TimeNom120)<0")$hypothesis,

hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupInsRCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Male.LFD:WeeksNom39+GroupInsRCre.Male.LFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom120)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom15)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom30)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom60)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom90)<0")$hypothesis,
hypothesis(fit_GTT_LFD,"(Intercept+GroupHetCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupHetCre.Male.LFD:WeeksNom39+GroupHetCre.Male.LFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Male.LFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Male.LFD:WeeksNom39+GroupWTCre.Male.LFD:WeeksNom39:TimeNom120)<0")$hypothesis)

######HFD#####
HFDtab<-rbind(
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD)-(Intercept+GroupWTCre.Female.HFD)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD)-(Intercept)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept)-(Intercept+GroupWTCre.Female.HFD)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom15+GroupWTCre.Female.HFD:WeeksNom9:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom15)-
           (Intercept+WeeksNom9:TimeNom15)>0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom15+GroupWTCre.Female.HFD:WeeksNom9:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom30+GroupWTCre.Female.HFD:WeeksNom9:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom30)-
           (Intercept+WeeksNom9:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom30+GroupWTCre.Female.HFD:WeeksNom9:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom60+GroupWTCre.Female.HFD:WeeksNom9:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom60)-
           (Intercept+WeeksNom9:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom60+GroupWTCre.Female.HFD:WeeksNom9:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom90+GroupWTCre.Female.HFD:WeeksNom9:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom90)-
           (Intercept+WeeksNom9:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom90+GroupWTCre.Female.HFD:WeeksNom9:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom120+GroupWTCre.Female.HFD:WeeksNom9:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom9:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom9:TimeNom120)-
           (Intercept+WeeksNom9:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom9:TimeNom120+GroupWTCre.Female.HFD:WeeksNom9:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21)-(Intercept+WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom21)-(Intercept+GroupWTCre.Female.HFD+WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom15)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom30)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom60)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom90)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom21+GroupInsRCre.Female.HFD:WeeksNom21:TimeNom120)-
           (Intercept+WeeksNom21+WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom21+WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Female.HFD:WeeksNom21+GroupWTCre.Female.HFD:WeeksNom21:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39)-(Intercept+WeeksNom39)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom39)-(Intercept+GroupWTCre.Female.HFD+WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom15)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom30)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom60)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom90)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom39+GroupInsRCre.Female.HFD:WeeksNom39:TimeNom120)-
           (Intercept+WeeksNom39+WeeksNom39:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom39+WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Female.HFD:WeeksNom39+GroupWTCre.Female.HFD:WeeksNom39:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54)-(Intercept+WeeksNom54)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom54)-(Intercept+GroupWTCre.Female.HFD+WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom15)-
           (Intercept+WeeksNom54+WeeksNom54:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom54+WeeksNom54:TimeNom15)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom30)-
           (Intercept+WeeksNom54+WeeksNom54:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom54+WeeksNom54:TimeNom30)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom60)-
           (Intercept+WeeksNom54+WeeksNom54:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom54+WeeksNom54:TimeNom60)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom90)-
           (Intercept+WeeksNom54+WeeksNom54:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom54+WeeksNom54:TimeNom90)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupInsRCre.Female.HFD:WeeksNom54+GroupInsRCre.Female.HFD:WeeksNom54:TimeNom120)-
           (Intercept+WeeksNom54+WeeksNom54:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+WeeksNom54+WeeksNom54:TimeNom120)-
           (Intercept+GroupWTCre.Female.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupWTCre.Female.HFD:WeeksNom54+GroupWTCre.Female.HFD:WeeksNom54:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD)-(Intercept+GroupWTCre.Male.HFD)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom15+GroupWTCre.Male.HFD:WeeksNom9:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom30+GroupWTCre.Male.HFD:WeeksNom9:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom60+GroupWTCre.Male.HFD:WeeksNom9:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom90+GroupWTCre.Male.HFD:WeeksNom9:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom120+GroupWTCre.Male.HFD:WeeksNom9:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD)-(Intercept+GroupHetCre.Male.HFD)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom15+GroupHetCre.Male.HFD:WeeksNom9:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom30+GroupHetCre.Male.HFD:WeeksNom9:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom60+GroupHetCre.Male.HFD:WeeksNom9:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom90+GroupHetCre.Male.HFD:WeeksNom9:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom9:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom120+GroupHetCre.Male.HFD:WeeksNom9:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD)-(Intercept+GroupWTCre.Male.HFD)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom15+GroupHetCre.Male.HFD:WeeksNom9:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom15+GroupWTCre.Male.HFD:WeeksNom9:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom30+GroupHetCre.Male.HFD:WeeksNom9:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom30+GroupWTCre.Male.HFD:WeeksNom9:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom60+GroupHetCre.Male.HFD:WeeksNom9:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom60+GroupWTCre.Male.HFD:WeeksNom9:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom90+GroupHetCre.Male.HFD:WeeksNom9:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom90+GroupWTCre.Male.HFD:WeeksNom9:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom9:TimeNom120+GroupHetCre.Male.HFD:WeeksNom9:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom9:TimeNom120+GroupWTCre.Male.HFD:WeeksNom9:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom21+GroupInsRCre.Male.HFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom15+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom30+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom60+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom90+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupHetCre.Male.HFD:WeeksNom21+GroupHetCre.Male.HFD:WeeksNom21:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom21+WeeksNom21:TimeNom120+GroupWTCre.Male.HFD:WeeksNom21+GroupWTCre.Male.HFD:WeeksNom21:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom39+GroupInsRCre.Male.HFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom15+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom30+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom60+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom90+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupHetCre.Male.HFD:WeeksNom39+GroupHetCre.Male.HFD:WeeksNom39:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom39+WeeksNom39:TimeNom120+GroupWTCre.Male.HFD:WeeksNom39+GroupWTCre.Male.HFD:WeeksNom39:TimeNom120)<0")$hypothesis,
  
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom15)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom30)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom60)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom90)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupInsRCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupInsRCre.Male.HFD:WeeksNom54+GroupInsRCre.Male.HFD:WeeksNom54:TimeNom120)-
           (Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom15)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom15+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom15)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom30)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom30+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom60)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom60+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom90)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom90+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GTT_HFD,"(Intercept+GroupHetCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupHetCre.Male.HFD:WeeksNom54+GroupHetCre.Male.HFD:WeeksNom54:TimeNom120)-
           (Intercept+GroupWTCre.Male.HFD+WeeksNom54+WeeksNom54:TimeNom120+GroupWTCre.Male.HFD:WeeksNom54+GroupWTCre.Male.HFD:WeeksNom54:TimeNom120)<0")$hypothesis)

#### Results Table ####
LFDtab_sub<-filter(LFDtab,Star=="*")
HFDtab_sub<-filter(HFDtab,Star=="*")

compsTab<-rbind(LFDtab_sub,HFDtab_sub)

h<-str_replace_all(compsTab$Hypothesis,"[(]exp[(]","")|>str_replace("[)]{2}","")|>str_replace("Intercept","")
  
compsTab$Sex<-ifelse(grepl("Female",h),"Female","Male") 
compsTab$Diet<-ifelse(grepl("LFD",h),"LFD","HFD") 
compsTab$Diet<-factor(compsTab$Diet,levels=c("LFD","HFD"))
compsTab$ContrGeno<-str_split_fixed(h, boundary("word"),18)[,1]
compsTab$ContrGeno<-str_replace(compsTab$ContrGeno,".Female.LFD","")|>
  str_replace(".Male.LFD","")|>str_replace(".Female.HFD","")|>
  str_replace(".Male.HFD","")|>str_replace("Group","")
compsTab$ContrGeno<-str_replace(compsTab$ContrGeno,"WeeksNom[0-9]{1,2}","HetCre")
compsTab$RefGeno<-str_split_fixed(h, "[-(]Intercept",2)[,2]
compsTab$RefGeno<-str_split_fixed(compsTab$RefGeno,boundary("word"),3)[,1]
compsTab$RefGeno<-str_replace(compsTab$RefGeno,".Female.LFD","")|>
  str_replace(".Male.LFD","")|>str_replace(".Female.HFD","")|>
  str_replace(".Male.HFD","")|>str_replace("Group","")
compsTab$RefGeno<-case_when(grepl("Weeks",compsTab$RefGeno)~"HetCre",
                              TRUE~compsTab$RefGeno)

GenoPal2<-GenoPal[c(3:4,1)]

compsTab$Age<-str_split_fixed(h, boundary("word"),18)[,2]|>str_replace("exp","WeeksNom4")|>
  str_replace("TimeNom15","WeeksNom4")
compsTab$Age<-str_replace_all(compsTab$Age,"WeeksNom","")
compsTab$Age<-factor(compsTab$Age,levels=unique(compsTab$Age))

compsTab$Time<-str_extract(compsTab$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
compsTab$Time[is.na(compsTab$Time)]<-0

compsTab<-filter(compsTab,RefGeno!="HetCre")
compsTab$RefGeno<-factor(compsTab$RefGeno,levels=c("InsRCre","WTCre","HetCre"))
compsTab$ContrGeno<-factor(compsTab$ContrGeno,levels=c("InsRCre","WTCre","HetCre"))


# htmlLabs<-c(html("<em>Insr</em><sup>f/f</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup> vs.<br><em>Insr</em><sup>wt/wt</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup>"),
# html("<em>Insr</em><sup>f/wt</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup> vs.<br><em>Insr</em><sup>wt/wt</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup>"),
# html("<em>Insr</em><sup>f/f</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup> vs.<br><em>Insr</em><sup>f/wt</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup>"))
# 
# # # 
# htmlLabs<-c("Insrf/fIns1cre/wtnTnG+/wt vs. Insrwt/wtIns1cre/wtnTnG+/wt",
#             "Insrf/wtIns1cre/wtnTnG+/wt vs.Insrwt/wtIns1cre/wtnTnG+/wt",
#             "Insrf/fIns1cre/wtnTnG+/wt vs.Insrf/wt~Ins1cre/wtnTnG+/wt")
# 
# labsGeno2<-c("HetCre"="<em>Insr</em><sup>f/wt</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup>",
#             "InsRCre"="<em>Insr</em><sup>f/f</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup>",
#             "WTCre"="<em>Insr</em><sup>wt/wt</sup><em>Ins1</em><sup>cre/wt</sup>nTnG<sup>+/wt</sup>")

tab1<-compsTab |> 
  mutate(Evidence=case_when(compsTab$Evid.Ratio<=10^(3/2)~"Strong",
                            compsTab$Evid.Ratio<=10^(2)~"Very strong",
                            compsTab$Evid.Ratio>10^(2)~"Decisive")) |>
  arrange(Diet,Sex,RefGeno,ContrGeno,Age,Time)|>
  group_by(Diet,Sex) |> 
  gt() |>
  cols_hide(columns=c(Hypothesis,Est.Error,Post.Prob,Star))|>
  tab_header(title = md("Table. Comparisons Between Genotype, Within Sex and Diet")) |>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  tab_spanner(
    label = "Genotype Comparisons",
    columns = c(RefGeno,ContrGeno)
  )|>
  cols_move_to_start(
    columns = c(RefGeno,ContrGeno,Age,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Bayes Factor",
             RefGeno = "Reference",
             ContrGeno = "Contrast",
             Estimate = "Difference in\nBlood Glucose (mM)",
             Age = "Age (weeks)",
             Time="Time (minutes)")|>
fmt_number(
  columns = c(Estimate,CI.Lower,CI.Upper),
  decimals = 1,
  use_seps = TRUE
) |>
  fmt_number(
    columns = c(Evid.Ratio),
    decimals = 0,
    use_seps = TRUE
  ) |>
  # data_color(
  #   columns = c(RefGeno),
  #   colors = scales::col_factor(
  #     palette = GenoPal[4],
  #     domain = levels(compsTab$RefGeno)
  #   )) |>
  # data_color(
  #   columns = c(ContrGeno),
  #   colors = scales::col_factor(
  #     palette = GenoPal[c(1,3,4)],
  #     domain = levels(compsTab$ContrGeno)
  #   )) |>
  tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
tab1

gtsave(tab1,"resultsTab_orig.pdf",path="./")
gtsave(tab1,"resultsTab_orig.png",path="./")

compsTab2<-compsTab |> 
  mutate(Evidence=case_when(compsTab$Evid.Ratio<=10^(3/2)~"Strong",
                            compsTab$Evid.Ratio<=10^(2)~"Very strong",
                            compsTab$Evid.Ratio>10^(2)~"Decisive")) |>
  arrange(Diet,Sex,RefGeno,ContrGeno,Age,Time)|>
  mutate(across(c(Estimate,CI.Lower,CI.Upper,Evid.Ratio,Post.Prob), round,1))|>
  select(-Hypothesis)|>
  select(-Est.Error)|>
  group_by(Diet,Sex)

write_csv(compsTab2,file="comparisonsTab.csv")
#