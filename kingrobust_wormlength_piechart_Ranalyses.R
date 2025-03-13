#R script
#eddy dowle March 2025
#this includes the R code for all the plots in the paper (relatedness, worm size and pie charts on genome repeats) 
#raw input files are included in this git repository

###############################
#hairworm relatedness analyses#
###############################


#want to plot 
#R0 vs R1 and kingrobust vs R1
#following the ryan waples paper https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14954
#and Helen's paper has a really nice picture in the supplementary breaking down the different relationship guides

#table from angsd + ngsRelate
setwd("C:/Users/hrlexd/Dropbox/Otago_2023 (1)/Hairworm/GBS_2023/Hairworm_GBS_2024")

library(tidyverse)
library(ggrepel)
ngsrelate_out<-read.table('newres.cleanKING.txt', header =T, row.names=NULL,sep='\t')
samplelist<-read.table('Samplelist.txt',header=T, row.names=NULL,sep='\t')

#filtering comparisons with less than 400 SNPS
ngsrelate_out_filter<-ngsrelate_out %>% filter( nSites>= 400)

#general plot KING vs R1
ggplot(ngsrelate_out_filter, aes(x=R1, y=KING)) +
  geom_point()

#filtering out conspecific comparisons
ngsrelate_out_filter<-ngsrelate_out_filter %>% mutate(name=paste0(a,b))

ngsrelate_out_conspecific<-ngsrelate_out_filter %>% filter(a %in% 0:5 & b %in% 1:6|a==7 & b==8|a==11 &b==12|a==13&b==14|a %in% 18:19 & b %in% 19:20|a==35&b==36|a==37&b==38|a==41&b==42|a %in% 43:50 & b %in% 44:51)

#comparisons of interest in our dataset
#HW14 0:6
#HW17 7:8
#HW24 11:12
#HW29 13:14
#HW_ch11 18:20
#HW_w133 35:36
#HW_w134 37:38
#HW_w41:42
#HW_w58 43:51

#KING vs R1 conspecific worms
ggplot(ngsrelate_out_conspecific, aes(x=R1, y=KING)) +
  geom_point()

ngsrelate_out_filter<-ngsrelate_out_filter %>% mutate(colouring=case_when(name %in% ngsrelate_out_conspecific$name ~ 'firebrick2', !name %in% ngsrelate_out_conspecific$name ~'black' ))

ngsrelate_out_filter<-ngsrelate_out_filter %>% mutate(WormType=case_when(name %in% ngsrelate_out_conspecific$name ~ 'Co-occurring', !name %in% ngsrelate_out_conspecific$name ~'General population' ))

#KING vs R1 coloured by worm type
ngsrelate_out_filter<-ngsrelate_out_filter %>% arrange(colouring)
ggplot(ngsrelate_out_filter, aes(x=R1, y=KING)) +
  geom_point(size=4,colour=ngsrelate_out_filter$colouring)+
  theme_bw()+
  geom_hline(yintercept=0.354, linetype="dashed", 
             color = "magenta1", linewidth=1) +
  geom_hline(yintercept=0.177, linetype="dashed", 
             color = "blueviolet", linewidth=1)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "dodgerblue", linewidth=1)
#from helen taylor#
#The KING manual https://www.kingrelatedness.com/manual.shtml provides relationship ranges: >0.354, [0.177, 0.354], [0#.0884, 0.177] and [0.0442, 0.0884] corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree #relationships respectively. The lines on B) are these ranges which I understand to be derived from the standard #deviations of the estimates in the original paper. They are very much a "guide" which assists with visualising the #data and are certainly not fixed values. The delineation between categories gets even more vague at lower #relationships, e.g. I would advocate against interpreting 2nd vs 3rd degree relationships. 

#Just plotting KING values to keep things easy to understand for paper
ggplot(ngsrelate_out_filter,aes(y=KING,x=WormType))+
  geom_boxplot(aes(fill=WormType))+
  theme_bw()+
  geom_hline(yintercept=0.354, linetype="dashed", 
             color = "magenta1", size=1) +
  geom_hline(yintercept=0.177, linetype="dashed", 
             color = "blueviolet", size=1)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "dodgerblue", size=1)+
  scale_fill_manual(values=c("firebrick1", "cyan3"))+
  labs(x="Worm type",y='KING-robust',fill='Worm type')



####################
#plots on worm size#
####################

size_info_worms<-read.csv('comparisonsofinterest.csv', header =T, row.names=NULL)

#size range of mature worms
size_info_worms_mature<-size_info_worms %>% filter(WormMaturity=='Mature')

#size range of mature worms
size_info_worms_mature<-size_info_worms %>% filter(WormMaturity=='Mature')
size_info_worms_juvenille<-size_info_worms %>% filter(WormMaturity=='Juv')
size_info_worms<-size_info_worms %>% mutate(conspecific_name=case_when(dual=='N' ~ 'Single infection',dual=='Y' ~'Multiple infection',dual=='' ~ 'Free-living',is.na(dual) ~ 'Free-living'))

size_info_worms_mature_cass<-size_info_worms_mature %>% filter(!is.na(Size)) %>% filter(Worms_collected_for!='Jeff_freelivingcass') %>% filter(Worms_collected_for!='Jeff_eggcollection')

size_info_worms<-size_info_worms %>% mutate(mature_worm=case_when(Code %in% size_info_worms_mature$Code ~ 'Mature' ))
size_info_worms<-size_info_worms %>% mutate(mature_worm_cave=case_when(Code %in% size_info_worms_mature_cass$Code ~ 'Mature_cass' ))
size_info_worms<-size_info_worms %>% mutate(juv_worm_cave=case_when(Code %in% size_info_worms_juvenille$Code ~ 'Juv_cass' ))
size_info_worms$GeneralPopulation<-'Gordius paranensis'

#go long I guess
long <- size_info_worms %>% 
  pivot_longer(
    cols = c(mature_worm,mature_worm_cave,GeneralPopulation,juv_worm_cave), 
    names_to = "type",
    values_to = "value"
  )
unique(long$type)
long$type <- factor(long$type, levels = c("GeneralPopulation", "mature_worm", "mature_worm_cave", "juv_worm_cave"))
test<-long %>% filter(!is.na(value)) %>% filter(!is.na(Size))

#dot plot of various groups
ggplot(long %>% filter(!is.na(value)) %>% filter(!is.na(Size)) ,aes(y=Size,x=type))+
  # geom_boxplot(aes(fill=Sex))+
  geom_point(position=position_jitterdodge(jitter.width=0.5,dodge.width = 0.1),aes(colour=Sex,shape=conspecific_name),size=3)+
  theme_bw()+
  scale_colour_manual(values=c("firebrick1", "cyan3"))+
  labs(x="Type",y='Size (mm)')
# scale_fill_manual(values=c('black','purple','green'))
size_info_worms %>% filter(WormMaturity=='Juv') %>% nrow()

################################
#plots on genome repeat content#
################################
repeats<-read.table('TE_classification proportion.txt', header=T, row.names=NULL,sep='\t')

repeats_2 <- repeats %>% 
  arrange((proportion)) %>%
  mutate(prop =proportion *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop, 
  SubSegment = factor(TE_classification, levels=TE_classification[order(-(prop))], ordered=TRUE))

repeats_2[,'prop']=round(repeats_2[,'prop'],1)


gb2020label <- 
  repeats %>% 
  mutate(prop =proportion *100) %>%
  arrange((proportion)) %>% ## arrange in the order of the legend
  mutate(text_y = cumsum(prop) - prop/2,
         SubSegment = factor(TE_classification, levels=TE_classification[order(-(prop))], ordered=TRUE)) ### calculate 

#where to place the text labels
gb2020label[,'prop']=round(gb2020label[,'prop'],1)

#pie charts
ggplot(repeats_2, aes(x="", y=prop, fill=SubSegment)) +
  geom_bar(stat="identity", width=1,color='white',alpha=0.8) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_brewer(palette="Set1")# +

ggplot(gb2020label, aes(x="", y=prop, fill=SubSegment)) +
  geom_bar(stat="identity", width=1,color='white',alpha=0.8) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_brewer(palette="Set1") +
  #theme(legend.position="none") +
  geom_label_repel(aes(label = prop, y = text_y), 
                   force=2,nudge_x = 0.8, nudge_y = 0.8,
                   size = 7, show.legend = T) +
  guides(fill=guide_legend(title="TE"))
