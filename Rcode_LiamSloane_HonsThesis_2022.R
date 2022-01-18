#############################################################
##R code calculation for Liam Sloane Hons Thesis


#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

###########################################################
##packages
#########################################################

library ("ggplot2")
library("PMCMRplus")




##working directory

setwd ("  ")
getwd()

######################################################
###1. BSA_quantification of bacterial and fungal DNA
#######################################################

#Import data file all lines
q_Fungal <-(read.delim("ITS_copies.txt", sep = "\t"))


############################################################
###############Plotting  Log Fungal ITS copies per genotype
############################################################

q_Fungal$Genotype <-  ordered(q_Fungal$Genotype, levels=c( "Bulk","HID144","Barke", "Int-17","Int-52", "Int-19", "Int-56", "Morex", "GP", "RGT Planet", "Hockett"))


q <-ggplot(q_Fungal, aes(x=Genotype, y=Log_ITS_copies, fill=Genotype)) + geom_boxplot() + ggtitle(" Log Fungal ITS copies")+ylim(5,15)
q + geom_jitter( size=5,shape=21, position=position_jitter(0))+ scale_fill_viridis_d()


##################################################################
#Stats Fungal genotypes
shapiro.test(q_Fungal$Log_ITS_copies)


#ANOVA
Log_ITS_copies_aov<-aov(Log_ITS_copies~Genotype, data=q_Fungal)
summary(Log_ITS_copies_aov)



###################################################################
###############plotting  Log Fungal ITS copies BSA
##################################################################

##Subset the bulk soil
q_Fungal_rhizo<-subset(q_Fungal, Genotype!="Bulk")

q <-ggplot(q_Fungal_rhizo, aes(x=P_A, y=Log_ITS_copies, fill=P_A)) + geom_boxplot() + ggtitle(" Log Fungal ITS copies BSA")+ylim(5,15)
q + geom_jitter( size=5,shape=21, position=position_jitter(0))+ scale_fill_viridis_d()

#Stats Fungal 
shapiro.test(q_Fungal_rhizo$Log_ITS_copies)

#T-test

t.test(Log_ITS_copies~P_A, data=q_Fungal_rhizo)

##End

##############################################################
##############################################################
## 2. Dry weight 
##############################################################
##############################################################


#Import data file all lines
W_results <- read.delim("Dryweights.txt", sep = "\t")

W_results_Q3<-subset(W_results, Soil=="Q3")

W_results_Q4<-subset(W_results, Soil=="Q4")

#############################################################
##Dry weight per genotypes
#############################################################

shapiro.test(W_results_Q4$Dry_w)

kruskal.test(Dry_w~Genotype , data=W_results_Q4)

kwAllPairsDunnTest (x= W_results_Q4$Dry_w, g= W_results_Q4$Genotype, p.adjust.method="BH")


#Order the levels according to a defined order  in Q4
W_results_Q4$Genotype <- ordered(W_results_Q4$Genotype, levels=c("HID144","Barke","124_17","124_52","96_19", "96_56", "Morex", "Golden Promise", "RGT Planet", "Hockett", "FT11"))

p<-ggplot(W_results_Q4, aes(x=Genotype, y=Dry_w, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,0.75)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


p+scale_fill_viridis_d()+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

#Order the levels according to a defined order in Q3
W_results_Q3$Genotype <- ordered(W_results_Q3$Genotype, levels=c("HID144","Barke","124_17","124_52","96_19", "96_56"))

p<-ggplot(W_results_Q3, aes(x=Genotype, y=Dry_w, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,1)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))


p+scale_fill_viridis_d()+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

##########################################################
### BSA Dry_weight 
#########################################################

W_results_Q4_BSA<-subset(W_results_Q4, Genotype!="FT11")

q <-ggplot(W_results_Q4_BSA, aes(x=P_A, y=Dry_w, fill=P_A)) + geom_boxplot() + ggtitle("BSA Dry weight")+ylim(0,0.75)
q + geom_jitter( size=5,shape=21, position=position_jitter(0))+ scale_fill_viridis_d()

#Stats Fungal 
shapiro.test(W_results_Q4_BSA$Dry_w)


kruskal.test(Dry_w~P_A , data=W_results_Q4_BSA)


kwAllPairsDunnTest (x= W_results_Q4$Dry_w, g= W_results_Q4$P_A, p.adjust.method="BH")

##End

#####################################################
#####################################################
##3. Nitrogen fertilization
####################################################
####################################################

#Import data file all lines
Nitrogen_ex <- read.delim("Nitrogen.txt", sep = "\t")

p<-ggplot(Nitrogen_ex, aes(x=Nitrogen, y= Dry.weight, fill=Nitrogen))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,0.75)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

p+scale_fill_viridis_d()+  geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)
p+facet_wrap(~Genotype)


##End
