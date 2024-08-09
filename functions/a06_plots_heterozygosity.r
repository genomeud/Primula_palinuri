# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              mean expected and observed heterozygosity according to stacks
# To add:		     
# Suggestions: 
# Fixes:  

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character",
  default="populations.sumstats.tsv",help="input file 'populations.sumstats.tsv' generated with populations comamand from stacks [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)) {
  stop("WARNING: No populations.sumstats.tsv file specified with '-I' flag.")
} else {  cat ("input file ", opt$input, "\n")
  input <- opt$input  
  }



######################################
# heterozygosity according to stacks #
######################################
library(dplyr)
tbl <- read.table(input, header=F)
colnames(tbl)=c("Locus_ID","Chr","BP","Col","Pop_ID","P_Nuc","Q_Nuc","N","P","Obs_Het","Obs_Hom","Exp_Het","Exp_Hom","Pi","Smoothed_Pi","Smoothed_Pi_P-value","Fis","Smoothed_Fis","Smoothed_Fis_P-value","HWE_P-value","Private")
# exp heterozygosity
full <- as.data.frame(tbl[,c(12,5)])
full <- na.omit(full)
colnames(full)=c('EXP_HET', 'NAME')
# use same colors as those used for PCA
color_list <- c("CIM" = "steelblue","FIU" = "orchid2","ID" = "seagreen", "LAM" = "gray68","PC" = "tomato2","PPA" = "lightgoldenrod","SGP" = "hotpink4", "ALL"="white")
# output in pdf format
# fileOut=("boxplot_het_categories_populations_ALL_stacks_sumstats.pdf")
# pdf(fileOut, paper="A4")
# par(cex.axis = 0.9, cex.lab = 1.2)
# par(cex.main = 1.6)
# par(plt = c(0.1, 0.5, 0.1, 0.5), mar = c(6, 4, 4, 3) + 0.1)
# par(mfrow=c(2,1))

# output in jpeg format
fileOut=("boxplot_he_categories_populations_stacks_sumstats.jpeg")
jpeg(fileOut,width=16,height=15,units="cm",res=300, type="cairo")
par(cex.axis = 0.8, cex.lab = 1)
par(cex.main = 1.1)
par(mar = c(5, 6, 6, 4),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.6,0.8,0))

boxplot(full$EXP_HET ~ full$NAME, xlab = "Category",ylab = "Expected heterozygosity Stacks", col = color_list,range = 0.1,main = "Expected heterozygosity",frame = FALSE, top = TRUE, outline=FALSE,ylim=c(0,0.3))
dev.off()

# Mean expected heterozygosity #

library(ggplot2)

# mean exp heterozygosity
library(tidyverse)
EXP_HET_mean_tbl <- (full %>% 
   select(NAME, EXP_HET) %>% 
   group_by(NAME) %>% 
   summarise(EXP_HET_mean = mean(EXP_HET)))
   write.table(EXP_HET_mean_tbl,"mean_exp_HET_pop_stacks.tbl",row.names=F,sep="\t",quote=F)

# mean obs heterozygosity
full <- as.data.frame(tbl[,c(10,5)])
full <- na.omit(full)
colnames(full)=c('OBS_HET', 'NAME')
OBS_HET_mean_tbl <- (full %>% 
   select(NAME, OBS_HET) %>% 
   group_by(NAME) %>% 
   summarise(OBS_HET_mean = mean(OBS_HET)))
   write.table(OBS_HET_mean_tbl,"mean_obs_HET_pop_stacks.tbl",row.names=F,sep="\t",quote=F)

# Barplot
# output in jpeg. mean exp and obs heterozygosity according to
fileOut=("barplot_mean_exp_obs_het_populations_stacks.jpeg")
jpeg(fileOut,width=35,height=18,units="cm",res=300, type="cairo")
par(cex.axis = 1.1, cex.lab = 1.3)
par(cex.main = 1.4)
par(mar = c(4,4,4,5),oma=c(1.5,0.3,0.3,0.1), mgp=c(2.2,0.8,0), mfrow=c(1,2))

barplot(EXP_HET_mean_tbl$EXP_HET_mean , border=T ,names.arg=EXP_HET_mean_tbl$NAME, 
    xlab = "Category", 
    ylab = "Mean Expected heterozygosity",
    col=color_list,
    range = 0.1,
    main = "Mean expected heterozygosity",
    frame = FALSE, top = TRUE,ylim=c(0,0.3))

barplot(OBS_HET_mean_tbl$OBS_HET_mean , border=T ,names.arg=OBS_HET_mean_tbl$NAME, 
    xlab = "Category", 
    ylab = "Mean Observed heterozygosity",
    col=color_list,
    range = 0.1,
    main = "Mean observed heterozygosity",
    frame = FALSE, top = TRUE ,ylim=c(0,0.3))
dev.off() 

# boxplot obs heterozygosity
full <- as.data.frame(tbl[,c(10,5)])
full <- na.omit(full)
colnames(full)=c('OBS_HET', 'NAME')
# output in jpeg format
fileOut=("boxplot_ho_categories_populations_stacks_sumstats.jpeg")
jpeg(fileOut,width=16,height=15,units="cm",res=300, type="cairo")
par(cex.axis = 0.8, cex.lab = 1)
par(cex.main = 1.1)
par(mar = c(5, 6, 6, 4),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.6,0.8,0))

boxplot(full$OBS_HET ~ full$NAME, xlab = "Category",ylab = "Observed heterozygosity", col = color_list,range = 0.1,main = "Observed heterozygosity",frame = FALSE, top = TRUE, outline=FALSE,ylim=c(0,0.3))
dev.off()
