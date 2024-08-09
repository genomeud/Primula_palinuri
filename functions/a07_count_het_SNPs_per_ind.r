# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              count heterozygous SNPs, informative SNPs and calculate ratio for each individual present in the vcf file. Do a scatterplot HE_SNPs/informative SNPs vs informative SNPs/ total SNPs and boxplot with heterozygosity per population
# To add:		     
# Suggestions: 
# Fixes: 

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-V", "--vcffile"), type="character",
  default="populations.snps.filtered.recode_MIS_filt_header.vcf",help="vcf file including the header [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character",
  default="genotype_counts.tbl",help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vcffile)) {
  stop("WARNING: No vcf file specified with '-V' flag.")
} else {  cat ("vcf file ", opt$vcffile, "\n")
  vcffile <- opt$vcffile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


library(vcfR)

# Read vcf file
invcf <- vcffile
vcf <- read.vcfR(invcf)

# extract genotypes
gt <- extract.gt(vcf)

# Function to count heterozygous genotypes
count_het <- function(genotypes) {
  sum(grepl("^0[/]1$|^1[/]0$", genotypes))
}

# Function to count informative SNPs (non-missing)
count_informative <- function(genotypes) {
  sum(grepl("^0[/]1$|^1[/]0$|^0[/]0$|^1[/]1$", genotypes))
}

# Function to count ALL SNPs (even missing)
count_all <- function(genotypes) {
  length(genotypes)
}

# Function to extract the third argument from the individual name (population name)
extract_third_argument <- function(name) {
  parts <- strsplit(name, "-")[[1]]
  if (length(parts) >= 3) {
    return(parts[3])
  } else {
    return(NA)
  }
}

# Apply the functions to each individual
het_counts <- apply(gt, 2, count_het)
informative_counts <- apply(gt, 2, count_informative)
total_counts <- apply(gt, 2, count_all)
individual_names <- colnames(gt)
third_arguments <- sapply(individual_names, extract_third_argument)

# average het SNPs
ratio <- het_counts / informative_counts

result <- data.frame(
  Individual = colnames(gt),
  Population = third_arguments,
  Heterozygous_Count = het_counts,
  Informative_SNP_Count = informative_counts,
  Total_SNP_Count = total_counts,
  Heterozygous_Ratio = ratio
)

# sort by ID
# result <- result[order(Individual),]
write.table(result, outfile, row.names = FALSE, quote=FALSE,sep="\t")

##########################################################################################
# scatterplot with ratio He_count / informative count vs informative count / Total count #
##########################################################################################

#define different colors for the Populations (same colors used before for the other plots)

col4="gray68"
col5="orchid2"
col6="seagreen"
col7="tomato2"
col8="hotpink4"
col9="steelblue"
col10="lightgoldenrod"
result$color<-"black"

result$color[result$Population=="LAM"]<-col4
result$color[result$Population=="FIU"]<-col5
result$color[result$Population=="ID"]<-col6
result$color[result$Population=="PC"]<-col7
result$color[result$Population=="SGP"]<-col8
result$color[result$Population=="CIM"]<-col9
result$color[result$Population=="PPA"]<-col10

# create a column with the ratios. Lo pongo en escala logaritmica que sino queda feo. first ratio HE/informative SNPs and second ratio informative SNPs /total count SNPs
result$He_ratio<-log10(result$Heterozygous_Count/result$Informative_SNP_Count)
result$snp_ratio<-log10(result$Informative_SNP_Count/result$Total_SNP_Count)
#Regression line
z <- lm(result$He_ratio ~ result$snp_ratio)
# rsquare
rsquare<-round(summary(z)$adj.r.squared,3)
# calculate p-value
pvalue<-cor.test(result$He_ratio,result$snp_ratio)$p.value
fileOut=("scatterplot_ratio_heterozygosity.jpeg")
jpeg(outplot,width=10,height=10,units="cm",res=300,type="cairo")
par(cex.axis = 0.7, cex.lab = 0.7)
par(cex.main = 1)
par(mar = c(3,4,2,2),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.6,0.6,0))
plot(result$He_ratio,result$snp_ratio,xlab="Log10(Heterozygous SNPs / Informative SNPs)",ylab="Log10(Informative SNPs / Total SNPs)",pch=21,las=1,cex=0.3, col=result$color,bg=result$color, main="Heterozygosity Ratio") 
legend ("bottomright",  ncol = 1, legend = c("LAM","FIU","ID","PC","SGP","CIM","PPA"), col=c("gray68","orchid2","seagreen","tomato2","hotpink4","steelblue","lightgoldenrod"), text.col=c("gray68","orchid2","seagreen","tomato2","hotpink4","steelblue","lightgoldenrod"), x.intersp = 0.5 , y.intersp = 1 , pt.cex = 1 , cex = 0.6, bty = "n" )
text(x=-1.11,y=-0.33,labels=bquote(R^2==~.(rsquare)),cex=0.5)
text(x=-1.11,y=-0.35,labels=bquote(p==~.(pvalue)),cex=0.5)
dev.off()


################################################################
# boxplot with heterozygosity distribution for each population #
################################################################

library(ggplot2)
color_list=c("steelblue","orchid2","seagreen","gray68","tomato2","lightgoldenrod","hotpink4")
fileOut=("boxplot_he_ratio.jpeg")
jpeg(fileOut,width=8,height=8,units="cm",res=300, type="cairo")
par(cex.axis = 0.8, cex.lab = 1)
par(cex.main = 1.1)
par(mar = c(5, 6, 6, 4),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.6,0.8,0))
ggplot(result, aes(x=Population, y=Heterozygous_Ratio, fill=Population)) + 
    geom_boxplot(outlier.shape = NA) +
    ylim(0, 0.3) +
    theme(text = element_text(size = 6)) +
    xlab("Populations") +
    ylab("Heterozygous positions / Informative positions") +
    ggtitle("Heterozygosity Ratio") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=8),
    axis.text.x = element_text(angle=90, hjust=1),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
    scale_fill_manual(values = color_list)
dev.off()

# Boxplot including wilcoxon test. This boxplot include letters depending on the distribution observed. Same letters means similar distribution according to wilcoxon test
library(dplyr)
library(multcompView)

full <- result[,c(6,2)]

tri.to.squ<-function(x)
{
rn<-row.names(x)
cn<-colnames(x)
an<-unique(c(cn,rn))
myval<-x[!is.na(x)]
mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
for(ext in 1:length(cn))
{
 for(int in 1:length(rn))
 {
 if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
 mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 }
}
return(mymat)
}

library(dplyr)

pp<-pairwise.wilcox.test(full[,1], full[,2], p.adjust.method = "none", paired = FALSE)
mymat<-tri.to.squ(pp$p.value)
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)

fileOut=("boxplot_he_ratio_wilcx.jpeg")
jpeg(fileOut,width=16,height=16,units="cm",res=300, type="cairo")
par(cex.axis = 0.9, cex.lab = 1.1)
par(cex.main = 1.4)
par(mar = c(5, 6, 4, 4),oma=c(1.5,0.3,0.3,0.1), mgp=c(2.2,0.8,0))
color_list=c("steelblue","orchid2","seagreen","gray68","tomato2","lightgoldenrod","hotpink4")

boxplot(full[,1] ~ full[,2],ylab="Heterozygous positions / Informative positions",xlab="Population",col=color_list,names=c("LAM","FIU","ID","PC","SGP","CIM","PPA"),main=c("Heterozygosity Ratio"),ylim=c(min(full[,1]),0.1+max(full[,1])))
text(c(1,2,3,4,5,6,7,8),0.1+max(full[,1]),c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3],myletters$Letters[4],myletters$Letters[5],myletters$Letters[6],myletters$Letters[7],myletters$Letters[8]))
dev.off()

