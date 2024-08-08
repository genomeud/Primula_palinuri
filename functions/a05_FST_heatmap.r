# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              Heatmap plot with FST from stacks
# To add:		     
# Suggestions: 
# Fixes:  

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

# Required input files are missing fileterd vcf and population file. Last option "to_do" needs to be filled with (pca_analysis,ibd_analysis,FST) depending on the analysis to perform
option_list = list(
  make_option(c("-I", "--infile"), type="character",
  default="NULL",help="stacks vcf file [default= %default]", metavar="character"),
  make_option(c("-O", "--outpath"), type="character",
  default=NULL,help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
  
if (is.null(opt$infile)) {
  stop("WARNING: No input file specified with '-I' flag.")
} else {  cat ("input file ", opt$infile, "\n")
  infile <- opt$infile  
  }
  
if (is.null(opt$outpath)) {
  stop("WARNING: No outpath specified with '-O' flag.")
} else {  cat ("outpath ", opt$outpath, "\n")
  outpath <- opt$outpath  
  }


# --- #
# FST # STACKS
# --- #

fst <- read.table(infile, row.names=1, header=T, check.names = FALSE)

# we need to change triangular matrix into square matrix
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

newfst=tri.to.squ(fst)

# plot heatmap
library(RColorBrewer)
library(gplots)
fileOut="heatmap_Fst_pops_stacks.jpeg"
jpeg(fileOut,width=16,height=15,units="cm",res=300, type="cairo")
par(cex.axis = 0.8, cex.lab = 1)
par(cex.main = 1.1)
par(mar=c(0,2,2.6,0.1), mgp=c(0.5,0.5,0),tck=-0.03,oma=c(1,0.1,1,0.1))
heatmap(as.matrix(newfst), scale="none", xlab="", ylab="", main="FST", col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
legend(x="topleft", legend=c("0","0.5", "1"), fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
dev.off()
