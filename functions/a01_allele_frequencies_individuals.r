# Copyright:	    Fabio Marroni & Paloma Perez 2024
# Aim:              Count how many times we have 
# To add:		     
# Suggestions: 
# Fixes:  

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
# here include the sample IDs to be analysed within the vcf
  make_option(c("-S", "--samples"), type="character", default="NULL", 
              help="Comma separated list of samples to keep for analysis [default= %default]", metavar="character"),
# put the path and file with vcf for individuals approach
  make_option(c("-V", "--vcffile"), type="character", default="NULL", 
              help="path to vcf input file [default= %default]", metavar="character"),
# put path and name for the output file
  make_option(c("-O", "--outfile"), type="character", default="allele_frequency.tbl", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


  if (is.null(opt$samples)) {
  stop("WARNING: No samples specified with '-S' flag.")
} else {  cat ("samples is ", opt$samples, "\n")
  samples <- opt$samples  
  }

  if (is.null(opt$vcffile)) {
  stop("WARNING: No vcffile specified with '-V' flag.")
} else {  cat ("vcffile is ", opt$vcffile, "\n")
  vcffile <- opt$vcffile  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }
  

filter_vcf<-function()
{
library(data.table)
#We have a parameter so that you do not hardcode the file name in the function
# remove all lines that contains ##
vcf<-fread(cmd=paste("grep -v ^##", vcffile),data.table=F)
cat("Selecting samples...\n")
# each sample is sepparated by commas
mysample<-unlist(strsplit(samples,","))
# first column header
# firstcol<-c("#CHROM","POS","REF","ALT","QUAL","INFO","FORMAT")
firstcol<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

allcol<-c(firstcol,mysample)
myvcf<-vcf[,names(vcf)%in%allcol]
# orig<-myvcf
# myvcf<-myvcf[1:100,]
# create 3 columns one for ind REF, one for ind_ALT and one for ind missing data 
myvcf$N_REF<-0
myvcf$N_ALT<-0
myvcf$N_MIS<-0

# for each sample in my sample do the following loop
for(aaa in 1:length(mysample))
{
# get the haplotype info which could be 0/0 when both haplotypes are like ref; 1/1 when both are like alt and 0/1 when we have 1 ref and 1 alt. Then we will sum how many times we have one or the other (ref or alt) and how many times is missing ./.
tmp<-unlist(lapply(strsplit(myvcf[,mysample[aaa]],":"),"[",1))
all1<-unlist(lapply(strsplit(tmp,"/"),"[",1))
all2<-unlist(lapply(strsplit(tmp,"/"),"[",2))
myvcf$N_REF<-apply(cbind(myvcf$N_REF,all1=="0"),1,sum)
myvcf$N_REF<-apply(cbind(myvcf$N_REF,all2=="0"),1,sum)
myvcf$N_ALT<-apply(cbind(myvcf$N_ALT,all1=="1"),1,sum)
myvcf$N_ALT<-apply(cbind(myvcf$N_ALT,all2=="1"),1,sum)
myvcf$N_MIS<-apply(cbind(myvcf$N_MIS,all1=="."),1,sum)
myvcf$N_MIS<-apply(cbind(myvcf$N_MIS,all2=="."),1,sum)
}
myvcf<-myvcf[,!names(myvcf)%in%mysample]

write.table(myvcf,outfile,quote=F,sep="\t",row.names=F)
cat("Done\n")
}

filter_vcf()
