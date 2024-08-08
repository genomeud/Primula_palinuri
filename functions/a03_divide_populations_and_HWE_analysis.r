# Copyright:	    Fabio Marroni & Paloma Perez 2024
# Aim:              Divide vcf by populations + compute MAF and heterozygosity for each population + histogram with MAF for each population
# To add:		     
# Suggestions: 
# Fixes:  


suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
# here include the sample files to be analysed within the vcf
  make_option(c("-S", "--samples"), type="character", default=NULL, 
              help="Comma separated list of samples to keep for analysis [default= %default]", metavar="character"),
# put the path and file with vcf filtered my missing data
  make_option(c("-V", "--vcffile"), type="character", default="populations.snps.filtered.recode_MIS_filt.vcf", 
              help="input vcf file [default= %default]", metavar="character"),
# In case that we want to randomly take a certain number of individuals for the analysis
  make_option(c("-N", "--sampnumb"), type="numeric", default=0, 
              help="Number of samples to randomly select for each population, 0 for no subsampling [default= %default]", metavar="character"),
# put path and name for the pdf output file including histogram with MAF for each population
  make_option(c("-F", "--pdfplot"), type="character", default="histogram_MAF_pops_ALL.pdf",
              help="output file name for the plot in pdf format [default= %default]", metavar="character")              
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

  if (is.null(opt$sampnumb)) {
  stop("WARNING: No sampnumb specified with '-N' flag.")
} else {  cat ("subsample is ", opt$sampnumb, "\n")
  sampnumb <- opt$sampnumb  
  }
  
    if (is.null(opt$pdfplot)) {
  stop("WARNING: No pdfplot specified with '-F' flag.")
} else {  cat ("pdfplot is ", opt$pdfplot, "\n")
  pdfplot <- opt$pdfplot  
  }

filter_vcf<-function()

{
library(data.table)
#We have a parameter so thar you do not hardcode the file name in the function
# remove all lines that contains ##
cat("Reading vcf...\n")
vcf<-fread(cmd=paste("grep -v ^##", vcffile),data.table=F)
cat("Selecting samples...\n")
# each sample is sepparated by commas
mysample<-unlist(strsplit(samples,","))
# first column header
# firstcol<-c("#CHROM","POS","REF","ALT","QUAL","INFO","FORMAT")
firstcol<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

allcol<-c(firstcol,mysample)
myvcf<-vcf[,names(vcf)%in%allcol]
# myvcf<-myvcf[1:100,]
#Create population table and codes
poptab<-data.frame(Sample=mysample,Pop=unlist(lapply(strsplit(mysample,"-"),"[",3)),stringsAsFactors=F)
#Subsampling if we want to compare samples only with same number
#In the future we should run a series of random sampling
if(sampnumb>0)
{
    #Count number of samples per population
    popnum<-table(poptab$Pop)
    abthr<-names(popnum)[popnum>=sampnumb]
    #Only keep populations with number of samples greater than subsample size
    poptab<-poptab[poptab$Pop%in%abthr,]
    poptab$keep<-0
    #Convoluted way to select random samples for each population
    for(ccc in 1:length(abthr))
    {
        selnames<-sample(poptab$Sample[poptab$Pop==abthr[ccc]],6)
        poptab$keep[poptab$Sample%in%selnames]<-1
    }
    poptab<-poptab[poptab$keep>0,]
    poptab$keep<-NULL
}

fileOut=(pdfplot)
# open pdf output file. 
pdf(fileOut, paper="A4")
popcodes<-unique(poptab$Pop)
toplot<-data.frame(Pop=popcodes,He=NA,stringsAsFactors=F)
for(bbb in 1:length(popcodes))
{
    cat("Processing population",bbb,"out of", length(popcodes),"\n")

    smallsample<-poptab$Sample[poptab$Pop%in%popcodes[bbb]]
    myvcf[,paste0(popcodes[bbb],"_REF")]<-0
    myvcf[,paste0(popcodes[bbb],"_ALT")]<-0
    myvcf[,paste0(popcodes[bbb],"_MIS")]<-0

    for(aaa in 1:length(smallsample))
    {
        tmp<-unlist(lapply(strsplit(myvcf[,smallsample[aaa]],":"),"[",1))
        all1<-unlist(lapply(strsplit(tmp,"/"),"[",1))
        all2<-unlist(lapply(strsplit(tmp,"/"),"[",2))
        myvcf[,paste0(popcodes[bbb],"_REF")]<-apply(cbind(myvcf[,paste0(popcodes[bbb],"_REF")],all1=="0"),1,sum)
        myvcf[,paste0(popcodes[bbb],"_REF")]<-apply(cbind(myvcf[,paste0(popcodes[bbb],"_REF")],all2=="0"),1,sum)
        myvcf[,paste0(popcodes[bbb],"_ALT")]<-apply(cbind(myvcf[,paste0(popcodes[bbb],"_ALT")],all1=="1"),1,sum)
        myvcf[,paste0(popcodes[bbb],"_ALT")]<-apply(cbind(myvcf[,paste0(popcodes[bbb],"_ALT")],all2=="1"),1,sum)
        myvcf[,paste0(popcodes[bbb],"_MIS")]<-apply(cbind(myvcf[,paste0(popcodes[bbb],"_MIS")],all1=="."),1,sum)
        myvcf[,paste0(popcodes[bbb],"_MIS")]<-apply(cbind(myvcf[,paste0(popcodes[bbb],"_MIS")],all2=="."),1,sum)
    }
    myvcf[,paste0(popcodes[bbb],"_N")]<-myvcf[,paste0(popcodes[bbb],"_REF")]+myvcf[,paste0(popcodes[bbb],"_ALT")]+myvcf[,paste0(popcodes[bbb],"_MIS")]
    myvcf[,paste0(popcodes[bbb],"_He")]<-2*(myvcf[,paste0(popcodes[bbb],"_REF")]/myvcf[,paste0(popcodes[bbb],"_N")])*(myvcf[,paste0(popcodes[bbb],"_ALT")]/myvcf[,paste0(popcodes[bbb],"_N")])
    # minor allele column
    myvcf[,paste0(popcodes[bbb],"_minALL")]<-ifelse(myvcf[,paste0(popcodes[bbb],"_REF")] <= myvcf[,paste0(popcodes[bbb],"_ALT")],myvcf[,paste0(popcodes[bbb],"_REF")],myvcf[,paste0(popcodes[bbb],"_ALT")])
    # alt allele freq
    myvcf[,paste0(popcodes[bbb],"_minALLfreq")]<-myvcf[,paste0(popcodes[bbb],"_minALL")]/(myvcf[,paste0(popcodes[bbb],"_REF")]+myvcf[,paste0(popcodes[bbb],"_ALT")])
    myvcf[,paste0(popcodes[bbb],"_minALLfreq")][myvcf[,paste0(popcodes[bbb],"_MIS")]>=0.5*myvcf[,paste0(popcodes[bbb],"_N")]]<-NA
    #For each population, if we have 50% missing data or more, we set He to NA
    myvcf[,paste0(popcodes[bbb],"_He")][myvcf[,paste0(popcodes[bbb],"_MIS")]>=0.5*myvcf[,paste0(popcodes[bbb],"_N")]]<-NA
    toplot$He[bbb]<-mean(myvcf[,paste0(popcodes[bbb],"_He")],na.rm=T)
    # histogram 
    # colors <- c("seagreen","tomato2","hotpink4","steelblue","lightgoldenrod","grey","white")
    hist(myvcf[,paste0(popcodes[bbb],"_minALLfreq")], breaks=10,
    main=(paste(popcodes[bbb]," Alelle frequency", sep="")),
    xlab="Alelle frequency",
    ylab="Total SNPs",
    col="steelblue"
    # freq=FALSE
    )  
    #Select the names of the population specific allele stats (we exploit the fact that we search for pop name and "_")
    popcol<-names(myvcf)[grep(paste0(popcodes[bbb],"_"),names(myvcf))]
    #Create the small vcf
    smallvcf<-myvcf[c(firstcol,popcol)]
    smallvcf$Popname<-popcodes[bbb]
    # output file for each population: the number of reference alleles (POP1_REF); the number of alternative alleles (POP1_ALT); the number of missing alleles (POP1_MIS); the total number of alleles for that population (POP1_N); expected heterozygosity (POP1_He); count for the minor allele (POP1_minALL) and the minor allele frequency (POP1_minALLfreq)
    write.table(smallvcf,file = paste("allele_info_",popcodes[bbb],".tbl", sep=""),quote=F,sep="\t",row.names=F)
}
dev.off()  
# write.table(myvcf,outfile,quote=F,sep="\t",row.names=F)

}

filter_vcf()
