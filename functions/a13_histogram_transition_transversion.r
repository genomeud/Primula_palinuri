# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              Histogram with distribution of transition/transposition
# To add:		     
# Suggestions: 
# Fixes:  

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character",
  default="SNP_type_final.txt",help="Table with SNP classification [default= %default]", metavar="character"),
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
  

tbl <- read.table(infile, header=T)
fileout <- file.path(outpath,"histogram_transition_transposition.jpeg")
jpeg(fileout,width=35,height=20,units="cm",res=300, type="cairo")
par(cex.axis = 1, cex.lab = 1.3)
par(cex.main = 1.4)
par(mar = c(5,5,4,5),oma=c(1.5,0.3,0.3,0.1), mgp=c(3.5,1,0))

barplot(tbl$SNPs , border=T ,names.arg=tbl$type, 
    xlab = "Category", 
    ylab = "Number of SNPs",
    col="skyblue1",
    las=1)
dev.off()

     