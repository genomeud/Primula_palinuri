# Copyright:	    Fabio Marroni & Paloma Perez 2024
# Aim:              Compute allele frequency and heterozygosity considering ALL the individuals + histogram with minor allele frequency
# To add:		     
# Suggestions: 
# Fixes:  


suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
# put the path and tblfile with allele info
  make_option(c("-F", "--tblfile"), type="character", default="/projects/marroni/seedforce/primula_palinuri/github/tests/allele_info_primula_palinuri_MinAll.tbl", 
              help="path input tblfile with minor allele previously generated [default= %default]", metavar="character"),
# put path and name for the output tblfile
  make_option(c("-O", "--outfile"), type="character", default="/projects/marroni/seedforce/primula_palinuri/github/tests/allele_freq_and_exp_het.tbl", 
              help="output tblfile name [default= %default]", metavar="character"),   
# put path and name for the output histogram in jpeg format
  make_option(c("-P", "--outplot"), type="character", default="/projects/marroni/seedforce/primula_palinuri/github/tests/histogram_MAF.jpeg", 
              help="output plot name in jpeg format [default= %default]", metavar="character")              
              
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


  if (is.null(opt$tblfile)) {
  stop("WARNING: No tblfile specified with '-F' flag.")
} else {  cat ("tblfile is ", opt$tblfile, "\n")
  tblfile <- opt$tblfile  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

  if (is.null(opt$outplot)) {
  stop("WARNING: No outplot specified with '-O' flag.")
} else {  cat ("outplot is ", opt$outplot, "\n")
  outplot <- opt$outplot  
  }
  
# Hardy-Weinberg-Equilibrium (HWE)
# p^2 +2pq + q^2 = 1
# p=frequency of allele homozygous dominant
# q=frequency of allele homozygous recesive (p-1)
# 2pq=frequency of heterozygous 

# Input table
# #CHROM   POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT          N_REF  N_ALT  N_MIS   MinAll
# 1        298  .   A    C    .     PASS    .     GT:DP:AD:GQ:GL  10     2      48      2
# 1        299  .   A    T    .     PASS    .     GT:DP:AD:GQ:GL  45     5      10      5
# 2        266  .   C    A    .     PASS    .     GT:DP:AD:GQ:GL  48     2      10      2
# 3        109  .   T    C    .     PASS    .     GT:DP:AD:GQ:GL  35     7      23      7
# 3        265  .   T    G    .     PASS    .     GT:DP:AD:GQ:GL  10     50      0      10


# N_REF = nº of alleles like REF. 
# N_ALT = nº of alleles like ALT. 
# N_MIS = nº of alleles MISSING. 
# MinAll = minor allele


##################################################################
# # # # ALLELE FREQUENCY AND HETEROZYGOSITY IN INDIVIDUALS # # # # ALL
##################################################################
# Use minor allele frequency!! From 0 to 0.5

HWE<-function(tblfile,outfile,outplot)
{

# browser()
allele <- tblfile

tbl = read.table(allele, header=F)
# add header 
names(tbl)=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "N_REF", "N_ALT", "N_MIS", "MinAll")
# filter inds with too many missing data
tbl<-tbl[tbl$N_MIS<(tbl$N_REF+tbl$N_ALT),]
# # minor allele frequency
min_all_FREQ = (tbl$MinAll/(tbl$N_ALT+tbl$N_REF))
# paste new column in the previous table and create a new tbl2
tbl2 <- cbind(tbl, min_all_FREQ)
# calculate Major allele frequency (calculate q knowing p) p + q = 1; q = 1 - p
# p=min_all_FREQ
# q=1-min_all_FREQ
max_all_FREQ = (1 - tbl2$min_all_FREQ)
# paste new column in the previous table and creates a new table 3
tbl3 <- cbind(tbl2, max_all_FREQ)
# calculate expected heterozygosity. 2pq
N_HET_EXP = (2 * tbl3$min_all_FREQ * tbl3$max_all_FREQ)
# paste new column in the previous table and creates a new table 4
tbl4 <- cbind(tbl3, N_HET_EXP)
# save the final table in tblfile called 09_mixed_01.txt
write.table(tbl4, file = outfile, quote=FALSE, row.names=FALSE, sep="\t")

# otput format
# CHROM   POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT          N_REF  N_ALT  N_MIS  MinAll  min_all_FREQ         max_all_FREQ       N_HET_EXP
# 2       30   .   G    T    .     PASS    .     GT:DP:AD:GQ:GL  163    1      18     1       0.00609756097560976  0.99390243902439   0.012120761451517
# 2       96   .   C    T    .     PASS    .     GT:DP:AD:GQ:GL  146    18     18     18      0.109756097560976    0.890243902439024  0.195419393218322

# min_all_FREQ = minor allele frequency
# max_all_FREQ = max allele frequency
# N_HET_EXP = expected heterozygosity under Hardy-Weinberg equilibrium

######################
# histogram with MAF #
######################

library(gplots)

fileOut=outplot
jpeg(fileOut,width=16,height=15,units="cm",res=300, type="cairo")
par(cex.axis = 1.2, cex.lab = 1.5)
par(cex.main = 1.6)
par(mar = c(2.9, 3.2, 2, 1),oma=c(1.5,0.2,0.3,0.1), mgp=c(1.8,0.5,0))
# histogram with added parameters
hist(tbl4$min_all_FREQ, breaks=10,
main=c("Minor alelle frequency"),
xlab="Minor alelle frequency",
ylab="Total SNPs",
col="seagreen"
# freq=FALSE
)
dev.off()

}

HWE(tblfile=tblfile,outfile=outfile,outplot=outplot)  
  