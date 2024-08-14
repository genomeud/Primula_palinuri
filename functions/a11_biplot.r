# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              Biplot analysis with all SNPs and top 40,50,60,80 and 100 SNPs  responsible for differentiating populations. Extract Individual genotypes of the top SNPs with the most extreme coordinates on the rotated system
# To add:		     
# Suggestions: 
# Fixes: 

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

# input files
option_list = list(
  make_option(c("-I", "--allelefreqtbl"), type="character",
  default="/de_novo_pipeline/populations_stacks/biplot/transposed_ind_SNPs_no_NA.tbl",help="allele frequency table file [default= %default]", metavar="character"),
  make_option(c("-M", "--metadata"), type="character",
  default="/raw_reads/sample_ID_and_populations.txt",help="metadata with individual IDs and pop names [default= %default]", metavar="character"),
  make_option(c("-O", "--outpath"), type="character",
  default="/de_novo_pipeline/populations_stacks/biplot/",help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

  
if (is.null(opt$allelefreqtbl)) {
  stop("WARNING: No allele frequency table file specified with '-I' flag.")
} else {  cat ("allelefreqtbl file ", opt$allelefreqtbl, "\n")
  allelefreqtbl <- opt$allelefreqtbl  
  }

if (is.null(opt$metadata)) {
  stop("WARNING: No metadata file specified with '-M' flag.")
} else {  cat ("metadata file ", opt$metadata, "\n")
  metadata <- opt$metadata  
  }
  
if (is.null(opt$outpath)) {
  stop("WARNING: No outpath specified with '-O' flag.")
} else {  cat ("outpath ", opt$outpath, "\n")
  outpath <- opt$outpath  
  }
  

##########
# Biplot #
##########

library(data.table)

infile=allelefreqtbl

mysnp<-fread(infile,data.table=F,fill=TRUE)
row.names(mysnp)<-mysnp[,1]
mysnp[,1]<-NULL
#Remove columns containing NAs. (is already done in bash because otherwise it takes very long)
myNA<-apply(apply(mysnp,2,is.na),2,sum)
mysnp<-mysnp[,myNA<1]
geno<-apply(mysnp,2,sum)
mysnp<-mysnp[,geno>0]

ori.pca <- prcomp(mysnp, center = T, scale = T)

# metadata with pop info for each ind
metadata <- read.table(metadata, header=T)
colnames(metadata) <- NULL
colnames(metadata) <- c("Sample_name","POP")
# assign colors to each population
col4="gray68"
col5="orchid2"
col6="seagreen"
col7="tomato2"
col8="hotpink4"
col9="steelblue"
col10="lightgoldenrod"
metadata$color_group[metadata$POP=="LAM"]<-col4
metadata$color_group[metadata$POP=="FIU"]<-col5
metadata$color_group[metadata$POP=="ID"]<-col6
metadata$color_group[metadata$POP=="PC"]<-col7
metadata$color_group[metadata$POP=="SGP"]<-col8
metadata$color_group[metadata$POP=="CIM"]<-col9
metadata$color_group[metadata$POP=="PPA"]<-col10
########################
# Biplot with ALL SNPs #
########################
pca_data <- data.frame(Sample_name = rownames(ori.pca$x), ori.pca$x[,1:2])
merged_data <- merge(pca_data, metadata, by = "Sample_name")
top.ori.one<-sort(abs(ori.pca$rotation[,1]),decreasing=T)
top.ori.two<-sort(abs(ori.pca$rotation[,2]),decreasing=T)
sel.snp<-unique(c(names(top.ori.one),names(top.ori.two)))
top.rotation<-ori.pca$rotation[sel.snp,]
    outplot <- file.path(outpath,"biplot.jpeg")
    jpeg(outplot,width=11,height=11,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(2.5,2.5,2,2),oma=c(0.6,0.1,0.2,0.1), mgp=c(1.6,0.6,0))
scaling_factor <- 9
plot(merged_data$PC1, merged_data$PC2, 
     xlab = "PC1", ylab = "PC2", 
     main = "Biplot",
     pch = 16,cex=0.7, col = merged_data$color_group)
arrows(0, 0, scaling_factor * top.rotation[,1] * max(merged_data$PC1), 
       scaling_factor * top.rotation[,2] * max(merged_data$PC2), 
       col = "red", length = 0.1)
text(top.rotation[,1] * max(merged_data$PC1), 
     top.rotation[,2] * max(merged_data$PC2), 
     labels = NULL, 
     col = "red", pos = 4)
     # Extract unique populations and their corresponding colors
unique_pops <- unique(merged_data$POP)
unique_colors <- unique(merged_data$color_group)
# Add a legend
legend("bottomright", legend = unique_pops, col = unique_colors, pch = 16,cex=0.6, title = "Populations")
dev.off()

# additional plot with individual names and SNPs_ID
outplot <- file.path(outpath,"biplot_names.pdf")
pdf(outplot)
biplot(ori.pca,
       col=c('palegreen4', 'tomato2'),
       cex=c(0.6, 0.5),
       xlim=c(-.2, .2),
       main="PCA Biplot",
       xlab='PC1',
       ylab='PC2'
       )    
dev.off()

###############
# top 40 SNPs #
###############
top.ori.one<-sort(abs(ori.pca$rotation[,1]),decreasing=T)[1:20]
top.ori.two<-sort(abs(ori.pca$rotation[,2]),decreasing=T)[1:20]
sel.snp<-unique(c(names(top.ori.one),names(top.ori.two)))
top.rotation<-ori.pca$rotation[sel.snp,]

# metadata with pop info for each ind
pca_data <- data.frame(Sample_name = rownames(ori.pca$x), ori.pca$x[,1:2])
merged_data <- merge(pca_data, metadata, by = "Sample_name")
    outplot <- file.path(outpath,"biplot_top40SNPs.jpeg")
    jpeg(outplot,width=11,height=11,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(2.5,2.5,2,2),oma=c(0.6,0.1,0.2,0.1), mgp=c(1.6,0.6,0))
scaling_factor <- 9
plot(merged_data$PC1, merged_data$PC2, 
     xlab = "PC1", ylab = "PC2", 
     main = "Biplot Top 40 SNPs",
     pch = 16,cex=0.7, col = merged_data$color_group)
arrows(0, 0, scaling_factor * top.rotation[,1] * max(merged_data$PC1), 
       scaling_factor * top.rotation[,2] * max(merged_data$PC2), 
       col = "red", length = 0.1)
text(top.rotation[,1] * max(merged_data$PC1), 
     top.rotation[,2] * max(merged_data$PC2), 
     labels = NULL, 
     col = "red", pos = 4)
     # Extract unique populations and their corresponding colors
unique_pops <- unique(merged_data$POP)
unique_colors <- unique(merged_data$color_group)
# Add a legend
legend("bottomright", legend = unique_pops, col = unique_colors, pch = 16,cex=0.6, title = "Populations")
dev.off()

top.rotation <- cbind(SNP_ID = rownames(top.rotation), top.rotation)
output_file <- file.path(outpath, "top_40SNPs.tbl")
write.table(top.rotation,output_file,quote=F,row.names=F,col.names=T,sep="\t") 

###############
# top 50 SNPs #
###############
top.ori.one<-sort(abs(ori.pca$rotation[,1]),decreasing=T)[1:25]
top.ori.two<-sort(abs(ori.pca$rotation[,2]),decreasing=T)[1:25]
sel.snp<-unique(c(names(top.ori.one),names(top.ori.two)))
top.rotation<-ori.pca$rotation[sel.snp,]

# metadata with pop info for each ind
pca_data <- data.frame(Sample_name = rownames(ori.pca$x), ori.pca$x[,1:2])
merged_data <- merge(pca_data, metadata, by = "Sample_name")
    outplot <- file.path(outpath,"biplot_top50SNPs.jpeg")
    jpeg(outplot,width=11,height=11,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(2.5,2.5,2,2),oma=c(0.6,0.1,0.2,0.1), mgp=c(1.6,0.6,0))
scaling_factor <- 9
plot(merged_data$PC1, merged_data$PC2, 
     xlab = "PC1", ylab = "PC2", 
     main = "Biplot Top 50 SNPs",
     pch = 16,cex=0.7, col = merged_data$color_group)
arrows(0, 0, scaling_factor * top.rotation[,1] * max(merged_data$PC1), 
       scaling_factor * top.rotation[,2] * max(merged_data$PC2), 
       col = "red", length = 0.1)
text(top.rotation[,1] * max(merged_data$PC1), 
     top.rotation[,2] * max(merged_data$PC2), 
     labels = NULL, 
     col = "red", pos = 4)
     # Extract unique populations and their corresponding colors
unique_pops <- unique(merged_data$POP)
unique_colors <- unique(merged_data$color_group)
# Add a legend
legend("bottomright", legend = unique_pops, col = unique_colors, pch = 16,cex=0.6, title = "Populations")
dev.off()

top.rotation <- cbind(SNP_ID = rownames(top.rotation), top.rotation)
output_file <- file.path(outpath, "top_50SNPs.tbl")
write.table(top.rotation,output_file,quote=F,row.names=F,col.names=T,sep="\t") 

###############
# top 60 SNPs #
###############
top.ori.one<-sort(abs(ori.pca$rotation[,1]),decreasing=T)[1:30]
top.ori.two<-sort(abs(ori.pca$rotation[,2]),decreasing=T)[1:30]
sel.snp<-unique(c(names(top.ori.one),names(top.ori.two)))
top.rotation<-ori.pca$rotation[sel.snp,]

# metadata with pop info for each ind
pca_data <- data.frame(Sample_name = rownames(ori.pca$x), ori.pca$x[,1:2])
merged_data <- merge(pca_data, metadata, by = "Sample_name")
    outplot <- file.path(outpath,"biplot_top60SNPs.jpeg")
    jpeg(outplot,width=11,height=11,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(2.5,2.5,2,2),oma=c(0.6,0.1,0.2,0.1), mgp=c(1.6,0.6,0))
scaling_factor <- 9
plot(merged_data$PC1, merged_data$PC2, 
     xlab = "PC1", ylab = "PC2", 
     main = "Biplot Top 60 SNPs",
     pch = 16,cex=0.7, col = merged_data$color_group)
arrows(0, 0, scaling_factor * top.rotation[,1] * max(merged_data$PC1), 
       scaling_factor * top.rotation[,2] * max(merged_data$PC2), 
       col = "red", length = 0.1)
text(top.rotation[,1] * max(merged_data$PC1), 
     top.rotation[,2] * max(merged_data$PC2), 
     labels = NULL, 
     col = "red", pos = 4)
     # Extract unique populations and their corresponding colors
unique_pops <- unique(merged_data$POP)
unique_colors <- unique(merged_data$color_group)
# Add a legend
legend("bottomright", legend = unique_pops, col = unique_colors, pch = 16,cex=0.6, title = "Populations")
dev.off()

top.rotation <- cbind(SNP_ID = rownames(top.rotation), top.rotation)
output_file <- file.path(outpath, "top_60SNPs.tbl")
write.table(top.rotation,output_file,quote=F,row.names=F,col.names=T,sep="\t") 

###############
# top 80 SNPs #
###############
top.ori.one<-sort(abs(ori.pca$rotation[,1]),decreasing=T)[1:40]
top.ori.two<-sort(abs(ori.pca$rotation[,2]),decreasing=T)[1:40]
sel.snp<-unique(c(names(top.ori.one),names(top.ori.two)))
top.rotation<-ori.pca$rotation[sel.snp,]

# metadata with pop info for each ind
pca_data <- data.frame(Sample_name = rownames(ori.pca$x), ori.pca$x[,1:2])
merged_data <- merge(pca_data, metadata, by = "Sample_name")
    outplot <- file.path(outpath,"biplot_top80SNPs.jpeg")
    jpeg(outplot,width=11,height=11,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(2.5,2.5,2,2),oma=c(0.6,0.1,0.2,0.1), mgp=c(1.6,0.6,0))
scaling_factor <- 9
plot(merged_data$PC1, merged_data$PC2, 
     xlab = "PC1", ylab = "PC2", 
     main = "Biplot Top 80 SNPs",
     pch = 16,cex=0.7, col = merged_data$color_group)
arrows(0, 0, scaling_factor * top.rotation[,1] * max(merged_data$PC1), 
       scaling_factor * top.rotation[,2] * max(merged_data$PC2), 
       col = "red", length = 0.1)
text(top.rotation[,1] * max(merged_data$PC1), 
     top.rotation[,2] * max(merged_data$PC2), 
     labels = NULL, 
     col = "red", pos = 4)
     # Extract unique populations and their corresponding colors
unique_pops <- unique(merged_data$POP)
unique_colors <- unique(merged_data$color_group)
# Add a legend
legend("bottomright", legend = unique_pops, col = unique_colors, pch = 16,cex=0.6, title = "Populations")
dev.off()

top.rotation <- cbind(SNP_ID = rownames(top.rotation), top.rotation)
output_file <- file.path(outpath, "top_80SNPs.tbl")
write.table(top.rotation,output_file,quote=F,row.names=F,col.names=T,sep="\t") 

################
# top 100 SNPs #
################
top.ori.one<-sort(abs(ori.pca$rotation[,1]),decreasing=T)[1:50]
top.ori.two<-sort(abs(ori.pca$rotation[,2]),decreasing=T)[1:50]
sel.snp<-unique(c(names(top.ori.one),names(top.ori.two)))
top.rotation<-ori.pca$rotation[sel.snp,]

# metadata with pop info for each ind
pca_data <- data.frame(Sample_name = rownames(ori.pca$x), ori.pca$x[,1:2])
merged_data <- merge(pca_data, metadata, by = "Sample_name")
    outplot <- file.path(outpath,"biplot_top100SNPs.jpeg")
    jpeg(outplot,width=11,height=11,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(2.5,2.5,2,2),oma=c(0.6,0.1,0.2,0.1), mgp=c(1.6,0.6,0))
scaling_factor <- 9
plot(merged_data$PC1, merged_data$PC2, 
     xlab = "PC1", ylab = "PC2", 
     main = "Biplot Top 100 SNPs",
     pch = 16,cex=0.7, col = merged_data$color_group)
arrows(0, 0, scaling_factor * top.rotation[,1] * max(merged_data$PC1), 
       scaling_factor * top.rotation[,2] * max(merged_data$PC2), 
       col = "red", length = 0.1)
text(top.rotation[,1] * max(merged_data$PC1), 
     top.rotation[,2] * max(merged_data$PC2), 
     labels = NULL, 
     col = "red", pos = 4)
     # Extract unique populations and their corresponding colors
unique_pops <- unique(merged_data$POP)
unique_colors <- unique(merged_data$color_group)
# Add a legend
legend("bottomright", legend = unique_pops, col = unique_colors, pch = 16,cex=0.6, title = "Populations")
dev.off()

top.rotation <- cbind(SNP_ID = rownames(top.rotation), top.rotation)
output_file <- file.path(outpath, "top_100SNPs.tbl")
write.table(top.rotation,output_file,quote=F,row.names=F,col.names=T,sep="\t") 

