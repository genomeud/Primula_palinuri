# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              Procrustes analysis was used to compute the similarity score, between the first two PCs of genetic variation and the geographic coordinates of the populations
# To add:		     
# Suggestions: 
# Fixes: 

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-C", "--geogcoords"), type="character",
  default="coordinates_individuals.txt",help="geographic coordinates for each individual [default= %default]", metavar="character"),
  make_option(c("-G", "--genetcoords"), type="character",
  default="genetic_coordinates.txt",help="genetic coordinates including PC1 and PC2 generated with SNPrelate [default= %default]", metavar="character"),
  make_option(c("-M", "--metadata"), type="character",
  default="raw_reads/sample_ID_and_populations.txt",help="metadata with individual IDs and pop names [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$geogcoords)) {
  stop("WARNING: No geographic coords file specified with '-C' flag.")
} else {  cat ("geographic coordinates for each individual ", opt$geogcoords, "\n")
  geogcoords <- opt$geogcoords  
  }

if (is.null(opt$genetcoords)) {
  stop("WARNING: No genetic coordinates specified with '-G' flag.")
} else {  cat ("genetcoords ", opt$genetcoords, "\n")
  genetcoords <- opt$genetcoords  
  }

if (is.null(opt$metadata)) {
  stop("WARNING: No metadata file specified with '-M' flag.")
} else {  cat ("metadata file ", opt$metadata, "\n")
  metadata <- opt$metadata  
  }
  
library(vegan)

coords <- read.table(geogcoords, header = TRUE)
gen_coords <- read.table(genetcoords, header = TRUE)

# Ensure both datasets have the same order of individuals
coords <- coords[order(coords$ID), ]
gen_coords <- gen_coords[order(gen_coords$ID), ]

# Remove the ID column for analysis
coords_matrix <- as.matrix(coords[, -1])
gen_coords_matrix <- as.matrix(gen_coords[, -1])

# Create a distance matrix for both coordinates and genetic data
coords_dist <- dist(coords_matrix)
gen_coords_dist <- dist(gen_coords_matrix)

# Perform Procrustes analysis
# proc_result <- procrustes(coords_matrix, gen_coords_matrix)
# tres=proc_result

protest_result <- protest(coords_matrix, gen_coords_matrix)

# plot(protest_result)

tres=protest_result

# metadata with pop info for each ind
metadata <- read.table(metadata, header=T)
colnames(metadata) <- NULL
colnames(metadata) <- c("Sample_name","POP")

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

matful<-data.frame(tres$X)
    matred<-data.frame(tres$Yrot)
    matful$name<-coords$ID

    matful$colore<-metadata$color_group[match(matful$name,metadata$Sample_name)]
    myxlim<-c(min(matful$Latitude,matred$X1),max(matful$Latitude,matred$X1))
    myylim<-c(min(matful$Longitude,matred$X2),max(matful$Longitude,matred$X2))
    
    outplot="procrustes.jpeg"
    jpeg(outplot,width=10,height=10,units="cm",res=300,type="cairo")
    par(cex.axis = 0.6, cex.lab = 0.75)
    par(cex.main = 1)
    par(mar = c(3,4,2,2),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.6,0.6,0))
    
    plot(matful$Latitude,matful$Longitude,type="n",xlab="DIM1",ylab="DIM2",xlim=myxlim,
        ylim=myylim)
    arrows(matful$Latitude,matful$Longitude,matred$X1,matred$X2,col=matful$colore,code=2,length=0.04)
    text(mean(myxlim),0.9*myylim[2],label=bquote(italic(r) == .(tres$t0)),cex=0.6)
    legend("topright",border=NA,bty="n",fill=c("gray68","orchid2","seagreen","tomato2","hotpink4","steelblue","lightgoldenrod"),legend=c("LAM","FIU","ID","PC","SGP","CIM","PPA"), col=c("gray68","orchid2","seagreen","tomato2","hotpink4","steelblue","lightgoldenrod"), text.col=c("gray68","orchid2","seagreen","tomato2","hotpink4","steelblue","lightgoldenrod"),cex=0.7)
    dev.off()
