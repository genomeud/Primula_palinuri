# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              mantel test between genetic distances (measured by fst) and geographic coordinates and create scatterplot with pairwise comparisons
# To add:		     
# Suggestions: 
# Fixes: 

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-F", "--genfile"), type="character",
  default="populations.fst_summary_for_heatmap.tsv",help="FST matrix with FST between populations [default= %default]", metavar="character"),
  make_option(c("-G", "--geofile"), type="character",
  default="pop_geographic_coordinates.txt",help="Geographic coordinates for each population[default= %default]", metavar="character"),
  make_option(c("-O", "--outpath"), type="character",
  default="NULL",help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genfile)) {
  stop("WARNING: No FST file specified with '-F' flag.")
} else {  cat ("FST matrix file ", opt$genfile, "\n")
  genfile <- opt$genfile  
  }

if (is.null(opt$geofile)) {
  stop("WARNING: No geofile with geographic coordinates fopr each population specified with '-G' flag. Desired input file contains Latitude\tLongitude\tpopulation")
} else {  cat ("geofile ", opt$geofile, "\n")
  geofile <- opt$geofile  
  }

if (is.null(opt$outpath)) {
  stop("WARNING: No outpath specified with '-O' flag.")
} else {  cat ("outpath ", opt$outpath, "\n")
  outpath <- opt$outpath  
  }
  
  
library(vegan)
library(geosphere)

fst_data  <- read.table(genfile, header=T, row.names=1, sep="\t")
coords  <- read.table(geofile, header=T, sep="\t")
colnames(coords) <- c("Longitude", "Latitude", "population")
fst_matrix <- as.matrix(fst_data)

# Copy upper triangle values to the lower triangle to symmetrize the matrix
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]

# The diagonal should remain zeros, so we do not need to change it
diag(fst_matrix) <- 0

# Ensure that the populations are in the same order in both the FST matrix and geographic coordinates
coords <- coords[match(rownames(fst_matrix), coords$population), ]

# Calculate geographic distance matrix
geo_dist_matrix <- distm(coords[, c("Longitude", "Latitude")], fun = distHaversine)

rownames(geo_dist_matrix) <- coords$population
colnames(geo_dist_matrix) <- coords$population

mantel_test <- mantel(fst_matrix, geo_dist_matrix)
print(mantel_test)

# scatterplot

# Flatten the distance matrices into vectors
genetic_dist <- as.vector(fst_matrix)
geographic_dist <- as.vector(geo_dist_matrix)
# Filter out the self-comparison points (where both distances are 0)
valid_points <- geographic_dist != 0 & genetic_dist != 0

# Create scatterplot
fileOut=file.path(outpath,"scatterplot_FST_geographic_populations.jpeg")
jpeg(fileOut,width=16,height=15,units="cm",res=300, type="cairo")
par(cex.axis = 0.8, cex.lab = 1)
par(cex.main = 1.1)
par(mar = c(5, 6, 4, 4),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.7,0.8,0))
plot(geographic_dist[valid_points], genetic_dist[valid_points],
     xlab = "Geographic Distance",
     ylab = "Genetic Distance (FST)",
     main = "Genetic Distance vs Geographic Distance",
     pch = 19, col = "blue")

# Add a regression line
abline(lm(genetic_dist[valid_points] ~ geographic_dist[valid_points]), col = "red")

# Add Mantel statistic and significance to the plot
text(x = 43000, y = 0.1,
     labels = paste("r =", round(mantel_test$statistic, 4), "\np =", round(mantel_test$signif, 3)),
     cex = 0.9, col = "black")
dev.off()
