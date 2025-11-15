library(vegan)
library(dplyr)
library(permute)

#read metadata, with site as numeric
meta <- read.table("InvasMetaMantel.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

# Convert the column to numeric (in case it's a character or factor)
meta$CopiesPerGSoil <- as.numeric(meta$CopiesPerGSoil)

# Create Euclidean distance from that numeric variable
copies_dist <- dist(meta$CopiesPerGSoil, method = "euclidean")

# Verify its a distance matrix 
inherits(copies_dist, "dist")   # returns TRUE if it's a distance matrix

#import bc distance matrix
bray <- read.table("bray-distance-matrix.tsv", 
                   header = TRUE, 
                   row.names = 1, 
                   sep = "\t", 
                   check.names = FALSE)


#control for site with partial mantel, make a "site" distance matrix
meta$site <- as.numeric(meta$site)
site_dist <- dist(meta$site, method = "euclidean")
inherits(site_dist, "dist") 


#run partial mantel
mantel.partial(bray, copies_dist, site_dist, permutations = 9999)
