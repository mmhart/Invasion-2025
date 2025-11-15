library(vegan)
library(dplyr)
library(permute)

getwd()

#import dm into R from Qiime as .tsv
bray <- read.table("bray-distance-matrix.tsv", 
                    header = TRUE, 
                    row.names = 1, 
                    sep = "\t", 
                    check.names = FALSE)

                   
jaccard <- read.table("jaccard-distance-matrix.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)   

unweighted<- read.table("unweighted-distance-matrix.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)

meta <- read.table("InvasMeta.tsv",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

all(rownames(meta) %in% rownames(bray))
all(rownames(bray) %in% rownames(meta))





# Option 1 one test by terms
op <- options(contrasts = c("contr.sum","contr.poly"))

set.seed(42)
adonis2(bray ~ site * Inoculant, data = meta, permutations = 9999, by = "terms")

# restrict permuations to site to remove effect of site
set.seed(42)
permanova_fixed <- adonis2(
  +   bray ~ Inoculant,
  +   data = meta,
  +   permutations = ctrl, 
  +   by = "margin"
  + )

# Define restricted permutation scheme: permutations occur within levels of 'site'
ctrl <- how(nperm = 9999)
setBlocks(ctrl) <- meta$site

set.seed(42)

# Test the effect of Inoculant while restricting permutations within site
permanova_fixed <- adonis2(
  bray ~ Inoculant,
  data        = meta,
  permutations = ctrl,
  by          = "margin"
)

permanova_fixed


