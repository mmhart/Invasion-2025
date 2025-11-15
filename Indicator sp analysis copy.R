library(indicspecies)
library(permute)
library(tidyverse)
library(dplyr)

## ---- read metadata ----
meta <- read.table("InvasMeta.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE,
                   comment.char = "",
                   stringsAsFactors = FALSE)

## ---- read feature table ----
ft_raw <- read.table("exported-collapsed-table-L7-m2.txt",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE,
                     comment.char = "")

## ---- ensure samples x taxa ----
clean_ids <- function(x) {
  x <- trimws(x)
  x <- gsub("\u00A0", " ", x)  # non-breaking space
  x <- gsub("\\s+", " ", x)
  x
}
rownames(meta)  <- clean_ids(rownames(meta))
rownames(ft_raw) <- clean_ids(rownames(ft_raw))
colnames(ft_raw) <- clean_ids(colnames(ft_raw))

# If samples are columns (QIIME export), transpose to samples x taxa
match_rows <- sum(rownames(ft_raw) %in% rownames(meta))
match_cols <- sum(colnames(ft_raw) %in% rownames(meta))
FT <- if (match_cols > match_rows) t(as.matrix(ft_raw)) else as.matrix(ft_raw)

mode(FT) <- "numeric"


## ---- align samples ----
keep_ids <- intersect(rownames(FT), rownames(meta))
if (length(keep_ids) < 2L) stop("No overlapping sample IDs between FT and metadata.")
FT   <- FT[keep_ids, , drop = FALSE]
meta <- meta[keep_ids, , drop = FALSE]

## ---- build groups from CopiesPerGSoil (Low/High) ----
# make numeric safely (strip commas if present)
meta$CopiesPerGSoil <- as.numeric(gsub(",", "", as.character(meta$CopiesPerGSoil)))
if (all(is.na(meta$CopiesPerGSoil))) stop("CopiesPerGSoil is all NA after numeric conversion.")

#set your groups, here 25000 copies was cut off between “High” and “Low”
group <- cut(meta$CopiesPerGSoil,
             breaks = c(-Inf, 25000, Inf),
             labels = c("Low", "High"),
             include.lowest = TRUE)


# drop samples with NA group (if any)
has_group <- !is.na(group)
FT   <- FT[has_group, , drop = FALSE]
meta <- meta[has_group, , drop = FALSE]
group <- droplevels(group[has_group])

if (nlevels(group) < 2L || any(table(group) < 2L)) {
  print(table(group))
  stop("Indicator analysis needs at least 2 groups with ≥2 samples each.")
}

#look at levels to make sure ok
nlevels(group)
table(group)



## ---- run IndVal.g ----
set.seed(123)
ctrl <- how(nperm = 9999)
indval <- multipatt(FT, cluster = group, func = "IndVal.g", control = ctrl)

#view results, including non-significant ones:
summary(indval, indvalcomp = TRUE, alpha = 1)
  





#############HOW TO SUMMARIZE IN A TABLE FOR ALL TAXA REGARDLESS OF P
# --- 1) Base table from FT: all taxa, no filtering ---
pseudocount <- 1
mean_high <- colMeans(FT[group == "High", , drop = FALSE])
mean_low  <- colMeans(FT[group == "Low",  , drop = FALSE])

all_taxa <- tibble(
  Taxon    = colnames(FT),
  meanHigh = as.numeric(mean_high[colnames(FT)]),
  meanLow  = as.numeric(mean_low[colnames(FT)]),
  log2FC   = log2((mean_high[colnames(FT)] + pseudocount) /
                    (mean_low[colnames(FT)]  + pseudocount))
) %>%
  mutate(Direction = case_when(
    log2FC > 0 ~ "Higher_in_High",
    log2FC < 0 ~ "Higher_in_Low",
    TRUE       ~ "No_change"
  ))

# --- 2) (Optional) add IndVal stats if available, but DON'T filter by them ---
sign_tbl <- as.data.frame(indval$sign)
sign_tbl$Taxon <- rownames(sign_tbl)

comb_df <- as.data.frame(indval$comb)
comb_df$Taxon <- rownames(comb_df)

grp_names <- setdiff(colnames(comb_df), "Taxon")
comb_df$GroupAssoc <- apply(comb_df[, grp_names, drop = FALSE], 1, function(x)
  paste(grp_names[which(x == 1)], collapse = if (sum(x) > 1) "+" else ""))

# keep only needed columns from IndVal
iv_stats <- sign_tbl %>% select(Taxon, stat, p.value)
iv_groups <- comb_df %>% select(Taxon, GroupAssoc)

# FULL join keeps every taxon from the base table, even if IndVal didn’t flag it
res_all <- all_taxa %>%
  left_join(iv_stats,  by = "Taxon") %>%
  left_join(iv_groups, by = "Taxon") %>%
  arrange(desc(log2FC))

# --- 3) Inspect everything (no significance filter applied) ---
print(res_all, n = Inf)

write.table(res_all, "indicator_taxa_results2.txt",
            sep = "\t",          # tab-separated
            quote = FALSE,       # don’t add quotes around text
            row.names = FALSE)   # omit row numbers


####plot FIGURE

library(ggplot2)
library(dplyr)
library(stringr)

# Create a shorter species label from taxonomy string
plot_df <- all_taxa %>%
  mutate(
    Species = ifelse(str_detect(Taxon, "s__"),
                     str_replace(Taxon, ".*s__", ""),  # take text after s__
                     Taxon),
    Species = str_replace_all(Species, "_", " "),     # prettier names
    Species = str_trim(Species)
  )


#order by unique spp

plot_df <- plot_df %>%
  mutate(Species = factor(Species, 
                          levels = unique(Species[order(log2FC)])))


#plot by species
ggplot(plot_df, aes(x = log2FC, y = Species, fill = Direction)) +
  geom_col(color = "black", width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    name = "Inoculant",  # legend title
    values = c(
      "Higher_in_High" = "#FC8D62",
      "Higher_in_Low"  = "#66C2A5"
    ),
    labels = c(
      "Higher_in_High" = "Inoculant present",
      "Higher_in_Low"  = "Inoculant absent"
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = expression(Log[2]*"(Fold Change: Inoculant present / Inoculant absent)"),
    y = "AMF Taxa"
  ) +
  scale_y_discrete(labels = c(
    "uncultured Glomeromycotina" = "uncultured Glomeraceae",
    "Diversispora spurca" = expression(italic("Diversispora spurca")),
    "Claroideoglomus etunicatum" = expression(italic("Claroideoglomus etunicatum")),
    "uncultured fungus" = "uncultured Glomeraceae",
    "Archaeospora schenckii" = expression(italic("Archaeospora schenckii")),
    "Dentiscutata reticulata" = expression(italic("Dentiscutata reticulata")),
    "Racocetra castanea" = expression(italic("Racocetra castanea")),
    "Acaulospora lacunosa" = expression(italic("Acaulospora lacunosa")),
    "Diversispora trimurales" = expression(italic("Diversispora trimurales")),
    "Paraglomus occultum" = expression(italic("Paraglomus occultum")),
    "uncultured Acaulospora" = "uncultured Acaulospora sp.",
    "uncultured Eimeriidae" = "uncultured Eimeriidae",
    "Acaulosporaceae sp." = expression(italic("Acaulosporaceae sp.")),
    "Archaeospora trappei" = expression(italic("Archaeospora trappei")),
    "Glomeromycotina sp." = "uncultured Glomeromycotina",
    "Glomus mycorrhizal" = "uncultured Glomus",
    "Cetraspora nodosa" = expression(italic("Cetraspora nodosa")),
    "Acaulospora laevis" = expression(italic("Acaulospora laevis")),
    "Paraglomus brasilianum" = expression(italic("Paraglomus brasilianum")),
    "Pacispora franciscana" = expression(italic("Pacispora franciscana"))
  )) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank()
  )
