library(devtools)
install_github("jokergoo/ComplexHeatmap")
install.packages("viridis")  # Install
library("viridis") 
library(ComplexHeatmap)
library(plyr)
library(readr)
library(stringr)
library(data.table)
library(tidyverse)

format_df <- function(filenames, family_linker){
  df <- ldply(filenames, read.csv, header=TRUE, sep='\t') 
  view(df)
  df <- merge(df, family_linker, by = ".id")
  df <- df[c("module", "family", "module_subcategory", "pathwise_module_completeness")]
  df <- unique(df)
  return(df)
}

family_mean <- function(df){
  # Take the mean of the values due to different completeness in genomes so 
  # That they can be used in a cross comparison.
  keys <- colnames(df)[!grepl('pathwise_module_completeness',colnames(df))]
  mean_vals <- as.data.table(df)
  
  df <- mean_vals[,list(pathwise_module_completeness = mean(pathwise_module_completeness)),keys]
  return(df)
}

meta_reshape <- function(df){
  # Reshape the dataframe
  reshaped <- spread(df, key = family, value = pathwise_module_completeness)
  reshaped <- reshaped %>% mutate_all(~replace(., is.na(.), 0))
  rownames(reshaped) <- paste0(paste0(paste(reshaped$module, '('), reshaped$module_subcategory), ')')
  reshaped <- as.matrix(reshaped[, -1][, -1])
  return(reshaped)
}

# Graph a species level heatmap for Lachnospiraceae family
family_linker <- data.frame(read.csv("/Users/kananen/Desktop/ImHere/bac120_metadata_r214.tsv",
                              sep='\t', header = TRUE))
family_linker$accession <- gsub('RS_', '', family_linker$accession)
family_linker$accession <- gsub('GB_', '', family_linker$accession)
family_linker <- family_linker[c("accession", "gtdb_taxonomy")]
family_linker$family <- gsub('f__', '', str_split_fixed(family_linker$gtdb_taxonomy, ";",7)[,5])
names(family_linker)[1] <- ".id"
family_linker <- family_linker[c(".id", "family")]

anvio <- list.files(list.files("/Users/kananen/Desktop/ImHere/set_36/anvio/stray/metabolism/", full.names = TRUE), full.names = TRUE)
microbeannotator <- list.files(list.files("/Users/kananen/Desktop/ImHere/set_36/microbeannotator/metabolism/refined/", full.names = TRUE), full.names = TRUE)
kofamscan <- list.files(list.files("/Users/kananen/Desktop/ImHere/set_36/kofamscan/raw/metabolism/", full.names = TRUE), full.names = TRUE)

names(anvio) <- paste(gsub('.txt','', basename(dirname(anvio))))
names(microbeannotator) <- paste(gsub('.txt','', basename(dirname(microbeannotator))))
names(kofamscan) <- paste(gsub('.txt','', basename(dirname(kofamscan))))


df_anvio <- format_df(anvio, family_linker)
df_microbeannotator <- format_df(microbeannotator, family_linker)
df_kofamscan <- format_df(kofamscan, family_linker)
ldply(df_kofamscan, read.csv, header=FALSE, sep='\t') 











df_anvio <- family_mean(df_anvio)
df_microbeannotator <- family_mean(df_microbeannotator)
df_kofamscan <- family_mean(df_kofamscan)


mx_anvio <- meta_reshape(df_anvio)
mx_microbeannotator <- meta_reshape(df_microbeannotator)
mx_kofamscan <- meta_reshape(df_kofamscan)

anvio_subsample <- mx_anvio[sample(nrow(mx_anvio), 20), ]
microbeannotator_subsample <- mx_microbeannotator[rownames(mx_microbeannotator) %in% rownames(anvio_subsample), ]

missing_rows <- rownames(anvio_subsample)[!(rownames(anvio_subsample) %in% rownames(microbeannotator_subsample))]

# Create a dataframe with missing rows filled with 0s
missing_df <- data.frame(matrix(0, nrow = length(missing_rows), ncol = ncol(anvio_subsample)))
rownames(missing_df) <- missing_rows
colnames(missing_df) <- colnames(anvio_subsample)

# Combine with missing_df
microbeannotator_subsample <- rbind(microbeannotator_subsample, missing_df)
colnames <- colnames(microbeannotator_subsample)

anvio_subsample <- anvio_subsample[order(rownames(anvio_subsample)), order(colnames(anvio_subsample))]
microbeannotator_subsample <- microbeannotator_subsample[order(rownames(microbeannotator_subsample)), order(colnames(microbeannotator_subsample))]
colnames(anvio_subsample) <- colnames
colnames(microbeannotator_subsample) <- colnames

h1 <- Heatmap(anvio_subsample, name = "Module Completeness",
              column_title = "Anvi'o",
              row_names_side = "left",
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_max_width = unit(10, "cm"),
              col=viridis(100))
# h2 <- Heatmap(microbeannotator_subsample, column_title = "Kofamscan",
#               show_heatmap_legend = FALSE,
#               show_row_names = FALSE,
#               cluster_rows = FALSE,
#               cluster_columns = FALSE,
#               row_names_max_width = unit(10, "cm"),
#               col=viridis(100))
h3 <- Heatmap(microbeannotator_subsample, column_title = "MicrobeAnnotator",
              show_heatmap_legend = FALSE, 
              show_row_names = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_max_width = unit(2, "cm"),
              col=viridis(100))
h_list = h1 + h3

draw(h_list, row_title = "KEGG Metabolic Modules", column_title = "Mean Pathway Completeness Per Family", column_title_gp = gpar(fontsize = 16))

# Noticed drug related pathways are more complete with stray ko recovery
# Checking if this is the case or if I am making patterns where they are not.
drugs_anvio <- mx_anvio[grep("Pathogenicity", row.names(mx_anvio)),]
drugs_microbeannotator <- mx_microbeannotator[grep("Pathogenicity", row.names(mx_microbeannotator)),]

missing_rows <- rownames(drugs_anvio)[!(rownames(drugs_anvio) %in% rownames(drugs_microbeannotator))]
missing_df <- data.frame(matrix(0, nrow = length(missing_rows), ncol = ncol(drugs_anvio)))
rownames(missing_df) <- missing_rows
colnames(missing_df) <- colnames(drugs_anvio)
drugs_microbeannotator <- rbind(drugs_microbeannotator, missing_df)

missing_rows <- rownames(drugs_microbeannotator)[!(rownames(drugs_microbeannotator) %in% rownames(drugs_anvio))]
missing_df <- data.frame(matrix(0, nrow = length(missing_rows), ncol = ncol(drugs_microbeannotator)))
rownames(missing_df) <- missing_rows
colnames(missing_df) <- colnames(drugs_microbeannotator)
drugs_anvio <- rbind(drugs_anvio, missing_df)

drugs_anvio <- drugs_anvio[order(rownames(drugs_anvio)), order(colnames(drugs_anvio))]
drugs_microbeannotator <- drugs_microbeannotator[order(rownames(drugs_microbeannotator)), order(colnames(drugs_microbeannotator))]
colnames(drugs_anvio) <- colnames
colnames(drugs_microbeannotator) <- colnames

h4 <- Heatmap(drugs_anvio, name = "Module Completeness",
              column_title = "Anvi'o",
              row_names_side = "left",
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_max_width = unit(10, "cm"),
              col=viridis(100))
# h5 <- Heatmap(drugs_anvio, column_title = "Anvi'o (Stray KO)",
#               show_heatmap_legend = FALSE,
#               show_row_names = FALSE,
#               cluster_rows = FALSE,
#               cluster_columns = FALSE,
#               row_names_max_width = unit(10, "cm"),
#               col=viridis(100))
h6 <- Heatmap(drugs_microbeannotator, column_title = "MicrobeAnnotator",
              show_heatmap_legend = FALSE, 
              show_row_names = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_max_width = unit(2, "cm"),
              col=viridis(100))
h_list = h4 + h6

draw(h_list, row_title = "KEGG Metabolic Modules", column_title = "Mean Drug Metabolic Pathway Completeness Per Family", column_title_gp = gpar(fontsize = 16))










# Graph a species level heatmap for Lachnospiraceae family
linker <- data.frame(read.csv("/Users/kananen/Desktop/ImHere/bac120_metadata_r214.tsv",
                              sep='\t', header = FALSE))

lachno_anvio <- ldply(anvio, read.csv, header=TRUE, sep='\t') 
lachno_anvio <- merge(lachno_anvio, family_linker, by = ".id")
tmp <- lachno_anvio[grep("Cyanobacteriaceae", lachno_anvio[["family"]]),][['.id']]
unique(tmp)

# lachno_anvio_subsample <- mx_anvio[sample(nrow(lachno_anvio), 20), ]
# 
# lachno_microbeannotator <- ldply(microbeannotator, read.csv, header=TRUE, sep='\t') 
# lachno_microbeannotator <- merge(lachno_microbeannotator, family_linker, by = ".id")
# lachno_microbeannotator <- lachno_microbeannotator[grep("Lachnospiraceae", lachno_microbeannotator[["family"]]),]
# lachno_microbeannotator["genome_name"] <- gsub('_microbeannotator', '', lachno_microbeannotator[["genome_name"]])
# 
# 
# Heatmap(lachno_microbeannotator)
