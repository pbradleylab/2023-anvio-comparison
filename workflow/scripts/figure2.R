library(tidyverse)
library(plyr)
library(data.table)
library(ComplexHeatmap)
library(viridis) 
library(patchwork)


read_family <- function(path){
  family_linker <- read.csv(path, sep='\t', header = TRUE)
  # Clean accession IDs
  family_linker$accession <- gsub('(RS_|GB_)', '', family_linker$accession)
  # Extract families and clean taxonomy
  family_linker$family <- gsub('f__', '', str_split_fixed(family_linker$gtdb_taxonomy, ";", 7)[, 5])
  colnames(family_linker)[1] <- ".id"
  # Clean taxonomy labels and Update specific accession IDs
  family_linker$family <- gsub("Microcystaceae_[AB]", "Microcystaceae", family_linker$family)
  family_linker$.id <- gsub("GCA_018500005\\.1", "GCA_018500005.2", family_linker$.id)
  return(family_linker)
}

format_df <- function(filenames, family_linker){
  names(filenames) <- paste(gsub('.txt','', basename(dirname(filenames))))
  df <- ldply(filenames, read.csv, header=TRUE, sep='\t') 
  df <- merge(df, family_linker, by = ".id")
  df <- df[c(".id","module", "family", "module_subcategory", "pathwise_module_completeness")]
  return(df)
}

family_mean <- function(df){
  # Take the mean of the values due to different completeness in genomes so 
  # That they can be used in a cross comparison.
  return(df[, .(mean_completeness = mean(pathwise_module_completeness)), by = .(module, family, module_subcategory)])
}

meta_reshape <- function(df){
  # Reshape the dataframe and replace NAs with 0
  reshaped <- spread(df, key = family, value = mean_completeness, fill = 0)
  rownames(reshaped) <- paste0(reshaped$module, ' (', reshaped$module_subcategory, ')')
  # Convert to matrix and remove extra columns
  return(as.matrix(reshaped[, -c(1, 2), drop = FALSE]))
}

add_missing_rows <- function(to_add_df, to_check_df){
  missing_rows <- setdiff(rownames(to_check_df), rownames(to_add_df))
  missing_df <- as.data.frame(matrix(0, nrow = length(missing_rows), ncol = ncol(to_check_df)))
  rownames(missing_df) <- missing_rows
  colnames(missing_df) <- colnames(to_check_df)
  return(rbind(to_add_df, missing_df))
}

rowwise_manhattan_distance <- function(matrix1, matrix2) {
  # Ensure row names match
  if (!identical(rownames(matrix1), rownames(matrix2))) {
    stop("Row names of the two matrices do not match.")
  }
  # Initialize vector to store distances
  distances <- numeric(nrow(matrix1))
  # Calculate Manhattan distance row by row
  for (i in 1:nrow(matrix1)) {
    distances[i] <- sum(abs(matrix1[i, ] - matrix2[i, ]))
  }
  return(distances)
}

extract_family_data <- function(family_name, anvio, kofamscan, microbeannotator) {
  # Extract data for the specified bacterial family from each dataframe
  data_anvio <- anvio[, family_name, drop = FALSE]  # Retain the specified columns, keeping all rows
  data_kofamscan <- kofamscan[, family_name, drop = FALSE]
  data_microbeannotator <- microbeannotator[, family_name, drop = FALSE]
  # Combine data for the bacterial family from all tools into a single dataframe
  merged_data <- data.frame(
    Anvio = data_anvio[, family_name],  # Access the specified column within each dataframe
    Kofamscan = data_kofamscan[, family_name],
    MicrobeAnnotator = data_microbeannotator[, family_name]
  )
  # Set the row names of the merged dataframe to match the original dataframe's row names
  rownames(merged_data) <- rownames(anvio)
  return(merged_data)
}



# Read in data -----------------------------------------------------------------
family_linker <- read_family("/Users/kananen/Desktop/ImHere/bac120_metadata_r214.tsv")
anvio <- list.files(list.files("/Users/kananen/Desktop/ImHere/unfiltered/anvio/metabolism/stray/", full.names = TRUE), full.names = TRUE)
microbeannotator <- list.files(list.files("/Users/kananen/Desktop/ImHere/unfiltered/microbeAnnotator/metabolism/refined/", full.names = TRUE), full.names = TRUE)
kofamscan <- list.files(list.files("/Users/kananen/Desktop/ImHere/unfiltered/kofamscan/metabolism/default/", full.names = TRUE), full.names = TRUE)
# Format and calculate the mean for the pathway completness for each family
anvio <- format_df(anvio, family_linker)
anvio <- subset(anvio, pathwise_module_completeness >= .80)
microbeannotator <- format_df(microbeannotator, family_linker)
microbeannotator <- subset(microbeannotator, pathwise_module_completeness >= .80)
kofamscan <- format_df(kofamscan, family_linker)
kofamscan <- subset(kofamscan, pathwise_module_completeness >= .80)

completness <- aggregate(pathwise_module_completeness ~.id, data = anvio, FUN = function(x) sum(x))
colnames(completness) <- c(".id", "anvio_completness")
completness <- merge(completness, aggregate(pathwise_module_completeness ~.id, data = microbeannotator, FUN = function(x) sum(x)), all = TRUE)
colnames(completness) <- c(".id", "anvio_completness", "microbeannotator_completness")
completness <- merge(completness, aggregate(pathwise_module_completeness ~.id, data = kofamscan, FUN = function(x) sum(x)), all = TRUE)
colnames(completness) <- c(".id", "anvio_completness", "microbeannotator_completness", "kofamscan_completness")
completness <- merge(completness, aggregate(pathwise_module_completeness ~.id, data = anvio, FUN = function(x) length(x)), all = TRUE)
colnames(completness) <- c(".id", "anvio_completness", "microbeannotator_completness", "kofamscan_completness", "anvio_completness_total")
completness <- merge(completness, aggregate(pathwise_module_completeness ~.id, data = microbeannotator, FUN = function(x) length(x)), all = TRUE)
colnames(completness) <- c(".id", "anvio_completness", "microbeannotator_completness", "kofamscan_completness", "anvio_completness_total", "microbeannotator_completness_total")
completness <- merge(completness, aggregate(pathwise_module_completeness ~.id, data = kofamscan, FUN = function(x) length(x)), all = TRUE)
colnames(completness) <- c(".id", "anvio_completness", "microbeannotator_completness", "kofamscan_completness", "anvio_completness_total", "microbeannotator_completness_total", "kofamscan_completness_total")
completness <- replace(completness, is.na(completness), 0) 


sum(completness$anvio_completness_total)/396
sum(completness$microbeannotator_completness_total)/396
sum(completness$kofamscan_completness_total)/396

microbeannotator_mean_80_comp <- sum(aggregate(pathwise_module_completeness ~.id, data = microbeannotator, FUN = function(x) sum(x)/length(x))$pathwise_module_completeness)/length(unique(microbeannotator$.id))
kofamscan_mean_80_comp <- sum(aggregate(pathwise_module_completeness ~.id, data = kofamscan, FUN = function(x) sum(x)/length(x))$pathwise_module_completeness)/length(unique(kofamscan$.id))

anvio <- meta_reshape(family_mean(setDT(format_df(anvio, family_linker))))
microbeannotator <- meta_reshape(family_mean(setDT(format_df(microbeannotator, family_linker))))
kofamscan <- meta_reshape(family_mean(setDT(format_df(kofamscan, family_linker))))

# Here to make sure we are being fair, we add in any missing modules found by the other
# metabolism pathway estimates for each method.
anvio <- as.matrix(add_missing_rows(add_missing_rows(anvio, microbeannotator), kofamscan))
microbeannotator <- as.matrix(add_missing_rows(add_missing_rows(microbeannotator, kofamscan), anvio))
kofamscan <- as.matrix(add_missing_rows(add_missing_rows(kofamscan, microbeannotator), anvio))
# Finally order them 
anvio <- anvio[order(rownames(anvio)), ]
microbeannotator <- microbeannotator[order(rownames(microbeannotator)), ]
kofamscan <- kofamscan[order(rownames(kofamscan)), ]

# Calculate the most dissimilar rows between the the methods
distances <- as.data.frame(rowwise_manhattan_distance(anvio, kofamscan))
distances$module <- rownames(anvio) 
top_20 <- head(distances[order(-distances$`rowwise_manhattan_distance(anvio, kofamscan)`), ], 20)
top_20_anvio <- anvio[rownames(anvio) %in% top_20$module, , drop = FALSE]
top_20_kofamscan <- kofamscan[rownames(kofamscan) %in% top_20$module, , drop = FALSE]
top_20_microbeannotator <- microbeannotator[rownames(microbeannotator) %in% top_20$module, , drop = FALSE]

an10 <- Heatmap(top_20_anvio, name = "Module Completeness",
                row_names_side = "left",
                show_row_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                column_title_rot = 45,
                row_names_max_width = unit(10, "cm"),
                col=viridis(100))
ko10 <- Heatmap(top_20_kofamscan, name = "Module Completeness",
                row_names_side = "left",
                show_row_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                column_title_rot = 45,
                row_names_max_width = unit(10, "cm"),
                col=viridis(100))
ma10 <- Heatmap(top_20_microbeannotator, name = "Module Completeness",
                row_names_side = "left",
                show_row_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                column_title_rot = 45,
                row_names_max_width = unit(10, "cm"),
                col=viridis(100))
an10 + ko10 + ma10

# Create separate data frames for each bacterial family
bacteroidaceae <- extract_family_data("Bacteroidaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
cyanobiaceae <- extract_family_data("Cyanobiaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
enterobacteriaceae <- extract_family_data("Enterobacteriaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
lachnospiraceae <- extract_family_data("Lachnospiraceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
nanosynbacteraceae <- extract_family_data("Nanosynbacteraceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
pelagibacteraceae <- extract_family_data("Pelagibacteraceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
rhizobiaceae <- extract_family_data("Rhizobiaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
rhodobacteraceae <- extract_family_data("Rhodobacteraceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
streptococcaceae <- extract_family_data("Streptococcaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
streptomycetaceae <- extract_family_data("Streptomycetaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)
synechococcaceae <- extract_family_data("Microcystaceae", top_20_anvio, top_20_kofamscan, top_20_microbeannotator)

h7 <- Heatmap(bacteroidaceae, name = "Module Completeness",
              column_title = "Bacteroidaceae",
              row_names_side = "left",
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_title_rot = 45,
              row_names_max_width = unit(10, "cm"),
              col=viridis(100))
h8 <- Heatmap(cyanobiaceae, column_title = "Cyanobiaceae",
              show_heatmap_legend = FALSE,
              show_row_names = FALSE,
              cluster_rows = FALSE,
              column_title_rot = 45,
              cluster_columns = FALSE,
              row_names_max_width = unit(2, "cm"),
              col=viridis(100))
h9 <- Heatmap(enterobacteriaceae, column_title = "Enterobacteriaceae",
              show_heatmap_legend = FALSE,
              show_row_names = FALSE,
              cluster_rows = FALSE,
              column_title_rot = 45,
              cluster_columns = FALSE,
              row_names_max_width = unit(2, "cm"),
              col=viridis(100))
h10 <- Heatmap(lachnospiraceae, column_title = "Lachnospiraceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               column_title_rot = 45,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h11 <- Heatmap(nanosynbacteraceae, column_title = "Nanosynbacteraceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 45,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h12 <- Heatmap(pelagibacteraceae, column_title = "Pelagibacteraceae",
               show_heatmap_legend = FALSE,
               column_title_rot = 45,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h13 <- Heatmap(rhizobiaceae, column_title = "Rhizobiaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               column_title_rot = 45,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h14 <- Heatmap(rhodobacteraceae, column_title = "Rhodobacteraceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 45,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h15 <- Heatmap(streptococcaceae, column_title = "Streptococcaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               column_title_rot = 45,
               col=viridis(100))
h16 <- Heatmap(streptomycetaceae, column_title = "Streptomycetaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 45,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h17 <- Heatmap(synechococcaceae, column_title = "Microcystaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               column_title_rot = 45,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h_list = h7 + h8 + h9 + h10 + h11 + h12 + h13 + h14 + h15 + h16 + h17

drugs_anvio <- anvio[grep("Drug", row.names(anvio)),]
drugs_microbeannotator <- microbeannotator[grep("Drug", row.names(microbeannotator)),]
drugs_kofamscan <- kofamscan[grep("Drug", row.names(kofamscan)),]

# Create separate data frames for each bacterial family
bacteroidaceae <- extract_family_data("Bacteroidaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
cyanobiaceae <- extract_family_data("Cyanobiaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
enterobacteriaceae <- extract_family_data("Enterobacteriaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
lachnospiraceae <- extract_family_data("Lachnospiraceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
nanosynbacteraceae <- extract_family_data("Nanosynbacteraceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
pelagibacteraceae <- extract_family_data("Pelagibacteraceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
rhizobiaceae <- extract_family_data("Rhizobiaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
rhodobacteraceae <- extract_family_data("Rhodobacteraceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
streptococcaceae <- extract_family_data("Streptococcaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
streptomycetaceae <- extract_family_data("Streptomycetaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)
synechococcaceae <- extract_family_data("Microcystaceae", drugs_anvio, drugs_kofamscan, drugs_microbeannotator)

h7 <- Heatmap(bacteroidaceae, name = "Module Completeness",
              column_title = "Bacteroidaceae",
              row_names_side = "left",
              show_row_dend = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_title_rot = 45,
              row_names_max_width = unit(10, "cm"),
              col=viridis(100))
h8 <- Heatmap(cyanobiaceae, column_title = "Cyanobiaceae",
              show_heatmap_legend = FALSE,
              show_row_names = FALSE,
              cluster_rows = FALSE,
              column_title_rot = 45,
              cluster_columns = FALSE,
              row_names_max_width = unit(2, "cm"),
              col=viridis(100))
h9 <- Heatmap(enterobacteriaceae, column_title = "Enterobacteriaceae",
              show_heatmap_legend = FALSE,
              show_row_names = FALSE,
              cluster_rows = FALSE,
              column_title_rot = 45,
              cluster_columns = FALSE,
              row_names_max_width = unit(2, "cm"),
              col=viridis(100))
h10 <- Heatmap(lachnospiraceae, column_title = "Lachnospiraceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               column_title_rot = 45,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h11 <- Heatmap(nanosynbacteraceae, column_title = "Nanosynbacteraceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 45,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h12 <- Heatmap(pelagibacteraceae, column_title = "Pelagibacteraceae",
               show_heatmap_legend = FALSE,
               column_title_rot = 45,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h13 <- Heatmap(rhizobiaceae, column_title = "Rhizobiaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               column_title_rot = 45,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h14 <- Heatmap(rhodobacteraceae, column_title = "Rhodobacteraceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 45,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h15 <- Heatmap(streptococcaceae, column_title = "Streptococcaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               column_title_rot = 45,
               col=viridis(100))
h16 <- Heatmap(streptomycetaceae, column_title = "Streptomycetaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 45,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
h17 <- Heatmap(synechococcaceae, column_title = "Microcystaceae",
               show_heatmap_legend = FALSE,
               show_row_names = FALSE,
               cluster_rows = FALSE,
               column_title_rot = 45,
               cluster_columns = FALSE,
               row_names_max_width = unit(2, "cm"),
               col=viridis(100))
d_list = h7 + h8 + h9 + h10 + h11 + h12 + h13 + h14 + h15 + h16 + h17

pdf(file="f2ab.pdf", width=15, height=20)
draw(h_list, row_title = "KEGG Metabolic Modules", column_title = "Top 20 Disimilar Mean Metabolic Method Pathways Completeness Per Family", column_title_gp = gpar(fontsize = 16))
draw(d_list, row_title = "KEGG Metabolic Modules", column_title = "Mean Drug Metabolic Method Pathways Completeness Per Family", column_title_gp = gpar(fontsize = 16))
dev.off()

for (i in unique(modules$category)) {
  filename <- paste0("heatmap_", gsub(" ", "_", str_to_title(i)), ".png")
  png(file = filename, width = 11.0, height = 9.25, units = "in", res = 1200)
  
  module_anvio <- mx_anvio[grep(i, rownames(mx_anvio)),]
  module_microbeannotator <- mx_microbeannotator[grep(i, rownames(mx_microbeannotator)),]
  module_kofamscan <- mx_kofamscan[grep(i, rownames(mx_kofamscan)),]
  
  if (!length(rownames(module_anvio)) <= 1) {
    module_anvio <- module_anvio[order(rownames(module_anvio)), order(colnames(module_anvio))]
    module_kofamscan <- module_kofamscan[order(rownames(module_kofamscan)), order(colnames(module_kofamscan))]
    module_microbeannotator <- module_microbeannotator[order(rownames(module_microbeannotator)), order(colnames(module_microbeannotator))]
    rownames(module_anvio) <- gsub(i, '', rownames(module_anvio))
    rownames(module_anvio) <- gsub("\\(|\\)", "", rownames(module_anvio))
    rownames(module_microbeannotator) <- gsub(i, '', rownames(module_microbeannotator))
    rownames(module_microbeannotator) <- gsub("\\(|\\)", "", rownames(module_microbeannotator))
    rownames(module_kofamscan) <- gsub(i, '', rownames(module_kofamscan))
    rownames(module_kofamscan) <- gsub("\\(|\\)", "", rownames(module_kofamscan))
  }
  common_min <- min(c(as.matrix(module_anvio), as.matrix(module_microbeannotator), as.matrix(module_kofamscan)))
  common_max <- max(c(as.matrix(module_anvio), as.matrix(module_microbeannotator), as.matrix(module_kofamscan)))
  
  colors <- viridis(100)
  col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 100), colors)
  cell_fun <- function(j, k, x, y, width, height, fill, ...) {
    if(k != 1 && j != 1){
      if (module_anvio[k, j] > 0 && module_kofamscan[k, j] == 0 && module_microbeannotator[k, j] == 0) {
        grid.text("KM", x, y, gp=gpar(fontsize = 10, col = "white"))
      } else if((module_anvio[k, j] > 0 && module_kofamscan[k, j] == 0 && module_microbeannotator[k, j] != 0)) {
        grid.text("K", x, y, gp=gpar(fontsize = 10, col = "white"))
      } else if((module_anvio[k, j] > 0 && module_kofamscan[k, j] != 0 && module_microbeannotator[k, j] == 0)) {
        grid.text("M", x, y, gp=gpar(fontsize = 10, col = "white"))
      }
    }
  }
  
  h1 <- Heatmap(as.matrix(module_anvio), name = "Module Completeness",
                column_title = "anvi'o",
                row_names_side = "left",
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                row_names_max_width = unit(10, "cm"),
                col = col_fun,
                cell_fun = cell_fun)
  h2 <- Heatmap(as.matrix(module_kofamscan), column_title = "Kofamscan",
                show_heatmap_legend = FALSE,
                show_row_names = FALSE,
                show_column_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                row_names_max_width = unit(10, "cm"),
                col = col_fun)
  h3 <- Heatmap(as.matrix(module_microbeannotator), column_title = "MicrobeAnnotator",
                show_heatmap_legend = FALSE, 
                show_row_names = FALSE,
                show_column_dend = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                row_names_max_width = unit(2, "cm"),
                col = col_fun)
  h_list = h1 + h2 + h3
  draw(h_list, row_title = "", column_title = paste(str_to_title(i), "Mean Pathway Completeness"), column_title_gp = gpar(fontsize = 16))
  dev.off()
}
