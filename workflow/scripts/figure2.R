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

manhattan_distance <- function(mat) {
  # Function to compute Manhattan distance per row
  apply(mat, 1, function(row) {
    sum(abs(row[1] - row[2]))
  })
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

# Calculates the differences between two dataframes between one column
get_aggregation <- function(colname1, colname2, df1, df2, name, output){
  # Creates a frequency table for the id column in the data data frame. Result 
  # is a table showing how many times a unique id value appears in the column.
  count1 <- table(df1$.id); count2 <- table(df2$.id)
  # Generate a dataframe from the count data.
  df <- data.frame(.id = unique(c(names(count1), names(count2))),
                    colname1 = count1[match(unique(c(names(count1), names(count2))), names(count1))],
                    colname2 = count2[match(unique(c(names(count1), names(count2))), names(count2))])
  df[is.na(df)] <- 0
  
  # Calculate the increase and add to output dataframe.
  if(nrow(output) == nrow(df)){
    # If the number of rows is the same for the output previously run then
    # we can skip the addition of the .id column.
    #
    # Take the percent increase
    output[[name]] <- ((df$colname1.Freq - df$colname2.Freq) / df$colname2.Freq)
    output <- merge(output, df[c("colname1.Freq", ".id")])
    output <- merge(output, df[c("colname2.Freq", ".id")])
  } else {
    # Add in the .id column to the frame first.
    df[[name]] <-((df$colname1.Freq - df$colname2.Freq) / df$colname2.Freq)
    output <- merge(output, df[c(".id", name)], all=TRUE)
    output <- merge(output, df[c("colname1.Freq", ".id")], all=TRUE)
    output <- merge(output, df[c("colname2.Freq", ".id")], all=TRUE)
  }
  names(output)[which(names(output) == "colname1")] <- colname1
  names(output)[which(names(output) == "colname2")] <- colname2
  return(output)
}


# Read in data -----------------------------------------------------------------
family_linker <- read_family("/Users/user/folder/ImHere/bac120_metadata_r214.tsv")
anvio <- format_df(list.files(list.files("/Users/user/folder/unfiltered/anvio/metabolism/default/", full.names = TRUE), full.names = TRUE), family_linker)
microbeannotator <- format_df(list.files(list.files("/Users/user/folder/unfiltered/microbeAnnotator/metabolism/default/", full.names = TRUE), full.names = TRUE), family_linker)
kofamscan <- format_df(list.files(list.files("/Users/user/folder/unfiltered/kofamscan/metabolism/default/", full.names = TRUE), full.names = TRUE), family_linker)

# Take a completeness above .80% for comparisons of more complete pathways
anvio_80 <- subset(anvio, pathwise_module_completeness >= .80)
microbeannotator_80 <- subset(microbeannotator, pathwise_module_completeness >= .80)
kofamscan_80 <- subset(kofamscan, pathwise_module_completeness >= .80)

# Make a combined dataframe for all of the modules for ease of use
all_modules <- merge(anvio, microbeannotator, all=TRUE, by=c(".id","module","family", "module_subcategory"))
names(all_modules) <- c(".id","module","family", "module_subcategory", "an_pwc", "ma_pwc")
all_modules <- merge(all_modules, kofamscan, all=TRUE, by=c(".id","module","family", "module_subcategory"))
names(all_modules) <- c(".id","module","family", "module_subcategory", "an_pwc", "ma_pwc", "ko_pwc")
# Set missing to be 0% complete
all_modules[is.na(all_modules)] <- 0.0

write_tsv(all_modules, file = "all_completeness_all_tools.csv")

# Generate Supplamental Figure 1
completeness <- read_tsv("all_completeness_all_tools.csv")
median_cmp_per_family <- completeness %>% 
  group_by(family, module) %>% 
  dplyr::summarize(anvio = median(an_pwc), ma = median(ma_pwc), kofamscan = median(ko_pwc))
median_cmp_per_family_long <- median_cmp_per_family %>%
  pivot_longer(cols=c(anvio, ma), names_to="cmp_method", values_to="cmp_completeness") %>%
  dplyr::rename(kfs_completeness = kofamscan)
sf1 <- median_cmp_per_family_long %>%
  ggplot(aes(x=kfs_completeness, y=cmp_completeness, color=cmp_method, label=module)) + 
  geom_point(alpha=0.8) + 
  geom_abline(slope=1,intercept=0,lty=2,col="#AAAAAA") + facet_wrap(~ family) + 
  theme_minimal() + scale_color_manual(values=c(anvio="#2A9D8F", ma="#6DB9E4"), 
                                       labels = c("ma" = "MicrobeAnnotator", "anvio" = "anvi'o")) +
  labs(title = "Module Completeness of Methods over Kofamscan by Family", 
       x = "Pathwise-Completeness (other Kofamscan)", 
       y = "Pathwise-Completeness (other method)",
       color = "Method")

# Calculate difference in modules found for anvio vs kofamscan and anvio vs
# microbeannotator. NOTE: This is NOT the percent increase. For this we can just
# take the counts with a table and then take the difference.
#
# Anvio vs Kofamscan
modules <- unique(anvio[c(".id")])
modules <- get_aggregation("modules","modules", anvio_80, kofamscan_80, "anvio_increase_kofamscan", modules)
modules$anvio_increase_kofamscan <- ((modules$colname1.Freq+1)-(modules$colname2.Freq+1))/(modules$colname2.Freq+1)
print(paste("Number of unique modules change percent in of anvio over kofamscan:", sum(modules$anvio_increase_kofamscan)/nrow(modules)))
# Anvio vs microbeannotator
modules <- unique(anvio[c(".id")])
modules <- get_aggregation("modules","modules", anvio_80, microbeannotator_80, "anvio_increase_microbeannotator", modules)
modules$anvio_increase_microbeannotator <- ((modules$colname1.Freq+1)-(modules$colname2.Freq+1))/(modules$colname2.Freq+1)
print(paste("Number of unique modules change percent in of anvio over microbeannotator:", sum(modules$anvio_increase_microbeannotator)/nrow(modules)))

# Calculate 20 most dissimilar rows accross the tools
#
# Here to make sure we are being fair, we add in any missing modules found by the other
# metabolism pathway estimates for each method prior to averaging by family. Then
# we shape average them by family.
anvio_fam <- anvio %>%
  group_by(family, module, module_subcategory) %>%
  dplyr::summarize(mean_completeness = mean(pathwise_module_completeness))
microbeannotator_fam <- kofamscan %>%
  group_by(family, module, module_subcategory) %>%
  dplyr::summarize(mean_completeness = mean(pathwise_module_completeness))
kofamscan_fam <- microbeannotator %>%
  group_by(family, module, module_subcategory) %>%
  dplyr::summarize(mean_completeness = mean(pathwise_module_completeness))
# Reshape them to have the format we want for graphing later on and take the mean
anvio_fam <- meta_reshape(family_mean(setDT(anvio)))
microbeannotator_fam <- meta_reshape(family_mean(setDT(microbeannotator)))
kofamscan_fam <- meta_reshape(family_mean(setDT(kofamscan, family_linker)))
# Here to make sure we are being fair, we add in any missing modules found by the other
# metabolism pathway estimates for each method.
anvio_mx <- as.matrix(add_missing_rows(add_missing_rows(anvio_fam, microbeannotator_fam), kofamscan_fam))
microbeannotator_mx <- as.matrix(add_missing_rows(add_missing_rows(microbeannotator_fam, kofamscan_fam), anvio_fam))
kofamscan_mx <- as.matrix(add_missing_rows(add_missing_rows(kofamscan_fam, microbeannotator_fam), anvio_fam))
# Finally order them 
anvio_mx <-  anvio_mx[order(rownames(anvio_mx)), ]
microbeannotator_mx <-  microbeannotator_mx[order(rownames(microbeannotator_mx)), ]
kofamscan_mx <-  kofamscan_mx[order(rownames(kofamscan_mx)), ]
# However we also don't want to compare this across families as then we are 
# not just looking at the tool completeness but are looking at family differences
# which could be biological.
an <- list(); ma <- list()
for (i in 1:ncol(anvio_mx)) {
  familyAn <- cbind(anvio_mx[, i], kofamscan_mx[, i])
  # Calculate Manhattan distance between rows
  an_manhattan <- cbind(familyAn, manhattan_distance(familyAn))
  an_10 <- head(an_manhattan[order(-an_manhattan[, 3]), ], 10)
  top_10_an <- anvio_mx[rownames(anvio_mx) %in% rownames(an_10), , drop = FALSE][, i]
  top_10_an <- as.data.frame(top_10_an)
  colnames(top_10_an) <- colnames(anvio_mx)[i]
  tmp <- top_10_an %>% rownames_to_column(var = "modules") %>% pivot_longer(cols = c(colnames(anvio_mx)[i]), names_to = "family", values_to = "completeness")
  tmp <- tmp %>%
    separate(modules, into = c("id", "description"), sep = " \\(")
  tmp$description <- gsub("\\)", "", tmp$description)
  an[[i]] <- tmp
  
  familyMa <- cbind(microbeannotator_mx[, i], kofamscan_mx[, i])
  # Calculate Manhattan distance between rows
  ma_manhattan <- cbind(familyMa, manhattan_distance(familyMa))
  ma_10 <- head(ma_manhattan[order(-ma_manhattan[, 3]), ], 10)
  top_10_ma <- microbeannotator_mx[rownames(microbeannotator_mx) %in% rownames(ma_10), , drop = FALSE][, i]
  top_10_ma <- as.data.frame(top_10_ma)
  colnames(top_10_ma) <- colnames(microbeannotator_mx)[i]
  tmp <- top_10_ma %>% rownames_to_column(var = "modules") %>% pivot_longer(cols = c(colnames(microbeannotator_mx)[i]), names_to = "family", values_to = "completeness")
  tmp <- tmp %>%
    separate(modules, into = c("id", "description"), sep = " \\(")
  tmp$description <- gsub("\\)", "", tmp$description)
  ma[[i]] <- tmp 
}
an_mx <- bind_rows(an)
an_mx$tool <- "anvio"
ma_mx <- bind_rows(ma)
ma_mx$tool <- "microbeannotator"

# Add in the average completeness socres per family
# Anvio
anvio_ave <- as.data.frame(anvio_mx) %>% 
  rownames_to_column(var = "modules") %>% 
  pivot_longer(cols = c(colnames(as.data.frame(anvio_mx))), names_to = "family", values_to = "completeness") %>% 
  separate(modules, into = c("id", "description"), sep = " \\(")
anvio_ave$description <- gsub("\\)", "", anvio_ave$description)
# MicrobeAnnotator
microbeannotator_ave <- as.data.frame(microbeannotator_mx) %>% 
  rownames_to_column(var = "modules") %>% 
  pivot_longer(cols = c(colnames(as.data.frame(microbeannotator_mx))), names_to = "family", values_to = "completeness") %>% 
  separate(modules, into = c("id", "description"), sep = " \\(")
microbeannotator_ave$description <- gsub("\\)", "", microbeannotator_ave$description)
#Kofamscan
kofamscan_ave <- as.data.frame(kofamscan_mx) %>% 
  rownames_to_column(var = "modules") %>% 
  pivot_longer(cols = c(colnames(as.data.frame(kofamscan_mx))), names_to = "family", values_to = "completeness") %>% 
  separate(modules, into = c("id", "description"), sep = " \\(")
kofamscan_ave$description <- gsub("\\)", "", kofamscan_ave$description)

sup3 <- rbind(an_mx, ma_mx) 
colnames(sup3) <- c("module","module_description", "family", "completeness", "tool")
kofamscan_10 <- kofamscan_ave %>%
  filter(id %in% sup3$module) %>%
  inner_join(sup3, by = c("id" = "module", "family" = "family"))
kofamscan_10$tool <- "kofamscan"
colnames(kofamscan_10) <- c("module","module_description", "family", "completeness", "dumm1","dummy","tool")
sup3 <- rbind(sup3, kofamscan_10[c("module","module_description", "family", "completeness", "tool")])

anvio_10 <- anvio_ave %>%
  filter(id %in% sup3$module) %>%
  inner_join(sup3, by = c("id" = "module", "family" = "family"))
anvio_10$tool <- "anvio"
colnames(anvio_10) <- c("module","module_description", "family", "completeness", "dumm1","dummy","tool")
sup3 <- rbind(sup3, anvio_10[c("module","module_description", "family", "completeness", "tool")])

microbeannotator_10 <- microbeannotator_ave %>%
  filter(id %in% sup3$module) %>%
  inner_join(sup3, by = c("id" = "module", "family" = "family"))
microbeannotator_10$tool <- "microbeannotator"
colnames(microbeannotator_10) <- c("module","module_description", "family", "completeness", "dumm1","dummy","tool")
sup3 <- rbind(sup3, microbeannotator_10[c("module","module_description", "family", "completeness", "tool")])
sup3 <- unique(sup3)

# Write out to make table 3
write_tsv(sup3, file = "supp_table_3.csv")

# Calculate number of modules only found by anvio and not kofamscan or microbeannotator
# both on raw and only 80% completeness. This is done at the sample level first.
anvio_not_in_kofamscan <- anti_join(anvio, kofamscan, by = c(".id","module"))
anvio_not_in_other <- anti_join(anvio_not_in_kofamscan, microbeannotator, by = c(".id","module"))
anvio_not_in_kofamscan_80 <- anti_join(anvio_80, kofamscan_80, by = c(".id","module"))
anvio_not_in_other_80 <- anti_join(anvio_not_in_kofamscan_80, microbeannotator_80, by = c(".id","module"))
# now do microbeannotator
microbeannotator_not_in_kofamscan <- anti_join(microbeannotator, kofamscan, by = c(".id","module"))
microbeannotator_not_in_other <- anti_join(microbeannotator_not_in_kofamscan, anvio, by = c(".id","module"))
microbeannotator_not_in_kofamscan_80 <- anti_join(microbeannotator_80, kofamscan_80, by = c(".id","module"))
microbeannotator_not_in_other_80 <- anti_join(microbeannotator_not_in_kofamscan_80, anvio_80, by = c(".id","module"))

# Report the numer of modules found in total and at 80% only
print(paste("Number of modules found in total anvio:", nrow(anvio_not_in_other)))
print(paste("Number of unique modules found in total anvio:", length(unique(anvio_not_in_other$module))))
print(paste("Number of modules found in total at 80% anvio:", nrow(anvio_not_in_other_80)))
print(paste("Number of unique modules found in total at 80% anvio:", length(unique(anvio_not_in_other_80$module))))
print(paste("Number of modules found in total microbeannotator:", nrow(microbeannotator_not_in_other)))
print(paste("Number of unique modules found in total microbeannotator:", length(unique(microbeannotator_not_in_other$module))))
print(paste("Number of modules found in total at 80% microbeannotator:", nrow(microbeannotator_not_in_other_80)))
print(paste("Number of unique modules found in total at 80% microbeannotator:", length(unique(microbeannotator_not_in_other_80$module))))

# Calculate the drug resistance 
drugs_anvio <- anvio_mx[grep("Drug", row.names(anvio_mx)),]
drugs_microbeannotator <- microbeannotator_mx[grep("Drug", row.names(microbeannotator_mx)),]
drugs_kofamscan <- kofamscan_mx[grep("Drug", row.names(kofamscan_mx)),]

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
draw(d_list, row_title = "KEGG Metabolic Modules", column_title = "Mean Drug Metabolic Method Pathways Completeness Per Family", column_title_gp = gpar(fontsize = 16))
dev.off()

modules <- unique(anvio$module_subcategory)
for (i in unique(modules)) {
  filename <- paste0("heatmap_", gsub(" ", "_", str_to_title(i)), ".png")
  png(file = filename, width = 11.0, height = 9.25, units = "in", res = 1200)
  
  module_anvio <- anvio_mx[grep(i, rownames(anvio_mx)),]
  module_microbeannotator <- microbeannotator_mx[grep(i, rownames(microbeannotator_mx)),]
  module_kofamscan <- kofamscan_mx[grep(i, rownames(kofamscan_mx)),]
  
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
