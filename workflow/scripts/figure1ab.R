library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(plyr)


make_linker <- function(path){
  linker <- data.frame(read.csv(path, sep='\t', header = TRUE))
  linker$accession <- gsub('RS_', '', linker$accession)
  linker$accession <- gsub('GB_', '', linker$accession)
  family_linker <- linker[c("accession", "gtdb_taxonomy")]
  family_linker$family <- gsub('f__', '', str_split_fixed(linker$gtdb_taxonomy, ";",7)[,5])
  names(family_linker)[1] <- ".id"
  family_linker <- family_linker[c(".id", "family")]
  family_linker$.id <- sub("..$", "", family_linker$.id)
  return(family_linker)
}

# Function to split a string every nth space
split_every_nth_space <- function(string, n) {
  split_string <- unlist(strsplit(string, " "))
  split_indices <- seq(1, length(split_string), by = n)
  split_result <- lapply(split_indices, function(i) {
    end_index <- min(i + n - 1, length(split_string))
    paste(split_string[i:end_index], collapse = " ")
  })
  return(split_result)
}

read_in_files <- function(path, type){
  if(type == "kofamscan"){
    filenames <- list.files(path, full.names = TRUE)
    names(filenames) <- paste(gsub('.tsv','',basename(filenames)))
    return(ldply(filenames, read.csv, header=FALSE, sep='\t', comment.char='#'))
  } else if(type == "microbeannotator") {
    filenames <- list.files(list.files(list.files(path, full.names = TRUE), 
                                       full.names = TRUE, pattern = "kofam"), 
                            full.names = TRUE, pattern = "filt")
    names(filenames) <- paste(gsub('_cds.faa.kofam.filt','',basename(filenames)))
    ma <- ldply(filenames, read.csv, header=FALSE, sep='\t')
    return(ma[!str_detect(ma$V1, "#"),])
  } else {
    filenames <- list.files(path, full.names = TRUE)
    names(filenames) <- paste(gsub('.tsv','',basename(filenames)))
    an <- ldply(filenames, read.csv, header=FALSE, sep='\t', skip = 1)
    return(an)
  }
}

make_stray_kofam <- function(path){
  stray_kofam_hmm <- data.frame(gsub("NAME  ", "", grep("^NAME", readLines(path), value = TRUE)))
  names(stray_kofam_hmm) <- c("V3")
  return(stray_kofam_hmm)
}

clean_an <- function(an){
  an <- an[!str_detect(an$V2, "KEGG_BRITE"),]
  an <- an[!str_detect(an$V2, "KEGG_Class"),]
  an <- an[!str_detect(an$V2, "KEGG_Module"),] 
  an$.id <- sub("..$", "", an$.id)
  return(an)
}

clean_ko <- function(ko){
  ko <- ko[!str_detect(ko$V2, "KEGG_BRITE"),]
  ko <- ko[!str_detect(ko$V2, "KEGG_Class"),]
  ko <- ko[!str_detect(ko$V2, "KEGG_Module"),]  
  ko <- ko[ko[,2] == "*", ]
  ko$.id <- sub("..$", "", ko$.id)
  return(ko)
}

clean_ma <- function(ma){
  ma$.id <- sub(".faa.kofam.filt", "", ma$.id)
  ma$.id <- sub("..$", "", ma$.id)
  return(ma)
}

make_ortholog_per_tool_plot <- function(annr_num, anra_num, anst_num, kora_num, 
                                        kore_num, mare_num, mara_num, mare_false_num,
                                        mare_unknown_num, kore_false_num, 
                                        kore_unknown_num){
  df <- data.frame(Orthologs_retrieved=c(annr_num, anra_num, anst_num, kora_num,
                                         (kore_num - (kore_false_num+kore_unknown_num)),
                                         (mare_num - (mare_false_num+mare_unknown_num)), 
                                         (mare_num - (mare_false_num+mare_unknown_num))),
                   False=c(0, 0, 0, 0, kore_false_num, mare_false_num, mare_false_num),
                   Unverified=c(0, 0, 0, 0, kore_unknown_num, mare_unknown_num, mare_unknown_num),
                   Method=c("Anvio (No Threshold)", "Anvio (Threshold)", "Anvio (Threshold+Strays)", "Kofamscan (Default)",
                            "Kofamscan (Relaxed)", "MicrobeAnnotator (Refined Mode)", "MicrobeAnnotator (Default)"))
  
  # Assigning Group based on Method
  df$Group <- ifelse(df$Method %in% c("Anvio (No Threshold)", "Anvio (Threshold)", "Anvio (Threshold+Strays)"), "Anvio", 
                     ifelse(df$Method %in% c("MicrobeAnnotator (Refined Mode)", "MicrobeAnnotator (Default)"), "MicrobeAnnotator", "Kofamscan"))
  
  # Create the plot
  max_y <- max(df$Orthologs_retrieved) + 300
  a <- ggplot(df) +
    geom_bar(aes(x = Method, y = Orthologs_retrieved, fill = Group), 
             stat = "identity", position = "stack") +
    geom_bar(aes(x = Method, y = False, fill = "False Positive"), 
             stat = "identity", position = "stack", alpha = 1) +  
    geom_bar(aes(x = Method, y = Unverified, fill = "Unverified"), 
             stat = "identity", position = "stack", alpha = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c("Anvio"="#2a9d8f", "Kofamscan"="#e9c46a", "MicrobeAnnotator"="#6DB9E4", 
                                 "False Positive" = "#708090", "Unverified"="#9e9e9e")) +
    labs(x = "Method", y = "Orthologs Recovered", title = "Total Unique Orthologs Recovered by Method") +
    coord_cartesian(ylim = c(0, max_y + 500)) +
    geom_text(aes(x = Method, y = Orthologs_retrieved - 250, label = Orthologs_retrieved), 
              vjust = 1, size = 3.5, color = "white",
              fontface = 'bold') +
    guides(fill = guide_legend(title = "Software", override.aes = list(alpha = 1))) +
    geom_rect(aes(fill = Group, xmin = 5.55, xmax = 6.45, 
                  ymin = 0, ymax = 1000), alpha = 1)  +
    geom_rect(aes(fill = Group, xmin = 6.55, xmax = 7.45, 
                  ymin = 0, ymax = 1000), alpha = 1) + 
    annotate('rect', xmin = 4.555, xmax = 5.445, 
             ymin = 0, ymax = 2500,
             fill = "#e9c46a", col = "#e9c46a") +
    annotate('rect', xmin = 4.555, xmax = 5.445, 
             ymin = (kore_num - kore_false_num), ymax = (kore_num),
             fill = "#708090", col = "#708090") +
    annotate('rect', xmin = 4.555, xmax = 5.445, 
             ymin = (kore_num - (kore_false_num+kore_unknown_num)), 
             ymax = (kore_num-kore_false_num), fill = "#9e9e9e", col = "#9e9e9e") +
    annotate('rect', xmin = 6.555, xmax = 7.445, 
             ymin = (mare_num - mare_false_num), ymax = (mare_num),
             fill = "#708090", col = "#708090") +
    annotate('rect', xmin = 6.555, xmax = 7.445, 
             ymin = (mare_num - (mare_false_num+mare_unknown_num)), 
             ymax = (mare_num-mare_false_num), fill = "#9e9e9e", col = "#9e9e9e") +
    annotate('text', x = 5,
             y = (kore_num - (kore_false_num+kore_unknown_num)) + 300,
             label = kore_unknown_num,
             fontface = 'bold',  # font style
             size = 3.5,  # font size
             color = "white"  # font color
    ) + 
    annotate('text', x = 5,
             y = (kore_num - kore_false_num) + 250,
             label = kore_false_num,
             fontface = 'bold',  # font style
             size = 3.5,  # font size
             color = "white"  # font color
    ) +
    annotate('text', x = 7,
             y = (mare_num - (mare_false_num+mare_unknown_num)) + 300,
             label = mare_unknown_num,
             fontface = 'bold',  # font style
             size = 3.5,  # font size
             color = "white"  # font color
    ) + 
    annotate('text', x = 7,
             y = (mare_num - mare_false_num) + 250,
             label = mare_false_num,
             fontface = 'bold',  # font style
             size = 3.5,  # font size
             color = "white"  # font color
    ) +
    annotate('rect', xmin = 5.555, xmax = 6.445, 
             ymin = (mare_num - mare_false_num), ymax = (mare_num),
             fill = "#708090", col = "#708090") +
    annotate('rect', xmin = 5.555, xmax = 6.445, 
             ymin = (mare_num - (mare_false_num+mare_unknown_num)), 
             ymax = (mare_num-mare_false_num), fill = "#9e9e9e", col = "#9e9e9e") +
    annotate('text', x = 6,
             y = (mare_num - (mare_false_num+mare_unknown_num)) + 300,
             label = mare_unknown_num,
             fontface = 'bold',  # font style
             size = 3.5,  # font size
             color = "white"  # font color
    ) + 
    annotate('text', x = 6,
             y = (mare_num - mare_false_num) + 250,
             label = mare_false_num,
             fontface = 'bold',  # font style
             size = 3.5,  # font size
             color = "white"  # font color
    )
  
  return(a)
}

uniq <- function(ko, an, ma, maf, family){
  ko_family <- split(ko, ifelse(ko$family %in% family, family[match(ko$family, family)], "Other"))
  an_family <- split(an, ifelse(an$family %in% family, family[match(an$family, family)], "Other"))
  ma_family <- split(ma, ifelse(ma$family %in% family, family[match(ma$family, family)], "Other"))
  maf_family <- split(maf, ifelse(maf$family %in% family, family[match(maf$family, family)], "Other"))
  
  out <- list()
  for(i in seq_along(ko_family)){
    ko_uniq_an <- unique(ko_family[[i]]$V3[!(ko_family[[i]]$V3 %in% an_family[[i]]$V3)])
    ko_uniq_ma <- unique(ko_family[[i]]$V3[!(ko_family[[i]]$V3 %in% ma_family[[i]]$V3)])
    an_uniq_ko <- unique(an_family[[i]]$V3[!(an_family[[i]]$V3 %in% ko_family[[i]]$V3)])
    an_uniq_ma <- unique(an_family[[i]]$V3[!(an_family[[i]]$V3 %in% ma_family[[i]]$V3)])
    ma_uniq_ko <- unique(ma_family[[i]]$V3[!(ma_family[[i]]$V3 %in% ko_family[[i]]$V3)])
    ma_uniq_an <- unique(ma_family[[i]]$V3[!(ma_family[[i]]$V3 %in% an_family[[i]]$V3)])
    maf_out <- length(unique(maf_family[[i]]$V3))
                      
    if (length(setdiff(ko_uniq_an, an_uniq_ko)) == 0 || length(setdiff(ko_uniq_ma, ma_uniq_ko)) == 0){
      ko_uniq <- 0
    } else {
      ko_uniq <- length((unique(c(setdiff(ko_uniq_an, an_uniq_ko), setdiff(ko_uniq_ma, ma_uniq_ko)))))
    }
    if (length(setdiff(an_uniq_ma, ma_uniq_an)) == 0 || length(setdiff(an_uniq_ko, ko_uniq_an)) == 0){
      an_uniq <- 0
    } else {
      an_uniq <- length((unique(c(setdiff(an_uniq_ko, ko_uniq_an), setdiff(an_uniq_ma, ma_uniq_an)))))
    }
    if (length(setdiff(ma_uniq_an, an_uniq_ma)) == 0 || length(setdiff(ma_uniq_ko, ko_uniq_ma)) == 0){
      ma_uniq=0
    } else {
      ma_uniq <- length((unique(c(setdiff(ma_uniq_an, an_uniq_ma), setdiff(ma_uniq_ko, ko_uniq_ma)))))
    }
    out <- c(out,list(ko_uniq,an_uniq,(ma_uniq-maf_out),maf_out))
  }
  return(out)
}


make_stacked_bar_data <- function(ko_num, family, ma_false){
  family <-  c(rep("Bacteroidaceae", 4), rep("Nanosynbacteraceae", 4), 
               rep("Streptococcaceae", 4),rep("Lachnospiraceae", 4), 
               rep("Enterobacteriaceae", 4), rep("Synechococcaceae", 4),
               rep("Rhodobacteraceae", 4), rep("Cyanobiaceae", 4), 
               rep("Pelagibacteraceae", 4), rep("Streptomycetaceae", 4), 
               rep("Rhizobiaceae", 4))
  condition <- rep(c("Only Kofamscan", "Only Anvi'o", "Only MicrobeAnnotator", 
                     "Only MicrobeAnnotator Strays"), 11)
  value <- c(ko_num[[ 1 ]],ko_num[[ 2 ]],ko_num[[ 3 ]],ko_num[[ 4 ]],
             ko_num[[ 5 ]],ko_num[[ 6 ]],ko_num[[ 7 ]],ko_num[[ 8 ]],
             ko_num[[ 9 ]],ko_num[[ 10 ]],ko_num[[ 11 ]],ko_num[[ 12 ]],
             ko_num[[ 13 ]],ko_num[[ 14 ]],ko_num[[ 15 ]],ko_num[[ 16 ]],
             ko_num[[ 17 ]],ko_num[[ 18 ]],ko_num[[ 19 ]],ko_num[[ 20 ]],
             ko_num[[ 21 ]],ko_num[[ 22 ]],ko_num[[ 23 ]],ko_num[[ 24 ]],
             ko_num[[ 25 ]],ko_num[[ 26 ]],ko_num[[ 27 ]],ko_num[[ 28 ]],
             ko_num[[ 29 ]],ko_num[[ 30 ]],ko_num[[ 31 ]],ko_num[[ 32 ]],
             ko_num[[ 33 ]],ko_num[[ 34 ]],ko_num[[ 35 ]],ko_num[[ 36 ]],
             ko_num[[ 37 ]],ko_num[[ 38 ]],ko_num[[ 39 ]],ko_num[[ 40 ]],
             ko_num[[ 41 ]],ko_num[[ 42 ]],ko_num[[ 43 ]],ko_num[[ 44 ]])
  environment=c("Gut","Gut","Gut","Gut",
                "Oral","Oral","Oral","Oral",
                "Oral","Oral","Oral","Oral",
                "Gut","Gut","Gut","Gut",
                "Control","Control","Control","Control",
                "Coastal", "Coastal", "Coastal", "Coastal", 
                "Coastal", "Coastal", "Coastal", "Coastal",
                "Open Ocean", "Open Ocean", "Open Ocean", "Open Ocean", 
                "Open Ocean","Open Ocean","Open Ocean","Open Ocean",
                "Soil","Soil", "Soil","Soil",
                "Soil","Soil","Soil","Soil")
  return(data.frame(environment,family,condition,value))
}

make_stacked_bar <- function(ko_num, family, ma_false) {
  custom_colors <- c("#2a9d8f","#e9c46a","#6DB9E4","#708090")
  # Get data and calculate counts for each family
  df <- make_stacked_bar_data(ko_num, family, ma_false)
  df$condition <- factor(df$condition, levels=c("Only Anvi'o", "Only Kofamscan", 
                                                "Only MicrobeAnnotator", 
                                                "Only MicrobeAnnotator Strays"))
  
  family_counts <- table(df$family)
  # Reorder levels of the family factor based on counts (largest to smallest)
  df$family <- factor(df$family, levels = names(sort(family_counts, decreasing = TRUE)))
  
  b <- ggplot() +
    geom_bar(data = df, aes(x = family, y = value, fill = condition), 
             stat = "identity", position = position_dodge(width = 0.3)) +
    theme_bw() + 
    facet_wrap(~ environment, scales = "free_x", nrow = 1, strip.position = "top", switch = "x") +
    ggtitle("Annotations Found for Methods per Family") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ylab("Number of Kegg Orthologues") +
    xlab("Family") +
    scale_fill_manual(values = custom_colors) +
    guides(fill = guide_legend(title = "Methods"))
  
  return(b)
}

make_stacked_bar(ko_num, family, strays_mare)

# Set 11 families we are interested in
family <- c("Bacteroidaceae","Nanosynbacteraceae", "Streptococcaceae",
            "Lachnospiraceae", "Enterobacteriaceae", "Synechococcaceae",
            "Rhodobacteraceae", "Cyanobiaceae", "Pelagibacteraceae",
            "Streptomycetaceae", "Rhizobiaceae")

# Read in metadata files needed to generate graphs
family_linker <- make_linker("/Users/kananen/Desktop/ImHere/bac120_metadata_r214.tsv")
stray_kofam_hmm <- make_stray_kofam("/Users/kananen/Desktop/ImHere/hmm_profiles_with_kofams_with_no_threshold.hmm")

# Read in Kofamscan raw and refined files
kora <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/kofamscan/raw/functions/", "kofamscan")
kore <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/kofamscan/refined/functions/", "kofamscan")
# Read in MicrobeAnnotator raw and refined files
mara <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/microbeannotator/functions/raw", "microbeannotator")
mare <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/microbeannotator/functions/refined", "microbeannotator")
# Read in Anvio raw, stray, and no recovery files
anra <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/anvio/default/functions", "anvio")
anst <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/anvio/stray/functions", "anvio")
annr <- read_in_files("/Users/kananen/Desktop/ImHere/set_36/anvio/no_rescue/functions", "anvio")

# Clean up anvi'o and kofamscan output
anra <- clean_an(anra); anst <- clean_an(anst); annr <- clean_an(annr)
kora <- clean_ko(kora); kore <- clean_ko(kore)
mara <- clean_ma(mara); mare <- clean_ma(mare)
# Match the species to the genome
kora <- merge(kora, family_linker, by = ".id")
kore <- merge(kore, family_linker, by = ".id")
annr <- merge(annr, family_linker, by = ".id")
anra <- merge(anra, family_linker, by = ".id")
anst <- merge(anst, family_linker, by = ".id")
mara <- merge(mara, family_linker, by = ".id")
mare <- merge(mare, family_linker, by = ".id")

# Identify Strays
strays_anst <- merge(anst, stray_kofam_hmm, by = "V3")
strays_mare <- merge(mare, stray_kofam_hmm, by="V3")

# Get unique number of ko found between the the tools. We use the modes that
# Return the most found hits.
annr_num <- length(unique(annr$V3))
anra_num <- length(unique(anra$V3))
anst_num <- length(unique(anst$V3))
kora_num <- length(unique(kora$V3))
kore_num <- length(unique(kore$V3))
mare_num <- length(unique(mare$V3))
mara_num <- length(unique(mara$V3))

# Make dataframe of all the strays that would meet filtering criteria of eval=0.01
# The default used by Kofamkoala online.
not_standard_kora <- length(unique(subset(kora, as.numeric(V6) > 0.01)))
not_standard_kore <- length(unique(subset(kore, as.numeric(V6) > 0.01)))
not_standard_mare <- length(unique(subset(mare, as.numeric(V6) > 0.01)))
not_standard_anra <- length(unique(subset(anra, as.numeric(V5) > 0.01)))
not_standard_mara <- length(unique(subset(mara, as.numeric(V6) > 0.01)))
not_standard_anst <- length(unique(subset(anst, as.numeric(V5) > 0.01)))
not_standard_annr <- length(unique(subset(annr, as.numeric(V5) > 0.01)))

# We need to test kofamscan's relaxed params explicitly 
not_standard_kore <- subset(kore, !(kore$V3 %in% kora$V3))[c(".id","V2", "V3")]
not_standard_kore <- unique(not_standard_kore)


# Pause here to see who has e-values above 0.01 and run them in eggnog mapper
# if they do.
eggnog_mare <- read.csv("/Users/kananen/Documents/upset/upset_R/MM__px08q43.emapper.annotations.tsv", sep='\t', comment.char = "#",
                        header=FALSE)
eggnog_kore <- read.csv("/Users/kananen/Documents/upset/upset_R/MM_4pqcene3.emapper.annotations.tsv", sep='\t', comment.char = "#",
header=FALSE)

eggnog_mare <- eggnog_mare[c("V1","V3","V12")]
eggnog_kore <- eggnog_kore[c("V1","V3","V12")]
eggnog_kore$V12 <- gsub("ko:",'',eggnog_kore$V12)
eggnog_mare$V12 <- gsub("ko:",'',eggnog_mare$V12)

split_values <- strsplit(as.character(eggnog_kore$V1), "_")
eggnog_kore$Initial <- sapply(split_values, function(x) x[1])
eggnog_kore$.id <- sapply(split_values, function(x) paste(x[-1], collapse = "_"))

split_values <- strsplit(as.character(eggnog_mare$V1), "_")
eggnog_mare$Initial <- sapply(split_values, function(x) x[1])
eggnog_mare$.id <- sapply(split_values, function(x) paste(x[-1], collapse = "_"))

# Remove the old 'V1' column
eggnog_mare <- subset(eggnog_mare, select = -V1)
eggnog_mare$V1 <- eggnog_mare$Initial
eggnog_mare$V6 <- eggnog_mare$V3
eggnog_mare$V3 <- eggnog_mare$V12
eggnog_kore <- subset(eggnog_kore, select = -V1)
eggnog_kore$V1 <- eggnog_kore$Initial
eggnog_kore$V6 <- eggnog_kore$V3
eggnog_kore$V3 <- eggnog_kore$V12

mare_sub <- subset(mare, select = c(.id, V1, V3, V6))
kore_sub <- subset(not_standard_kore, select = c(.id, V2, V3))
kore_sub$V1 <- kore_sub$V2

eggnog_mare <- merge(subset(mare_sub, as.numeric(V6) > 0.01), eggnog_mare[c("V1", ".id", "V3", "V6")], by=c("V1", ".id"))
eggnog_kore <- merge(kore_sub, eggnog_kore[c("V1", ".id", "V3", "V6")], by=c("V1", ".id"))


# Filter the dataframe based on whether values from V3 are present in rows of V2
mare_unknown <- eggnog_mare %>% 
  filter(str_detect(V3.y, "-"))
mare_known <- eggnog_mare %>% 
  filter(!str_detect(V3.y, "-"))
kore_unknown <- eggnog_kore %>% 
  filter(str_detect(V3.y, "-"))
kore_known <- eggnog_kore %>% 
  filter(!str_detect(V3.y, "-"))


kore_false <- kore_known %>%
  filter(!str_detect(V3.y, pull(., V3.x))) %>%
  select(V1, .id, V3.x, V6, V3.y)
kore_false <- merge(kore_false, family_linker, by=".id")
mare_false <- mare_known %>%
  filter(!str_detect(V3.y, pull(., V3.x))) %>%
  select(V1, .id, V3.x, V6.x, V3.y, V6.y)
mare_false <- merge(mare_false, family_linker, by=".id")


kore_unknown_num  <- length(unique(kore_unknown$V3.x))
kore_false_num <- length(unique(kore_false$V3.x))
mare_unknown_num  <- length(unique(mare_unknown$V3.x))
mare_false_num <- length(unique(mare_false$V3.x))

# Make F1a Barplot
a <- make_ortholog_per_tool_plot(annr_num, anra_num, anst_num, 
                                 kora_num, kore_num,
                                 mare_num, mara_num,
                                 mare_false_num, mare_unknown_num,
                                 kore_false_num, kore_unknown_num)

# Generate the difference in the sets for a cross tool comparison
ko_num <- uniq(kora, anst, mare, strays_mare, family)

# Graph stacked barplot (F1c) of the total unique values per tool compared to the 
# other tools / methods.
b <- make_stacked_bar(ko_num, family, strays_mare)

# Write files out that are needed for subsequent analysis
write.csv(not_standard_mare, file = "not_standard_mare_strays.csv", row.names = FALSE, col.names = FALSE)
write.csv(not_standard_kore[c(".id", "V2")], file = "not_standard_kore_strays.csv", row.names = FALSE)
ggsave(filename = "plot1.png", plot = a, width = 6, height = 4)
ggsave(filename = "plot1.png", plot = b, width = 6, height = 4)
