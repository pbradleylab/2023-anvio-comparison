library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(plyr)
library(tidyverse)
library(rstatix)
library(jsonlite)
library(tidyverse)
library(ggpubr)


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

read_in_files <- function(path, type){
  filenames <- list.files(path, full.names = TRUE)
  if(type == "kofamscan"){
    names(filenames) <- paste(gsub('.tsv','',basename(filenames)))
    return(ldply(filenames, read.csv, header=FALSE, sep='\t', comment.char='#'))
  } else if(type == "microbeannotator") {
    filenames <- list.files(list.files(list.files(path, full.names = TRUE), 
                                       full.names = TRUE, pattern = "kofam"), 
                            full.names = TRUE, pattern = "filt")
    names(filenames) <- paste(gsub('_cds.faa.kofam.filt','',basename(filenames)))
    ma <- ldply(filenames, read.csv, header=FALSE, sep='\t')
    return(ma[!str_detect(ma$V1, "#"),])
  } else if(type == "eggnog"){
    filenames <- list.files(path, full.names = TRUE, pattern = ".emapper.annotations")
    names(filenames) <- paste(gsub('.emapper.annotations','',basename(filenames)))
    en <- ldply(filenames, read.csv, header=FALSE, sep='\t', comment.char='#')
    en <- en[c(".id","V1","V3","V11","V12","V17")]
    en$V12 <- gsub("ko:",'', en$V12)
    en$.id <- sub("..$", "", en$.id)
    return(en)
  } else {
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
  an <- an %>% 
    group_by(.id, V1) %>% 
    slice_min(V5)
  return(an)
}

clean_ko <- function(ko){
  ko <- ko[ko[,2] == "*", ]
  ko$.id <- sub("..$", "", ko$.id)
  ko <- ko %>% 
    group_by(.id, V2) %>% 
    slice_min(V6, n=1) %>%
    slice_max(V5, n=1) %>%
    slice_max(V4, n=1) %>%
    arrange(V1) %>%
    filter(row_number() == 1) %>%
    ungroup()
  return(ko)
}

clean_ma <- function(ma){
  ma$.id <- sub(".faa.kofam.filt", "", ma$.id)
  ma$.id <- sub("..$", "", ma$.id)
  ma <- ma %>% 
    group_by(.id, V1) %>% 
    slice_min(V6)
  return(ma)
}

check_not_ko <- function(text, ko_needed) {
  all(!str_detect(text, ko_needed))
}

check_ko <- function(text, ko_needed) {
  all(str_detect(text, ko_needed))
}

count_function <- function(data) {
  table(data$.id)
}

get_aggregation <- function(colname1, colname2, df1, df2, name, stats){
  # Create a data frame for the counts
  count1 <- count_function(df1)
  count2 <- count_function(df2)
  tmp <- data.frame(.id = unique(c(names(count1), names(count2))),
                    colname1 = count1[match(unique(c(names(count1), names(count2))), names(count1))],
                    colname2 = count2[match(unique(c(names(count1), names(count2))), names(count2))])
  tmp[is.na(tmp)] <- 0
  
  # Calculate the increase and add to stats
  stats[[name]] <- ((tmp$colname1.Freq - tmp$colname2.Freq) / tmp$colname2.Freq)
  stats <- merge(stats, tmp[c("colname1.Freq", ".id")])
  stats <- merge(stats, tmp[c("colname2.Freq", ".id")])
  names(stats)[which(names(stats) == "colname1")] <- colname1
  names(stats)[which(names(stats) == "colname2")] <- colname2
  return(stats)
}

create_threshold_plot <- function(data_main, data_stray, title, show_below = TRUE) {
  ggplot() +
    geom_jitter(data = data_main, aes(x = name, y = value, color = color_group, shape = "nt-KO"),
                position = position_jitter(width = 0.2), alpha = 1, size = 3) +
    geom_point(data = data_stray, aes(x = name, y = value, color = color_group, shape = "nt-KO"),
               size = 3, position = position_jitter(width = 0.2), alpha = ifelse(show_below, 1, 0.5)) +
    scale_color_manual(values = custom_colors, 
                       labels = ifelse(show_below, c("Above", "Below"), "Above")) +
    scale_shape_manual(values = c("nt-KO" = 15)) +
    labs(title = title,
         x = "Category",
         y = "log10(E-Value)",
         color = "Threshold",
         shape = "KO Type") +
    theme_light() +
    theme(text = element_text(size = 20)) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue")
}


merge_eggnog_no_ko_match <- function(df, eggnog, type) {
  df <- unique(merge(df[c(".id", "V1", type)], eggnog, by = c(".id", "V1"))) #684790
  df <- df[!str_detect(df$V12, df[[type]]), ] #576587
  df <- df[complete.cases(df$.id), ] #576587
  return(df)
}

merge_eggnog <- function(df, eggnog, type) {
  df <- unique(merge(df[c(".id", "V1", type)], eggnog, by = c(".id", "V1")))
  df <- df[str_detect(df$V12, df[[type]]), ]
  df <- df[complete.cases(df$.id), ]
  return(df)
}

# Set 11 families we are interested in
family <- c("Bacteroidaceae","Nanosynbacteraceae", "Streptococcaceae",
            "Lachnospiraceae", "Enterobacteriaceae", "Microcystaceae",
            "Rhodobacteraceae", "Cyanobiaceae", "Pelagibacteraceae",
            "Streptomycetaceae", "Rhizobiaceae")

# Read in metadata files needed to generate graphs
family_linker <- make_linker("/Users/kananen/Desktop/Keep_until_anvio_published/bac120_metadata_r214.tsv")
stray_kofam_hmm <- make_stray_kofam("/Users/kananen/Desktop/Keep_until_anvio_published/hmm_profiles_with_kofams_with_no_threshold.hmm")

# Read in Kofamscan raw and refined files
kora <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/kofamscan/functions/default/", "kofamscan")
kore <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/kofamscan/functions/refined/", "kofamscan")
# Read in MicrobeAnnotator raw and refined files
mara <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/microbeAnnotator/functions/default/", "microbeannotator")
mare <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/microbeAnnotator/functions/refined/", "microbeannotator")
# Read in Anvio raw, stray, and no recovery files
anra <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/anvio/functions/default/", "anvio")
anst <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/anvio/functions/stray/", "anvio")
annr <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/anvio/functions/no_heuristic/", "anvio")

# Pause here to see who has e-values above 0.01 and run them in eggnog mapper
# if they do. Read in all the files in a given directory
eggnog <- read_in_files("/Users/kananen/Desktop/Keep_until_anvio_published/unfiltered/eggnog/", "eggnog")
# Get the EC and BRITE 
eggnog_ecbrite <- eggnog[c(".id","V1","V11","V12","V17")]
colnames(eggnog_ecbrite) <- c(".id", "V1", "eg_ec", "eg_ko", "eg_brite")

# Calculated ORFS from worflow in table 1 based off of the number
orfs <- read.csv("/Users/kananen/Documents/Upset/upset_R/f1_table.tsv", sep='\t')
orfs$accession <- sub("..$", "", orfs$accession)

# Here we need to make three versions of kofam in the cleaning process for later comparisons
kora <- clean_ko(kora); kore <- clean_ko(kore)
mara <- clean_ma(mara); mare <- clean_ma(mare)
anra <- clean_an(anra); anst <- clean_an(anst); annr <- clean_an(annr)

# Match the species to the genome
kora <- merge(kora, family_linker, by = ".id")
kore <- merge(kore, family_linker, by = ".id")
mara <- merge(mara, family_linker, by = ".id")
mare <- merge(mare, family_linker, by = ".id")
annr <- merge(annr, family_linker, by = ".id")
anra <- merge(anra, family_linker, by = ".id")
anst <- merge(anst, family_linker, by = ".id")

# Only take the columns we need and use the same naming convention
kora_sub <- kora[c(".id", "V2", "V3", "V6")]; colnames(kora_sub) <- c(".id", "V1", "kora", "kora_eval")
kore_sub <- kore[c(".id", "V2", "V3", "V6")]; colnames(kore_sub) <- c(".id", "V1", "kore", "kore_eval")
anra_sub <- anra[c(".id", "V1", "V3", "V5")]; colnames(anra_sub) <- c(".id", "V1", "anra", "anra_eval")
annr_sub <- annr[c(".id", "V1", "V3", "V5")]; colnames(annr_sub) <- c(".id", "V1", "annr", "annr_eval")
anst_sub <- anst[c(".id", "V1", "V3", "V5")]; colnames(anst_sub) <- c(".id", "V1", "anst", "anst_eval")
mara_sub <- mara[c(".id", "V1", "V3", "V6")]; colnames(mara_sub) <- c(".id", "V1", "mara", "mara_eval")
mare_sub <- mare[c(".id", "V1", "V3", "V6")]; colnames(mare_sub) <- c(".id", "V1", "mare", "mare_eval")

# Calculate the eggnog hits that match
eggnog_kora <- merge_eggnog(kora_sub[c(".id", "V1", "kora")], eggnog, "kora")
eggnog_kore <- merge_eggnog(kore_sub[c(".id", "V1", "kore")], eggnog, "kore")
eggnog_mara <- merge_eggnog(mara_sub[c(".id", "V1", "mara")], eggnog, "mara")
eggnog_mare <- merge_eggnog(mare_sub[c(".id", "V1", "mare")], eggnog, "mare")
eggnog_anra <- merge_eggnog(anra_sub[c(".id", "V1", "anra")], eggnog, "anra")
eggnog_annr <- merge_eggnog(annr_sub[c(".id", "V1", "annr")], eggnog, "annr")
eggnog_anst <- merge_eggnog(anst_sub[c(".id", "V1", "anst")], eggnog, "anst")

# Calculate statistics needed for paper
stats <- orfs[c("accession", "number_genes")]; colnames(stats) <- c(".id", "num_ORFs")
# Percentage change for kofamscan refined and kofamscan default
stats <- merge(stats, get_aggregation("kore","kora", kore_sub, kora_sub, "kore_increase_kora", stats), by = c(".id", "num_ORFs"), all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'kore'
names(stats)[names(stats) == 'colname2.Freq'] <- 'kora'
stats <- merge(stats, get_aggregation("kore_valid","kora_valid", eggnog_kore, eggnog_kora, "kore_increase_kora_valid", stats), all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'kore_valid'
names(stats)[names(stats) == 'colname2.Freq'] <- 'kora_valid'
# Percentage change for MicrobeAnnotator refined and MicrobeAnnotator default
stats <- merge(stats, get_aggregation("mare","mara", mare_sub, mara_sub, "mare_increase_mara", stats), all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'mare'
names(stats)[names(stats) == 'colname2.Freq'] <- 'mara'
stats <- merge(stats, get_aggregation("mare_valid","mara_valid", eggnog_mare, eggnog_mara, "mare_increase_mara_valid", stats), all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'mare_valid'
names(stats)[names(stats) == 'colname2.Freq'] <- 'mara_valid'
# Percentage change for anvi'o default and anvi'o no heuristic
stats <- merge(stats, get_aggregation("anra","annr", anra_sub, annr_sub, "anra_increase_annr", stats), all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'anra'
names(stats)[names(stats) == 'colname2.Freq'] <- 'annr'
stats <- merge(stats, get_aggregation("anra_valid","annr_valid", eggnog_anra, eggnog_annr, "anra_increase_annr_valid", stats), all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'anra_valid'
names(stats)[names(stats) == 'colname2.Freq'] <- 'annr_valid'
# Percentage change for anvi'o stray and anvi'o default
stats <- merge(stats, get_aggregation("anst","anra", anst_sub, anra_sub, "anst_increase_anra", stats)[c(".id","colname1.Freq")], all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'anst'
stats <- merge(stats, get_aggregation("anst_valid","anra_valid", eggnog_anst, eggnog_anra, "anst_increase_anra_valid", stats)[c(".id","colname1.Freq")], all=TRUE)
names(stats)[names(stats) == 'colname1.Freq'] <- 'anst_valid'

kora_stat <- nrow(eggnog_kora)/nrow(kora_sub)
anra_stat <- (nrow(eggnog_anra) - nrow(eggnog_kora)) / (nrow(anra_sub) - nrow(kora_sub))
mara_stat <- (nrow(eggnog_mara) - nrow(eggnog_kora)) / (nrow(mara_sub) - nrow(kora_sub))

# Calculate fraction matching to the EC
mara_ec <- mara %>%
  mutate(EC_mara = sub(".*\\[EC:([0-9\\.]+)]", "\\1", V20)) %>%
  mutate(EC_mara = sub(".*\\[EC:([^]]+)]", "\\1", V20),
         EC_mara = gsub(" ", ",", EC_mara))
mara_ec <- mara_ec[!(grepl("[a-zA-Z]", mara_ec$EC)), ]
mara_ec$.id <- gsub("\\.faa\\.kofam\\.filt", "", mara_ec$.id)

anra_brite <- anra[anra$V2 == "KEGG_BRITE", ]
anra_brite$V3 <- gsub("!!!", ",", anra_brite$V3)
anra_brite <- anra_brite[c(".id","V1","V3")]
anra_ec <- anra[anra$V2 == "KOfam", ]
anra_ec <- anra_ec %>%
  mutate(EC_anra = sub(".*\\[EC:([0-9\\.]+)]", "\\1", V4)) %>%
  mutate(EC_anra = sub(".*\\[EC:([^]]+)]", "\\1", V4),
         EC_anra = gsub(" ", ",", EC_anra))
anra_ec <- anra_ec[c(".id","V1","V3","EC_anra")]

eggnog_anra_anno <- unique(merge(merge_eggnog_no_ko_match(anra_sub[c(".id", "V1", "anra")], eggnog, "anra"), 
                                 eggnog_ecbrite, by = c(".id", "V1"), all.x = TRUE))
eggnog_anra_anno  <- merge(eggnog_anra_anno, anra_ec[c(".id","V1","EC_anra")], by=c(".id", "V1"), all.x=TRUE)
eggnog_mara_anno  <- unique(merge(merge_eggnog_no_ko_match(mara_sub[c(".id", "V1", "mara")], eggnog, "mara"), 
                                  eggnog_ecbrite, by=c(".id", "V1"), , all.x = TRUE))
eggnog_mara_anno  <- merge(eggnog_mara_anno, mara_ec[c(".id","V1","EC_mara")], by=c(".id", "V1"), all.x=TRUE)

no_eggnog_anra_anno <- eggnog_anra_anno %>%
  rowwise() %>%
  filter(!anra %in% strsplit(eg_ko, ",")) %>%
  filter((grepl("-", eg_ec)))
no_eggnog_mara_anno <- eggnog_mara_anno %>%
  rowwise() %>%
  filter(!mara %in% strsplit(eg_ko, ",")) %>%
  filter((grepl("-", eg_ec)))

eggnog_mara_anno <- eggnog_mara_anno %>%
  rowwise() %>%
  filter(!mara %in% strsplit(eg_ko, ",")) %>%
  filter(!(grepl("-", eg_ec)))
eggnog_anra_anno <- eggnog_anra_anno %>%
  rowwise() %>%
  filter(!anra %in% strsplit(eg_ko, ",")) %>%
  filter(!(grepl("-", eg_ec)))

anra_no_ec <- eggnog_anra_anno %>% 
  filter(is.na(EC_anra))
# If any value matches 100% i.e 1.2.3.4,2.3.4.5 and 2.3.4.5 would be a match.
anra_match_100 <- eggnog_anra_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Loop-like action on each value
      e1 <- eg_ec
      stre1 <- strsplit(e1, ",")[[1]]
      e2 <- EC_anra
      stre2 <- strsplit(e2, ",")[[1]]
      any(stre1 %in% stre2)
    }
  ) %>%
  filter(`match`)
# Match the first digit and only count if there are no matches. Cases such as 
# 1.2.3.4,2.3.4.5 compared with 2.1.1.1 should still count as a match under this 
# logic and thus not be included here.
anra_match_1_not_2 <- eggnog_anra_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_anra, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      any(is_digit & is_true_1 & !(is_true_2))
    }
  ) %>%
  filter(`match`) %>%
  ungroup()
# Now if it matches the first two but not the third one
anra_match_12_not_3 <- eggnog_anra_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_anra, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb1c <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2c <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))

      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      is_true_3 <- length(intersect(sb1c, sb2c)) > 0
      
      any(is_digit & is_true_1 & is_true_2 & !(is_true_3))
    }
  ) %>%
  filter(`match`) %>%
  ungroup()
# Now if it matches the first three but not the fourth one
anra_match_123_not_4 <- eggnog_anra_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_anra, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb1c <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb1d <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2c <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb2d <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      is_true_3 <- length(intersect(sb1c, sb2c)) > 0
      is_true_4 <- length(intersect(sb1d, sb2d)) > 0
      
      any(is_digit & is_true_1 & is_true_2 & is_true_3 & !(is_true_4))
          }
  ) %>%
  filter(`match`) %>%
  ungroup()
# No match
anra_no_match <- eggnog_anra_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_anra, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb1c <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb1d <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2c <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb2d <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      is_true_3 <- length(intersect(sb1c, sb2c)) > 0
      is_true_4 <- length(intersect(sb1d, sb2d)) > 0
      
      any(is_digit & !(is_true_1) & !(is_true_2) & !(is_true_3) & !(is_true_4))
    }
  ) %>%
  filter(`match`) %>%
  ungroup()

mara_no_ec <- eggnog_mara_anno %>% 
  filter(is.na(EC_mara))
# If any value matches 100% i.e 1.2.3.4,2.3.4.5 and 2.3.4.5 would be a match.
mara_match_100 <- eggnog_mara_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Loop-like action on each value
      e1 <- eg_ec
      stre1 <- strsplit(e1, ",")[[1]]
      e2 <- EC_mara
      stre2 <- strsplit(e2, ",")[[1]]
      any(stre1 %in% stre2)
    }
  ) %>%
  filter(`match`)

# Match the first digit and only count if there are no matches. Cases such as 
# 1.2.3.4,2.3.4.5 compared with 2.1.1.1 should still count as a match under this 
# logic and thus not be included here.
mara_match_1_not_2 <- eggnog_mara_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_mara, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      any(is_digit & is_true_1 & !(is_true_2))
    }
  ) %>%
  filter(`match`)
  
# Now if it matches the first two but not the third one
mara_match_12_not_3 <- eggnog_mara_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_mara, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb1c <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2c <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      is_true_3 <- length(intersect(sb1c, sb2c)) > 0
      
      any(is_digit & is_true_1 & is_true_2 & !(is_true_3))
    }
  ) %>%
  filter(`match`) %>%
  ungroup()

# Now if it matches the first three but not the fourth one
mara_match_123_not_4 <- eggnog_mara_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_mara, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb1c <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb1d <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2c <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb2d <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      is_true_3 <- length(intersect(sb1c, sb2c)) > 0
      is_true_4 <- length(intersect(sb1d, sb2d)) > 0
      
      any(is_digit & is_true_1 & is_true_2 & is_true_3 & !(is_true_4))
    }
  ) %>%
  filter(`match`) %>%
  ungroup()
# No match
mara_no_match <- eggnog_mara_anno %>%
  rowwise() %>%
  mutate(
    `match` = {
      # Split both strings into vectors
      stre1 <- strsplit(eg_ec, ",")[[1]]
      stre2 <- strsplit(EC_mara, ",")[[1]]
      
      # Extract the first and second components of EC numbers
      sb1a <- sapply(stre1, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x)) 
      sb1b <- sapply(stre1, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb1c <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb1d <- sapply(stre1, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      sb2a <- sapply(stre2, function(x) sub("(\\d+)\\.\\d+.*", "\\1", x))
      sb2b <- sapply(stre2, function(x) sub("\\d+\\.(\\d+).*", "\\1", x))
      sb2c <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\.(\\d+).*", "\\1", x))
      sb2d <- sapply(stre2, function(x) sub("\\d+\\.\\d+\\..\\d+\\.(\\d+).*", "\\1", x))
      
      is_digit <- grepl("^\\d+$", sb2a[[1]])
      is_true_2 <- length(intersect(sb1b, sb2b)) > 0
      is_true_1 <- length(intersect(sb1a, sb2a)) > 0
      is_true_3 <- length(intersect(sb1c, sb2c)) > 0
      is_true_4 <- length(intersect(sb1d, sb2d)) > 0
      
      any(is_digit & !(is_true_1) & !(is_true_2) & !(is_true_3) & !(is_true_4))
    }
  ) %>%
  filter(`match`) %>%
  ungroup()

print(paste("mara only (no kora) validated with eggnog:", mara_stat))
print(paste("anra only (no kora) validated with eggnog:", anra_stat))

# Calculate the average for each increase or decrease
print(paste("anra increase from annr:", sum(stats$anra_increase_annr)/nrow(stats)))
print(paste("anra min increase from annr:", min(stats$anra_increase_annr)))
print(paste("anra max increase from annr:", max(stats$anra_increase_annr)))

# Calculate the refined difference
print(paste("kore increase from annr:", sum(stats$anre_increase_annr)/nrow(stats)))
print(paste("kore increase from kora:", sum(stats$kore_increase_kora)/nrow(stats)))
print(paste("kore increase from kora valid:", sum(stats$kore_increase_kora_valid)/nrow(stats)))

# Calculate the number of KO found per tool accross all genomes (not unique set)
# that are validated with eggnog
print(paste("total number of ORFs:", sum(stats$num_ORFs)))
print(paste("total number of KO found with kofamscan:", sum(stats$kora_valid)))
print(paste("Percent KO annotated:", sum(stats$kora)/sum(stats$num_ORFs)))
print(paste("Percent KO validated:", sum(stats$kora_valid)/sum(stats$kora)))
print(paste("total number of KO found with microbeannotator:", sum(stats$mara_valid)))
print(paste("Percent KO annotated:", sum(stats$mara)/sum(stats$num_ORFs)))
print(paste("Percent KO validated:", sum(stats$mara_valid)/sum(stats$mara)))
print(paste("total number of KO found with anvio:", sum(stats$anra_valid)))
print(paste("Percent KO annotated:", sum(stats$anra)/sum(stats$num_ORFs)))
print(paste("Percent KO validated:", sum(stats$anra_valid)/sum(stats$anra)))

# Calculate the percent increase or decrease between tools
print(paste("Percent KO annotated in MicrobeAnnotator over Kofamscan:", sum((stats$mara-stats$kora)/stats$kora)/nrow(stats)))
print(paste("Percent KO annotated in MicrobeAnnotator over anvio:", sum((stats$mara-stats$anra)/stats$anra)/nrow(stats)))
print(paste("Percent KO annotated in microbeannotator over refined", sum((stats$mare-stats$mara)/stats$mara)/nrow(stats)))
print(paste("Percent KO annotated in anvi'o over Kofamscan:", sum((stats$anra-stats$kora)/stats$kora)/nrow(stats)))
print(paste("Percent KO annotated in anvi'o over microbeannotator:", sum((stats$anra-stats$mara)/stats$mara)/nrow(stats)))
print(paste("Percent KO annotated in kofamscan over microbeannotator:", sum((stats$kora-stats$mara)/stats$mara)/nrow(stats)))
print(paste("Percent KO annotated over kofamscan over anvi'o", sum((stats$kora-stats$anra)/stats$anra)/nrow(stats)))

# Calculate the percent increase or decrease between tools
print(paste("Percent valid KO annotated in MicrobeAnnotator over Kofamscan:", sum((stats$mara_valid-stats$kora_valid)/stats$kora_valid)/nrow(stats)))
print(paste("Percent valid KO annotated in MicrobeAnnotator over anvio:", sum((stats$mara_valid-stats$anra_valid)/stats$anra_valid)/nrow(stats)))
print(paste("Percent valid KO annotated in microbeannotator over refined", sum((stats$mare_valid-stats$mara_valid)/stats$mara_valid)/nrow(stats)))
print(paste("Percent valid KO annotated in anvi'o over Kofamscan:", sum((stats$anra_valid-stats$kora_valid)/stats$kora_valid)/nrow(stats)))
print(paste("Percent valid KO annotated in anvi'o over microbeannotator:", sum((stats$anra_valid-stats$mara_valid)/stats$mara_valid)/nrow(stats)))
print(paste("Percent valid KO annotated in kofamscan over microbeannotator:", sum((stats$kora_valid-stats$mara_valid)/stats$mara_valid)/nrow(stats)))
print(paste("Percent valid KO annotated in kofamscan over anvi'o", sum((stats$kora_valid-stats$anra_valid)/stats$anra_valid)/nrow(stats)))

# Calculate min and max
print(paste("Min valid KO annotated in anvi'o over microbeannotator:", min((stats$anra_valid-stats$mara_valid)/stats$mara_valid)))
print(paste("Max valid KO annotated in anvi'o over microbeannotator:", max((stats$anra_valid-stats$mara_valid)/stats$mara_valid)))

# Calculate Wilcoxon test
wilcox.test(stats$annr, stats$anra, paired = TRUE, alternative = "two.sided")
wilcox.test(stats$anra_valid, stats$mara_valid, paired = TRUE, alternative = "two.sided")
wilcox.test(stats$anra, stats$anst, paired = TRUE, alternative = "two.sided")

# Calculate percent increase of decrease between anvio nt-ko and default
print(paste("Percent valid KO annotated in stray anvi'o over default:", sum((stats$anst-stats$anra)/stats$anra)/nrow(stats)))

# Calculate average ORFs annotated per genome per tool
print(paste("Average number of ORFs annotated anvio: ", sum(stats$anst_valid/nrow(stats))))
print(paste("Average number of ORFs annotated microbeannotator: ", sum(stats$mara_valid/nrow(stats))))

names(eggnog_anra)[names(eggnog_anra) == 'anra_valid'] <- "anra"
names(eggnog_annr)[names(eggnog_annr) == 'annr_valid'] <- "annr"
names(eggnog_anst)[names(eggnog_anst) == 'anst_valid'] <- "anst"
names(eggnog_mare)[names(eggnog_mare) == 'mare_valid'] <- "mare"
names(eggnog_mara)[names(eggnog_mara) == 'mara_valid'] <- "mara"
names(eggnog_kore)[names(eggnog_kore) == 'kore_valid'] <- "kore"
names(eggnog_kora)[names(eggnog_kora) == 'kora_valid'] <- "kora"

# Plot the number of functions compared to ORFS per genome regardless of duplication
ko_per_tool <- aggregate(kora ~ .id, data = kora_sub[c(".id", "V1", "kora")], FUN = function(x) length(x))
ko_per_tool$V4 <- aggregate(mara ~ .id, data = mara_sub[c(".id", "V1", "mara")], FUN = function(x) length(x))$mara
ko_per_tool$V5 <- aggregate(anra ~ .id, data = anra_sub[c(".id", "V1", "anra")], FUN = function(x) length(x))$anra
ko_per_tool$V6 <- aggregate(kora ~ .id, data = eggnog_kora[c(".id", "V1", "kora")], FUN = function(x) length(x))$kora
ko_per_tool$V7 <- aggregate(mara ~ .id, data = eggnog_mara[c(".id", "V1", "mara")], FUN = function(x) length(x))$mara
ko_per_tool$V8 <- aggregate(anra ~ .id, data = eggnog_anra[c(".id", "V1", "anra")], FUN = function(x) length(x))$anra
colnames(ko_per_tool) <- c("accession","number_hits_kofamscan", "number_hits_microbeannotator", 
                           "number_hits_anvio", "kofamscan_valid", "microbeannotator_valid", "anvio_valid")
# Make the data into long format
ko_per_tool <- merge(ko_per_tool, orfs[c("gtdb_family", "accession", "number_genes")], by="accession")
ko_per_tool_long <- reshape2::melt(ko_per_tool, 
                                   id.vars = c("accession", "gtdb_family", "number_genes", "kofamscan_valid", "microbeannotator_valid", "anvio_valid"), 
                                   measure.vars = c("number_hits_kofamscan", "number_hits_microbeannotator", "number_hits_anvio"),
                                   variable.name = "method", value.name = "number_hits")
# Plot F1a
scatter <- ggplot(ko_per_tool_long, aes(x = number_genes, y = number_hits, color = method, size = ifelse(method == "number_hits_kofamscan", kofamscan_valid, ifelse(method == "number_hits_microbeannotator", microbeannotator_valid, anvio_valid)))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("#e8b437", "#3093CB", "#038575")) +
  labs(x = "Number of ORFs per Species", y = "Number of Orthologs Annotated", color = "Method", size = "Eggnog Matches") +
  ggtitle("Orthologs Found per Method for Species") +
  theme_light() + 
  theme(text=element_text(size=20)) +
  scale_color_manual(values = c("#e8b437", "#3093CB", "#038575"),
                     labels = c("Kofamscan", "MicrobeAnnotator", "anvi'o"))

# Calculate the numbers per genome ---------------------------------------------
# Microbeannotator not in anvi'o and kofamscan
all_tools <- mara_sub
all_tools <- merge(all_tools, kora_sub, by=c(".id","V1"), all=TRUE)
all_tools <- all_tools %>%
  group_by_at(vars(".id", "V1")) %>%
  filter((is.na(kora) & !is.na(mara)) | (!is.na(kora) & check_not_ko(mara, kora))) %>%
  ungroup()
all_tools <- merge(all_tools, anra_sub, all=TRUE)
all_tools <- all_tools %>%
  group_by_at(vars(".id", "V1")) %>%
  filter((is.na(anra) & !is.na(mara)) | (!is.na(anra) & check_not_ko(mara, anra))) %>%
  ungroup()
tool_comparison <- aggregate(mara ~ .id, data = all_tools, FUN = function(x) length(x))

df_mara <- data.frame(
  title = c(
    "Specific Enzyme (Full)",
    "No Match", "No EC", "Enzyme Class (First)",
    "Substrate Group (Second)", "Specific Reaction Type (Third)"),
  value = unlist(list(
    merge(all_tools, mara_match_100, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, mara_no_match, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, mara_no_ec, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, mara_match_1_not_2, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, mara_match_12_not_3, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, mara_match_123_not_4, by=c(".id","V1")) %>% nrow())
  )
)
df_mara <- df_mara %>%
  select(title, value) %>%
  mutate(Method = "MicrobeAnnotator")  # Add a column to indicate Tool1
df_mara <- df_mara %>%
  group_by(Method) %>%
  mutate(percentage = (value / sum(value)) * 100)

# anvi'o not in Microbeannotator and kofamscan
all_tools <- anra_sub
all_tools <- merge(all_tools, kora_sub, by=c(".id","V1"), all=TRUE)
all_tools <- all_tools %>%
  group_by_at(vars(".id", "V1")) %>%
  filter((is.na(kora) & !is.na(anra)) | (!is.na(anra) & check_not_ko(anra, kora))) %>%
  ungroup()
all_tools <- merge(all_tools, mara_sub, all=TRUE)
all_tools <- all_tools %>%
  group_by_at(vars(".id", "V1")) %>%
  filter(is.na(mara) | (!is.na(anra) & check_not_ko(anra, mara))) %>%
  ungroup()
tool_comparison <- merge(tool_comparison, aggregate(anra ~ .id, data = all_tools, FUN = function(x) length(x)), by=".id")

df_anra <- data.frame(
  title = c(
    "Specific Enzyme (Full)",
    "No Match", "No EC", "Enzyme Class (First)",
    "Substrate Group (Second)", "Specific Reaction Type (Third)"),
  value = unlist(list(
    merge(all_tools, anra_match_100, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, anra_no_match, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, anra_no_ec, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, anra_match_1_not_2, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, anra_match_12_not_3, by=c(".id","V1")) %>% nrow(),
    merge(all_tools, anra_match_123_not_4, by=c(".id","V1")) %>% nrow())
  )
)
df_anra <- df_anra %>%
  select(title, value) %>%
  mutate(Method = "anvi'o")  # Add a column to indicate Tool1
df_anra <- df_anra %>%
  group_by(Method) %>%
  mutate(percentage = (value / sum(value)) * 100)

# Combine the datasets
combined <- bind_rows(df_mara, df_anra)
eg1 <- ggplot(combined, aes(x = title, y = percentage, fill = Method)) + 
  geom_col(position = position_dodge(width = 0.9)) +  
  scale_fill_manual(values = c("#038575", "#3093CB")) + 
  labs(
    title = "Unique KO per Tool Comparison by EC", 
    x = "EC Number Match", 
    y = "Percentage of Unique KOs"
  ) + 
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percentage, value)),  
            position = position_dodge(width = 0.8),  
            vjust = -0.5, size = 3) +  # Add percentage and value as label
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.text.x = element_text(angle = 45, hjust = 1)  # Ensure the x-axis labels are rotated
  )
ggsave("eggnog_ec_match.svg", plot = eg1, device = "svg", width = 20, height = 20)

# Perform binomial tests
filtered <- combined %>%
  select(title, Method, value) %>%  # Selecting the necessary columns
  pivot_wider(
    names_from = Method, 
    values_from = value
  )
fisher_test_results <- filtered %>%
  rowwise() %>%
  mutate(
    fisher_p_value = ifelse(
      MicrobeAnnotator == 0 | `anvi'o` == 0, 
      NA,  # Skip rows with zero counts for either method
      fisher.test(matrix(c(MicrobeAnnotator, sum(filtered$MicrobeAnnotator) - MicrobeAnnotator, 
                           `anvi'o`, sum(filtered$`anvi'o`) - `anvi'o`), nrow = 2))$p.value
    )
  ) %>%
  ungroup()
z_test_results <- filtered %>%
  rowwise() %>%
  mutate(
    z_test_p_value = prop.test(
      c(MicrobeAnnotator, `anvi'o`), 
      c(sum(filtered$MicrobeAnnotator), sum(filtered$`anvi'o`))
    )$p.value
  ) %>%
  ungroup()
write.csv(fisher_test_results, "fisher_test_results_ec_matches.csv", row.names = FALSE)
write.csv(z_test_results, "z_test_results_ec_matches.csv", row.names = FALSE)

# Shared Microbeannotator and anvi'o not in kofamscan
all_tools <- mara_sub
all_tools <- merge(all_tools, anra_sub, by=c(".id","V1"), all=TRUE)
all_tools <- all_tools %>%
  group_by_at(vars(".id", "V1")) %>%
  filter(check_ko(mara, anra)) %>%
  ungroup()
all_tools <- merge(all_tools, kora_sub, all=TRUE)
all_tools <- all_tools %>%
  group_by_at(vars(".id", "V1")) %>%
  filter(check_not_ko(mara, kora)) %>%
  ungroup()

tool_comparison <- merge(tool_comparison, aggregate(mara ~ .id, data = all_tools, FUN = function(x) length(x)), by=".id", all=TRUE)
tool_comparison[is.na(tool_comparison)] <- 0

colnames(tool_comparison) <- c(".id", "Only MicrobeAnnotator", "Only anvi'o", "MicrobeAnnotator + anvi'o")
tool_comparison <- merge(tool_comparison, family_linker)

# Add back in family
tool_comparison <- tool_comparison %>%
  mutate(environment = if_else(family %in% c("Rhizobiaceae", "Streptomycetaceae"), "Soil",
                               if_else(family %in% c("Pelagibacteraceae", "Cyanobiaceae"), "Open Ocean",
                                       if_else(family %in% c("Rhodobacteraceae", "Microcystaceae"), "Coastal",
                                               if_else(family %in% c("Lachnospiraceae", "Enterobacteriaceae", "Bacteroidaceae"), "Gut",
                                                       if_else(family %in% c("Nanosynbacteraceae", "Streptococcaceae"), "Oral", NA_character_))))))

# Make longer pivot
tool_comparison_long <- tool_comparison %>%
  pivot_longer(cols = c(`Only MicrobeAnnotator`, `MicrobeAnnotator + anvi'o`, `Only anvi'o`),
               names_to = "Category",
               values_to = "Count")
tool_comparison_long$Category <- factor(tool_comparison_long$Category, 
                                        levels=c("Only MicrobeAnnotator", 
                                                 "MicrobeAnnotator + anvi'o", 
                                                 "Only anvi'o"))

# Reorder levels of the family factor based on counts (largest to smallest)
family_counts <- table(tool_comparison_long$family)
tool_comparison_long$family <- factor(tool_comparison_long$family, 
                                      levels = names(sort(family_counts, decreasing = TRUE)))
tool_comparison_long <- tool_comparison_long[tool_comparison_long$Category != "MicrobeAnnotator + anvi'o", ]

# Plot the families and orthologs found for F1b.
ortholog_cnt <- ggplot(tool_comparison_long, aes(x = family, y = Count, fill = Category)) +
  geom_boxplot(width=1) +  # Changed to geom_boxplot
  facet_wrap(~ environment, scales = "free_x", nrow = 1, strip.position = "top", switch = "x") +
  labs(x = "Family", y = "Number of Orthologs Annotated", fill = "Method Set") +
  ggtitle("Orthologs Found for Family per Environment") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text=element_text(size=20)) +
  scale_fill_manual(values = c("#3093CB", "#038575", "#05d0b7", "#e8b437"),
                    labels = c("MicrobeAnnotator (No anvi'o)", "anvi'o (no Kofamscan)", "anvi'o (no Kofamscan)", "Kofamscan (no anvi'o)"))

# Plot F1c
strays_anst <- merge(anst, stray_kofam_hmm, by = "V3")
strays_mare <- merge(mare, stray_kofam_hmm, by="V3")
# to plot values that are exactly the same.
strays_anst_unq <- strays_anst %>% distinct(V5, .keep_all = TRUE)
not_standard_mare_unq <- unique(subset(mare, as.numeric(V6) > 0.01)) %>% distinct(V6, .keep_all = TRUE)
strays_mare_unq <- strays_mare %>% distinct(V6, .keep_all = TRUE)
anst_unq <- anst %>% distinct(V5, .keep_all = TRUE)
anra_unq <- anra %>% distinct(V5, .keep_all = TRUE)
mare_unq <- mare %>% distinct(V6, .keep_all = TRUE)
# Generate evals dataframes to compare strays
colnames(anst_unq)[colnames(anst_unq) == "V5"] <- "anvi'o (Heuristic + nt-KO)"
colnames(anra_unq)[colnames(anra_unq) == "V5"] <- "anvi'o (Heuristic)"
colnames(mare_unq)[colnames(mare_unq) == "V6"] <- "MicrobeAnnotator"
colnames(strays_anst_unq)[colnames(strays_anst_unq) == "V5"] <- "anvi'o (Heuristic + nt-KO)"
colnames(strays_mare_unq)[colnames(strays_mare_unq) == "V6"] <- "MicrobeAnnotator"
colnames(not_standard_mare_unq)[colnames(not_standard_mare_unq) == "V6"] <- "MicrobeAnnotator"

mare_unq$MicrobeAnnotator <- as.numeric(mare_unq$MicrobeAnnotator)
strays_mare_unq$MicrobeAnnotator <- as.numeric(strays_mare_unq$MicrobeAnnotator)
not_standard_mare_unq$MicrobeAnnotator <- as.numeric(not_standard_mare_unq$MicrobeAnnotator)

# Create the value columns
evals <- bind_rows(dplyr::select(anst_unq, `anvi'o (Heuristic + nt-KO)`, `anvi'o (Heuristic + nt-KO)`),
                   dplyr::select(anra_unq, `anvi'o (Heuristic)`, `anvi'o (Heuristic)`),
                   dplyr::select(mare_unq, MicrobeAnnotator, MicrobeAnnotator),
                   .id = "Source") %>% 
  pivot_longer(cols = c(`anvi'o (Heuristic + nt-KO)`, `anvi'o (Heuristic)`, MicrobeAnnotator)) %>%
  drop_na()
evals_stray <- bind_rows(dplyr::select(strays_anst_unq, `anvi'o (Heuristic + nt-KO)`, `anvi'o (Heuristic + nt-KO)`),
                         dplyr::select(not_standard_mare_unq, MicrobeAnnotator, MicrobeAnnotator), 
                         .id = "Source") %>% 
  pivot_longer(cols = c(`anvi'o (Heuristic + nt-KO)`, MicrobeAnnotator)) %>%
  drop_na()

# Define constants and colors
constant <- min(evals$value[evals$value > 0], na.rm = TRUE)
threshold <- 0.01
threshold_log <- log10(0.01 + constant)
custom_colors <- c(above_threshold = "red", below_threshold = "gray")

# Transform the values with log10 and add columns for graphing
evals_log <- evals %>%
  mutate(
    value = log10(value + constant),
    color_group = ifelse(value >= threshold_log, "above_threshold", "below_threshold"),
    type = "Other")
evals_stray_log <- evals_stray %>%
  mutate(
    value = log10(value + constant),
    color_group = ifelse(value >= threshold_log, "above_threshold", "below_threshold"),
    type = "nt-KO")
evals <- evals %>%
  mutate(
    value = value,
    color_group = ifelse(value >= threshold, "above_threshold", "below_threshold"),
    type = "Other")
evals_stray <- evals_stray %>%
  mutate(
    value = value,
    color_group = ifelse(value >= threshold, "above_threshold", "below_threshold"),
    type = "nt-KO")

# Make a combined dataset
evals_filtered_log <- evals_log %>%
  anti_join(evals_stray_log, by = c("Source", "name", "value"))
evals_combined_log <- bind_rows(evals_filtered_log, evals_stray_log)
evals_combined_log$color_group <- ifelse(evals_combined_log$type == "nt-KO", "nt-KO", 
                                     ifelse(evals_combined_log$value > threshold, "above_threshold", "below_threshold"))
evals_filtered <- evals %>%
  anti_join(evals_stray, by = c("Source", "name", "value"))
evals_combined <- bind_rows(evals_filtered, evals_stray)
evals_combined$color_group <- ifelse(evals_combined$type == "nt-KO", "nt-KO", 
                                     ifelse(evals_combined$value > threshold, "above_threshold", "below_threshold"))

if(signif(nrow(evals_combined_log %>% filter(color_group == "above_threshold") %>% filter(name == "anvi'o (Heuristic)"))/nrow(evals_combined_log %>% filter(name == "anvi'o (Heuristic)")),1) != 0){
  stop("anvio values above threshold I")
}
if(signif(nrow(evals_combined_log %>% filter(color_group == "above_threshold") %>% filter(name == "anvi'o (Heuristic + nt-KO)"))/nrow(evals_combined_log %>% filter(name == "anvi'o (Heuristic + nt-KO)")),1) != 0){
  stop("anvio values above threshold II")
}
if(signif(nrow(evals_combined %>% filter(color_group == "above_threshold") %>% filter(name == "anvi'o (Heuristic)"))/nrow(evals_combined %>% filter(name == "anvi'o (Heuristic)")),1) != 0){
  stop("anvio values above threshold I")
}
if(signif(nrow(evals_combined %>% filter(color_group == "above_threshold") %>% filter(name == "anvi'o (Heuristic + nt-KO)"))/nrow(evals_combined %>% filter(name == "anvi'o (Heuristic + nt-KO)")),1) != 0){
  stop("anvio values above threshold II")
}

# This is only set with the paste0("0.0") because it passes the check above.
above_threshold_ma <- paste0((signif(nrow(evals_combined %>% filter(color_group == "above_threshold") %>% filter(name == "MicrobeAnnotator"))/nrow(evals_combined %>% filter(name == "MicrobeAnnotator")),1)*100), "%")
above_threshold_an <- paste0((signif(nrow(evals_combined %>% filter(color_group == "above_threshold") %>% filter(name == "anvi'o (Heuristic)"))/nrow(evals_combined %>% filter(name == "anvi'o (Heuristic)")),1)*100), "%")
above_threshold_anst <- paste0((signif(nrow(evals_combined %>% filter(color_group == "above_threshold") %>% filter(name == "anvi'o (Heuristic + nt-KO)"))/nrow(evals_combined %>% filter(name == "anvi'o (Heuristic + nt-KO)")),1)*100), "%")

p1 <- ggplot(evals_combined, aes(x = name, y = value, color = color_group, shape = type)) + 
  geom_jitter(size = 2, width = 0.2) +   
  scale_shape_manual(values = c("Other" = 16, "nt-KO" = 15)) +   
  scale_color_manual(values = c("above_threshold" = "red", 
                                "below_threshold" = "gray", 
                                "nt-KO" = "black"),
                     labels = c("above_threshold" = "Above Threshold", 
                                "below_threshold" = "Below Threshold", 
                                "nt-KO" = "nt-KO")) +  
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue") + 
  theme_light() +  
  ggtitle("nt-KO Recovery Methods E-Values") +
  theme(
    text = element_text(size = 20),  
    axis.text.x = element_text(angle = 0)  
  ) + 
  labs(
    x = "Category",  
    y = "E-Value",  
    color = "Threshold",  
    shape = "Type"
  ) 

p2 <- ggplot(evals_combined_log, aes(x = name, y = value, color = color_group, shape = type)) + 
  geom_jitter(size = 2, width = 0.2) +   
  scale_shape_manual(values = c("Other" = 16, "nt-KO" = 15)) +   
  scale_color_manual(values = c("above_threshold" = "red", 
                                "below_threshold" = "gray", 
                                "nt-KO" = "black"),
                     labels = c("above_threshold" = "Above Threshold", 
                                "below_threshold" = "Below Threshold", 
                                "nt-KO" = "nt-KO")) +  
  geom_hline(yintercept = threshold_log, linetype = "dashed", color = "blue") + 
  theme_light() +  
  theme(
    text = element_text(size = 20),  
    axis.text.x = element_text(angle = 0),  
  ) + 
  labs(
    x = "Category",  
    y = "log10(E-Value)",  
    color = "Threshold",  
    shape = "Type"
  ) + 
  coord_cartesian(ylim = c(threshold_log - 1.5, threshold_log + 3.5)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(color = "black", fill = NA, size = 1, linetype = "dotted"), # Dotted border around entire plot
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    text = element_text(size = 18)
  )

bitscores <- p1 +
  annotation_custom(ggplotGrob(p2), xmin = 0.5, xmax = 2.7, 
                    ymin = 2.5, ymax = 10.45) +
  geom_path(aes(x, y, group = grp), 
            data=data.frame(x = c(1,1.1,2,1.7,2.8,2.3), y=c(0.3,2.5,0.3,2.5,0.3,2.5),grp=c(1,1,2,2,3,3)),
            linetype = 'dashed', 
            inherit.aes = FALSE)
  

# plot in patchwork 
combined <- (scatter / bitscores / ortholog_cnt) + plot_layout(height = c(4, 4, 4)) + plot_annotation(tag_levels = 'A')
ggsave("f1abc.svg", plot = combined, device = "svg", width = 20, height = 20)

pdf("f1abc.pdf", width=20, height=20)
combined
dev.off()
