library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(plyr)
library(tidyverse)


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
    en <- en[c(".id","V1","V3","V12")]
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
    group_by(.id, V2) %>% 
    slice_min(V5)
  return(an)
}

clean_ko <- function(ko){
  ko <- ko[!str_detect(ko$V2, "KEGG_BRITE|KEGG_Class|KEGG_Module"), ]
  ko$.id <- sub("..$", "", ko$.id)
  ko <- ko[ko[,2] == "*", ]
  ko <- ko %>% 
    group_by(.id, V2) %>% 
    slice_min(V6)
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


# Set 11 families we are interested in
family <- c("Bacteroidaceae","Nanosynbacteraceae", "Streptococcaceae",
            "Lachnospiraceae", "Enterobacteriaceae", "Microcystaceae",
            "Rhodobacteraceae", "Cyanobiaceae", "Pelagibacteraceae",
            "Streptomycetaceae", "Rhizobiaceae")

# Read in metadata files needed to generate graphs
family_linker <- make_linker("/Users/kananen/Desktop/ImHere/bac120_metadata_r214.tsv")
stray_kofam_hmm <- make_stray_kofam("/Users/kananen/Desktop/ImHere/hmm_profiles_with_kofams_with_no_threshold.hmm")

# Pause here to see who has e-values above 0.01 and run them in eggnog mapper
# if they do. Read in all the files in a given directory
eggnog <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/eggnog/", "eggnog")
# Calculated ORFS from worflow in table 1 based off of the number
orfs <- read.csv("/Users/kananen/Documents/Upset/upset_R/f1_table.tsv", sep='\t')
orfs$accession <- sub("..$", "", orfs$accession)

# Read in Kofamscan raw and refined files
kora <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/kofamscan/functions/default/", "kofamscan")
kore <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/kofamscan/functions/refined/", "kofamscan")
# Read in MicrobeAnnotator raw and refined files
mara <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/microbeAnnotator/functions/default/", "microbeannotator")
mare <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/microbeAnnotator/functions/refined/", "microbeannotator")
# Read in Anvio raw, stray, and no recovery files
anra <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/anvio/functions/default/", "anvio")
anst <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/anvio/functions/stray/", "anvio")
annr <- read_in_files("/Users/kananen/Desktop/ImHere/unfiltered/anvio/functions/no_hueristic/", "anvio")

# Clean up anvi'o and kofamscan output
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
# Make one dataframe we will work with for simplicity
all <- merge(merge(merge(merge(merge(merge(kora_sub, mara_sub, by=c(".id", "V1"), all=TRUE), 
                                     anra_sub, by=c(".id", "V1"), all=TRUE), 
                               kore_sub, by=c(".id", "V1"), all=TRUE), 
                         mare_sub, by=c(".id", "V1"), all=TRUE),
                   annr_sub, by=c(".id", "V1"), all=TRUE),
             anst_sub, by=c(".id", "V1"), all=TRUE)

# Calculate the eggnog hits that match
# Kofamscan default
eggnog_kora <- unique(merge(all[c(".id", "V1", "kora")], eggnog, by=c(".id", "V1")))
eggnog_kora <- eggnog_kora[str_detect(eggnog_kora$V12, eggnog_kora$kora), ]
# Kofamscan refined
eggnog_kore <- unique(merge(all[c(".id", "V1", "kore")], eggnog, by=c(".id", "V1")))
eggnog_kore <- eggnog_kore[str_detect(eggnog_kore$V12, eggnog_kore$kore), ]
# MicrobeAnnotator default
eggnog_mara <- unique(merge(all[c(".id", "V1", "mara")], eggnog, by=c(".id", "V1")))
eggnog_mara <- eggnog_mara[str_detect(eggnog_mara$V12, eggnog_mara$mara), ]
# MicrobeAnnotator refined
eggnog_mare <- unique(merge(all[c(".id", "V1", "mare")], eggnog, by=c(".id", "V1")))
eggnog_mare <- eggnog_mare[str_detect(eggnog_mare$V12, eggnog_mare$mare), ]
# Anvio default
eggnog_anra <- unique(merge(all[c(".id", "V1", "anra")], eggnog, by=c(".id", "V1")))
eggnog_anra <- eggnog_anra[str_detect(eggnog_anra$V12, eggnog_anra$anra), ]
# Anvio no Hueristic
eggnog_annr <- unique(merge(all[c(".id", "V1", "annr")], eggnog, by=c(".id", "V1")))
eggnog_annr <- eggnog_annr[str_detect(eggnog_annr$V12, eggnog_annr$annr), ]
# Anvio stray
eggnog_anst <- unique(merge(all[c(".id", "V1", "anst")], eggnog, by=c(".id", "V1")))
eggnog_anst <- eggnog_anst[str_detect(eggnog_anst$V12, eggnog_anst$anst), ]

eggnog_kora <- eggnog_kora[complete.cases(eggnog_kora$.id), ]
eggnog_kore <- eggnog_kore[complete.cases(eggnog_kore$.id), ]
eggnog_mara <- eggnog_mara[complete.cases(eggnog_mara$.id), ]
eggnog_mare <- eggnog_mare[complete.cases(eggnog_mare$.id), ]
eggnog_anra <- eggnog_anra[complete.cases(eggnog_anra$.id), ]
eggnog_annr <- eggnog_annr[complete.cases(eggnog_annr$.id), ]
eggnog_anst <- eggnog_anst[complete.cases(eggnog_anst$.id), ]

# Calculate statistics needed for paper
stats <- orfs[c("accession", "number_genes")]; colnames(stats) <- c(".id", "num_ORFs")
# Percentage change for kofamscan refined and kofamscan default
stats$kore_increase_kora <- ((aggregate(kore ~ .id, data = all[c(".id", "V1", "kore")][complete.cases(all[c(".id", "V1", "kore")]),], FUN = function(x) length(x))$kore -
                               aggregate(kora ~ .id, data = all[c(".id", "V1", "kora")][complete.cases(all[c(".id", "V1", "kora")]),], FUN = function(x) length(x))$kora) /
                             aggregate(kora ~ .id, data = all[c(".id", "V1", "kora")][complete.cases(all[c(".id", "V1", "kora")]),], FUN = function(x) length(x))$kora)
stats$kore_increase_kora_valid <- ((aggregate(kore ~ .id, data = eggnog_kore[c(".id", "V1", "kore")][complete.cases(eggnog_kore[c(".id", "V1", "kore")]),], FUN = function(x) length(x))$kore - 
                                      aggregate(kora ~ .id, data = eggnog_kora[c(".id", "V1", "kora")][complete.cases(eggnog_kora[c(".id", "V1", "kora")]),], FUN = function(x) length(x))$kora) / 
                                     aggregate(kora ~ .id, data = eggnog_kora[c(".id", "V1", "kora")][complete.cases(eggnog_kora[c(".id", "V1", "kora")]),], FUN = function(x) length(x))$kora)
# Percentage change for MicrobeAnnotator refined and MicrobeAnnotator default
stats$mare_increase_mara <- ((aggregate(mare ~ .id, data = all[c(".id", "V1", "mare")][complete.cases(all[c(".id", "V1", "mare")]),], FUN = function(x) length(x))$mare -
                                    aggregate(mara ~ .id, data = all[c(".id", "V1", "mara")][complete.cases(all[c(".id", "V1", "mara")]),], FUN = function(x) length(x))$mara) /
                                   aggregate(mara ~ .id, data = all[c(".id", "V1", "mara")][complete.cases(all[c(".id", "V1", "mara")]),], FUN = function(x) length(x))$mara)
stats$mare_increase_mara_valid <- ((aggregate(mare ~ .id, data = eggnog_mare[c(".id", "V1", "mare")][complete.cases(eggnog_mare[c(".id", "V1", "mare")]),], FUN = function(x) length(x))$mare - 
                                      aggregate(mara ~ .id, data = eggnog_mara[c(".id", "V1", "mara")][complete.cases(eggnog_mara[c(".id", "V1", "mara")]),], FUN = function(x) length(x))$mara) /
                                     aggregate(mara ~ .id, data = eggnog_mara[c(".id", "V1", "mara")][complete.cases(eggnog_mara[c(".id", "V1", "mara")]),], FUN = function(x) length(x))$mara)
# Percentage change for anvi'o default and anvi'o no Hueristic
stats$anra_increase_annr <- ((aggregate(anra ~ .id, data = all[c(".id", "V1", "anra")][complete.cases(all[c(".id", "V1", "anra")]),], FUN = function(x) length(x))$anra -
                                aggregate(annr ~ .id, data = all[c(".id", "V1", "annr")][complete.cases(all[c(".id", "V1", "annr")]),], FUN = function(x) length(x))$annr) /
                               aggregate(annr ~ .id, data = all[c(".id", "V1", "annr")][complete.cases(all[c(".id", "V1", "annr")]),], FUN = function(x) length(x))$annr)
stats$anra_increase_annr_valid <- ((aggregate(anra ~ .id, data = eggnog_anra[c(".id", "V1", "anra")][complete.cases(eggnog_anra[c(".id", "V1", "anra")]),], FUN = function(x) length(x))$anra - 
                                      aggregate(annr ~ .id, data = eggnog_annr[c(".id", "V1", "annr")][complete.cases(eggnog_annr[c(".id", "V1", "annr")]),], FUN = function(x) length(x))$annr) /
                                     aggregate(annr ~ .id, data = eggnog_annr[c(".id", "V1", "annr")][complete.cases(eggnog_annr[c(".id", "V1", "annr")]),], FUN = function(x) length(x))$annr)
# Percentage change for anvi'o stray and anvi'o default
stats$anst_increase_anra <- ((aggregate(anst ~ .id, data = all[c(".id", "V1", "anst")][complete.cases(all[c(".id", "V1", "anst")]),], FUN = function(x) length(x))$anst -
                                aggregate(anra ~ .id, data = all[c(".id", "V1", "anra")][complete.cases(all[c(".id", "V1", "anra")]),], FUN = function(x) length(x))$anra) /
                               aggregate(anra ~ .id, data = all[c(".id", "V1", "anra")][complete.cases(all[c(".id", "V1", "anra")]),], FUN = function(x) length(x))$anra)
stats$anst_increase_anra_valid <- ((aggregate(anst ~ .id, data = eggnog_anst[c(".id", "V1", "anst")][complete.cases(eggnog_anst[c(".id", "V1", "anst")]),], FUN = function(x) length(x))$anst - 
                                      aggregate(anra ~ .id, data = eggnog_anra[c(".id", "V1", "anra")][complete.cases(eggnog_anra[c(".id", "V1", "anra")]),], FUN = function(x) length(x))$anra) /
                                     aggregate(anra ~ .id, data = eggnog_anra[c(".id", "V1", "anra")][complete.cases(eggnog_anra[c(".id", "V1", "anra")]),], FUN = function(x) length(x))$anra)

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
# Plot F1b
scatter <- ggplot(ko_per_tool_long, aes(x = number_genes, y = number_hits, color = method, size = ifelse(method == "number_hits_kofamscan", kofamscan_valid, ifelse(method == "number_hits_microbeannotator", microbeannotator_valid, anvio_valid)))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("#e8b437", "#3093CB", "#038575")) +
  labs(x = "Number of ORFs per Species", y = "Number of Orthologs Annotated", color = "Method", size = "Eggnog Matches") +
  ggtitle("Orthologs Found per Method for Species") +
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
if (nrow(all_tools) != 0) {
  tool_comparison <- merge(tool_comparison, aggregate(mara ~ .id, data = all_tools, FUN = function(x) length(x)), by=".id")
} else {
  tool_comparison$mara_anra <- 0
}
colnames(tool_comparison) <- c(".id", "Only MicrobeAnnotator", "Only anvi'o", "MicrobeAnnotator + anvi'o")

# Add back in family
tool_comparison <- tool_comparison %>%
  mutate(environment = if_else(family %in% c("Rhizobiaceae", "Streptomycetaceae"), "Soil",
                               if_else(family %in% c("Pelagibacteraceae", "Cyanobiaceae"), "Open Ocean",
                                       if_else(family %in% c("Rhodobacteraceae", "Microcystaceae"), "Coastal",
                                               if_else(family %in% c("Lachnospiraceae", "Enterobacteriaceae", "Bacteroidaceae"), "Gut",
                                                       if_else(family %in% c("Nanosynbacteraceae", "Streptococcaceae"), "Oral", NA_character_))))))
tool_comparison <- merge(tool_comparison, family_linker)

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

# Plot the families and orthologs found for F1b.
ortholog_cnt <- ggplot(tool_comparison_long, aes(x = family, y = Count, fill = Category)) +
  geom_boxplot(width=1) +  # Changed to geom_boxplot
  facet_wrap(~ environment, scales = "free_x", nrow = 1, strip.position = "top", switch = "x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Family", y = "Number of Orthologs Predicted log10(y+1)", fill = "Method Set") +
  ggtitle("Orthologs Found for Family per Environment") +
  scale_fill_manual(values = c("#3093CB", "#038575", "#05d0b7", "#e8b437"),
                    labels = c("MicrobeAnnotator (No anvi'o)", "anvi'o (no MicrobeAnnotator)", "anvi'o (no Kofamscan)", "Kofamscan (no anvi'o)"))

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
colnames(anst_unq)[colnames(anst_unq) == "V5"] <- "anvi'o (Hueristic + nt-KO)"
colnames(anra_unq)[colnames(anra_unq) == "V5"] <- "anvi'o (Hueristic)"
colnames(mare_unq)[colnames(mare_unq) == "V6"] <- "MicrobeAnnotator"
colnames(strays_anst_unq)[colnames(strays_anst_unq) == "V5"] <- "anvi'o (Hueristic + nt-KO)"
colnames(strays_mare_unq)[colnames(strays_mare_unq) == "V6"] <- "MicrobeAnnotator"
colnames(not_standard_mare_unq)[colnames(not_standard_mare_unq) == "V6"] <- "MicrobeAnnotator"

mare_unq$MicrobeAnnotator <- as.numeric(mare_unq$MicrobeAnnotator)
strays_mare_unq$MicrobeAnnotator <- as.numeric(strays_mare_unq$MicrobeAnnotator)
not_standard_mare_unq$MicrobeAnnotator <- as.numeric(not_standard_mare_unq$MicrobeAnnotator)

evals <- bind_rows(dplyr::select(anst_unq, `anvi'o (Hueristic + nt-KO)`, `anvi'o (Hueristic + nt-KO)`),
                   dplyr::select(anra_unq, `anvi'o (Hueristic)`, `anvi'o (Hueristic)`),
                   dplyr::select(mare_unq, MicrobeAnnotator, MicrobeAnnotator),
                   .id = "Source") %>% 
  pivot_longer(cols = c(`anvi'o (Hueristic + nt-KO)`, `anvi'o (Hueristic)`, MicrobeAnnotator)) %>%
  drop_na()

evals_stray <- bind_rows(dplyr::select(strays_anst_unq, `anvi'o (Hueristic + nt-KO)`, `anvi'o (Hueristic + nt-KO)`),
                         dplyr::select(not_standard_mare_unq, MicrobeAnnotator, MicrobeAnnotator), 
                         .id = "Source") %>% 
  pivot_longer(cols = c(`anvi'o (Hueristic + nt-KO)`, MicrobeAnnotator)) %>%
  drop_na()

evals_subset <- evals %>% 
  filter(name %in% c("anvi'o (Hueristic + nt-KO)", "anvi'o (Hueristic)"))
evals_stray_subset <- evals_stray %>% 
  filter(name %in% c("anvi'o (Hueristic + nt-KO)", "anvi'o (Hueristic)"))

# Combine the datasets and add a column to distinguish between them
evals_combined <- rbind(
  mutate(evals_subset, set = "Subset"),
  mutate(evals_stray_subset, set = "nt-KO")
)

# Define custom colors
custom_colors <- c("#2a9d8f", "#2a9d8f", "#6DB9E4", "black")

# Plotting  
p1 <- ggplot() +
  geom_jitter(data = evals, aes(x = name, y = value, color = name),
              position = position_jitter(width = 0.2), alpha = 1) +
  geom_jitter(data = evals_stray, aes(x = name, y = value, color = "nt-KO"),
              position = position_jitter(width = 0.2), alpha = 1) +
  scale_color_manual(values = custom_colors, 
                     labels = c("Anvio", "Anvio", "MicrobeAnnotator", "nt-KO")) + 
  labs(title = "nt-KO Recovery Methods E-Values",
       x = "Category",
       y = "E-Value",
       color = "Method")

# Create the plot
p2 <- ggplot(data = evals_combined, aes(x = name, y = value, color = set)) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6) +
  scale_color_manual(values = c("Subset" = custom_colors[1], "nt-KO" = "black")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(color = "black", fill = NA, size = 1, linetype = "dotted"), # Dotted border around entire plot
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Category",
       y = "E-Value") 

bitscores <- p1 + 
  annotation_custom(ggplotGrob(p2), xmin = .5, xmax = 2.5, 
                    ymin = 2.5, ymax = 10) +
  geom_rect(aes(xmin = .5, xmax = 2.5, ymin = 2.5, ymax = 10), 
            color='black', linetype='dashed', alpha=0) +
  geom_path(aes(x,y,group=grp), 
            data=data.frame(x = c(1,.5,2,2.5), y=c(0.3,2.5,0.3,2.5),grp=c(1,1,2,2)),
            linetype='dashed') +
  geom_hline(yintercept=0.01, linetype="dashed", color = "red")

# plot in patchwork 
(scatter / bitscores / ortholog_cnt) + plot_layout(height = c(4, 4, 4)) + plot_annotation(tag_levels = 'A')

pdf(file="f1abc.pdf", width=15, height=20)
(scatter / bitscores / ortholog_cnt) + plot_layout(height = c(4, 4, 4)) + plot_annotation(tag_levels = 'A')
dev.off()
