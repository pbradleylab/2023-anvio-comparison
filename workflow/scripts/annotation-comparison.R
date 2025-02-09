# NB: Most of these are actually needed for plotting later on

# Note: ko.tsv is required in the current directory to get KOfam names and descriptions (can be downloaded from KEGG).
# The GTDB taxonomy r202 also needs to be downloaded (bac120_taxonomy.tsv.gz).

library(impute)
library(arrow)
library(Biostrings)
library(MSA2dist)
library(ggplot2)
library(ggtree)
library(ggmsa)
library(ggtreeExtra)
library(ggnewscale)
library(ape)
library(phytools)
library(cowplot)
library(patchwork)
library(ggeasy)
library(gggenes)
library(tidyverse)
library(ape)
library(phytools)
map <- purrr::map

### MISC FUNCTIONS/DEFINITIONS

tbl_to_mtx <- function(x) {
  y <- data.matrix(x[,-1])
  rownames(y) <- x[[1]]
  y
}

chart_cols <- c(salmon = "#f9977b",      # 1
                oceanblue = "#4d779e",   # 2
                green = "#007600",       # 3
                orange2 = "#ee9a00",     # 4
                tomato = "#ff6347",      # 1
                steelblue2 = "#4f94cd",  # 2
                lime = "#28e158",        # 3
                peach = "#e3bb65",       # 4
                skyblue = "#18c9ef",     # 5
                violetred = "#e889b4",   # 6
                tan = "#a8552c",         # 7
                puce = "#9e466e",        # 6
                umber = "#6e2d0d",       # 7
                darkgray = "#444444",    # 8
                pewter = "#a4a4a4",      # 8
                lightgray = "#cccccc",   # 8
                black = "#000000")       # 8
colors_in_order <- c("lime", "oceanblue", "peach", "violetred", "skyblue",
                     "salmon", "puce", "orange2", "steelblue2", "tomato",
                     "tan")


### READ DATA

# Read in KO descriptions
ko_tsv <- read_tsv("ko.tsv",
                   col_names = c("accession", "desc")) %>%
  mutate(short_desc = gsub("([^;]+);.*", "\\1", desc))

anvio_dir <- "../../results/annotations/anvio/default/"
ma_dir <- "../../results/annotations/microbeannotator/default/"
kfs_dir <- "../../results//annotations/kofamscan/default/"

kfs_tsvs <- list.files(kfs_dir, "*tsv")

if (file.exists("kfs_fxns_all.parquet")) {
  kfs_fxns_all <- open_dataset("kfs_fxns_all.parquet")
} else {
  # Read Kofamscan data; transiently uses a lot of memory
  # Note that the full Kofamscan data reports bitscores and thresholds for
  # non-significant hits as well. We will use those later to determine the
  # bitscore and threshold for tools that do not automatically report these
  kfs_fxns_all <- open_dataset(sources=file.path(kfs_dir, kfs_tsvs),
                               format = "tsv",
                               schema = schema(best=utf8(),
                                               gene_callers_id=utf8(),
                                               accession=utf8(),
                                               threshold=double(),
                                               bitscore=double(),
                                               e_value=double(),
                                               `function`=utf8()),
                               skip_rows=2) %>%
    mutate(fn = add_filename()) %>%
    compute() %>%
    mutate(genome = gsub(".*(GC.*).tsv", "\\1", fn)) %>%
    select(-fn) %>%
    compute()
  
  write_dataset(kfs_fxns_all, "kfs_fxns_all.parquet")
}
#

kfs_fxns <- kfs_fxns_all %>%
  filter(best == "*") %>%
  collect() %>%
  group_by(genome, gene_callers_id) %>%
  slice_min(e_value, n=1) %>%
  slice_max(bitscore, n=1) %>%
  slice_max(threshold, n=1) %>%
  arrange(accession) %>%
  filter(row_number() == 1) %>%
  ungroup()

all_kfs_thresholds <- kfs_fxns_all %>%
  select(accession, threshold) %>% distinct() %>% collect()

# Read anvi'o results
anvio_tsvs <- list.files(anvio_dir, "*tsv")
anvio_list_fxn <- map(anvio_tsvs, ~ read_tsv(file.path(anvio_dir, .x),
                                             col_types="ccccd"))
anvio_fxns <- map2(anvio_list_fxn, anvio_tsvs, ~ {
  .x %>% 
    filter(source=="KOfam") %>%
    mutate(genome = gsub("\\.tsv", "", .y))
  }, .progress = TRUE) %>%
  Reduce(bind_rows, .) %>%
  as_arrow_table %>% 
  left_join(.,
            select(kfs_fxns_all,
                   gene_callers_id,
                   accession,
                   genome,
                   threshold,
                   bitscore),
            by=c("gene_callers_id", "accession", "genome")) %>%
  collect() %>%
  group_by(genome, gene_callers_id) %>%
  slice_min(e_value, n=1) %>%
  slice_max(bitscore, n=1) %>%
  slice_max(threshold, n=1) %>%
  arrange(accession) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Note: Don't collapse for this one because sometimes the "best" e-value is
# something that does not pass that KOfam's threshold Read MicrobeAnnotator data
ma_subdirs <- list.files(ma_dir, "*", include.dirs = TRUE)
ma_list_fxn <- map(ma_subdirs, ~ 
                      read_tsv(
                        file.path(ma_dir, .x, "kofam_results",
                                  paste0(.x, ".faa.kofam.filt")),
                               col_names=str_split("gene_callers_id	q_orig_accession	accession	t_orig_accession	score_type	e_value	score	bias	E-value_b1d	score_b1d	bias_b1d	exp	reg	clu	ov	env	dom	rep	inc	description_of_target", "\t")[[1]],
                               col_types="cccccddddddddddddddc",
                               skip=3))
ma_fxns <- map2(ma_list_fxn, ma_subdirs, ~ {
  .x %>% # note, should already be filtered for best hit
    mutate(genome = .y)
  }) %>%
  Reduce(bind_rows, .) %>%
  as_arrow_table %>% 
  left_join(.,
            select(kfs_fxns_all,
                   gene_callers_id,
                   accession,
                   genome,
                   threshold,
                   bitscore),
            by=c("gene_callers_id", "accession", "genome")) %>%
  collect() %>%
  group_by(genome, gene_callers_id) %>%
  slice_min(e_value, n=1) %>%
  slice_max(bitscore, n=1) %>%
  slice_max(threshold, n=1) %>%
  arrange(accession) %>%
  filter(row_number() == 1) %>%
  ungroup()
rm(ma_list_fxn)
gc()

# Combine best hits into one table and write to disk
if (file.exists("combined_fxns.csv")) {
  combined_fxns <- read_csv("combined_fxns.csv")
} else {
  combined_fxns <- bind_rows(mutate(kfs_fxns, tool="kfs"),
                             mutate(anvio_fxns, tool="anvio"),
                             mutate(ma_fxns, tool="ma"))
  write_csv(combined_fxns, "combined_fxns.csv")
}


# Add GTDB information (especially family)
genomes <- read_tsv("gtdb_sample_table.tsv")
gtdb <- read_tsv("bac120_taxonomy.tsv.gz", col_names=c("genome", "taxonomy"))
gtdb_align <- gtdb %>%
  mutate(genome = gsub("RS_", "", genome)) %>%
  mutate(genome=gsub("GB_", "", genome))

# we refer to these a lot so this is useful
lachnos <- genomes %>% filter(gtdb_family=="Lachnospiraceae") %>% .[[1]]
enteros <- genomes %>% filter(gtdb_family=="Enterobacteriaceae") %>% .[[1]]

if (file.exists("combined_fxns_fams.csv.xz")) {
  combined_fxns_fams <- read_csv("combined_fxns_fams.csv.xz")
} else {
  combined_fxns_fams <- combined_fxns %>%
    left_join(., genomes, by=c("genome" = "sample_name"))
  write_csv(combined_fxns_fams, "combined_fxns_fams.csv.xz")
}

# Get basic statistics about the number of distinct KOs and the number of total
# annotated genes, either per genome or overall:

total_distinct_KO <- combined_fxns %>%
  group_by(tool) %>%
  select(accession) %>%
  distinct() %>%
  count()

n_distinct_KO <- combined_fxns %>%
  group_by(genome, tool) %>%
  select(accession) %>%
  distinct() %>% 
  summarize(n_annotated = n(), .groups="keep")

total_distinct_annot <- combined_fxns %>%
  group_by(tool) %>%
  select(genome, gene_callers_id) %>%
  distinct() %>%
  count()

n_distinct_annots <- combined_fxns %>%
  group_by(genome, tool) %>%
  select(gene_callers_id) %>%
  distinct() %>%
  summarize(n_annotated = n(), .groups="keep")

total_distinct_KO_per_genome <- n_distinct_KO %>%
  ungroup() %>%
  group_by(tool) %>%
  summarize(sum(n_annotated))

print(total_distinct_annot)
print(total_distinct_KO)
print(total_distinct_KO_per_genome)

nd_summarize <- function(tbl) {
  tbl %>%
    pivot_wider(names_from=tool,values_from=n_annotated) %>%
    mutate(anvio_extra = anvio/kfs - 1, ma_extra = anvio/ma - 1) %>%
    ungroup() %>%
    reframe(anvio_q = quantile(anvio_extra, seq(0,1,0.1)) %>%
              enframe(name="anvio_q", value="anvio_val"),
            ma_q = quantile(ma_extra, seq(0,1,0.1)) %>%
              enframe(name="ma_q", value="ma_val"))
}

nd_annot_summary <- n_distinct_annots %>% nd_summarize
nd_KO_summary <- n_distinct_KO %>% nd_summarize
print(nd_annot_summary)

### Assess differences between families and generate Fig 1c

# Figure 1c

unique_gene_count <- combined_fxns_fams %>%
  group_by(accession, tool, genome, gtdb_family) %>%
  count() %>%
  pivot_wider(names_from=tool, values_from=n, values_fill = 0) %>%
  group_by(genome, gtdb_family) %>%
  count(a=anvio>0, k=kfs>0, m=ma>0)

Fig1c_data <- unique_gene_count %>%
  mutate(comparison = pmap_chr(list(a,k,m), function(a,k,m) {
    ifelse(a & !k & !m, "anvio_only",
           ifelse(!a & !k & m, "ma_only",
                  ifelse(a & !k & m, "both_not_kfs", "other")))
    })) %>%
  mutate(exact_type = paste0(ifelse(a,"A","a"),
                             ifelse(k,"K","k"),
                             ifelse(m,"M","m"))) %>%
  select(-a, -k, -m) %>%
  ungroup()

n_anvio_ma_agree <- Fig1c_data %>%
  filter(comparison=="both_not_kfs") %>%
  ungroup() %>%
  (\(.) sum(.$n))
print(n_anvio_ma_agree) # only 75 genes

Fig1c <- Fig1c_data %>%
  filter(comparison != "other") %>%
  filter(comparison != "both_not_kfs") %>%
  ggplot(aes(x=gtdb_family, y=n, fill=comparison)) + geom_boxplot()

print(Fig1c)
pdf("fig1c.pdf", width=8,height=5); print(Fig1c); dev.off()


# top 20 families with biggest differences in anvi'o or ma

fxn_diffs_all <- combined_fxns_fams %>%
  group_by(gtdb_family, accession, tool) %>%
  count() %>%
  pivot_wider(names_from=tool, values_from=n, values_fill = 0) %>%
  mutate(anvio_diff = anvio-kfs, ma_diff = ma-kfs)

# This version doesn't count multiple annotations to the same KO per genome
fxn_detect_diffs_all <- combined_fxns_fams %>%
  select(gtdb_family, accession, tool, genome) %>%
  distinct() %>%
  group_by(gtdb_family, accession, tool) %>%
  count() %>%
  pivot_wider(names_from=tool, values_from=n, values_fill = 0) 

top_20_anvio_diffs_all <- fxn_diffs_all %>%
  select(gtdb_family, accession, anvio_diff) %>%
  group_by(gtdb_family) %>%
  slice_max(anvio_diff, n=20) %>%
  distinct() 

top_20_ma_diffs_all <- fxn_diffs_all %>%
  select(gtdb_family, accession, ma_diff) %>%
  group_by(gtdb_family) %>%
  slice_max(ma_diff, n=20) %>%
  distinct() 

# now get the ones with the biggest lachno-*specific* differences (i.e. we see
# more annotations in lachnos, but not other taxa)

top_20_anvio_lachnos <- top_20_anvio_diffs_all %>%
  filter(gtdb_family == "Lachnospiraceae") %>%
  .$accession

# NB, this only plots the ones retained by kfs
lachno_plot_kfs <- combined_fxns_fams %>%
  mutate(family = ifelse(gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
  filter(tool=="kfs") %>%
  mutate(bit_norm = map2_dbl(bitscore, threshold, ~ {(.x/.y)})) %>%
  filter(accession %in% filter(top_20_anvio_diffs_all,
                               gtdb_family=="Lachnospiraceae")$accession) %>%
  ggplot(aes(x=bit_norm, fill=family, color=family)) +
  geom_density(bw="sj", adjust=5, alpha=0.5, bounds=c(1,Inf)) +
  facet_wrap(~ accession)  +
  xlim(.75,2) +
  geom_vline(xintercept=1, lty=2) +
  theme_minimal() +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1")
pdf("lachno_distro_plot_kfs.pdf"); lachno_plot_kfs; dev.off()

### Butyrate pathway analysis
kegg_butyrate <- read_tsv("kegg-butyrate-4.txt", col_types="ccccc")
custom_op <- c("K00074", "K17865", "K01715", "K03522", "K03521", "K00248")
kegg_butyrate_found <- inner_join(kegg_butyrate,
                                  filter(fxn_detect_diffs_all,
                                         gtdb_family=="Lachnospiraceae"),
                                  by="accession") %>%
  print(n=nrow(.))
write_csv(kegg_butyrate_found, "kegg_butyrate_found.csv")

### Completeness analysis

completeness <- read_tsv("all_completeness_all_tools.csv")
median_cmp_per_family <- completeness %>%
  group_by(family, module) %>%
  summarize(anvio = median(an_pwc),
            ma = median(ma_pwc),
            kofamscan = median(ko_pwc))
median_cmp_per_family_long <- median_cmp_per_family %>%
  pivot_longer(cols=c(anvio, ma),
               names_to="Method",
               values_to="cmp_completeness") %>%
  rename(kfs_completeness = kofamscan)
median_cmp_per_family_long %>%
  ggplot(aes(x=kfs_completeness,
             y=cmp_completeness,
             color=Method,
             label=module)) +
  geom_point(alpha=0.8)  +
  geom_abline(slope=1,intercept=0,lty=2,col="#AAAAAA") +
  facet_wrap(~ family) +
  theme_minimal() +
  xlab("Completeness (kofamscan)") +
  ylab("Completeness (other method)") +
  scale_color_manual(values=c(anvio="#2A9D8F", ma="#6DB9E4"))
# ggplotly()
  
# Generate some additional statistics

frac_higher_cmp_by_family_tbl <- completeness %>%
  mutate(anvio_higher = an_pwc > ko_pwc,
         ma_higher = ma_pwc > ko_pwc) %>%
  group_by(family) %>%
  summarize(frac_anvio_higher = mean(anvio_higher),
            frac_ma_higher = mean(ma_higher)) %>%
  arrange(-frac_anvio_higher)
write_tsv(frac_higher_cmp_by_family_tbl,
          "frac_higher_completeness_per_family.tsv")
frac_higher_avg_cmp_by_family_tbl <- median_cmp_per_family %>%
  mutate(anvio_higher = anvio > kofamscan, ma_higher = ma > kofamscan) %>%
  group_by(family) %>%
  summarize(frac_anvio_higher = mean(anvio_higher),
            frac_ma_higher = mean(ma_higher)) %>%
  arrange(-frac_anvio_higher)
write_tsv(frac_higher_avg_cmp_by_family_tbl,
          "frac_higher_median_completeness_per_family.tsv")

# slight misnomer since we now already have the bitscores included; here we're
# just returning the anvio scores along with whether they were "best" by kfs
anvio_vs_kfs_bitscore <- combined_fxns_fams %>%
  filter(tool == "anvio") %>%
  select(gene_callers_id,accession,threshold,bitscore,genome,gtdb_family) %>%
  distinct() %>%
  left_join(., select(filter(combined_fxns_fams, tool=="kfs"),
                      genome,
                      gene_callers_id,
                      accession,
                      best),
            by=c("genome", "gene_callers_id", "accession")) %>%
  select(gtdb_family, gene_callers_id, accession, genome, threshold,
         bitscore, best) %>%
  distinct()

# Note that sometimes the MA score and the kfs bitscore differ
ma_with_thresh <- combined_fxns_fams %>%
  filter(tool=="ma") %>%
  select(gtdb_family, gene_callers_id, accession,
         genome, bitscore, score, threshold) %>%
  distinct()

top_20_ma_diffs <- fxn_diffs_all %>%
  select(gtdb_family,accession,ma_diff) %>%
  inner_join(., all_kfs_thresholds) %>%
  filter(ma_diff > 0) %>%
  group_by(gtdb_family) %>%
  slice_max(ma_diff, n=20) %>%distinct() 

# Actually, the most MA-added hits have no threshold! So figure can't be made. 
print(top_20_ma_diffs %>% filter(gtdb_family=="Lachnospiraceae"), n=20)

# Consdering only those with thresholds:
top_20_ma_diffs_all_thresh <- fxn_diffs_all %>%
  select(gtdb_family,accession,ma_diff) %>%
  inner_join(., all_kfs_thresholds) %>%
  filter(!is.na(threshold)) %>%
  filter(ma_diff > 0) %>%
  group_by(gtdb_family) %>%
  slice_max(ma_diff, n=20) %>%distinct() 


# Get Lachno specific sequences...
write_csv(file="top_20_anvio_lachnos_accessions.csv",
          anvio_fxns %>%
            filter(accession %in% top_20_anvio_lachnos) %>%
            left_join(., genomes, by=c("genome"="sample_name")) %>%
            filter(gtdb_family=="Lachnospiraceae") %>%
            select(gene_callers_id, genome, accession) %>%
            arrange(accession,genome,gene_callers_id))

write_csv(file="top_20_MA_lachnos_accessions.csv",
          ma_fxns %>%
            filter(accession %in%
                     (filter(top_20_ma_diffs,
                             gtdb_family=="Lachnospiraceae")$accession)) %>%
            left_join(., genomes, by=c("genome"="sample_name")) %>%
            filter(gtdb_family=="Lachnospiraceae") %>%
            select(gene_callers_id, genome, accession) %>%
            arrange(accession,genome,gene_callers_id))

# Get all sequences
write_csv(file="top_20_anvio_lachnos_all_accessions.csv",
          anvio_fxns %>%
            filter(accession %in% top_20_anvio_lachnos) %>%
            left_join(., genomes, by=c("genome"="sample_name")) %>%
            select(gene_callers_id, genome, accession) %>%
            arrange(accession,genome,gene_callers_id))


# Also try this for the other citrate synthase related genes
write_csv(file="anvio_lachno_citr.csv",
          anvio_fxns %>%
            filter(accession %in% c("K05942", "K01643")) %>%
            left_join(., genomes, by=c("genome"="sample_name")) %>%
            filter(gtdb_family=="Lachnospiraceae") %>%
            select(gene_callers_id, genome, accession) %>%
            arrange(accession,genome,gene_callers_id))

# Plot Supplemental Figure showing bitscore distributions for Lachnos and others
anvio_kfs_figure <- anvio_vs_kfs_bitscore %>%
  mutate(family = ifelse(gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
  mutate(bit_norm = map2_dbl(bitscore, threshold, ~ {log2(.x/.y)})) %>%
  filter(accession %in% filter(top_20_anvio_diffs_all,
                               gtdb_family=="Lachnospiraceae")$accession) %>%
  ggplot(aes(x=bit_norm, fill=family, color=family)) +
  geom_density(bw="sj", adjust=5, alpha=0.7) +
  facet_wrap(~ accession) +
  scale_alpha_discrete(range=c(0.25,0.75)) +
  xlim(-2,2)+
  geom_vline(xintercept=0, lty=2) +
  theme_minimal() +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1")
print(anvio_kfs_figure)


## Not used: show what it looks like if we take all the hits
anvio_vs_all <- bind_rows(
  anvio_vs_kfs_bitscore %>%
    mutate(family = ifelse(
      gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
    mutate(bit_norm = map2_dbl(bitscore, threshold, ~ {log2(.x/.y)})) %>%
    filter(accession %in% filter(top_20_anvio_diffs_all,
                                 gtdb_family=="Lachnospiraceae")$accession) %>%
    select(bit_norm, family, accession, genome) %>%
    mutate(score_type="anvio_kfs"),
  kfs_fxns_all %>%
    filter(accession %in% top_20_anvio_lachnos) %>%
    collect() %>% 
    left_join(., genomes, by=c("genome" = "sample_name")) %>%
    mutate(bit_norm = map2_dbl(bitscore, threshold, ~ {log2(.x/.y)})) %>%
    mutate(family = ifelse(
      gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
    #filter(best == "") %>%
    select(bit_norm, family, accession, genome) %>%
    mutate(score_type="everything")
)
full_anvio_plot <- anvio_vs_all %>%
  mutate(score_type=factor(score_type,
                           levels=c("everything", "anvio_kfs"))) %>%
  ggplot(aes(x=bit_norm,
             fill=family,
             color=family,
             alpha=score_type,
             lty=desc(score_type))) +
  geom_density(bw="sj", adjust=5) +
  facet_wrap(~ accession) +
  scale_alpha_discrete(range=c(0.25,0.75)) +
  xlim(-6,2) +
  geom_vline(xintercept=0, lty=2) +
  theme_minimal() +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1") +
  scale_alpha_discrete(range=c(0.25,0.75)) +
  scale_linetype_binned()


# Barely any for MA, but keeping just in case
# most diff via anvio, considering only thresholded:
ma_with_thresh %>%
  mutate(family = ifelse(gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
  mutate(bit_norm = map2_dbl(score, threshold, ~ {log2(.x/.y)})) %>%
  filter(accession %in% filter(top_20_anvio_diffs_all,
                               gtdb_family=="Lachnospiraceae")$accession) %>%
  ggplot(aes(x=bit_norm, fill=family, color=family)) +
  geom_density(bw="sj", adjust=5) +
  facet_wrap(~ accession) +
  scale_alpha_discrete(range=c(0.25,0.75)) +
  xlim(-2,2)+
  geom_vline(xintercept=0, lty=2) +
  theme_minimal() +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1")

# most different via MA
ma_with_thresh %>%
  mutate(family = ifelse(gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
  mutate(bit_norm = map2_dbl(score, threshold, ~ {log2(.x/.y)})) %>%
  filter(accession %in% filter(top_20_ma_diffs_all_thresh,
                               gtdb_family=="Lachnospiraceae")$accession) %>%
  ggplot(aes(x=bit_norm, fill=family, color=family)) +
  geom_density(bw="sj", adjust=1) +
  facet_wrap(~ accession) +
  scale_alpha_discrete(range=c(0.25,0.75)) +
  xlim(-2,2)+
  geom_vline(xintercept=0, lty=2) +
  theme_minimal() +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1")

anvio_vs_kfs_bitscore %>%
  mutate(family = ifelse(gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
  mutate(bit_norm = map2_dbl(bitscore, threshold, ~ {log2(.x/.y)})) %>%
  filter(accession %in% filter(top_20_ma_diffs_all_thresh,
                               gtdb_family=="Lachnospiraceae")$accession) %>%
  ggplot(aes(x=bit_norm, fill=family, color=family)) +
  geom_density(bw="sj", adjust=1) +
  facet_wrap(~ accession) +
  scale_alpha_discrete(range=c(0.25,0.75)) +
  xlim(-2,2) +
  geom_vline(xintercept=0, lty=2) +
  theme_minimal() +
  scale_color_brewer(type="qual", palette="Set1") +
  scale_fill_brewer(type="qual", palette="Set1")

### Alignment/trees analysis

# Make annotations that we will use for the alignments
aln_annotations <- anvio_vs_kfs_bitscore  %>%
  filter(accession %in% top_20_anvio_lachnos) %>%
  select(gene_callers_id, genome, accession,
         gtdb_family, bitscore, threshold) %>%
  mutate(name = paste0(gene_callers_id, "_", genome)) %>%
  mutate(d_family = substr(gtdb_family, 1,9)) %>%
  mutate(d_bitscore = sprintf("b%05.0f", bitscore)) %>%
  mutate(d_above = ifelse(bitscore >= threshold, "Y", "N")) %>%
  mutate(d_combo = paste(d_family, d_above, d_bitscore, sep="-")) %>%
  relocate(name)

write_tsv(aln_annotations %>% select(name, annotation=d_family),
          "aln_annotation_family.tsv")
write_tsv(aln_annotations %>% select(name, annotation=d_bitscore),
          "aln_annotation_bitscore.tsv")
write_tsv(aln_annotations %>% select(name, annotation=d_above),
          "aln_annotation_above.tsv")
write_tsv(aln_annotations %>% select(name, annotation=d_combo),
          "aln_annotation_combo.tsv")
write_tsv(aln_annotations %>% select(name, d_family:d_combo),
          "aln_annotation_file.tsv")

# look at alignments

all_aln_files <- list.files("alignments/all_trim/", pattern = "aln$")
all_alns <- lapply(file.path("alignments/all_trim/", all_aln_files),
                   \(.) Biostrings::readAAMultipleAlignment(filepath=.,
                                                            format="fasta"))
names(all_alns) <- gsub("\\-trim.aln", "", all_aln_files)

lachno_aln_files <- list.files("alignments/top_fasta/", pattern = "aln$")
lachno_alns <- lapply(file.path("alignments/top_fasta//", lachno_aln_files),
                      \(.) Biostrings::readAAMultipleAlignment(filepath=.,
                                                               format="fasta"))
names(lachno_alns) <- gsub("\\.aln", "", lachno_aln_files)

all_tree_files <- list.files("trees/", pattern = ".tree$")
all_trees <- lapply(file.path("trees/", all_tree_files), read.tree)
names(all_trees) <- gsub("-trim.tree", "", all_tree_files)



# RPS-BLAST

blast6_header <- c(
  "qseqid",
  "sseqid",
  "pident",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore")
blast6_types = "ccdddddddddd"
cd_res_f <- list.files("./cdsearch/", "K......txt")
cd_res <- map(cd_res_f, ~ read_tsv(file.path("./cdsearch", .x),
                                   col_names = blast6_header,
                                   col_types=blast6_types))
names(cd_res) <- gsub("\\.txt", "", cd_res_f)

cdr_data <- map(cd_res, function(cdr) {
  top_each <- cdr %>% group_by(qseqid) %>%
    slice_max(order_by=bitscore, n=1) %>%
    ungroup()
})

# read in CDD database
cddids <- read_tsv("cddid_all.tbl.gz",
                   col_names = c("PSSM_id",
                                 "CD_accession",
                                 "CD_name",
                                 "CD_desc",
                                 "PSSM_length")) %>%
  mutate(CDD_id = paste0("CDD:", PSSM_id))


# parse out by best sseqid

K01647_which_cdds <- cdr_data$K01647 %>% count(sseqid) %>% arrange(-n)
print(K01647_which_cdds)
K01647_seqs <- read.FASTA("./alignments/top_fasta/K01647.faa", type="AA")
K01647_which_prk <- filter(cdr_data$K01647, sseqid=="CDD:184465")$qseqid
K01647_prks <- filter(cdr_data$K01647, sseqid=="CDD:184465")$qseqid
K01647_citls <- filter(cdr_data$K01647, sseqid=="CDD:99866")$qseqid
write.FASTA(file = "K01647_prk_184465.faa", K01647_seqs[K01647_prks])
write.FASTA(file = "K01647_cit_like_99866.faa", K01647_seqs[K01647_citls])

# Show that citls/prks are very close
print(cd_res$K01647 %>%
        filter(qseqid %in% K01647_citls) %>%
        arrange(qseqid, -bitscore) %>%
        group_by(qseqid) %>%
        slice_max(order_by=bitscore, n=2))

# same thing for HemNZ

K02495_which_cdds <- cdr_data$K02495 %>% count(sseqid) %>% arrange(-n)
print(K02495_which_cdds)
K02495_seqs <- read.FASTA("./alignments/top_fasta/K02495.faa", type="AA")
K02495_seq_ids <- map(K02495_which_cdds$sseqid, \(cdd) {
  return(dplyr::filter(cdr_data$K02495, sseqid==cdd)$qseqid)
})
names(K02495_seq_ids) <- K02495_which_cdds$sseqid
walk2(K02495_seq_ids, names(K02495_seq_ids), \(seqids, cdd) {
  fn = paste0("K02495_",
              gsub("CDD:", "", cdd),
              "_seqs.faa")
  print(fn)
  write.FASTA(file = fn,
              K02495_seqs[seqids])
})

# 440400 is not as close for these:
print(cd_res$K02495 %>%
        filter(qseqid %in% c(K02495_seq_ids$`CDD:181292`)) %>%
        arrange(qseqid, -bitscore) %>%
        group_by(qseqid) %>%
        slice_max(order_by=bitscore, n=5) %>%
        select(qseqid, sseqid, bitscore) %>%
        pivot_wider(names_from=sseqid, values_from=bitscore)) 

# but 274910 is equally as close for this one:
print(cd_res$K02495 %>%
        filter(qseqid %in% c(K02495_seq_ids$`CDD:236187`)) %>%
        arrange(qseqid, -bitscore) %>%
        group_by(qseqid) %>%
        slice_max(order_by=bitscore, n=5) %>%
        select(qseqid, sseqid, bitscore) %>%
        pivot_wider(names_from=sseqid, values_from=bitscore)) 

# Plotting

families <- sort(unique(genomes$gtdb_family))
fam_cols <- chart_cols[colors_in_order]
names(fam_cols) <- families

make_viz_cd_tree_2 <- function(phy, phy_n,
                               cdr = cdr_data,
                               annot = aln_annotations,
                               remove_guides=FALSE,
                               small_guides=FALSE,
                               tip_alpha=1) {
  viz_tree <- ggtree(midpoint_root(phy))
  viz_ylims <- ggplot_build(viz_tree)$layout$panel_scales_y[[1]]$range$range
  scalebar_ypos = max(viz_ylims) - (max(viz_ylims) / 20)
  viz_tree <- viz_tree + geom_treescale(y=scalebar_ypos,
                                        x=0,
                                        linesize = .8,
                                        offset=1)
  cd <- left_join(annot,  cdr[[phy_n]], by=c("name" = "qseqid")) %>%
    left_join(., cddids, by=c("sseqid" = "CDD_id")) %>%
    mutate(CD_label = paste0(CD_name, "(", sseqid, ")"))
  if (small_guides) cd <- mutate(cd, CD_label = gsub("\\(", "\n(", CD_label))
  cd_count <- cd %>% count(CD_label) %>% arrange(-n) %>% deframe
  cd <- cd %>%
    mutate(CD_label = map_chr(CD_label, ~ {
      if(.x %in% names(cd_count)[1:5]) return(.x) else return("other")
      }))
  man_cols <- fam_cols[1:6]
  names(man_cols) <- c(names(cd_count)[1:5], "other")
  if (small_guides) {
    man_cols["NA\n(NA)"] <- "#AAAAAA"
  } else {
    man_cols["NA(NA)"] <- "#AAAAAA"
  }
  man_cols["other"] <- "#444444"
  v2 <- viz_tree %<+% cd +
    geom_tippoint(aes(color=CD_label)) +
    scale_color_manual(values=man_cols) +
    geom_fruit(aes(fill=d_above), pwidth = 0.1, width=0.1, offset=0.2,
               geom=geom_tile, show.legend=FALSE) +
    scale_fill_manual(values=c(N="#000000FF", Y="#FFFFFF00")) +
    ggtitle(paste(phy_n, "—", ko_tsv %>%
                    filter(accession==phy_n) %>%
                    .$short_desc)) +
    theme(legend.key.spacing.y=unit(-5, "pt"),
          legend.spacing.y=unit(3, "pt"),
          legend.margin = margin(t=1,r=1,b=1,l=1, "pt"),
          legend.text = element_text(margin = margin(t = 1)), 
          legend.title = element_text(margin = margin(b = 2)))
  if (remove_guides) {
    v2 + theme(legend.position="none")
  } else if (small_guides) {
    v2 +
      guides(color = guide_legend(override.aes = list(size = 0.5))) +
      theme(legend.text = element_text(size=5,
                                       margin = margin(t = 1, unit = "pt")), 
            legend.spacing = unit(1, "pt"),
            legend.title = element_blank(),
            legend.position="bottom")
  } else {
    v2
  }
}

alt_viz_tree <- function(.x, .y, remove_guides=FALSE,
                         small_guides=FALSE, mask_frac=1) {
  v2 <- facet_plot(.x + new_scale_fill(),
             panel="alignment",
             data=tidy_msa(Biostrings::maskGaps(.y,
                                                min.fraction=mask_frac,
                                                min.block.width=1)),
             geom=geom_msa,
             font=NULL,
             border=NA) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(-2, "pt")
    )
  if (remove_guides) {
    v2 + theme(legend.position="none")
  } else {
    v2
  }
}

L_trees <- map(all_trees, ~ {
  keep.tip(.x, intersect(.x$tip.label,
                         (aln_annotations %>%
                            filter(gtdb_family=="Lachnospiraceae") %>%
                            .$name)))
})
viz_Ltrees <- map2(L_trees, names(L_trees), ~
                     make_viz_cd_tree_2(.x, .y, tip_alpha=1))
viz_Ltrees_rg <- map2(L_trees, names(L_trees), ~
                        make_viz_cd_tree_2(.x, .y, remove_guides=TRUE))
viz_Ltrees_sg <- map2(L_trees, names(L_trees), ~
                        make_viz_cd_tree_2(.x, .y, small_guides=TRUE))

viz_msas_L2 <- map2(viz_Ltrees, lachno_alns[names(viz_Ltrees)], ~
                      alt_viz_tree(.x, .y,
                                   remove_guides=FALSE,
                                   mask_frac = 0.75))
viz_msas_L2_rg <- map2(viz_Ltrees_rg, lachno_alns[names(viz_Ltrees)], ~
                         alt_viz_tree(.x, .y,
                                      remove_guides=FALSE,
                                      mask_frac = 0.75))
viz_msas_L2_sg <- map2(viz_Ltrees_sg, lachno_alns[names(viz_Ltrees)], ~
                         alt_viz_tree(.x, .y,
                                      mask_frac = 0.75))

mega_msas_2s <- Reduce(`+`, viz_msas_L2_sg) + plot_layout(ncol=5) &
  theme(legend.key.spacing.x=unit(0,"pt"),
        legend.key.spacing.y=unit(1.5,"pt"),
        legend.key.height=unit(1, "lines"),
        legend.spacing.x=unit(15,"pt"),
        legend.text=element_text(size=5),
        legend.key.size=unit(5,"pt"),
        legend.background = element_rect(color="#AAAAAA", linewidth = 0.5),
        legend.margin = margin(1,1,1,1))

pdf(width=25,height=15,file = "all_distros_lachno_cd.pdf")
anvio_kfs_figure
dev.off()

pdf(width=25,height=15,file = "all_tree_msas_lachno_cd.pdf")
#anvio_kfs_figure
mega_msas_2s
dev.off()

png(width=2500,height=1500,file = "all_distros.png")
anvio_kfs_figure & theme(text=element_text(size=15),
                         title=element_text(size=20),
                         axis.text=element_text(size=15),
                         axis.title =element_text(size=20)) &
  xlab("log2(bitscore/threshold)")
dev.off()
png(width=15,height=9, units = "in", res=300, file="all_tree_alns.png")
mega_msas_2s
dev.off()
#tiff(width=2500,height=1500, file="all_tree_alns.tiff", type="cairo", compression="webp")
#mega_msas_2s
#dev.off()

pdf(width=26,height=15.6, file="all_tree_alns.pdf")
mega_msas_2s
dev.off()

anvio_kfs_figure_indiv <- map(names(L_trees), ~ {
  anvio_vs_kfs_bitscore %>%
    mutate(family = ifelse(gtdb_family=="Lachnospiraceae", "lachno", "other")) %>%
    mutate(bit_norm = map2_dbl(bitscore, threshold, ~ {log2(.x/.y)})) %>%
    filter(accession == .x) %>%
    ggplot(aes(x=bit_norm, fill=family, color=family)) +
    geom_density(bw="sj", adjust=5, alpha=0.7) +
    scale_alpha_discrete(range=c(0.25,0.75)) +
    xlim(-2,2)+
    geom_vline(xintercept=0, lty=2) +
    theme_minimal() +
    scale_color_brewer(type="qual", palette="Set1") +
    scale_fill_brewer(type="qual", palette="Set1") +
    ggtitle(paste(.x, "—", ko_tsv %>% filter(accession==.x) %>% .$short_desc))
})
names(anvio_kfs_figure_indiv) <- names(L_trees)

# Generate new Fig2?
pdf(width=10,height=15,file="density_vs_msas.pdf")
subs <- c("K00615", "K02495", "K04771",  "K01647", "K02529", "K01992")
lh <- Reduce(`/`, map(anvio_kfs_figure_indiv[subs], ~
                  .x + xlab("log2(bitscore/threshold)"))) +
         plot_layout(guides="collect", axes="collect") &
  theme(legend.position = "bottom")
rh <- Reduce(`/`, map(viz_msas_L2[subs], ~
                        .x + guides(fill="none"))) &
  theme(legend.position="bottom",
        legend.spacing = unit(1, "pt"),
        legend.margin=margin(t=-10,0,0,0, unit="pt"),
        legend.key.spacing = unit(1, "pt"),
        plot.margin=margin(1,1,1,1, unit="pt"))
(lh | rh) + plot_layout(widths=c(1,2))
dev.off()

## Examples

# Write out butyrate annotations
write_csv(file="butyr_g.csv",
          anvio_vs_kfs_bitscore %>%
            filter(accession %in% kegg_butyrate_found$accession) %>%
            filter(gtdb_family=="Lachnospiraceae") %>%
            select(gene_callers_id, genome, accession) %>%
            arrange(accession,genome,gene_callers_id))

## After running alignments, etc., proceed:

# Butyrate example
buty_aln_files <- list.files("alignments/butyrate/aln/", pattern = "aln$")
buty_alns <- lapply(file.path("alignments/butyrate/aln/", buty_aln_files),
                    \(.) Biostrings::readAAMultipleAlignment(filepath=.,
                                                             format="fasta"))
names(buty_alns) <- gsub("\\.aln", "", buty_aln_files)

buty_cdd_files <- list.files("alignments/butyrate/cdd//", pattern = "txt$")
buty_cdd <- lapply(file.path("alignments/butyrate/cdd/", buty_cdd_files),
                    \(.) read_tsv(file.path(.), col_names = blast6_header,
                                  col_types=blast6_types))
names(buty_cdd) <- gsub("\\.txt", "", buty_cdd_files)
buty_cdr <- map(buty_cdd, function(cdr) {
  top_each <- cdr %>% group_by(qseqid) %>%
    slice_max(order_by=bitscore, n=1) %>%
    ungroup()
})

buty_tree_files <- list.files("alignments/butyrate/tree//", pattern = "tree$")
buty_trees <- lapply(file.path("alignments/butyrate/tree/", buty_tree_files),
                     read.tree)
names(buty_trees) <- gsub("\\.tree", "", buty_tree_files)

butyr_annotations <- anvio_vs_kfs_bitscore  %>%
  filter(accession %in% names(buty_trees)) %>%
  select(gene_callers_id, genome, accession,
         gtdb_family, bitscore, threshold) %>%
  mutate(name = paste0(gene_callers_id, "_", genome)) %>%
  mutate(d_family = substr(gtdb_family, 1,9)) %>%
  mutate(d_bitscore = sprintf("b%05.0f", bitscore)) %>%
  mutate(d_above = ifelse(bitscore >= threshold, "Y", "N")) %>%
  mutate(d_combo = paste(d_family, d_above, d_bitscore, sep="-")) %>%
  relocate(name)

viz_buty <- map2(buty_trees, names(buty_trees), ~
                   make_viz_cd_tree_2(.x, .y,
                                      cdr = buty_cdr,
                                      annot=butyr_annotations))
viz_msas_buty <- map2(viz_buty, buty_alns[names(viz_buty)], ~
                        alt_viz_tree(.x, .y,
                                     remove_guides=FALSE,
                                     mask_frac = 0.75))
viz_msas_L2_rg <- map2(viz_Ltrees_rg, lachno_alns[names(viz_Ltrees)], ~
                         alt_viz_tree(.x, .y,
                                      remove_guides=FALSE,
                                      mask_frac = 0.75))
viz_msas_L2_sg <- map2(viz_Ltrees_sg, lachno_alns[names(viz_Ltrees)], ~
                         alt_viz_tree(.x, .y, mask_frac = 0.75))

glutacon_plus <- Biostrings::readAAMultipleAlignment(
  filepath="alignments/butyrate_test/glutaconyl-K23351.aln", format="fasta")
glutacon_aln <- glutacon_plus
names(glutacon_aln@unmasked) <- map_chr(
  labels(glutacon_aln@unmasked), ~ gsub("([^ ]*) .*", "\\1", .x))
glutacon_tree <- read.tree("alignments/butyrate_test/glutaconyl-K23351.tree")

pdf("glutacon.pdf")
make_viz_cd_tree_2(midpoint.root(glutacon_tree),
                   "K23351",
                   cdr = buty_cdr,
                   annot=butyr_annotations) %>%
  alt_viz_tree(., glutacon_aln, mask_frac=.5) & theme(legend.position="bottom")
dev.off()

# Carbapenem example

carbapenem <- c(
  "K18768",
  "K18970",
  "K19316",
  "K22346",
  "K18794",
  "K19318",
  "K18971",
  "K18793",
  "K19319",
  "K19320",
  "K19321",
  "K19322",
  "K18972",
  "K19211",
  "K18976",
  "K21277",
  "K18782",
  "K18781",
  "K18780",
  "K19099",
  "K19216"
)

# Write out carbapenem annotations
write_csv(file="ma_carb_g.csv",
          ma_fxns %>%
            filter(accession %in% carbapenem) %>%
            left_join(., genomes, by=c("genome"="sample_name")) %>%
            select(gene_callers_id, genome, accession) %>%
            arrange(accession,genome,gene_callers_id))

# After running alignments, proceed:
carb_tree_files <- list.files("alignments/carbapenem/tree//", pattern = "tree$")
carb_trees <- lapply(file.path("alignments/carbapenem/tree/", carb_tree_files), read.tree)
names(carb_trees) <- gsub("\\.tree", "", carb_tree_files)

carbapenem_aln <- Biostrings::readAAMultipleAlignment(filepath="alignments/carbapenem/aln/K19099.aln", format="fasta")
carbapenem_tree <- read.tree("alignments/carbapenem/tree/K19099.tree")

carb_cdd_files <- list.files("alignments/carbapenem//cdd//", pattern = "txt$")
carb_cdd <- lapply(file.path("alignments/carbapenem//cdd/", carb_cdd_files),
                   \(.) read_tsv(file.path(.), col_names = blast6_header, col_types=blast6_types))
names(carb_cdd) <- gsub("\\.txt", "", carb_cdd_files)
carb_cdr <- map(carb_cdd, function(cdr) {
  top_each <- cdr %>% group_by(qseqid) %>% slice_max(order_by=bitscore, n=1) %>% ungroup()
})

carb_annotations <- ma_with_thresh %>%
  select(-bitscore) %>% # for this we use the score actually computed by MA
  rename(bitscore=score) %>%
  filter(accession %in% names(carb_trees)) %>%
  select(gene_callers_id, genome, accession, gtdb_family, bitscore, threshold) %>%
  mutate(name = paste0(gene_callers_id, "_", genome)) %>%
  mutate(d_family = substr(gtdb_family, 1,9)) %>%
  mutate(d_bitscore = sprintf("b%05.0f", bitscore)) %>%
  mutate(d_above = ifelse(bitscore >= threshold, "Y", "N")) %>%
  mutate(d_combo = paste(d_family, d_above, d_bitscore, sep="-")) %>%
  relocate(name)

make_viz_cd_tree_2(midpoint.root(carbapenem_tree), "K19099",
                   cdr = carb_cdr, annot=carb_annotations) %>%
  alt_viz_tree(., carbapenem_aln, mask_frac=1) & theme(legend.position="bottom")

gim_plus <- Biostrings::readAAMultipleAlignment(
  filepath="alignments/carb_test/GIM-K19099.aln",
  format="fasta")
gim_aln <- gim_plus
names(gim_aln@unmasked) <- map_chr(
  labels(gim_aln@unmasked), ~ gsub("([^ ]*) .*", "\\1", .x))
gim_tree <- read.tree("alignments/carb_test/GIM-K19099.tree")

pdf("gim.pdf")
(make_viz_cd_tree_2(midpoint.root(gim_tree),
                    "K19099", cdr = carb_cdr, annot=carb_annotations) %>%
    alt_viz_tree(., gim_aln, mask_frac=1) &
    theme(legend.position="bottom")) /
  (make_viz_cd_tree_2(midpoint.root(gim_tree),
                      "K19099",
                      cdr = carb_cdr, annot=carb_annotations) %>%
     alt_viz_tree(., gim_aln, mask_frac=.25) &
     theme(legend.position="bottom"))
dev.off()

### Glutaconate neighorhood analysis and plot

glutacon_tbl <- kfs_fxns_all %>%
  filter(genome %in% lachnos) %>%
  collect() %>%
  group_by(genome, gene_callers_id) %>%
  slice_max(bitscore, n=1) %>%
  filter(accession=="K23351") %>%
  arrange(-bitscore) %>%
  select(gene_callers_id, genome, bitscore, threshold, best, accession) %>%
  mutate(name = paste(gene_callers_id, genome, sep="_"))

glutacon_genomes <- glutacon_tbl$genome

glutacon_window <- glutacon_tbl %>%
  rename(central_gene = gene_callers_id) %>%
  mutate(gene_callers_id = map(as.numeric(central_gene), ~ (.x-10):(.x+10))) %>% 
  unnest() %>%
  mutate(gene_callers_id = as.character(gene_callers_id)) %>%
  select(central_gene, genome, gene_callers_id) %>%
  as_arrow_table() %>%
  left_join(kfs_fxns_all, by=c("genome", "gene_callers_id")) %>%
  collect() %>%
  ungroup()

glutacon_arranged <- glutacon_window %>%
  group_by(genome, gene_callers_id)  %>%
  slice_max(bitscore) %>%
  mutate(name = paste0(central_gene, "_", genome)) %>%
  mutate(position = as.numeric(gene_callers_id) - as.numeric(central_gene)) %>%
  select(name, accession, position, bitscore, threshold, best) %>%
  mutate(norm_bitscore = log2(bitscore/threshold)) %>%
  ungroup() %>%
  group_by(genome, name) %>%
  nest() %>%
  mutate(data = map(data, ~ {
    if (!("K23352" %in% .x$accession)) { return(.x) }
    which23352 = .x %>% filter(accession=="K23352") %>% .$position
    if (which23352 < 0) { .x <- mutate(.x, position = -position)}
    .x
  })) %>% 
  unnest() %>%
  ungroup()

# Parse GFF files -- this allows us to get coordinates etc. for plotting
lachno_gff_files <- file.path("anvio_make_gff", paste0(lachnos, ".gff"))
names(lachno_gff_files) <- lachnos
lachno_gff_files_found <- intersect(list.files("anvio_make_gff", 
                                               ".gff$",
                                               full.names = TRUE),
                                    lachno_gff_files)
lachno_gffs <- map(lachno_gff_files_found, ~
                     read_table(.x, col_names=c("seqname",
                                                "src",
                                                "feature_type",
                                                "start",
                                                "end",
                                                "score",
                                                "strand",
                                                "frame",
                                                "attribute"),
                                                  col_types="cccddcccc",
                                                  skip=1)) %>%
  Reduce(bind_rows, .) %>%
  mutate(gene_id = map_chr(attribute, ~  gsub("ID=(.*)_(\\d)___(.*)$",
                                              "\\3_\\1.\\2", .)))

anvio_L_names <- anvio_vs_kfs_bitscore %>%
  filter(gtdb_family=="Lachnospiraceae") %>%
  mutate(name = paste0(gene_callers_id, "_", genome)) %>%
  .$name

rename <- dplyr::rename
top12_ko <- glutacon_arranged %>%
  count(accession) %>% 
  filter(!is.na(accession)) %>%
  arrange(-n) %>%
  .[1:12,] %>%
  .$accession

ko_cols <- chart_cols[c("lime", "orange2", "steelblue2", "salmon", "skyblue",
                        "peach", "violetred", "oceanblue", "green", "tomato",
                        "puce", "tan")]
names(ko_cols) <- top12_ko
ko_cols["other"] <- "#888888"
ko_cols["NA"] <- "#cdcdcd"

glutacon_plot_data <- glutacon_arranged %>%
  rename(focal_gene = name) %>%
  mutate(indiv_gene = paste0(gene_callers_id, "_", genome)) %>%
  left_join(., lachno_gffs, by=c("indiv_gene"="gene_id")) %>%
  group_by(focal_gene) %>%
  filter(!is.na(strand)) %>%
  nest() %>%
  mutate(data=map(data, ~ {
   focal_strand <- filter(.x, accession=="K23351")$strand
    if (focal_strand == "-") {
      .x$xcoord1 <- -(.x$end)
      .x$xcoord2 <- -(.x$start)
      .x$indiv_fwd <- (.x$strand == focal_strand)
    } else {
      .x$xcoord1 <- (.x$start)
      .x$xcoord2 <- (.x$end)
      .x$indiv_fwd <- (.x$strand == focal_strand)
    }
    contig = filter(.x, accession=="K23351")$seqname
    filter(.x, seqname==contig) %>%
      mutate(focal_strand_fwd=(focal_strand=="+"))
  })) %>%
  unnest() %>%
  mutate(KO = map_chr(accession, ~ ifelse(.x %in% top12_ko, .x, "other"))) %>%
  filter(position >= -5, position <= 5)  %>%
  mutate(anvio_TF = map_lgl(indiv_gene, ~ .x %in% anvio_L_names)) %>%
  mutate(kfs_TF = map_lgl(best, ~ .x=="*")) %>%
  mutate(ntko_TF = map_lgl(threshold, is.na)) %>%
  mutate(gene_call = pmap_chr(list(a=anvio_TF, k=kfs_TF, n=ntko_TF),
                              function(a, k, n) {
    if (n) return("none")
    if (k) return("kfs")
    if (a) return("anvio")
    return("none")
  })) %>%
  ungroup()

# these are used to align the gene plots
dummies <- make_alignment_dummies(
  glutacon_plot_data,
  aes(xmin = xcoord1,
      xmax = xcoord2,
      y = genome,
      id = KO,
      forward=indiv_fwd,
      indiv_gene = indiv_gene),
  on = "K23351"
) %>%
  left_join(., select(glutacon_plot_data, -start, -end, -xcoord1, -xcoord2), by=c("genome", "KO"))

pdf(height=5, width=8, file="test_glutacon.pdf")
glutacon_plot_data %>%
  filter(position <= 4, position >= -4) %>%
  ggplot(aes(xmin=xcoord1,
             xmax=xcoord2,
             y=genome,
             fill=KO,
             size=gene_call,
             lty=gene_call, 
             color=gene_call,
             label=KO,
             forward=indiv_fwd
             )) +
  geom_gene_arrow() +
  geom_blank(data=dummies) +
  facet_wrap(~ genome, scales = "free", ncol = 2, strip.position="left") +
  scale_fill_manual(values=ko_cols) +
  scale_size_manual(values=c(kfs=.7, anvio=.5, none=0, ntKO=3)) +
  scale_linetype_manual(values=c(kfs=1, anvio=2, none=0)) +
  scale_color_manual(values=c(kfs="#000000", anvio="#888888", none="#dddddd")) +
  theme_void()
dev.off()

# citrate synthases in lachnos

pdf(width=6,height=6,"lachno_cit_synthases.pdf")
citrate_df <- combined_fxns_fams %>%
  filter(accession %in% c("K01643", "K01644", "K01647", "K05942"),
         genome %in% lachnos) %>%
  filter(tool=="anvio") %>%
  select(gene_callers_id, genome, accession) %>%
  distinct() %>%
  pivot_wider(names_from=accession, values_from=gene_callers_id,
              values_fn=length, values_fill=0) %>%
  pivot_longer(-genome) %>%
  mutate(name = forcats::fct_relevel(name, "K01647")) %>%
  pivot_wider(names_from=genome)
data.matrix(citrate_df[-1]) %>%
  (\(x) { rownames(x) <- citrate_df[[1]]; x }) %>%
  t() %>% 
  gheatmap(dist_method="canberra")
dev.off()


# correlations
tbl_to_mtx <- function(x) { m <- data.matrix(x[-1]); rownames(m) <- x[[1]]; m }
anvio_lachnos <- combined_fxns_fams %>% filter(tool=="anvio",
                                               genome %in% lachnos)
anvio_lachno_mtx <- anvio_lachnos %>%
  select(accession, genome) %>%
  group_by(genome) %>%
  count(accession)  %>%
  pivot_wider(names_from=accession,values_from=n, values_fill=0) %>%
  tbl_to_mtx
# we don't care about KOfams that are almost always at a single value 
anvio_testable <- anvio_lachno_mtx %>% (\(x) {
  n_not_at_modal_value <- apply(x, 2, \(.) {
    value_counts <- sort(c(table(.)), dec=T)
    modal_value <- names(value_counts[1])
    sum(. != modal_value)
  })
  names_testable <- names(which(n_not_at_modal_value > 2))
  x[, names_testable]
})
cor_anvio <- cor(anvio_testable)

# compute p-values for all correlations, so we get a better idea of the null
anvio_cor_pvals_tested <- matrix(
  nr=ncol(anvio_testable)
  nc=ncol(anvio_testable),
  NA, 
  dimnames=list(colnames(anvio_testable),
                colnames(anvio_testable)))

pb <- progress::progress_bar$new(
  total = (nrow(anvio_cor_pvals_tested) *
             (nrow(anvio_cor_pvals_tested) - 1)) / 2)

for (.x in 1:(ncol(anvio_testable)-1)) {
  for (.y in (.x+1):ncol(anvio_testable)) {
    pb$tick()
    anvio_cor_pvals_tested[.x, .y] <- cor.test(anvio_testable[,.x],
                                               anvio_testable[,.y])$p.value
  }
}
# correct p-values for multiple testing, then make symmetric
anvio_cor_pvals_df <- anvio_cor_pvals_tested %>%
  as_tibble(rownames="correlated_to") %>%
  pivot_longer(-correlated_to, names_to="accession", values_to="p.value") %>%
  filter(!is.na(p.value))
anvio_cor_qvals_df <- anvio_cor_pvals_df %>%
  mutate(q.value = qvalue::qvalue(p.value)$qvalues)
anvio_cor_full_qvals_df <- anvio_cor_qvals_df %>%
  (\(.) bind_rows(., rename(., accession2=correlated_to) %>%
                     rename(., correlated_to=accession, accession=accession2)))

top_btm_corrs <- list( # note, the top correlation is always that gene,
                       # which we do want to show the description of
  K01647_pos = cor_anvio["K01647",] %>% sort(dec=T) %>% .[1:11],
  K01647_neg = cor_anvio["K01647",] %>% sort(dec=F) %>% .[1:10],
  K01643_pos = cor_anvio["K01643",] %>% sort(dec=T) %>% .[1:11],
  K01643_neg = cor_anvio["K01643",] %>% sort(dec=F) %>% .[1:10],
  K05942_pos = cor_anvio["K05942",] %>% sort(dec=T) %>% .[1:11],
  K05942_neg = cor_anvio["K05942",] %>% sort(dec=F) %>% .[1:10]
)

top_btm_corrs_tbl <- imap(top_btm_corrs, \(corrs, name) {
  corrs %>%
    enframe(name="accession") %>%
    left_join(ko_tsv, by="accession") %>%
    mutate(Name=name)
}) %>%
  bind_rows() %>%
  relocate(Name, accession, value) %>%
  select(-short_desc) %>%
  separate_wider_delim(Name, delim="_",
                       names=c("correlated_to", "direction")) %>%
  left_join(., anvio_cor_full_qvals_df)

write_csv(top_btm_corrs_tbl, "top_btm_corrs_tbl.tsv")

# HemNZ analysis

hemnz_hits <- read_tsv("allProkHemNZ_hitdata.txt", skip=6) %>%
  separate_wider_delim(Query, delim=" - ", names=c("qnum", "Sequence_ID"))
hemnz_tsv <- read_tsv("hemn_seqs.tsv") %>% distinct()
hemnz_overall <- left_join(hemnz_tsv, hemnz_hits)
write_tsv(hemnz_overall, "hemnz_cdsearch_table.tsv")

hemnz_cdsearch_summary <- hemnz_overall %>%
  group_by(Type) %>%
  count(Accession) %>%
  mutate(Accession=replace_na(Accession, "none")) %>%
  group_by(Type) %>%
  (\(.) bind_rows(., summarize(., n = sum(n)) %>%
                    mutate(Accession="subfam_total"))) %>%
  ungroup() %>%
  group_by(Accession) %>%
  (\(.) bind_rows(., summarize(., n = sum(n)) %>%
                    mutate(Type="type_total"))) %>%
  pivot_wider(names_from=Type, values_from=n, values_fill=0) %>%
  left_join(., select(cddids, CD_accession, CD_name, CD_desc),
            by=c("Accession"="CD_accession"))  %>%
  arrange(-type_total)
write_tsv(hemnz_cdsearch_summary, "hemnz_cdsearch_summary.tsv")
