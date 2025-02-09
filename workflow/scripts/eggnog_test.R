# Run after annotation-comparison.R
# Note: ko00001.json.gz needs to be downloaded into the current directory from the KEGG BRITE website

library(jsonlite)
library(tidyverse)

read_eggnog <- function(file) {
  eggnog_coltypes <- "ccddccccccccccccccccc"
  read_tsv(file,
           col_types=eggnog_coltypes,
           comment="##") %>%
  rename(gene_callers_id = `#query`) %>%
  mutate(genome = gsub("\\.emapper\\.annotations$", "", basename(file))) %>%
  relocate(genome)
}

combine_eggnog <- function(eggnog, w = combined_fxns_fams_wide) {
  left_join(eggnog, w, by=c("genome", "gene_callers_id")) %>%
    select(genome, gene_callers_id, KEGG_ko, kfs, anvio, ma) %>%
    separate_longer_delim(KEGG_ko, ",") %>%
    mutate(KEGG_ko = gsub("ko:", "", KEGG_ko))
}

eval_rows <- function(tbl) {
  tbl %>%
    pivot_longer(c(kfs, anvio, ma),
                 names_to="tool",
                 values_to="accession") %>%
    mutate(match = ifelse(is.na(accession), NA, 
                          ifelse(KEGG_ko=="-", "no_eggnog",
                                 ifelse(accession==KEGG_ko,
                                        "match",
                                        "mismatch"))))
}

eval_test_cases <- function(tbl) {
  tbl %>%
    group_by(genome, gene_callers_id, tool) %>%
    summarize(#KEGG_kos = list(unique(KEGG_ko)),
              n_match = sum(match=="match"),
              n_mismatch = sum(match=="mismatch"),
              n_no_eggnog = sum(match=="no_eggnog"),
              n_NA = sum(is.na(match)),
              .groups="keep") %>%
    ungroup()
}

eval_summarize_per_row <- function(tbl) {
  tbl %>%
    mutate(summary = pmap_chr(
      .l = list(nm=n_match,
                nmm=n_mismatch,
                ne=n_no_eggnog,
                nna=n_NA),
      .f = function(nm, nmm, ne, nna) {
        if (nna > 0) { return("not_annotated") }
        if (ne > 0) { return("annot_no_eggnog") }
        if (nm > 0) {
          if (nmm > 0) return("annot_match_one") else return("annot_match_all")
        }
        return("annot_no_match")
      })) %>%
    select(genome, gene_callers_id, tool, summary)
}

eval_add_which_annot <- function(tbl) {
  tbl %>%
    mutate(annotated = map_lgl(summary, ~ .x != "not_annotated")) %>%
    select(-summary) %>%
    pivot_wider(names_from=tool,
                values_from=annotated,
                names_prefix = "annotated_") %>%
    left_join(tbl, ., by=c("genome", "gene_callers_id"))
}

make_annot_tbl <- function(tbl) {
  tbl %>% 
    pivot_wider(names_from=summary, values_from=n, values_fill = 0) %>%
    mutate(total_annot = annot_match_all + annot_match_one + annot_no_eggnog + annot_no_match) %>%
    mutate(
      frac_annot_any_match = (annot_match_all + annot_match_one) / total_annot,
      frac_annot_all_match = (annot_match_all) / total_annot,
      frac_annot_no_eggnog = annot_no_eggnog / total_annot,
      frac_annot_no_match = annot_no_match / total_annot) %>%
    relocate(c(total_annot, frac_annot_any_match:frac_annot_no_match), .after = "tool")
}

combined_fxns_fams <- read_csv("combined_fxns_fams.csv.xz") %>%
  mutate(gene_callers_id = as.character(gene_callers_id))
ko_tsv <- read_tsv("ko.tsv", col_names = c("accession", "desc")) %>%
  mutate(short_desc = gsub("([^;]+);.*", "\\1", desc))

## Reformat hits

combined_fxns_fams_wide <- combined_fxns_fams %>%
  select(genome, gene_callers_id, tool, accession) %>%
  pivot_wider(names_from=tool, values_from=accession)

## Read in EggNOG-mapper data and summarize

all_emapper_files <- list.files("./eggnog", ".*\\.emapper\\.annotations",
                                full.names = TRUE)
eggnog_tbl <- map(all_emapper_files, read_eggnog) %>%
  bind_rows() %>%
  mutate(gene_callers_id = as.character(gene_callers_id))
joined_eggnog <- combine_eggnog(eggnog_tbl, combined_fxns_fams_wide) %>% distinct()
joined_eval_rows <- eval_rows(joined_eggnog)
joined_eval_cases <- eval_test_cases(joined_eval_rows)
joined_eval_scases <- eval_summarize_per_row(joined_eval_cases)
joined_which_annot <- eval_add_which_annot(joined_eval_scases)
write_csv(joined_which_annot, file="joined_which_annot.csv")

## Calculate fractions

overall_annot_fracs <- joined_which_annot %>%
  count(tool, summary) %>%
  make_annot_tbl

genome_annot_fracs <- joined_which_annot %>%
  count(genome, tool, summary) %>%
  make_annot_tbl

overall_additional_annot_fracs <- joined_which_annot %>%
  filter(!annotated_kfs) %>%
  count(tool, summary) %>%
  make_annot_tbl %>%
  filter(tool != "kfs")

genome_additional_annot_fracs <- joined_which_annot %>%
  filter(!annotated_kfs) %>%
  count(genome, tool, summary) %>%
  make_annot_tbl %>%
  filter(tool != "kfs")

write_csv(overall_annot_fracs, "overall_annot_fracs.csv")
write_csv(genome_annot_fracs, "genome_annot_fracs.csv")
write_csv(overall_additional_annot_fracs, "overall_additional_annot_fracs.csv")
write_csv(genome_additional_annot_fracs, "genome_additional_annot_fracs.csv")

# Which families are the "worst", on average (median)?

af_overall_by_family <- genome_annot_fracs %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>%
  relocate(gtdb_family) %>%
  group_by(gtdb_family, tool) %>%
  summarize(across(total_annot:frac_annot_no_match, median)) %>%
  arrange(-frac_annot_no_match)

af_addl_by_family <- genome_additional_annot_fracs %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>%
  relocate(gtdb_family) %>%
  group_by(gtdb_family, tool) %>%
  summarize(across(total_annot:frac_annot_no_match, median)) %>%
  arrange(-frac_annot_no_match)

write_csv(af_by_family, "af_by_family.csv")
write_csv(af_addl_by_family, "af_addl_by_family.csv")

# Interestingly, both decline as total_annotations increases...
af_addl_by_family %>%
  ggplot(aes(x=total_annot, y=frac_annot_no_match, color=tool, shape=tool)) +
  geom_point()

## For anvio, Nanosynbacteraceae has the "worst" agreement...

anvio_nanosyn_tbl <- joined_which_annot %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>% # add gtdb info
  filter(gtdb_family=="Nanosynbacteraceae") %>%
  filter(tool=="anvio") %>%
  filter(annotated_anvio, !annotated_kfs, summary != "annot_no_eggnog") %>% # remove cases with no eggnog annotations
  left_join(., distinct(select(joined_eggnog, -KEGG_ko))) %>% # add which KOfams anvio (or ma) annotated; get rid of this KEGG_ko column
  left_join(., eggnog_tbl)  # add the rest of the eggnog information

anvio_lachno_tbl <- joined_which_annot %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>% # add gtdb info
  filter(gtdb_family=="Lachnospiraceae") %>%
  filter(tool=="anvio") %>%
  filter(annotated_anvio, !annotated_kfs, summary != "annot_no_eggnog") %>% # remove cases with no eggnog annotations
  left_join(., distinct(select(joined_eggnog, -KEGG_ko))) %>% # add which KOfams anvio (or ma) annotated; get rid of this KEGG_ko column
  left_join(., eggnog_tbl)  # add the rest of the eggnog information

ma_nanosyn_tbl <- joined_which_annot %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>% # add gtdb info
  filter(gtdb_family=="Nanosynbacteraceae") %>%
  filter(tool=="ma") %>%
  filter(annotated_anvio, !annotated_kfs, summary != "annot_no_eggnog") %>% # remove cases with no eggnog annotations
  left_join(., distinct(select(joined_eggnog, -KEGG_ko))) %>% # add which KOfams anvio (or ma) annotated; get rid of this KEGG_ko column
  left_join(., eggnog_tbl)  # add the rest of the eggnog information

ma_lachno_tbl <- joined_which_annot %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>% # add gtdb info
  filter(gtdb_family=="Lachnospiraceae") %>%
  filter(tool=="ma") %>%
  filter(annotated_anvio, !annotated_kfs, summary != "annot_no_eggnog") %>% # remove cases with no eggnog annotations
  left_join(., distinct(select(joined_eggnog, -KEGG_ko))) %>% # add which KOfams anvio (or ma) annotated; get rid of this KEGG_ko column
  left_join(., eggnog_tbl)  # add the rest of the eggnog information

write_csv(anvio_nanosyn_tbl, "anvio_nanosyn_eggnog.csv")
write_csv(anvio_lachno_tbl, "anvio_lachno_eggnog.csv")

# Load in BRITE hierarchy

brite_json <- jsonlite::read_json(path = "ko00001.json.gz") 

recursive_json <- function(x, lvl=0) {
  fields <- names(x)
  if ("children" %in% fields) {
    this_name = paste0("name", lvl)
    sub_tbl <- enframe(x) %>%
      pivot_wider() %>%
      unnest(c(name, children))
    colnames(sub_tbl)[which(colnames(sub_tbl)=="name")] <- this_name
    sub_tbl <- sub_tbl %>%
      mutate(children = map(children, ~ recursive_json(.x, lvl=lvl+1))) %>%
      mutate(ctbl = map_lgl(children, is_tibble))
    if (all(sub_tbl$ctbl)) {
      return(unnest(sub_tbl %>% select(-ctbl), children))
    }
    if (all(!sub_tbl$ctbl)) {
      return(unnest(sub_tbl %>% select(-ctbl), children))
    }
    dropped <- filter(sub_tbl, !ctbl)
    warning(paste0(dropped, " entries dropped"))
    return(filter(sub_tbl, ctbl) %>% select(-ctbl) %>% unnest(children))
  } else {
    return(x$name)
  }
}

brite_tbl <- recursive_json(brite_json) %>%
  separate(children, "  ",
           into=c("accession", "description"),
           extra = "merge") 

brite_jaccard <- function(cmp1, cmp2, brt=brite_tbl, lvl="name3") {
  if (is.na(cmp1)) return(NA)
  if (is.na(cmp2)) return(NA)
  br1 <- brt %>% filter(accession==cmp1) %>% pluck(lvl)
  br2 <- brt %>% filter(accession==cmp2) %>% pluck(lvl)
  length(intersect(br1, br2)) / length(union(br1, br2))
}

# note: 27.2K KOs in BRITE; we can filter out less-informative terms that include >5% of the total KOs
n_brite_accessions <- brite_tbl$accession %>% unique %>% length
rarer_brite_tbl <- brite_tbl %>%
  group_by(name3) %>%
  nest() %>%
  mutate(n_annots = map_dbl(data, nrow)) %>%
  filter(n_annots <= (n_brite_accessions * 0.05)) %>%
  unnest()

genomes <- read_tsv("gtdb_sample_table.tsv")

mismatches <- joined_eggnog %>%
  pivot_longer(c(kfs,anvio,ma),
               names_to="tool",
               values_to="accession") %>%
  left_join(.,
            joined_which_annot,
            by=c("genome", "gene_callers_id", "tool")) %>%
  filter(summary=="annot_no_match") %>%
  left_join(., genomes, by=c("genome"="sample_name"))

# at least one KO from eggNOG and one from the annotation tool must share all their KEGG BRITE annotations
# To restrict to terms with <=5% of the total accessions, replace `brite_tbl` with `rarer_brite_tbl`
mismatch_j <- mismatches %>%
  left_join(., select(brite_tbl, tool_name3=name3, accession), relationship = "many-to-many") %>%
  left_join(., select(brite_tbl, eggnog_name3=name3, accession), relationship = "many-to-many", by=c("KEGG_ko"="accession")) %>%
  filter(!is.na(tool_name3)) %>%
  filter(!is.na(eggnog_name3)) %>%
  group_by(genome, gene_callers_id, tool) %>%
  summarize(brite_jaccard = (length(intersect(tool_name3, eggnog_name3)) /
                               length(union(tool_name3, eggnog_name3))),
            .groups = "keep") %>%
  mutate(overlap = ifelse(brite_jaccard == 1, "Y", "N")) %>%
  mutate(overlap2 = ifelse(brite_jaccard >= 0.5, "Y", "N")) %>%
  group_by(genome, gene_callers_id, tool) %>%
  summarize(any_overlap = any(overlap=="Y"),
            mean_overlap = mean(overlap=="Y"),
            any_overlap2 = any(overlap2=="Y"),
            mean_overlap2 = mean(overlap2=="Y"),
            .groups="keep") %>%
  ungroup()

brite_name3_match_overall <- mismatch_j %>%
  left_join(., joined_which_annot) %>%
  group_by(tool) %>%
  summarize(n_overlap = sum(any_overlap), n_no_overlap=sum(!any_overlap)) %>%
  mutate(pct_overlap = n_overlap/(n_overlap+n_no_overlap))

brite_name3_match_additional <- mismatch_j %>%
  left_join(., joined_which_annot) %>%
  filter(!annotated_kfs) %>%
  group_by(tool) %>%
  summarize(n_overlap = sum(any_overlap), n_no_overlap=sum(!any_overlap)) %>%
  mutate(pct_overlap = n_overlap/(n_overlap+n_no_overlap))

ko_EC <- ko_tsv %>%
  mutate(EC = stringr::str_extract(desc, "\\[EC:(.*)\\]", group=1)) %>%
  separate_longer_delim(EC, delim=" ") %>%
  mutate(EC=gsub("\\.$", "", EC))


ko_EC_pos <- ko_EC %>%
  separate_wider_delim(EC, names = c("EC1","EC2","EC3","EC4"),
                       delim = ".",
                       too_few="align_start",
                       cols_remove=FALSE)

#, ~ str_split_1(.x, " "))) %>% select(accession,EC,desc) %>% unnest(EC)
ko_EC_spec <- ko_EC %>% filter(!grepl("-", EC))

eggnog_ECs_spec <- eggnog_tbl %>% filter(!grepl("-", EC)) %>% mutate(EC = map(EC, ~ str_split_1(.x, ","))) %>% unnest(EC)
eggnog_EC <- eggnog_tbl %>% separate_longer_delim(EC, delim=",") %>% mutate(EC = na_if(EC, "-"))
eggnog_EC_pos <- eggnog_EC %>%
  separate_wider_delim(EC, names = c("EC1","EC2","EC3","EC4"),
                       delim = ".",
                       too_few="align_start",
                       cols_remove = FALSE)

mm_eggnog_ECs <- mismatches %>%
  left_join(., (eggnog_ECs_spec %>%
                  select(genome, gene_callers_id, EC) %>%
                  distinct()),
            by=c("genome", "gene_callers_id"),
            relationship="many-to-many") %>%
  group_by(genome, gene_callers_id, annotated_anvio, annotated_kfs, annotated_ma) %>%
  select(EC) %>%
  distinct() %>%
  nest(.key = "eggnog_EC") 

mm_tool_ECs <- mismatches %>%
  left_join(., select(ko_EC_spec, accession, EC_KO=EC),
            by="accession",
            relationship="many-to-many") %>%
  group_by(genome, gene_callers_id, tool, annotated_anvio, annotated_kfs, annotated_ma) %>%
  select(EC_KO) %>%
  distinct() %>%
  nest(.key = "tool_EC")

mm_en_vs_tool_ECs <- left_join(mm_tool_ECs, mm_eggnog_ECs, by=c("genome", "gene_callers_id", "annotated_anvio", "annotated_kfs", "annotated_ma"))
mm_ec_jaccard <- mm_en_vs_tool_ECs %>%
  mutate(jaccard = map2_dbl(tool_EC, eggnog_EC, ~ {
    length(intersect(.y$EC, .x$EC_KO)) / length(union(.y$EC, .x$EC_KO))
    }, .progress = TRUE))
  
mm_J_tbl <- mm_ec_jaccard %>%
  group_by(tool) %>%
  mutate(J=ifelse(jaccard>0, "Y", "N")) %>%
  count(J) %>%
  pivot_wider(names_from=J, values_from=n) %>%
  mutate(pct = Y/(Y+N))

mm_J_tbl_2 <- mm_ec_jaccard %>%
  filter(!annotated_kfs) %>%
  group_by(tool) %>%
  mutate(J=ifelse(jaccard>0, "Y", "N")) %>%
  count(J) %>%
  pivot_wider(names_from=J, values_from=n) %>%
  mutate(pct = Y/(Y+N))

# Compare total fractions matching for each tool

genome_annot_fracs %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>%
  select(gtdb_family, genome, tool, frac_annot_no_match) %>%
  pivot_wider(names_from=tool, values_from=frac_annot_no_match) %>% # necessary to compare all to kfs
  pivot_longer(c(ma,anvio), names_to = "tool", values_to = "frac_annot_no_match") %>% # ditto
  left_join(., select(genome_additional_annot_fracs, genome, tool, addl_annot=total_annot)) %>%
  left_join(., select(genome_annot_fracs, genome, tool, total_annot)) %>%
  mutate(addl_annot_rate = (addl_annot / total_annot)) %>%
  ggplot(aes(x=kfs, y=frac_annot_no_match, color=tool, shape=tool, size=addl_annot_rate)) +
  geom_point() +
  scale_size(range = c(0.25, 2.5)) +
  geom_abline(slope=1,intercept=0) +
  theme_minimal() +
  theme(text=element_text(size=14),
        axis.text=element_text(size=12)) +
  xlab("Frac. Kofamscan annotations not matching EggNOG") +
  ylab("Frac. other annotations not matching EggNOG")


genome_annot_fracs %>%
  left_join(., genomes, by=c("genome"="sample_name")) %>%
  select(gtdb_family, genome, tool, frac_annot_no_match) %>%
  filter(tool=="kfs") %>%
  pivot_wider(names_from=tool, values_from=frac_annot_no_match) %>% # necessary to compare all to kfs
  left_join(., select(genome_additional_annot_fracs, frac_annot_any_match, frac_annot_all_match, frac_annot_no_match, genome, tool, addl_annot=total_annot)) %>%
  left_join(., select(genome_annot_fracs, genome, tool, total_annot)) %>%
  mutate(addl_annot_rate = (addl_annot / total_annot)) %>%
  ggplot(aes(x=total_annot, y=frac_annot_no_match, color=tool, shape=tool, size=addl_annot_rate)) +
  geom_point() +
  scale_size(range = c(0.25, 2.5)) +
  geom_abline(slope=1,intercept=0) +
  geom_hline(yintercept=0.1) +
  theme_minimal() +
  theme(text=element_text(size=14),
        axis.text=element_text(size=12)) +
  xlab("Frac. Kofamscan annotations not matching EggNOG") +
  ylab("Frac. additional annotations not matching EggNOG")


anvio_not_kfs <- combined_fxns_fams_wide %>% filter(is.na(kfs), !is.na(anvio))
ma_not_kfs <- combined_fxns_fams_wide %>% filter(is.na(kfs), !is.na(ma))

eggnog_plus_brite <- left_join(joined_eggnog,
                               select(brite_tbl, name0:accession),
                               by=c("KEGG_ko"="accession"),
                               relationship='many-to-many') %>%
  pivot_longer(kfs:ma, names_to = "tool", values_to="ann_ko") %>%
  left_join(.,
            select(brite_tbl, name0:accession),
            by=c("ann_ko"="accession"),
            relationship='many-to-many',
            suffix=c(".egg", ".ann")) %>%
  distinct()
write_tsv(eggnog_plus_brite, "eggnog_plus_brite.tsv.xz")

eggnog_brite_match <- eggnog_plus_brite %>%
  filter(!is.na(ann_ko)) %>%
  mutate(ko_match = KEGG_ko == ann_ko) %>%
  mutate(brite0_match = name0.ann==name0.egg) %>%
  mutate(brite1_match = name1.ann==name1.egg) %>%
  mutate(brite2_match = name2.ann==name2.egg) %>%
  mutate(brite3_match = name3.ann==name3.egg) %>%
  mutate(brite_match = brite0_match & brite1_match & brite2_match & brite3_match) %>%
  distinct() %>%
  group_by(genome, gene_callers_id, tool) %>%
  summarize(any_brite_match = any(na.omit(brite_match)),
            any_same_KO = any(na.omit(ko_match)),
            all_eggbrite_NA = all(is.na(name0.egg)),
            n_egg_KOs = length(unique(na.omit(setdiff(KEGG_ko, "-")))),
            n_annot_KOs = length(unique(na.omit(ann_ko))),
            .groups="keep") %>%
  ungroup()

match_with_eggnog <- function(input_annots, which_tool="anvio", brite=eggnog_brite_match) {
  this_brite <- brite %>%
    filter(!is.na(tool)) %>%
    filter(tool==which_tool) %>%
    select(-tool)
  this_by <- set_names("accession", quo_name(enquo(which_tool)))
  egg_pos <- left_join(input_annots,
                       eggnog_EC_pos,
                       by=c("genome", "gene_callers_id")) %>%
    select(genome, gene_callers_id, {{ which_tool }}, EC, EC1:EC4) %>%
    left_join(., ko_EC_pos, by=this_by,
              suffix = c(".egg", ".ann"),
              relationship="many-to-many")
  egg_match <- egg_pos %>%
    distinct() %>%
    mutate(EC1.match = EC1.egg==EC1.ann,
           EC2.match = (EC1.match & (EC2.egg==EC2.ann)),
           EC3.match = (EC2.match & (EC3.egg==EC3.ann)),
           EC4.match = (EC3.match & (EC4.egg==EC4.ann))) %>%
    group_by(genome, gene_callers_id) %>%
    summarize(
      n_tot = n(),
      n_egg_ECs = length(unique(na.omit(EC.egg))),
      n_ann_ECs = length(unique(na.omit(EC.ann))),
      n_EC1 = sum(na.omit(EC1.match)),
      n_EC2 = sum(na.omit(EC2.match)),
      n_EC3 = sum(na.omit(EC3.match)),
      n_EC4 = sum(na.omit(EC4.match)),
      allNA_egg = all(is.na(EC1.egg)),
      allNA_ann = all(is.na(EC1.ann)),
      .groups="keep") %>%
    ungroup()
  egg_match_brite <- egg_match %>%
    left_join(., this_brite, by=c("genome", "gene_callers_id")) %>%
    (\(x) left_join(input_annots, x, by=c("genome", "gene_callers_id")))
  return(list(pos=egg_pos, match=egg_match_brite))
}


anvio_not_kfs <- combined_fxns_fams_wide %>% filter(is.na(kfs), !is.na(anvio))
ma_not_kfs <- combined_fxns_fams_wide %>% filter(is.na(kfs), !is.na(ma))

anvio_ec_matches <- match_with_eggnog(anvio_not_kfs, "anvio")
ma_ec_matches <- match_with_eggnog(ma_not_kfs, "ma")
kfs_ec_matches <- match_with_eggnog(combined_fxns_fams_wide %>% filter(!is.na(kfs)),
                                    which_tool="kfs")


anvio_notkfs_ecbrite_match <- anvio_ec_matches$match %>%
  filter(n_egg_ECs > 0 | n_egg_KOs > 0) %>% # have to have some eggnog annotation
  mutate(match = any_brite_match | any_same_KO | (n_EC4>0)) %>% 
  count(match)

anvio_overall_rate <- anvio_notkfs_ecbrite_match %>% deframe %>% (\(.) .["TRUE"]/(.["TRUE"]+.["FALSE"]))

ma_notkfs_ecbrite_match <- ma_ec_matches$match %>%
  filter(n_egg_ECs > 0 | n_egg_KOs > 0) %>% # have to have some eggnog annotation
  mutate(match = any_brite_match | any_same_KO | (n_EC4>0)) %>% 
  count(match)

ma_overall_rate <- ma_notkfs_ecbrite_match %>% deframe %>% (\(.) .["TRUE"]/(.["TRUE"]+.["FALSE"]))

kfs_full_ecbrite_match <- kfs_ec_matches$match %>%
  filter(n_egg_ECs > 0 | n_egg_KOs > 0) %>% # have to have some eggnog annotation
  mutate(match = any_brite_match | any_same_KO | (n_EC4>0)) %>% 
  count(match)

kfs_overall_rate <- kfs_full_ecbrite_match %>% deframe %>% (\(.) .["TRUE"]/(.["TRUE"]+.["FALSE"]))

