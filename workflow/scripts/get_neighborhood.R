# Run after annotation-comparison.R
library(tidyverse)
get_neighborhood <- function(KO, fxns=kfs_fxns_all, genome_set=lachnos, nsize=5, matches_per_gene=3) {
  tbl <- fxns %>%
    filter(genome %in% genome_set) %>%
    collect() %>%
    group_by(genome, gene_callers_id) %>%
    slice_max(bitscore, n=1) %>%
    filter(accession==KO) %>%
    arrange(-bitscore) %>%
    select(gene_callers_id, genome, bitscore, threshold, best, accession) %>%
    mutate(name = paste(gene_callers_id, genome, sep="_")) 
  these_genomes <- tbl$genome
  window <- tbl %>%
    rename(central_gene = gene_callers_id) %>%
    mutate(gene_callers_id = map(as.numeric(central_gene), ~ (.x-nsize):(.x+nsize))) %>% 
    unnest(gene_callers_id) %>%
    mutate(gene_callers_id = as.character(gene_callers_id)) %>%
    select(central_gene, genome, gene_callers_id) %>%
    as_arrow_table() %>%
    left_join(filter(fxns, genome %in% these_genomes), by=c("genome", "gene_callers_id")) %>%
    collect() %>%
    ungroup() %>%
    group_by(genome, gene_callers_id) %>%
    slice_max(bitscore, n=matches_per_gene) %>%
    ungroup()
  rearrange_pivot <- window %>%
    filter(as.numeric(gene_callers_id) - as.numeric(central_gene) > 1) %>%
    count(accession) %>%
    arrange(-n) %>%
    filter(!is.na(accession)) %>%
    slice_max(n, n=1) %>%
    .$accession %>%
    .[1]
  message(paste0("pivoting by accession ", rearrange_pivot))
  arranged <- window %>%
    group_by(genome, gene_callers_id)  %>%
    slice_max(bitscore) %>%
    mutate(name = paste0(central_gene, "_", genome)) %>%
    mutate(position = as.numeric(gene_callers_id) - as.numeric(central_gene)) %>%
    select(name, accession, position, bitscore, threshold, best) %>%
    mutate(orientation="+") %>%
    mutate(norm_bitscore = log2(bitscore/threshold)) %>%
    ungroup() %>%
    group_by(genome, name) %>%
    nest() %>%
    mutate(data = map(data, ~ {
      if (!(rearrange_pivot %in% .x$accession)) { return(.x) }
      pivot_pos <- .x %>% filter(accession==rearrange_pivot) %>% .$position
      if (pivot_pos < 0) {
        .x <- mutate(.x, position = -position) %>%
          mutate(orientation="-")
        }
      .x
    })) %>% 
    unnest(c(data)) %>%
    ungroup()
  wide <- arranged %>%
    group_by(genome, gene_callers_id, name) %>%
    slice_max(bitscore, n=1) %>%
    slice_max(threshold, n=1) %>%
    ungroup() %>%
    mutate(pos = sprintf("%02+d", position)) %>%
    arrange(position) %>%
    select(genome, name, accession, pos) %>%
    pivot_wider(values_from=accession, names_from=pos)
  list(
    tbl=tbl,
    arranged=arranged,
    wide=wide
  )
}

