library(seqinr)

al_accessions <- read_csv("results/top_20_anvio_lachnos_accessions.csv")

l_genomes <- al_accessions$genome %>% unique
accession_tbls <- al_accessions %>% group_by(accession) %>% nest() %>% deframe()
genome_dir <- "results/annotation/anvio/get_sequences_for_gene_calls/"
genome_files <- file.path(genome_dir, paste0(l_genomes, ".faa"))
names(genome_files) <- l_genomes
seqs <- lapply(genome_files, \(f) seqinr::read.fasta(f, seqtype="AA", as.string=TRUE))


dir.create("results/top_fasta")
seqs_per_acc <- lapply(accession_tbls, \(acc) {
  acc_groups <- acc %>% group_by(genome) %>% nest() %>% deframe() %>% map(., ~ deframe(.x))
  all_seqs <- map2(names(acc_groups), acc_groups, \(n, g) {
    s <- seqs[[n]][as.character(g)]
    names(s) <- paste0(names(s), "_", n)
  }) %>% Reduce(c, .)
  write.fasta(all_seqs, names(all_seqs), file.out = file.path("results/top_fasta/", paste0(acc, ".faa")))
  all_seqs
})