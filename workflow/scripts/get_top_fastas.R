library(seqinr, quietly=TRUE)
library(argparser, quietly=TRUE)

p <- arg_parser("Get fasta sequences listed in a .csv file")
p <- add_argument(p, "--input", help="input file", type="character", short="-i", default="results/top_20_anvio_lachnos_accessions.csv")
p <- add_argument(p, "--outdir", help="output directory", type="character", short="-o", default="results/top_fasta")
p <- add_argument(p, "--seqdir", help="directory to search for protein sequences", default="results/annotation/anvio/get_sequences_for_gene_calls/", short="-s")

# Parse the command line arguments
argv <- parse_args(p)
al_accessions <- read_csv(argv$input)

l_genomes <- al_accessions$genome %>% unique
accession_tbls <- al_accessions %>% group_by(accession) %>% nest() %>% deframe()
genome_dir <- argv$seqdir
genome_files <- file.path(genome_dir, paste0(l_genomes, ".faa"))
names(genome_files) <- l_genomes
seqs <- lapply(genome_files, \(f) seqinr::read.fasta(f, seqtype="AA", as.string=TRUE))

dir.create(argv$outdir)
seqs_per_acc <- map2(accession_tbls, names(accession_tbls), \(acc, accN) {
  acc_groups <- acc %>% group_by(genome) %>% nest() %>% deframe() %>% map(., ~ deframe(.x))
  all_seqs <- map2(names(acc_groups), acc_groups, \(n, g) {
                 s <- seqs[[n]][as.character(g)]
                 names(s) <- paste0(names(s), "_", n)
                 s
               }) %>% Reduce(c, .)
  all_seqs
  output_f <- file.path(argv$outdir, paste0(accN, ".faa"))
  print(output_f)
  write.fasta(all_seqs, names(all_seqs), file.out = output_f)
})
