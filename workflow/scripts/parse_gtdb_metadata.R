library("tidyverse")
library("optparse")


option_list = list(
        make_option(c("-i", "--input"),
                        type="character", default=NULL),
        make_option(c("-l", "--low_completeness"),
                        type="character", default=NULL),
        make_option(c("-o", "--other_completeness"),
                        type="character", default=NULL),
        make_option(c("-c", "--family_count"),
                        type="character", default=NULL))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


genome_metadata <- read_tsv(opt$input)
metadata_L <- genome_metadata %>%
  separate_wider_delim(gtdb_taxonomy, delim = ";",
                       names = c("domain", "phylum", "class",
                                 "order", "family", "genus",
                                 "species"))
metadata_complete_count <- metadata_L %>%
  group_by(domain, phylum, class, order, family, genus, species) %>%
  summarize(n_complete = sum(contig_count <= 10),
            n_total = length(contig_count)) %>%
  ungroup()

metadata_L %>%
  ggplot(aes(x=log10(contig_count))) + geom_histogram()

family_count <- metadata_complete_count %>%
  group_by(domain, phylum, class, order, family) %>%
  summarize(n_complete = sum(n_complete), n_total = sum(n_total)) %>%
  ungroup() %>%
  mutate(frac_complete = n_complete / n_total) %>%
  arrange(-n_complete) 


most_genomes_low_completeness <- family_count %>%
  filter(frac_complete <= 0.05) %>%
  arrange(-n_total) %>%
  .[1:20, ]


most_genomes_other <- family_count %>%
  filter(frac_complete > 0.05) %>%
  arrange(-n_total) %>%
  .[1:20, ]

write_csv(family_count, opt$family_count)
write_csv(most_genomes_low_completeness, opt$low_completeness)
write_csv(most_genomes_other, opt$other_completeness)
