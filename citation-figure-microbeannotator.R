#install.packages("dplyr")
#install.packages("networkD3")

#library(networkD3)
#library(dplyr)

# This was a pain in the butt to make as there's not an easy api and
# no where has a citation tree that generates based on a single root
# paper. This had to be made by hand so may contain errors. As of Jan 23, 2024
# these are all the citations resulting from MicrobeAnnotator in google
# scholar.
src <- c("MicrobeAnnotator", "MicrobeAnnotator", 
         "MicrobeAnnotator","MicrobeAnnotator",
         "MicrobeAnnotator", "10.1016/j.envint.2023.108374",
         "MicrobeAnnotator", "MicrobeAnnotator", 
         "10.20517/mrr.2022.21", "10.1038/s41467-023-39931-2",
         "10.1038/s41467-023-39931-2", "10.1038/s41467-023-39931-2",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "other_book_1",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "10.1016/j.heliyon.2023.e17652", "10.1016/j.heliyon.2023.e17652",
         "10.1016/j.heliyon.2023.e17652", "MicrobeAnnotator",
         "other_book_3", "MicrobeAnnotator",
         "10.1186/s12934-023-02105-2", "10.1186/s12934-023-02105-2",
         "10.1186/s12934-023-02105-2", "10.1002/lol2.10337",
         "10.1002/lol2.10337", "10.1002/lob.10606",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "MicrobeAnnotator", "MicrobeAnnotator",
         "10.3389/fmicb.2023.1086021", "10.3389/fmicb.2023.1086021",
         "10.3389/fmicb.2023.1086021", "MicrobeAnnotator",
         "MicrobeAnnotator", "10.3390/ijms24043950",
         "10.3390/ijms24043950", "10.3390/ijms24043950", 
         "10.3390/ijms24043950", "10.3390/ijms24043950",
         "10.3390/ijms24043950", "10.3390/ijms24043950",
         "10.3390/ijms24043950", "10.3390/ijms24043950",
         "10.3390/ijms24043950", "10.3390/ijms24043950",
         "10.3390/ijms24043950", "10.3390/ijms24043950",
         "10.3390/ijms24043950", "10.3390/ijms24076253",
         "10.3390/ijms24076253", "10.3390/ijms24076253",
         "10.3390/v15091901", "10.3390/ijms24043950", 
         "10.3390/pharmaceutics15082078", "10.3390/pharmaceutics15082078",
         "10.3390/ijms24043950", "10.3390/app131810540",
         "10.3390/ijms24043950", "10.3390/microorganisms11081899",
         "10.3390/microorganisms11081899", "10.3390/microorganisms11081899",
         "MicrobeAnnotator", "10.1016/j.scitotenv.2023.162194",
         "10.1016/j.scitotenv.2023.162194", "10.1016/j.scitotenv.2023.162194",
         "10.1016/j.scitotenv.2023.162194", "10.1016/j.scitotenv.2023.162194",
         "10.1016/j.scitotenv.2023.162194", "10.1016/j.scitotenv.2023.162194",
         "10.1016/j.scitotenv.2023.162194", "10.1016/j.scitotenv.2023.162194",
         "10.1016/j.envint.2023.108374", "10.1016/j.scitotenv.2023.162194",
         "10.1016/j.scitotenv.2023.163933", "10.1016/j.scitotenv.2023.163933",
         "10.1016/j.scitotenv.2023.163933", "10.1016/j.scitotenv.2023.162194",
         "10.3389/fcimb.2023.1203663", "10.3389/fcimb.2023.1203663",
         "10.1016/j.scitotenv.2023.162194", "10.1016/j.jhazmat.2023.131823",
         "10.1016/j.jhazmat.2023.131823", "10.1016/j.jhazmat.2023.131823",
         "10.1016/j.jhazmat.2023.131823", "10.1016/j.jhazmat.2023.131823",
         "10.1016/j.jhazmat.2023.132250", "10.1021/acs.analchem.3c02973", 
         "10.1039/D3MH01588B"
         )
target <- c("10.1101/2022.11.28.518280", "10.1038/s41586-023-06893-w", 
            "10.1021/acs.est.3c07840", "10.1021/acs.est.3c07737", 
            "10.1016/j.envint.2023.108374", "10.1016/j.jhazmat.2024.133555",
            "10.1038/s41598-023-47897-w", "10.20517/mrr.2022.21", 
            "10.1038/s41467-023-39931-2", "10.3748/wjg.v29.i45.5945",
            "10.3389/fnut.2023.1342142", "10.5527/wjn.v5.i5.429", 
            "10.1186/s43141-023-00598-3", "10.1080/19490976.2023.2281010",
            "10.20944/preprints202311.1377.v1", "other_thesis", 
            "10.1021/acs.est.3c03114", "other_publication",
            "10.1089/omi.2023.0140", "conference_abstract",
            "10.12688/f1000research.139488.1", "10.1128/aem.01034-23",
            "10.1016/j.chemosphere.2023.140743", "10.1093/bib/bbad439",
            "other_thesis_1", "other_book", "other_book_1",
            "other_misc", "other_book_2", "10.1016/j.syapm.2023.126440",
            "other_thesis_1", "10.1016/j.ecolind.2023.110369",
            "other_thesis_2", "10.1016/j.heliyon.2023.e17652",
            "other_thesis_3", "onference_abstract_2", "10.3390/microbiolres14030089",
            "other_book_3", "10.1099/ijsem.0.006224", 
            "10.1186/s12934-023-02105-2", "10.1094/PHYTOFR-09-23-0124-A",
            "10.3389/fcimb.2023.1277176", "10.1002/lol2.10337",
            "other_misc_1", "10.1002/lob.10606", "10.1002/lob.10608",
            "10.21203/rs.3.rs-2617055/v1", "10.1101/2023.03.02.530749",
            "10.21203/rs.3.rs-2602810/v1", "10.3389/fmicb.2023.1086021",
            "10.1002/imt2.157", "10.1109/BIBM58861.2023.10386063",
            "10.21203/rs.3.rs-3318640/v1", "other_misc_2",
            "10.3390/ijms24043950", "other_misc_3",
            "other_blog", "10.21272/eumj.2023;11(3):241-259",
            "10.4239/wjd.v14.i11.1585", "10.1101/2023.12.19.569832",
            "10.1111/jocd.15927", "10.1016/j.jtv.2023.11.001",
            "other_book_4", "10.3390/biology12091187",
            "other_thesis_4", "other_blog_2", "other_misc_4",
            "other_blog_3", "10.3390/ijms24076253",
            "10.3390/biomedicines12010060", "10.4274/balkanmedj.galenos.2023.2023-10-25",
            "10.3390/v15091901", "10.3390/v15122452",
            "10.3390/pharmaceutics15082078", "10.3390/ma16247550",
            "10.3390/antibiotics12121698", "10.3390/app131810540",
            "10.20944/preprints202310.1337.v1", "10.3390/microorganisms11081899",
            "10.3390/cimb46010020", "10.3389/fmicb.2023.1306192",
            "other_blog_4", "10.1016/j.scitotenv.2023.162194",
            "10.3390/antibiotics12101509", "10.3390/agriculture13112143",
            "10.1016/j.envpol.2024.123293", "10.1016/j.biortech.2023.129557",
            "10.1016/j.jhazmat.2023.132161", "10.1016/j.jwpe.2023.104489",
            "10.1016/j.scitotenv.2024.170080", "10.1016/j.envpol.2023.122643",
            "10.1016/j.envint.2023.108374", "10.1016/j.jhazmat.2024.133555",
            "10.1016/j.scitotenv.2023.163933", "10.3390/vetsci10110630",
            "10.1016/j.scitotenv.2023.169794", "10.1016/j.jhazmat.2024.133555",
            "10.3389/fcimb.2023.1203663", "10.26655/JMCHEMSCI.2024.2.1",
            "10.1155/2023/4702866", "10.1016/j.jhazmat.2023.131823",
            "10.1016/j.aca.2024.342240", "10.1016/j.snb.2023.135251",
            "10.1016/j.jhazmat.2023.132250", "10.1021/acs.analchem.3c02973",
            "10.1039/D3MH01588B", "10.1016/j.jhazmat.2023.133161", 
            "10.1016/j.jhazmat.2024.133494", "10.1016/j.microc.2023.109874"
            )
networkData <- data.frame(src, target, stringsAsFactors = FALSE)

nodes <- data.frame(name = unique(c(src, target)), stringsAsFactors = FALSE)
nodes$id <- 0:(nrow(nodes) - 1)


# create a data frame of the edges that uses id 0:9 instead of their names
edges <- networkData %>%
  left_join(nodes, by = c("src" = "name")) %>%
  select(-src) %>%
  rename(source = id) %>%
  left_join(nodes, by = c("target" = "name")) %>%
  select(-target) %>%
  rename(target = id)

edges$width <- 1

# make a grouping variable that will match to colours
nodes <- nodes %>%
  mutate(group = ifelse(name=="MicrobeAnnotator", "origin", "citation"))

ColourScale <- 'd3.scaleOrdinal()
            .domain(["origin", "citation"])
           .range(["#FF6900", "#694489"]);'

forceNetwork(Links = edges, Nodes = nodes, 
             Source = "source",
             Target = "target",
             NodeID ="name",
             Group = "group",
             Value = "width",
             opacity = 0.9,
             zoom = TRUE,
             colourScale = JS(ColourScale))

