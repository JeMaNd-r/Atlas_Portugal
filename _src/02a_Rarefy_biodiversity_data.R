#- - - - - - - - - - - - - - - - -#
#        Bioinformatics           #
#                                 #
#     author: Romy Zeiss          #
#       date: 2025-04-24          #
#- - - - - - - - - - - - - - - - -#

library(tidyverse)
library(phyloseq)

library("pr2database") # to assign Protist taxa

set.seed(48293)

#- - - - - - - - - - - - - - - - -
## 16s data ####
#- - - - - - - - - - - - - - - - -
Taxon_name <- "Bacteria"

# load raw OTU tables
raw_16s <- read_delim(file="_data/raw_otutable_16s.txt")
raw_18s <- read_delim(file="_data/raw_otutable_euk.txt")

# split taxonomy into multiple columns
raw_16s <- raw_16s %>% separate_wider_regex( col = taxonomy,
                                                 patterns = c(D = "d_+\\w+",
                                                              "; ",
                                                              P = "p_+\\w+",
                                                              "; ",
                                                              C = "c_+\\w+",
                                                              "; ",
                                                              O = "o_+\\w+",
                                                              "; ",
                                                              F = "f_+\\w+",
                                                              "; ",
                                                              G = "g_+\\w+",
                                                              "; ",
                                                              S = "s_+\\w+"), 
                                                 too_few = "align_start")  %>%  
  mutate(D = gsub(x = D, pattern = "d__", replacement = ""),
         P = gsub(x = P, pattern = "p__", replacement = ""),
         C = gsub(x = C, pattern = "c__", replacement = ""),
         O = gsub(x = O, pattern = "o__", replacement = ""),
         F = gsub(x = F, pattern = "f__", replacement = ""),
         G = gsub(x = G, pattern = "g__", replacement = ""),
         S = gsub(x = S, pattern = "s__", replacement = ""))
unique(raw_16s %>% dplyr::select(D:S))

# exclude Archaea
raw_16s <- raw_16s %>% #only 1 OTU
  filter(D != "Archaea")

# define objects for phyloseq functions
otu_16s <- raw_16s %>% as.data.frame()
rownames(otu_16s) <- otu_16s$`#OTU ID`
otu_16s <- otu_16s %>% dplyr::select(!(c(D:G,S, `#OTU ID`))) %>% as.matrix()

tax_16s <- raw_16s %>% as.data.frame()
rownames(tax_16s) <- tax_16s$`#OTU ID`
tax_16s <- tax_16s %>% dplyr::select(c(D:G,S)) %>% as.matrix()

OTU <- otu_table(otu_16s, taxa_are_rows = TRUE)
TAX <- tax_table(tax_16s)

physeq_16s <- phyloseq(OTU, TAX)
physeq_16s # 33970 taxa and 407 samples

# # rarefaction on all OTUs
# physeq_16s_rar <- rarefy_even_depth(physeq_16s, sample.size = 5000, trimOTUs = T)
# # sample.size=1000:  7 samples, 12778 OTUs
# # sample.size=2000: 10 samples,  9310 OTUs
# # sample.size=3000: 13 samples,  7498 OTUs
# # sample.size=4000: 14 samples,  6410 OTUs
# # sample.size=5000: 14 samples,  5561 OTUs
# # sample.size=6000: 14 samples,  4964 OTUs

# rarefaction at genus-level OTUs
physeq_16s <- tax_glom(physeq_16s, taxrank = "G")
physeq_16s_rar <- rarefy_even_depth(physeq_16s, sample.size = 1000, trimOTUs = T)
# sample.size= 500:  6 samples, 89 OTUs
# sample.size=1000:  7 samples, 44 OTUs
# sample.size=2000: 12 samples, 20 OTUs
# sample.size=5000: 14 samples,  7 OTUs

# Extract OTU table
rar_otu <- as(otu_table(physeq_16s_rar), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Extract taxonomy table
rar_tax <- tax_table(physeq_16s_rar) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Merge by OTU_ID
rar_all <- left_join(rar_otu, rar_tax, by = "OTU_ID")

# save
write_csv(rar_all, paste0("_data/OTUtable_rarefied_16s_", Taxon_name, ".csv"))


#- - - - - - - - - - - - - - - - -
## 18s data (Eukaryotes) ####
#- - - - - - - - - - - - - - - - -

# load raw OTU tables
raw_18s <- read_delim(file="_data/raw_otutable_euk.txt")

# split taxonomy into multiple columns
raw_18s <- raw_18s %>% separate_wider_regex( col = taxonomy,
                                             patterns = c(D = "d_+\\w+",
                                                          "; ",
                                                          P = "p_+\\w+",
                                                          "; ",
                                                          C = "c_+\\w+",
                                                          "; ",
                                                          O = "o_+\\w+",
                                                          "; ",
                                                          F = "f_+\\w+",
                                                          "; ",
                                                          G = "g_+\\w+",
                                                          "; ",
                                                          S = "s_+\\w+"), 
                                             too_few = "align_start")  %>%  
  mutate(D = gsub(x = D, pattern = "d__", replacement = ""),
         P = gsub(x = P, pattern = "p__", replacement = ""),
         C = gsub(x = C, pattern = "c__", replacement = ""),
         O = gsub(x = O, pattern = "o__", replacement = ""),
         F = gsub(x = F, pattern = "f__", replacement = ""),
         G = gsub(x = G, pattern = "g__", replacement = ""),
         S = gsub(x = S, pattern = "s__", replacement = ""))
unique(raw_18s %>% dplyr::select(D:S))

#- - - - - - - - - - - - - - - - -
### Fungi ####
Taxon_name <- "Fungi"

clean_18s <- raw_18s %>% filter(str_detect(P, "mycota"))

# define objects for phyloseq functions
otu_16s <- clean_18s %>% as.data.frame()
rownames(otu_16s) <- otu_16s$`#OTU ID`
otu_16s <- otu_16s %>% dplyr::select(!(c(D:G,S, `#OTU ID`))) %>% as.matrix()

tax_16s <- clean_18s %>% as.data.frame()
rownames(tax_16s) <- tax_16s$`#OTU ID`
tax_16s <- tax_16s %>% dplyr::select(c(D:G,S)) %>% as.matrix()

OTU <- otu_table(otu_16s, taxa_are_rows = TRUE)
TAX <- tax_table(tax_16s)

physeq_18s <- phyloseq(OTU, TAX) 
physeq_18s # 2358 taxa and 407 samples

# rarefaction at genus-level OTUs
physeq_18s <- tax_glom(physeq_18s, taxrank = "G")
physeq_18s_rar <- rarefy_even_depth(physeq_18s, sample.size = 1000, trimOTUs = T)
# sample.size= 500:  30 samples,  4 OTUs
## sample.size=1000:  40 samples,  2 OTUs
# sample.size=2000: 100 samples,  2 OTUs
# sample.size=5000: 343 samples, 13 OTUs

# Extract OTU table
rar_otu <- as(otu_table(physeq_18s_rar), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") #note: one of the OTUs is choosen to represent merged species...

# Extract taxonomy table
rar_tax <- tax_table(physeq_18s_rar) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Merge by OTU_ID
rar_all <- left_join(rar_otu, rar_tax, by = "OTU_ID")
rar_all

# save
write_csv(rar_all, paste0("_data/OTUtable_rarefied_18s_", Taxon_name, ".csv"))


#- - - - - - - - - - - - - - - - -
### Protists ####
# filter per taxon group
Taxon_name <- "Protists"

# get list of protists using PR2Database
pr_data <- pr2database::pr2_database()

# Note: Eukaroyta_X = Eukaryotes that could not be assigned to any taxon group

# Vector of protist supergroups: 
# SAR groups: Stramenopiles, Alveolata, Rhizaria
# no Obazoa (Ophistokonta) as they include fungi & animals: 
# Excavata maybe excluding Algea
# Amoebozoa, but also some not-protist taxa
protist_supergroups <- c("Alveolata", "Amoebozoa", "Excavata", "Obazoa", "Rhizaria", "Stramenopiles", "TSAR")

# Keep only rows where the domain is "Eukaryota" and supergroup matches Protist groups
protist_data <- pr_data %>%
  filter((domain == "Eukaryota" & supergroup %in% protist_supergroups) | domain == "Eukaryota:apic") %>%
  filter(domain != "TSAR:apic")

# Additionally handle Opisthokonta by excluding known non-protists (Fungi and Metazoa)
protist_data <- protist_data %>%
  filter(!(supergroup == "Obazoa" & subdivision %in% c("Fungi", "Metazoa")))

protist_data <- protist_data %>%
  dplyr::select(domain, supergroup, division, subdivision, class, order, family, genus, species, sequence) %>%
  unique()
#View(unique(protist_data %>% dplyr::select(domain:species)))

# remove "wrong" genus and species
taxa_to_remove <- c(
  "metagenome",
  "uncultured"
)

clean_18s <- raw_18s %>%
  filter(!(G %in% taxa_to_remove)) %>% #12,829
  # match with PR2 and account for different taxonomy orders
  filter(P %in% c(unique(protist_data$domain),
                    unique(protist_data$supergroup),
                    unique(protist_data$division),
                    unique(protist_data$subdivision),
                    unique(protist_data$class),
                    unique(protist_data$order),
                  "SAR")) #325
unique(clean_18s %>% dplyr::select(D:S))

# # check removed ones #-> FINE
# unique(raw_18s %>% 
#          filter(!(P %in% unique(clean_18s$P))) %>%
#          dplyr::select(D:S))

# fill missing taxa names with next larger one
clean_18s$P[is.na(clean_18s$P)] <- clean_18s$D[is.na(clean_18s$P)]
clean_18s$C[is.na(clean_18s$C)] <- clean_18s$P[is.na(clean_18s$C)]
clean_18s$O[is.na(clean_18s$O)] <- clean_18s$C[is.na(clean_18s$O)]
clean_18s$F[is.na(clean_18s$F)] <- clean_18s$O[is.na(clean_18s$F)]
clean_18s$G[is.na(clean_18s$G)] <- clean_18s$F[is.na(clean_18s$G)]

# define objects for phyloseq functions
otu_18s <- clean_18s %>% as.data.frame()
rownames(otu_18s) <- otu_18s$`#OTU ID`
otu_18s <- otu_18s %>% dplyr::select(!(c(D:S, `#OTU ID`))) %>% as.matrix()

tax_18s <- clean_18s %>% as.data.frame()
rownames(tax_18s) <- tax_18s$`#OTU ID`
tax_18s <- tax_18s %>% dplyr::select(c(D:S)) %>% as.matrix()

OTU <- otu_table(otu_18s, taxa_are_rows = TRUE)
TAX <- tax_table(tax_18s)

physeq_18s <- phyloseq(OTU, TAX)
physeq_18s # 3184 taxa and 407 samples

# rarefaction at genus-level OTUs
physeq_18s <- tax_glom(physeq_18s, taxrank = "G") #159 taxa and 407 samples
physeq_18s_rar <- rarefy_even_depth(physeq_18s, sample.size = 1000, trimOTUs = T, rngseed = 831)
# sample.size= 500:  33 samples,  2 OTUs
## sample.size=1000:  92 samples,  1 OTUs
# sample.size=2000: 179 samples,  0 OTUs
# sample.size=5000: 357 samples, 11 OTUs

# Extract OTU table
rar_otu <- as(otu_table(physeq_18s_rar), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Extract taxonomy table
rar_tax <- tax_table(physeq_18s_rar) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Merge by OTU_ID
rar_all <- left_join(rar_otu, rar_tax, by = "OTU_ID")

# save
write_csv(rar_all, paste0("_data/OTUtable_rarefied_18s_", Taxon_name, ".csv"))
rm(rar_all)

#- - - - - - - - - - - - - - - - -
### Other Eukaryotes ####
#- - - - - - - - - - - - - - - - -
# filter per taxon group
Taxon_name <- "Eukaryotes"

protist_list <- read_csv(paste0("_data/OTUtable_rarefied_18s_Protists.csv"))

clean_18s <- raw_18s %>% 
  filter(!str_detect(P, "mycota")) %>% # no fungi
  filter(!(P %in% unique(protist_list$P))) %>% # no protists
  filter(P != "Nematozoa" & P != "Annelida") # no nematodes (earthworms anyway not present)

# check taxonomic level
#View(clean_18s %>% dplyr::select(D:S) %>% unique() %>% arrange(D, P, C, O, F, G, S))

# remove metagenome and uncultured taxa
clean_18s <- clean_18s %>%
  filter(!(str_detect(G, "metagenome"))) %>%
  filter(!(str_detect(G, "uncultured")))

# define objects for phyloseq functions
otu_18s <- clean_18s %>% as.data.frame()
rownames(otu_18s) <- otu_18s$`#OTU ID`
otu_18s <- otu_18s %>% dplyr::select(!(c(D:S, `#OTU ID`))) %>% as.matrix()

tax_18s <- clean_18s %>% as.data.frame()
rownames(tax_18s) <- tax_18s$`#OTU ID`
tax_18s <- tax_18s %>% dplyr::select(c(D:S)) %>% as.matrix()

OTU <- otu_table(otu_18s, taxa_are_rows = TRUE)
TAX <- tax_table(tax_18s)

physeq_18s <- phyloseq(OTU, TAX)
physeq_18s #1392 taxa and 407 samples

# rarefaction at genus-level OTUs
physeq_18s <- tax_glom(physeq_18s, taxrank = "G") #107 taxa and 407 samples
physeq_18s_rar <- rarefy_even_depth(physeq_18s, sample.size = 1000, trimOTUs = T, rngseed = 831)
# sample.size= 500:  27 samples, 254 OTUs
## sample.size=1000:  57 samples, 197 OTUs
# sample.size=2000: 147 samples, 204 OTUs
# sample.size=5000: 313 samples, 441 OTUs

# Extract OTU table
rar_otu <- as(otu_table(physeq_18s_rar), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Extract taxonomy table
rar_tax <- tax_table(physeq_18s_rar) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID")

# Merge by OTU_ID
rar_all <- left_join(rar_otu, rar_tax, by = "OTU_ID")

# save
write_csv(rar_all, paste0("_data/OTUtable_rarefied_18s_", Taxon_name, ".csv"))

