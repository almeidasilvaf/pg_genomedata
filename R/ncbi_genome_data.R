
#----Setup----------------------------------------------------------------------
library(tidyverse)
library(rvest)
library(taxize)
library(here)
source(here("R", "utils.R"))

# Extract from NCBI genomes

#----Get NCBI summary table-----------------------------------------------------
viridiplantae <- ncbi_datasets_taxon(taxon = "Viridiplantae")$assemblies

#----Create variables-----------------------------------------------------------
# Species names
species <- data.frame(
    Name = viridiplantae$assembly$org$sci_name,
    Species = viridiplantae$assembly$org$tax_id
)

# Taxonomy for each species
taxize_options(ncbi_sleep = 1)
tax <- get_family(viridiplantae$assembly$org$tax_id)

# Species and family
taxonomy <- full_join(species, tax) %>%
    select(Name, Family, Order, Class, Phylum) %>%
    rename(Species = Name) %>%
    drop_na()

# Technical details
buscos <- viridiplantae$assembly$annotation_metadata$busco$complete
ngenomes <- viridiplantae$assembly$org$assembly_counts$subtree
genome_size <- as.numeric(viridiplantae$assembly$estimated_size) / 10^6

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
assembly_list <- chunk2(viridiplantae$assembly$assembly_accession, 5)
tech1 <- get_seq_tech(assembly_list[[1]], verbose = TRUE)
tech2 <- get_seq_tech(assembly_list[[2]], verbose = TRUE)
tech3 <- get_seq_tech(assembly_list[[3]], verbose = TRUE)
tech4 <- get_seq_tech(assembly_list[[4]], verbose = TRUE)
tech5 <- get_seq_tech(assembly_list[[5]], verbose = TRUE)
technology <- c(tech1, tech2, tech3, tech4, tech5)

technical_df <- data.frame(
    Species = species$Name,
    BUSCOs = buscos,
    N_genomes = ngenomes,
    genome_size = round(genome_size, 2),
    technology = technology
)

# Ploidy level
ploidy <- readr::read_tsv(here::here("data", "ploidy_cvalues_db.tsv"),
                          show_col_types = FALSE) %>%
    mutate(species = paste(Genus, Species, Subspecies, sep = " ")) %>%
    select(species, `Ploidy Level (x)`) %>%
    rename(Species = species, Ploidy = `Ploidy Level (x)`) %>%
    mutate(Species = str_replace_all(Species, " NA", "")) %>%
    filter(Ploidy != "-") %>%
    mutate(Ploidy = as.numeric(Ploidy))


#----Create final data frame of genome info-------------------------------------
final_df <- merge(taxonomy, technical_df)
final_df <- merge(final_df, ploidy, all.x=TRUE)

final_df <- final_df %>%
    dplyr::distinct(., .keep_all = TRUE)

final_dfl <- split(final_df, final_df$Species)
final_dfl <- lapply(final_dfl, function(x) {
    if(nrow(x) == 1) { # only one entry per species
        df <- x
    } else {
        idx_maxbusco <- which.max(x$BUSCOs)
        if(length(idx_maxbusco) == 0) { # none of the genomes have BUSCO info
            idx_maxsize <- which.max(x$genome_size)
            df <- x[idx_maxsize, ]
        } else {
            df <- x[idx_maxbusco, ]
        }
    }
    return(df)
})
final_df <- Reduce(rbind, final_dfl)

# If BUSCO info is not available here and it is in the table from
# the Nature Plants paper, get it
np <- readr::read_tsv(here("data", "nature_plants_data.txt"),
                      show_col_types = FALSE) %>%
    mutate(Species = str_c(Genus, species, sep = " "),
           NP_BUSCOs = `BUSCO % complete (combined)`) %>%
    select(Species, NP_BUSCOs) %>%
    mutate(NP_BUSCOs = str_replace_all(NP_BUSCOs, ",", ".")) %>%
    mutate(NP_BUSCOs = str_replace_all(NP_BUSCOs, "%", "")) %>%
    mutate(NP_BUSCOs = as.numeric(NP_BUSCOs)) %>%
    full_join(., final_df, by = "Species") %>%
    mutate(BUSCOs = BUSCOs * 100)

max_buscos <- pmax(np$NP_BUSCOs, np$BUSCOs, na.rm = TRUE)

final_genomes <- np %>%
    mutate(BUSCO = round(max_buscos, 2)) %>%
    select(Phylum, Class, Order, Family, Species, N_genomes, BUSCO,
           genome_size, technology, Ploidy)

# Remove entries that only contain species info and nothing else
delete <- which(rowSums(is.na(final_genomes)) >= 9)
final_genomes <- final_genomes[-delete, ]

readr::write_tsv(
    final_genomes,
    file = here("data", "genome_info_table.tsv")
)
