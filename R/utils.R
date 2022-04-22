
#--- 1) To be used in ncbi_genome_data.R---------------------------------------

# Search NCBI genomes by taxon
ncbi_datasets_taxon <- function(taxon = NULL) {

    outfile <- tempfile(fileext = ".json")
    args <- c("summary genome taxon", taxon, " > ", outfile)

    datasets <- system2("datasets", args = args)

    parsed_data <- jsonlite::fromJSON(outfile)
    unlink(outfile)
    return(parsed_data)
}

# Get BioSample location from the output of \code{ncbi_datasets_taxon()}
get_biosample_locations <- function(outlist = NULL) {
    biosamples <- outlist$assemblies$assembly$biosample
    biosamples <- biosamples[!duplicated(biosamples$accession), ]
    locations <- unlist(lapply(biosamples$attributes, function(x) {
        loc <- x$value[x$name == "geo_loc_name"]
        if(length(loc) == 0) {
            loc <- NA
        } else if(length(loc) > 1) {
            loc <- paste(loc, collapse = ", ")
        }
        return(loc)
    }))
    final_df <- data.frame(
        Biosamples = biosamples$accession,
        Location = locations
    )
    return(final_df)
}

# Get family for NCBI species IDS
get_family <- function(species_ids = NULL) {
    uspecies <- unique(species_ids)
    classification <- Reduce(rbind, lapply(uspecies, function(x) {
        fam <- taxize::classification(x, db = "ncbi")[[1]]
        family <- fam$name[fam$rank == "family"]
        order <- fam$name[fam$rank == "order"]
        class <- fam$name[fam$rank == "class"]
        phylum <- fam$name[fam$rank == "phylum"]

        if(length(family) == 0) { family <- NA }
        if(length(order) == 0) { order <- NA }
        if(length(class) == 0) { class <- NA }
        if(length(phylum) == 0) { phylum <- NA }

        fam_df <- data.frame(
            Species = x,
            Family = family,
            Order = order,
            Class = class,
            Phylum = phylum
        )
        return(fam_df)
    }))

    return(classification)
}

# Get sequencing technologies for assemblies
get_seq_tech <- function(assembly_acc = NULL, verbose = FALSE) {
    tech <- unlist(lapply(assembly_acc, function(x) {
        if(verbose) {
            message("Working on genome accession ", x)
        }
        url <- paste0("https://www.ncbi.nlm.nih.gov/assembly/", x, "/")
        page <- rvest::read_html(url)
        name <- rvest::html_text(rvest::html_elements(page, "#summary dl > dt"))
        value <- rvest::html_text(rvest::html_elements(page, "#summary dd"))
        diff <- length(value) - length(name)
        if(diff > 0) {
            remove <- seq(3, 3 + (diff-1))
            value <- value[-c(remove)]
        }
        df <- data.frame(name, value)

        stech <- df$value[grepl("Sequencing technology", df$name)]
        if(length(stech) == 0) {
            stech <- NA
        }
        return(stech)
    }))
    return(tech)
}
