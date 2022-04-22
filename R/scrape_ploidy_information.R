
#----Setup----------------------------------------------------------------------
library(tidyverse)
library(rvest)


# Get ploidy by scraping the Plant DNA C-Values Database
get_ploidy <- function(species = NULL, port = 4567L) {

    genus <- vapply(strsplit(species, " "), `[`, 1, FUN.VALUE = character(1))
    genus <- unique(genus)
    rD <- RSelenium::rsDriver(browser = "firefox", verbose = FALSE,
                              port = port)
    remDr <- rD[["client"]]

    ploidy <- lapply(genus, function(x) {
        message("Working on genus ", x)
        # Go to Plant DNA C-values Database
        remDr$navigate("https://cvalues.science.kew.org/search")

        # Click on "Ploidy level" button
        remDr$findElement(
            using = "id", value = "ploidyCheckbox"
        )$clickElement()

        # Fill in "Genus" text box with the content of object `genus`
        remDr$findElement(
            using = "id", value = "genusField"
        )$sendKeysToElement(list(x))

        # Click on "Search" button
        remDr$findElements("class name", "btn")[[1]]$clickElement()
        Sys.sleep(1.5)

        # Get HTML table
        html <- remDr$getPageSource()[[1]]
        table <- rvest::read_html(html) %>%
            rvest::html_table() %>%
            purrr::pluck(2)

        # How many pages (tabs) are there?
        pages <- rvest::read_html(html) %>%
            rvest::html_elements(".pagination") %>%
            pluck(1) %>%
            rvest::html_elements("li") %>%
            html_text2()

        if(nrow(table) == 0) {
            message("No ploidy information for genus ", x)
            final_table <- NULL
        } else if(length(pages) > 1) {
            tables <- lapply(seq_along(pages)[-1], function(p) {
                css_sel <- paste0(
                    "#tabResults > div:nth-child(2) > div:nth-child(1) > " ,
                    "nav:nth-child(1) > ul:nth-child(1) > ",
                    "li:nth-child(", p, ") > a:nth-child(1)"
                )
                remDr$findElement("css selector", css_sel)$clickElement()
                Sys.sleep(1.5)
                tablen <- rvest::read_html(html) %>%
                    rvest::html_table() %>%
                    purrr::pluck(2)

                return(tablen)
            })
            tables <- Reduce(rbind, tables)
            final_table <- rbind(table, tables)
        } else {
            final_table <- table
        }
        return(final_table)
    })
    ploidy_df <- Reduce(rbind, ploidy)
    return(ploidy_df)
}

#----Get NCBI summary table-----------------------------------------------------
viridiplantae <- ncbi_datasets_taxon(taxon = "Viridiplantae")$assemblies


species <- viridiplantae$assembly$org$sci_name

sp <- unique(species)
sp <- sp[order(sp)]

# Repeat process in groups of 1000 projects to avoid losing everything for
# connection issues
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

sps <- chunk2(sp, 10)

ploidy1 <- get_ploidy(sps[[1]], port = 4563L)
ploidy2 <- get_ploidy(sps[[2]], port = 4550L)
ploidy3 <- get_ploidy(sps[[3]], port = 4551L)
ploidy4 <- get_ploidy(sps[[4]], port = 4553L)
ploidy5 <- get_ploidy(sps[[5]], port = 4554L)
ploidy6 <- get_ploidy(sps[[6]], port = 4555L)
ploidy7 <- get_ploidy(sps[[7]], port = 4556L)
ploidy8 <- get_ploidy(sps[[8]], port = 4557L)
ploidy9 <- get_ploidy(sps[[9]], port = 4558L)
ploidy10 <- get_ploidy(sps[[10]], port = 4559L)

ploidy_cvalues <- bind_rows(
    ploidy1, ploidy2, ploidy3, ploidy4, ploidy5,
    ploidy6, ploidy7, ploidy8, ploidy9, ploidy10
)

readr::write_tsv(
    ploidy_cvalues,
    file = here::here("data", "ploidy_cvalues_db.tsv")
)
