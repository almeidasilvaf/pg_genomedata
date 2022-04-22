
# Code description

## ncbi_genome_data.R

This file contains the code used to create the genome information table that is 
displayed in the **Overview** tab. 

**RUN:** To create the table from scratch and to update it.

## scrape_ploidy_information.R

This file contains code to scrape the Plant DNA C-values database using 
RSelenium and store the output in a .tsv file.

**RUN:** To create the table from scratch and after every update on 
the database. There is no need to update this table constantly, as the database
itself is now updated constantly.
