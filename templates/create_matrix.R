#!/usr/bin/env Rscript

## Required libraries
library(stringr)

## Get the list of taxon files
taxon_files <- system('find . -iname "*fasta.taxon"', intern = TRUE)

## Dynamically create vectors/list of taxons of each sample
the_list <- vector()
for(i in taxon_files){
    the_list <- c(the_list, paste(str_replace(basename(i), '_fasta.taxon', '')))
    assign(paste(str_replace(basename(i), '_fasta.taxon', '')), scan(i))
}

## Get the list of unique taxons for all datasets
taxons <- unique(get(the_list))

## Initiate a matrix to populate - taxons on rows and samples on columns
the_table <- as.data.frame(taxons)

## Populate the matrix - 0 is present in sample and 1 is absent
for(i in the_list){
    the_table[,paste0(i)] <- as.numeric(the_table$taxon %in% get(i))
}

## Use the taxonomy numbers as names
rownames(the_table) <- the_table$taxons
the_table <- the_table[,-1]
the_table[is.na(the_table)] <- 0

## Remove 0 (unknown), 1 (root), 9606 (homo sapiens), 374840 (phiX).
##the_table <- the_table[!(rownames(the_table) %in% c(0,1,9606,374840)),]

## ## Remove taxons that are found throughout the samples - not useful
the_table <- the_table[rowSums(the_table) != length(colnames(the_table)),]

## Order table by column names
the_table <- the_table[,order(names(the_table))]

## Add the names of the taxonomies
names <- read.table('names_table.dmp', header = FALSE, sep = '\t')
colnames(names) <- c("taxid","taxname")

## Add last column "Taxonomy" with the names of the taxids
the_table$Taxonomy <- names$taxname[match(rownames(the_table),names$taxid)]

## Write the matrix to a file
write.csv(the_table, file="upset_data.csv")
