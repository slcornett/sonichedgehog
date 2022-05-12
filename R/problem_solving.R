#preliminary packages to build {sonichedgehog}----
require(devtools) # to make package
require(usethis) # to make package
require(roxygen2) # to make package
require(withr) # to make package
# for our package
require(tidyverse)
require(ape) #phylo
require(ggplot2)
require(phytools) # plot tree
require(BiocManager) # to get packages from bioconductor (ie package not on CRAN)
require(ggtree) # extension of ggplot to make phylo trees
require(msa) # multiple sequence alignment

# for installing packages from bioconductor
BiocManager::install("ggtree") # install ggTree package. will give two prompts, answer "all" and "yes" to them respectively
BiocManager::install("msa") # install "msa" package.  will give two prompts, answer "all" and "yes" to them respectively


#sonic.f <- "https://raw.githubusercontent.com/slcornett/sonichedgehog/master/data/sonic.csv"

#Shh_fasta.f <- "https://raw.githubusercontent.com/slcornett/sonichedgehog/master/data/Shh_fasta.fasta"

# adding sonic to dataset
sonic.f <- "https://raw.githubusercontent.com/slcornett/sonichedgehog/master/data/sonic.csv"
usethis::use_data(sonic.f, overwrite = TRUE)
sonic <- readr::read.csv(sonic.f, col_names = TRUE)
usethis::use_data(sonic, overwrite = TRUE)
# adding Shh_fasta to dataset
Shh_fasta.f <- "https://raw.githubusercontent.com/slcornett/sonichedgehog/master/data/Shh_fasta.fasta"
usethis::use_data(Shh_fasta.f, overwrite = TRUE)
Shh_fasta <- Biostrings::readAAStringSet(Shh_fasta.f, format = "fasta", use.names = TRUE)
usethis::use_data(Shh_fasta, overwrite = TRUE)


sonic #prints the included sonic dataset

#ldata ready to load also included in the package
sonic.f # sonic the raw file git link (.f)

#use the function to read in the data and create the r dataset
sonic.csv <- read_file(sonic.f)
head(sonic.csv) # same as sonic




usethis::use_package("readr")

