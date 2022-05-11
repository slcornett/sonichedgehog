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


sonic.csv <- "https://raw.githubusercontent.com/slcornett/sonichedgehog/master/data/sonic.csv"

Shh_fast.f <- "https://raw.githubusercontent.com/slcornett/sonichedgehog/master/data/Shh_fasta.fasta"
