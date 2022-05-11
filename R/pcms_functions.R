# Phylogenetic Comparative Analysis
#
# This rscript contains the functions we created for our PCMs
# package.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#Create a function that reads in the file information----
#' @title read_file
#' @description
#' @param # one line for each parameter in our function
#' @keywords # to aid in searching functions
#' @export # makes the function available for others to use when your package is loaded
#' @examples # sample code

read_file <- function(f_name){
  data <-readr::read_csv(f_name, col_names = TRUE)
  return(data)
}

# fxn to read in multiple AA sequences----
#' @title f_file
#' @description
#' @param - one line for each parameter in our function
#' @keywords - to aid in searching functions
#' @export - makes the function available for others to use when your package is loaded
#' @examples - sample code
f_file <- function(fast_file){
  #for reading multiple AA sequences from msa package
  fast <- Biostrings::readAAStringSet(fast_file, format = "fasta", use.names = TRUE) # format: biostrings, AAString set
  return(fast)
}

# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file ----
#' @title read_file
#' @description
#' @param # one line for each parameter in our function
#' @keywords # to aid in searching functions
#' @export # makes the function available for others to use when your package is loaded
#' @examples # sample code
pcms_AlignAA <- function(file){
  #for reading multiple AA sequences from msa package
  fas <- Biostrings::readAAStringSet(file, format = "fasta", use.names = TRUE)
  # align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package
  fas_msa <- fas %>% msa::msa(method = c("Muscle"), type = "protein", order=c("aligned", "input"))
  return(fas_msa)
}

# Building a function to calculate the phylogenetic distance between AAs, then plot a tree ----
#' @title read_file
#' @description
#' @param # one line for each parameter in our function
#' @keywords # to aid in searching functions
#' @export # makes the function available for others to use when your package is loaded
#' @examples # sample code
pcms_treeAA <- function(x){#the input to this function is the out put of the above function
  class(x) == "MsaAAMultipleAlignment"
  #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
  dist_fas_msa <- x %>% ape::as.AAbin(show.aa = TRUE, check.names = TRUE) %>% # AAbin storage
    ape::dist.aa() # distance AA
  #neighbor joining method
  tree <- ape::nj(dist_fas_msa) # new tree
  ggt <- ggtree::ggtree(tree, cex = 1, aes(color=branch.length)) +
    scale_color_continuous(high='green',low='blue') +
    geom_tiplab(align = FALSE, size = 2) +
    geom_treescale(y = 0, color = "coral4", fontsize = 4)
  return(dist_fas_msa)
  return(tree)
  return(ggt)
}

#Plot the maximum likelihood----
#' @title read_file
#' @description
#' @param # one line for each parameter in our function
#' @keywords # to aid in searching functions
#' @export # makes the function available for others to use when your package is loaded
#' @examples # sample code
plot_Tree <- function(tree){
  phytools::plotTree(tree, fsize=0.8,lwd=1,offset=3) #Plot the ML, set font size, create space so nodes aren't on top of each other
}
