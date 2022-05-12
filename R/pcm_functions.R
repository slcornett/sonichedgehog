# Phylogenetic Comparative Methods Analysis
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
#' @description This function allows the user to read in the file data
#' @param file
#' @keywords file, read-in
#' @export makes the function available for others to use when your package is loaded
#' @examples # sample code

read_file <- function(sonic.f){
  data <-readr::read_csv(sonic.f, col_names = TRUE)
  return(data)
}

# fxn to read in multiple AA sequences----
#This function was created to read in the fasta data, however it became
# a bit redundant later on in our package.

#@title f_file
#@description
#@param # one line for each parameter in our function
#@keywords # to aid in searching functions
#@export # makes the function available for others to use when your package is loaded
#@examples - sample code
#f_file <- function(Shh_fasta){
  #for reading multiple AA sequences from msa package
  #fast <- Biostrings::readAAStringSet(Shh_fasta, format = "fasta", use.names = TRUE) # format: biostrings, AAString set
 # return(fast)
#}

# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file ----
#' @title AlignAA
#' @description Amino Acid multiple sequence alignment using the MUSCLE algorithm of fasta formatted data
#' @param file a fasta file assigned to a vector in your environment
#' @keywords multiple sequence alignment, fasta, peptide, AA alignment
#' @export # makes the function available for others to use when your package is loaded
#' @examples using the built-in Shh_fasta data of sonic hedgehog signal molecule orthologs
#' file <- "link to raw data on git hub"
#' pcms_AlignAA(file) # the output of this will be the sequences aligned by the MUSCLE algorithm
AlignAA <- function(Shh_fasta){
  class(Shh_fasta) == "character"
  #for reading multiple AA sequences from biostrings package
  fas <- Biostrings::readAAStringSet(Shh_fasta, format = "fasta", use.names = TRUE)
  # align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package
  fas_msa <- fas %>% msa::msa(method = c("Muscle"), type = "protein", order=c("aligned", "input"))
  fas_msa <- fas %>% msa::msa(method = c("Muscle"), type = "protein", order=c("aligned", "input"))
  print(fas_msa, showConsensus = TRUE, show = "complete")
  return(fas_msa)
}

# Building a function to calculate the phylogenetic distance between AAs, then plot a tree ----
#' @title phyloAA
#' @description This function conducts the multiple alignment and distance calculations for phylogenetic data
#' @param msa, formatted multiple alignment data
#' @keywords msa, distance
#' @export # makes the function available for others to use when your package is loaded
#' @examples fas_msa <- "MsaMultipleAlignment" phyloAA(fas_msa) #the output is phylo tree data
phyloAA <- function(fas_msa){ #the input to this function is the out put of the above function
  class(fas_msa) == "MsaAAMultipleAlignment" # the only class these functions work on
  #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
  dist_fas_msa <- fas_msa %>% ape::as.AAbin(show.aa = TRUE, check.names = TRUE) %>% # AAbin storage
    ape::dist.aa(pairwise.deletion = TRUE, scaled = TRUE) # distance AA
  print(dist_fas_msa) # prints this but not the other two
  #neighbor joining method
  tree <- ape::nj(dist_fas_msa) # new tree by neighbor joining method
  return(tree) #  need to be able to extract this from the vector so can modify outgroup
}

# Building a function to plot a tree, using a function ggtree----
#' @title treeAA
#' @description This function plots a phlogenetic tree
#' @param ptree, phylo formatted data
#' @keywords phylogenetric tree
#' @export # makes the function available for others to use when your package is loaded
#' @examples ptree <- "phylo" #prinnts out tree
treeAA <- function(ptree){
  class(ptree) == "phylo" # the only class these functions work on
  ggt <- ggtree::ggtree(ptree, #the new tree as data
    cex = 1, # branch line thickness
    aes(color = branch.length), # the color factor is branch length
    layout = "roundrect") + # style
    ggplot2::scale_color_continuous(high = 'green', low = 'darkgreen') + # the colors used for indicating branch length
    ggplot2::theme(legend.position = "left") + # position of the legend
    ggtree::geom_tiplab(align = FALSE, size = 3) + # tips not aligned, font size = 3
    ggtree::geom_treescale(x= 0, y = 0, color = "black", fontsize = 3) # the tree scale (number of substitutions)
  return(ggt)
}


#Plot the maximum likelihood----
#' @title plot_Tree
#' @description This function uses phyloTools to print out a maximum likelihood using phylogenetic data
#' @param tree
#' @keywords tree, maximum likelihood
#' @export # makes the function available for others to use when your package is loaded
#' @examples plot_Tree <- tree
plot_Tree <- function(tree){
  phytools::plotTree(tree, fsize=0.8,lwd=1,offset=3) #Plot the ML, set font size, create space so nodes aren't on top of each other
}
