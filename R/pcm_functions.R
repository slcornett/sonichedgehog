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
#' @param # one line for each parameter in our function
#' @keywords # to aid in searching functions
#' @export # makes the function available for others to use when your package is loaded
#' @examples - sample code
f_file <- function(file){
  #for reading multiple AA sequences from msa package
  fast <- Biostrings::readAAStringSet(file, format = "fasta", use.names = TRUE) # format: biostrings, AAString set
  return(fast)
}

# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file ----
#' @title AlignAA
#' @description Amino Acid multiple sequence alignment using the MUSCLE algorithm of fasta formatted data
#' @param file a fasta file assigned to a vector in your environment
#' @keywords multiple sequence alignment, fasta, peptide, AA alignment
#' @export # makes the function available for others to use when your package is loaded
#' @examples using the built-in Shh_fasta data of sonic hedgehog signal molecule orthologs
#' file <- "link to raw data on git hub"
#' pcms_AlignAA(file) # the output of this will be the sequences aligned by the MUSCLE algorithm
AlignAA <- function(file){
  class(file) == "character"
  #for reading multiple AA sequences from biostrings package
  fas <- Biostrings::readAAStringSet(file, format = "fasta", use.names = TRUE)
  # align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package
  fas_msa <- fas %>% msa::msa(method = c("Muscle"), type = "protein", order=c("aligned", "input"))
  print(fas_msa, showConsensus = TRUE, show = "complete")
  return(fas_msa)
}

# Building a function to calculate the phylogenetic distance between AAs, then plot a tree ----
#' @title treeAA
#' @description
#' @param # one line for each parameter in our function
#' @keywords # to aid in searching functions
#' @export # makes the function available for others to use when your package is loaded
#' @examples # sample code
treeAA <- function(x){ #the input to this function is the out put of the above function
  class(x) == "MsaAAMultipleAlignment" # the only class these functions work on
  #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
  dist_fas_msa <- x %>% ape::as.AAbin(show.aa = TRUE, check.names = TRUE) %>% # AAbin storage
    ape::dist.aa(pairwise.deletion = TRUE, scaled = TRUE) # distance AA
  print(dist_fas_msa) # prints this but not the other two
  #neighbor joining method
  tree <- ape::nj(dist_fas_msa) # new tree by neighbor joining method
  return(tree) #  need to be able to extract this from the vector so can modify outgroup
  ggt <- ggtree::ggtree(tree, #the new tree as data
                        cex = 1, # branch line thickness
                        aes(color = branch.length), # the color factor is branch length
                        layout = "roundrect") + # style
    scale_color_continuous(high = 'green', low = 'darkgreen') + # the colors used for indicating branch length
    theme(legend.position = "left") + # position of the legend
    geom_tiplab(align = FALSE, size = 3) + # tips not aligned, font size = 3
    geom_treescale(x= 0, y = 0, color = "black", fontsize = 3) # the tree scale (number of substitutions)
  print(ggt) # not printed
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
