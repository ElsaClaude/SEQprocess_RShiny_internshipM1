#Change lib path if necessary with : .libPaths()
.libPaths("/net/travail/elclaude/R/x86_64-pc-linux-gnu-library/3.5")

#source codes
source("https://bioconductor.org/biocLite.R")

#GenomeInfoDbData package
biocLite("GenomeInfoDbData")

#DelayedArray package
biocLite("DelayedArray")

#devtools packages
install.packages("devtools")

#call devtools library
library(devtools)

#install packages through github ?
install_github("omicsCore/SEQprocess")

#load packages and documentation
library(SEQprocess) # Loads the package
library(help="SEQprocess") # Lists package info

#install libraries
install.packages("sequenza")
install.packages("data.table")
install.packages("fastqcr")
install.packages("pander", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
install.packages("knitr")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("reshape2")