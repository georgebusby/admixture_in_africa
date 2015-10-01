#####################################################################
## WRITING FUNCTIONS ##
## THIS IS A LITTLE HELPER SCRIPT TO HELP WITH COMPILING FUNCTIONS ##
## TO BE EVENTUALLY USED AS A PACKAGE 

## SOURCE THIS FILE IF WRITING FUNCTIONS ##
## see http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
## THESE TWO LIBRARIES ALLOW THE FUNCTIONS DETAILS TO 
## BE WRITTEN WITHIN THE FUNCTION
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

## write some functions ....
## now document and install ...
setwd("~/repos/popgen/packages_dev/copyselection")
document()
install()
