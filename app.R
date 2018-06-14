source("https://bioconductor.org/biocLite.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, DT, limma, 
               psych, ggplot2,
               genefilter,reshape2,
               stringr,missForest)
runApp(".")