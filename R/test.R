args <- commandArgs(TRUE)

values <- NULL
namez <- NULL
for(a in args){
  values <- c(values, as.numeric(strsplit(a,"=")[[1]][2]))
  namez <- c(namez, strsplit(a,"=")[[1]][1])
}
names(values) <- namez

library(base)
library(tmvtnorm)
library(glasso)
library(simone)
library(Matrix)
library(snow)

setwd("C:/Users/Pariya/Desktop/crucialCommands/Codes/Code")
#setwd("~/Github/PariyaCode")
data <- read.csv("tests/test_data.txt",sep="\t")


setwd("C:/Users/Pariya/Desktop/clean")
source("codes/gnetMultiCores.R")
source("codes/calculate.cutoffs.R")
source("codes/calculate.lower.upper.R")
source("codes/calculate.R.R")

results <- gnet(data, 10, 100, values["iter"], 100, n.cores=3)
