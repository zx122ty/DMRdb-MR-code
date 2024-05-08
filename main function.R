#  main function
setwd('...')

library(data.table)
library(TwoSampleMR)
library(foreach)
library(doParallel)
library(RadialMR)
library(ggplot2)
library(MRPRESSO)

source('MR_analysis.R')
source('blank_plot.R')
source('plot_all.R')

outcome_file_adress = '...'
outcome_id = '...'
exposure_sample = '...'

outcome_file <- fread(outcome_file_adress, select =  c("SNP", "beta", "eaf","se","pval","effect_allele","other_allele"))
exposure_file_address = '...'   
exposure_file = list.files(exposure_file_address)

parallel::mclapply(exposure_file, function(x) MR_analysis(x, outcome_id, outcome_file, exposure_file_address), mc.cores=4)
