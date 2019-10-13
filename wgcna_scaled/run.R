#!/usr/bin/env Rscript
library(logging)
args = commandArgs(trailingOnly=TRUE)

basicConfig()
addHandler(writeToFile, logger="", file = args[3], sep = "")

source("helpers/run_module_pipeline_MA.R")
source("helpers/run_module_pipeline_RNASeq.R")

if(args[4] == "rna"){
	print("Processing RNASeq datasets")
	run_module_pipeline_RNASeq(GSE_id = args[1], module_file_path = "", gse_info_file_path = "", logs_path = args[2])
}
if(args[4] == "ma"){
	print("Processing Microarray datasets")
        run_module_pipeline_MA(GSE_id = args[1], module_file_path = "", gse_info_file_path = "", logs_path = args[2])
}
