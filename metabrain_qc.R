# Script to QC all metabrain tissues.


# Load all libraries needed
.libPaths("/data/kronos/kronos/acarrasco/R_libs/")
library(data.table)
library(purrr)
library(tidyverse)
library(colochelpR)
library(data.table)
library(DT)


# Create util to get the total sample size for each SNP
sumSampleSize = function(string) {
  unlist(regmatches(string, gregexpr("[[:digit:]]+", string))) %>%
    as.numeric %>% 
    sum
}

processing_function = function(tissue_chrfiles, 
                               tissue_name,
                               raw_data_path = rawpath,
                               qc_data_path = qcpath) {
  
  # MAF file processing
  maf_index = grep(pattern="MAF", x = tissue_chrfiles)
  maf = fread(paste(raw_data_path, tissue_name, tissue_chrfiles[maf_index], sep = "/" ))
  tissue_chrfiles = tissue_chrfiles[-maf_index]
  
  # map each tissue chunk to the function that does the actual processing
  map(tissue_chrfiles, function(x,
                                maf. = maf,
                                tissuename = tissue_name, 
                                raw_path = raw_data_path,
                                qc_path = qc_data_path) {
    
    tissuechr_name = x
    cat("\nProcessing", x, "\n")
    tissuechr_chunk = fread(paste(raw_path, tissuename, tissuechr_name, sep = "/" ))
    cat("Dims before any QC", dim(tissuechr_chunk), "\n")
  
    # Inner joining the maf file to the tissuechr_chunk and trim the ENDEMBL ID
    tissuechr_chunk_proc = tissuechr_chunk %>%
      dplyr::inner_join(maf. %>% dplyr::select(SNP, MAF), by = c("SNPName" = "SNP" ) ) %>%
      dplyr::mutate(ProbeName = gsub("\\..*", "", ProbeName))
    
    # Get rsID on a column and CHR:BP on another column
    snpname = as.data.frame(
      stringr::str_split_fixed(
        tissuechr_chunk_proc$SNPName, n = 4, pattern = ":"
      ))
    tissuechr_chunk_proc$SNP = snpname$V3
    tissuechr_chunk_proc$CHRPOS = paste(snpname$V1, snpname$V2, sep = ":")
    
    # Find the column number with the sample size and process it
    sampleSize_index = grep("DatasetsNrSamples", colnames(tissuechr_chunk_proc))
    NrSamples = apply(tissuechr_chunk_proc %>% select(all_of(sampleSize_index)), 1, sumSampleSize)
    tissuechr_chunk_proc = cbind(tissuechr_chunk_proc, NrSamples)
    
    # Add column stating eQTL
    tissuechr_chunk_proc$eQTL = paste0("eQTL_", tissuename)
    cat("Dims after QC was applied", dim(tissuechr_chunk_proc), "\n")
    
    # Write results in the QC path
    fwrite(tissuechr_chunk_proc, 
           paste(qc_path, tissuename, paste0("QCed_",tissuechr_name), sep = "/" ),
           quote = F, sep = "\t", row.names=F, col.names=T)
  })
}


setwd("/mnt/rreal/RDS/DATA/eQTLdata/METABRAIN")
rawpath = "/mnt/rreal/RDS/DATA/eQTLdata/METABRAIN/RAW_DATA/"
qcpath = "/mnt/rreal/RDS/DATA/eQTLdata/METABRAIN/METABRAIN_QC_V2/"

# Get all tissues
tissues = list.files(rawpath)

# Get all data for ea
tissues_chromosomes =  map(tissues, function(tissue) {
  allfiles = list.files(paste(rawpath, tissue, sep = "/"))
  allfiles
}) %>%
  setNames(tissues)

map2(tissues_chromosomes, names(tissues_chromosomes), processing_function)
