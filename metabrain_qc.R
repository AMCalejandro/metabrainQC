# Script to QC all metabrain tissues.


# Load all libraries needed
.libPaths("/data/kronos/kronos/acarrasco/R_libs/")
library(data.table)
library(purrr)
library(tidyverse)
library(colochelpR)
library(data.table)
library(DT)



processing_function = function(tissue_chrfiles, 
                               tissue_name,
                               .raw_data_path = raw_data_path,
                               .qc_data_path = qc_data_path) {
  
  # MAF file processing
  maf_index = grep(pattern="MAF", x = tissue_chrfiles)
  maf = fread(paste(.raw_data_path, tissue_name, tissue_chrfiles[maf_index], sep = "/" ))
  tissue_chrfiles = tissue_chrfiles[-maf_index]
  
  # map each tissue chunk to the function that does the actual processing
  map(tissue_chrfiles, function(x,
                                maf. = maf,
                                tissue_name. = tissue_name, 
                                raw_path = .raw_data_path,
                                qc_path = .qc_data_path) {
    
    tissuechr_name = x
    cat("\nProcessing", x, "\n")
    tissuechr_chunk = fread(paste(raw_path, tissue_name., tissuechr_name, sep = "/" ))
    cat("Dims before any QC", dim(tissuechr_chunk), "\n")
  
    # Inner joining the maf file to the tissuechr_chunk
    tissuechr_chunk_proc = tissuechr_chunk %>%
      inner_join(maf. %>% dplyr::select(SNP, MAF), by = c("SNPName" = "SNP" ) )
    
    
    # Another keyQC is to get the number of samples per SNP
    sumSampleSize = function(string) {
      unlist(regmatches(string, gregexpr("[[:digit:]]+", string))) %>% as.numeric %>% sum
    }
    
    # Find the column number with the sample size and process it
    sampleSize_index = grep("DatasetsNrSamples", colnames(tissuechr_chunk_proc))
    NrSamples = apply(tissuechr_chunk_proc %>% select(all_of(sampleSize_index)), 1, sumSampleSize)
    tissuechr_chunk_proc = cbind(tissuechr_chunk_proc, NrSamples)
    
    # Add column stating eQTL
    tissuechr_chunk_proc$eQTL = paste0("eQTL_", tissue_name.)
    cat("Dims after QC was applied", dim(tissuechr_chunk_proc), "\n")
    
    # Write results in the QC path
    fwrite(tissuechr_chunk_proc, 
           paste(qc_path, tissue_name., paste0("QCed_",tissuechr_name), sep = "/" ),
           quote = F, sep = "\t", row.names=F, col.names=T)
  }
    )

}


setwd("/mnt/rreal/RDS/DATA/eQTLdata/METABRAIN")
raw_data_path = "/mnt/rreal/RDS/DATA/eQTLdata/METABRAIN/RAW_DATA/"
qc_data_path = "/mnt/rreal/RDS/DATA/eQTLdata/METABRAIN/METABRAIN_QC/"

# Get all tissues
tissues = list.files(raw_data_path)

# Get all data for ea
tissues_chromosomes =  map(tissues, function(tissue) {
  allfiles = list.files(paste(raw_data_path, tissue, sep = "/"))
  allfiles
}) %>%
  setNames(tissues)

map2(tissues_chromosomes, names(tissues_chromosomes), processing_function)





