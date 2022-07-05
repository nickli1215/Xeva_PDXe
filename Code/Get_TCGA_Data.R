library(UCSCXenaTools)
studies <- availTCGA("ProjectID")
studies <- studies[studies != "PANCAN"]
studies <- studies[studies != "FPPP"]
temp_dir <- "/tmp/TCGA/"
output_dir <- "../Data/TCGA"

system("mkdir /tmp/TCGAdownload")


for (study in studies) {
  print(study)
  temp_dir <- "/tmp/TCGA/"
  system(paste('rm -r', temp_dir))
  system(paste('mkdir', temp_dir))
  
  downloadTCGA(
    project = study
    , data_type = 'Gene Expression RNASeq'
    , file_type = 'IlluminaHiSeq RNASeqV2'
    , destdir = temp_dir
    , force = TRUE
  )
  temp_dir <- sprintf('/tmp/TCGA/TCGA.%s.sampleMap/', study)
  # temp_dir <- paste0(temp_dir, "TCGA.", study, ".", "sampleMAP/")
  system(paste0('gzip -d ', temp_dir, "HiSeqV2.gz"))
  df <- read.table(paste0(temp_dir, "HiSeqV2"), sep = "\t", header = TRUE)
  system(paste('rm -r', temp_dir))
  
  system(paste0("mkdir ", output_dir,"/", study))
  print(paste0("mkdir ", output_dir,"/",study))
  saveRDS(df, sprintf("%s/%s/%s.rds", output_dir, study, "RNASeq"))
  
}