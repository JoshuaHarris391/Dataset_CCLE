
#' ## Setup
#----------------------------------------------
# Setting WD
setwd("/Users/joshua_harris/Dropbox/Research/PhD/Bioinformatics/Datasets/Dataset_CCLE")

# Packages
library(tidyverse)


#----------------------------------------------


#' # Reading in source data
#----------------------------------------------

# Setting WD 

# # Reading in source data files
# CCLE_SAMPLE_INFO_DF <- read.delim(file = "Source_Data/Legacy/CCLE_SAMPLE_INFO.txt", header = T, na.strings = c("", NA))
# class(CCLE_SAMPLE_INFO_DF)
# View(CCLE_SAMPLE_INFO_DF[1:5,1:5])
# 
# CCLE_RPKM_DF <- read.delim(file = "Source_Data/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", header = T, na.strings = c("", NA), skip = 2)
# class(CCLE_RPKM_DF)
# View(CCLE_RPKM_DF[1:5,1:5])
# 
# CCLE_READS_DF <- read.delim(file = "Source_Data/CCLE_DepMap_18q3_RNAseq_reads_20180718.gct", header = T, na.strings = c("", NA), skip = 2)
# class(CCLE_READS_DF)
# View(CCLE_READS_DF[1:5,1:10])
# 
# CCLE_MAF_DF <- read.delim(file = "Source_Data/CCLE_DepMap_18q3_maf_20180718.txt", header = T, na.strings = c("", NA))
# class(CCLE_MAF_DF)
# View(CCLE_MAF_DF[1:5,])
# 
# 
# CCLE_COPY_NUM_DF <- read.table(file = "Source_Data/Legacy/CCLE_COPYNUM.txt", header = T, na.strings = c("", NA))
# class(CCLE_COPY_NUM_DF)
# View(CCLE_COPY_NUM_DF[1:5,1:10])
# 
# # Saving an R data file for faster loading 
# save.image(file = "Source_Data/Source_Data.RData")

# Clearing Environment 
rm(list=ls())

# Reloading R Data file
load(file = "Source_Data/Source_Data.RData")

#----------------------------------------------


#' # Obtaining Sample names from multiple DFs
#----------------------------------------------

Copy_sam_nam <- names(CCLE_COPY_NUM_DF[, 6:ncol(CCLE_COPY_NUM_DF)])
MAF_sam_nam <- CCLE_MAF_DF$Tumor_Sample_Barcode %>% as.character() %>% unique()
Reads_sam_nam <- names(CCLE_READS_DF[, 3:ncol(CCLE_READS_DF)])
RPKM_sam_nam <- names(CCLE_RPKM_DF[, 3:ncol(CCLE_RPKM_DF)])
Sample_sam_nam <- CCLE_SAMPLE_INFO_DF$CCLE.name %>% as.character()

#----------------------------------------------


#' # Creating Gene annotation DF
#----------------------------------------------

CCLE_GENE_ANNOT_DF <- CCLE_RPKM_DF[, 1:2]

#----------------------------------------------



#' # Cleaning Reads and RPKM DF sample names 
#----------------------------------------------
# Finding names that start with 'X'
Reads_sam_nam[grep("^X", Reads_sam_nam )][1:10] %>% print()
RPKM_sam_nam[grep("^X", RPKM_sam_nam )][1:10] %>% print()

# Removing X from the start of names
Reads_sam_nam <- gsub("^X", "",Reads_sam_nam )
Reads_sam_nam[1:10] %>% print()
RPKM_sam_nam <- gsub("^X", "", RPKM_sam_nam )
RPKM_sam_nam[1:10] %>% print()

# Finding names with ..ACH
Reads_sam_nam[grep("\\.\\.ACH.*", Reads_sam_nam )][1:10] %>% print()
RPKM_sam_nam[grep("\\.\\.ACH.*", RPKM_sam_nam )][1:10] %>% print()

# Removing ..ACH from names
Reads_sam_nam <- gsub("\\.\\.ACH.*", "",Reads_sam_nam )
Reads_sam_nam[1:10] %>% print()
RPKM_sam_nam <- gsub("\\.\\.ACH.*", "", RPKM_sam_nam )
RPKM_sam_nam[1:10] %>% print()

#----------------------------------------------


#' # Renaming Reads and RPKM sample names 
#----------------------------------------------
colnames(CCLE_READS_DF) <- c("Name", "Description", Reads_sam_nam) 
print(names(CCLE_READS_DF[, 1:10]))
colnames(CCLE_RPKM_DF) <- c("Name", "Description", RPKM_sam_nam)
print(names(CCLE_RPKM_DF[, 1:10]))
#----------------------------------------------

#' # Making a TPM table
#----------------------------------------------
# defining function
RPKM_TPM <- function(x){
  (x/sum(x))*10^6
}
# Converting values to TPM
CCLE_TPM <- apply(CCLE_RPKM_DF[, 3:ncol(CCLE_RPKM_DF)], 2, RPKM_TPM)
# Converting back to DF
CCLE_TPM_DF <- as.data.frame(CCLE_TPM)
# Re-attaching name and description info
CCLE_TPM_DF <- data.frame(CCLE_RPKM_DF[, 1:2], CCLE_TPM_DF)


#' # Saving Data files
#----------------------------------------------
# Writing data frames
write.table(CCLE_COPY_NUM_DF, file = "Output_Data/CCLE_COPY_NUM_DF.txt", sep = "\t")
write.table(CCLE_MAF_DF, file = "Output_Data/CCLE_MAF_DF.txt", sep = "\t")
write.table(CCLE_READS_DF, file = "Output_Data/CCLE_READS_DF.txt", sep = "\t")
write.table(CCLE_RPKM_DF, file = "Output_Data/CCLE_RPKM_DF.txt", sep = "\t")
write.table(CCLE_SAMPLE_INFO_DF, file = "Output_Data/CCLE_SAMPLE_INFO_DF.txt", sep = "\t")
write.table(CCLE_TPM_DF, file = "Output_Data/CCLE_TPM_DF", sep = "\t")
write.table(CCLE_GENE_ANNOT_DF, file = "Output_Data/CCLE_GENE_ANNOT_DF", sep = "\t")

# Clearing Environment 
rm(list=ls())

# Loading Data frames
CCLE_COPY_NUM_DF <- read.table(file = "Output_Data/CCLE_COPY_NUM_DF.txt", sep = "\t", header = TRUE)
CCLE_MAF_DF <- read.table(file = "Output_Data/CCLE_MAF_DF.txt", sep = "\t", header = TRUE)
CCLE_READS_DF <- read.table(file = "Output_Data/CCLE_READS_DF.txt", sep = "\t", header = TRUE)
CCLE_RPKM_DF <- read.table(file = "Output_Data/CCLE_RPKM_DF.txt", sep = "\t", header = TRUE)
CCLE_SAMPLE_INFO_DF <- read.table(file = "Output_Data/CCLE_SAMPLE_INFO_DF.txt", sep = "\t", header = TRUE)
CCLE_TPM_DF <- read.table(file = "Output_Data/CCLE_TPM_DF", sep = "\t", header = TRUE)
CCLE_GENE_ANNOT_DF <- read.table(file = "Output_Data/CCLE_GENE_ANNOT_DF", sep = "\t", header = TRUE)

# Saving an R data file for faster loading 
save.image(file = "Output_Data/CCLE_Processed_20190517.RData")

#----------------------------------------------
