# Reading in files 
Exp <- read.table(file = "Source_Files/CCLE_RPKM.gct", skip = 2, header = T, na.strings = c("", NA))
CopNum <- read.table(file = "Source_Files/CCLE_COPYNUM.txt", header = T, na.strings = c("", NA))
SamInfo <- read.delim(file = "Source_Files/CCLE_SAMPLE_INFO.txt", header = T, na.strings = c("", NA))
Reads <- read.table(file = "Source_Files/CCLE_READS.gct", skip = 2, header = T, na.strings = c("", NA))

# Converting Exp.Data from RPKM to TPM 
RPKM_TPM <- function(x){
  (x/sum(x))*10^6
}


# Removing Spaces, Hyphens, and Periods in Cell Line Names 
SamInfo$Cell.line.primary.name <- gsub("-","", SamInfo$Cell.line.primary.name)
SamInfo$Cell.line.primary.name <- gsub(" ","", SamInfo$Cell.line.primary.name)
SamInfo$Cell.line.primary.name <- gsub("\\.","", SamInfo$Cell.line.primary.name)

# Making Annotation Table 
Anot <- Exp[, 1:2]


# Converting RPKM to TPM
Exp.TPM <- Exp
Exp.TPM[, 3:1021] <- apply(Exp.TPM[, 3:1021], 2, RPKM_TPM)


# Testing RPKM and TPM colsums 
colSums(Exp[, 3:1021])
colSums(Exp.TPM[, 3:1021])

# Writing Files
# dir.create("DATA")
write.table(Exp.TPM, file = "DATA/CCLE_TPM.txt", sep = "\t")
write.table(Anot, file = "DATA/CCLE_ANOT.txt", sep = "\t")
# write.table(SamInfo, file = "DATA/CCLE_SAMINFO.txt", sep = "\t")
# write.table(CopNum, file = "DATA/CCLE_COPNUM.txt", sep = "\t")
write.table(Exp, file = "DATA/CCLE_RPKM.txt", sep = "\t")
write.table(Reads, file = "DATA/CCLE_READS.txt", sep = "\t")

## Reading in Files and saving Global Environment 
# Clearing Global Environment 
rm(list=ls())
# Loading in files 
CCLE_READS <- read.delim(file = "DATA/CCLE_READS.txt", header = T, na.strings = c("", NA))
CCLE_RPKM  <- read.delim(file = "DATA/CCLE_RPKM.txt", header = T, na.strings = c("", NA))
CCLE_TPM   <- read.delim(file = "DATA/CCLE_TPM.txt", header = T, na.strings = c("", NA))
ANOT           <- read.delim(file = "DATA/CCLE_ANOT.txt", header = T, na.strings = c("", NA))
CCLE_SAMINFO        <- read.delim(file = "DATA/CCLE_SAMINFO.txt", header = T, na.strings = c("", NA))
CCLE_COPNUM         <- read.delim(file = "DATA/CCLE_COPNUM.txt", header = T, na.strings = c("", NA))
# Saving Global Env
save.image(file = "DATA/CCLE_GLOB_ENV.RData")
# Loading Glob Env 
load(file = "DATA/CCLE_GLOB_ENV.RData")
