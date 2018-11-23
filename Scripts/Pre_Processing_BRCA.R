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

# Subsetting Breast Cancer Expression DF in RPKM
BRCA.REF <-as.character(SamInfo[SamInfo$Site.Primary == "breast", 1])
Sub <- is.na(match(BRCA.REF, colnames(Exp)))
BRCA.REF <- BRCA.REF[!Sub]
BRCA.EXP.CCLE.RPKM <- Exp[ ,c(1,2, match(BRCA.REF, colnames(Exp)))]
# Writing BRCA.EXP.CCLE.RPKM file
write.table(BRCA.EXP.CCLE.RPKM, file = "Source_Files/BRCA.EXP.CCLE.RPKM.txt", sep = "\t")

# Subsetting Breast Cancer Reads DF
BRCA.REF <-as.character(SamInfo[SamInfo$Site.Primary == "breast", 1])
Sub.read <- is.na(match(BRCA.REF, colnames(Reads)))
BRCA.REF.read <- BRCA.REF[!Sub.read]
BRCA.EXP.CCLE.READS <- Reads[ ,c(1,2, match(BRCA.REF.read, colnames(Reads)))]

# Converting RPKM to TPM
BRCA.EXP.CCLE.TPM <- BRCA.EXP.CCLE.RPKM
BRCA.EXP.CCLE.TPM[, 3:59] <- apply(BRCA.EXP.CCLE.TPM[, 3:59], 2, RPKM_TPM)

# Loading Annotation 
Anot <- read.delim(file = "Annotation/Cell_Line_Annotation.txt", header = T, na.strings = c("", NA))
Anot <- Anot[1:43, ]
Anot$CELL_LINE <- gsub("-", "", Anot$CELL_LINE)

# Testing RPKM and TPM colsums 
colSums(BRCA.EXP.CCLE.RPKM[, 3:59])
colSums(BRCA.EXP.CCLE.TPM[, 3:59])

# Writing Files
dir.create("DATA")
write.table(BRCA.EXP.CCLE.TPM, file = "DATA/BRCA_EXP_CCLE_TPM.txt", sep = "\t")
write.table(Anot, file = "DATA/BRCA_ANOT.txt", sep = "\t")
write.table(SamInfo, file = "DATA/CCLE_SAMINFO.txt", sep = "\t")
write.table(CopNum, file = "DATA/CCLE_COPNUM.txt", sep = "\t")
write.table(BRCA.EXP.CCLE.RPKM, file = "DATA/BRCA_EXP_CCLE_RPKM.txt", sep = "\t")
write.table(BRCA.EXP.CCLE.READS, file = "DATA/BRCA_EXP_CCLE_READS.txt", sep = "\t")

## Reading in Files and saving Global Environment 
# Clearing Global Environment 
rm(list=ls())
# Loading in files 
BRCA_EXP_CCLE_READS <- read.delim(file = "DATA/BRCA_EXP_CCLE_READS.txt", header = T, na.strings = c("", NA))
BRCA_EXP_CCLE_RPKM  <- read.delim(file = "DATA/BRCA_EXP_CCLE_RPKM.txt", header = T, na.strings = c("", NA))
BRCA_EXP_CCLE_TPM   <- read.delim(file = "DATA/BRCA_EXP_CCLE_TPM.txt", header = T, na.strings = c("", NA))
BRCA_ANOT           <- read.delim(file = "DATA/BRCA_ANOT.txt", header = T, na.strings = c("", NA))
CCLE_SAMINFO        <- read.delim(file = "DATA/CCLE_SAMINFO.txt", header = T, na.strings = c("", NA))
CCLE_COPNUM         <- read.delim(file = "DATA/CCLE_COPNUM.txt", header = T, na.strings = c("", NA))
# Saving Global Env
save.image(file = "DATA/BRCA_GLOB_ENV.RData")
# Loading Glob Env 
load(file = "DATA/BRCA_GLOB_ENV.RData")
