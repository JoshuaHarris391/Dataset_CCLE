


# Function for RPKM to TPM 
# TPM = ( RPKM / sum of RPKM over all genes/transcripts) * 10^6

RPKM_TPM <- function(x){
  (x/sum(x))*10^6
}


# Making a subset to test 
Exp.RPKM <- Exp[, 1:10]

# Applying Function
Exp.TPM <- apply(Exp.RPKM, 2, RPKM_TPM)

# Testing 
colSums(Exp.RPKM)
colSums(Exp.TPM)

########################################################