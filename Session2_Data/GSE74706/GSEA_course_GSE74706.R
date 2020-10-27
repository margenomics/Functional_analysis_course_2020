##################################################################
##############Functional Analysis with GSEA#######################
#R commands in order to generate the functional analysis course###
##################################################################
#Load libraries
library(Biobase)
library(GEOquery)
library(limma)
library(AnalysisFunctions)
#Set paths
wd <- getwd()
resultsDir <- file.path(wd, "CursGSEA_BicoH/GSE74706")
##################################################################
############Load data from GEO####################################
gset <- getGEO("GSE74706", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13497", attr(gset, "names")) else idx <- 1
#ExpresionSet
gset <- gset[[idx]]
class(gset)

#Targets
Targets <- pData(phenoData(gset))
dim(Targets)#36 33
table(Targets$`tissue:ch1`)
# NSCLC Tumor-free lung 
# 18              18 
Targets$`tissue:ch1` <- gsub(" ", "_", Targets$`tissue:ch1`)
Targets$`tissue:ch1` <- gsub("-", "", Targets$`tissue:ch1`)
table(Targets$`tissue:ch1`)
# NSCLC Tumorfree_lung 
# 18             18 

#Expression matrix
data <- exprs(gset)
dim(data)#34183    36

#Change the column names
all.equal(colnames(data), rownames(Targets))#TRUE
colnames(data) <- Targets$title
colnames(data) <- gsub(" ", "_", colnames(data))

####################################################################
###############Output files for the GSEA analysis###################
#Make cls and gct files
generate_cls(Targets$`tissue:ch1`, cls=file.path(resultsDir, "NSCLCndNT.cls"))
generate_gct(data, gct=file.path(resultsDir, "NSCLCndNT.gct"))

#Write the expression matrix
write.table(data, file=file.path(resultsDir, "GSE74706_exprsMatrix.txt"), 
            sep="\t", quote=F, row.names = F)

###################################################################
###############DE analysis with R##################################


