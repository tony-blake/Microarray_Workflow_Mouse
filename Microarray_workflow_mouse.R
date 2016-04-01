# Workflow_for_parkisondata
# Version 3, Wednesday 29 July 2015
# Created by: Tony Blake 

# Description: This  Work Flow analyses microarray data (Human, Mouse and Rat)
# for Gene Expression in Parkinson Disease Cases and considers genes of interest that 
# relate to the neurotrophic sugnalling pathways involved in the () treatment of Parkinsons

# Genes of interest: 
# GDNF: GDNF Ret tyrosine kinase (receptor) and GFR alpha 1, Shc, Src, FRS2, DOK 4/5, IRS1/2, GAB1, GRB2, AKT
# Neurturin: Neurturin, Ret tyrosine kinase (receptor) and GFR alpha 2
# GDF5: GDF5, BMPR1b, BMPR2, ROR2, SMAD1, SMAD4, SMAD5, SMAD 8, SMAD2, SMAD 3, Akt
# CDNF
# MKP1 or DuSP1
# p38 alpha, beta, gamma and delta kinase
# JNK1 JNK2, JNK3 kinase
# ERK1, ERK2 and ERK5

# Workflow -> Microarray data normilisation, generation of differentially expressed genes, 
# selection of expression data for g.o.i,  prodcution of heatmap, GSEA Pathway Analysis using and GO

# Install Packages
source('http://www.bioconductor.org/biocLite.R')
biocLite('Sweave')
biocLite('arrayQualityMetrics')
biocLite('GEOquery')
biocLite('limma')
biocLite('marray')
biocLite('simpleaffy')
bioclite('affy')
biocLite('gplots')
#biocLite('biomaRt')
biocLite('hgu133plus2.db')
biocLite('hgu133a2.db')
biocLite('rgu34a.db')
biocLite('mgu74a.db')
biocLite('gage')
biocLite('gageData')
biocLite('pathview')
biocLite('stats')
#biocLite('multtest')
library(bioDist)

# Load librarys
library(arrayQualityMetrics)
library(GEOquery)
library(limma)
library(marray)
library(simpleaffy)
library(affy)
library(gplots)
#library(biomaRt)
library(hgu133plus2.db)
library(hgu133a2.db)
library(rgu34a.db)
library(mgu74a.db)
library(gage)
library(gageData)
library(pathview)
library(stats)
library(plyr)
library(GOstats)
library(multtest)
library(xlsx)
library(reshape2)
############################################################################################################################################################################
################################################################  Begin Paramter Block #####################################################################################
runName <- "GSE4788" # string that is applied to all output files.
runName2 <- "mGSE4788"
inputDirectory <- "/Users/tonyblake/Desktop/Bioinformatics/parkinson_project/FinalCorrectedGSE/GSE4788/data_GSE4788" # Path to directory that contains the .CEL files and phenodata text file.
outputDirectory <- "/Users/tonyblake/Desktop/Bioinformatics/parkinson_project/GSE4788/Output" # Path to directory that outputs heatmaps, tables, etc 

baseline <- c("control") # string to convey baseline condition from target column in phnodata text file
minFC <- 1.5 # This is the minimum log2 fold change for a gene to be differentially expressed.
ttestPVal <- 0.01 # This is the threshold p value for significance. 
hgCutoff <- 0.01 # This is the GOStats p value threshold.

################################################################ End Parameter Block #######################################################################################

# Download raw data file to a named data directory and uncompress files

getGEOSuppFiles("GSE4788")
untar("GSE4788/GSE4788_RAW.tar", exdir="data_GSE4788") # uncompress folder
cels <- list.files("data_GSE4788/", pattern = "[gz]") # list all files ending in ".gz" as character strings
sapply(paste("data_GSE4788", cels, sep="/"), gunzip) # unzip .CEL files in data directory

################################################# Begin  "phenodata.txt" file creation from commandline #############################################################################

# $ cd "/Users/tonyblake/Desktop/Bioinformatics/parkinson_project/WorkFlow/data_GSE17204"
# $ ls *.CEL > phenodata.txt
# ....open in spreadsheet and create extra columns of data, copy all text and use "paste and match style" command to paste over original text in "phenodata.txt" file
#  so that it looks like 
#
# Name  File_Name  Target
# GSM430339.CEL  GSM430339.CEL  antisense_DJI_B
# GSM430340.CEL  GSM430340.CEL  antisense_DJI_B
# GSM430341.CEL	GSM430341.CEL	antisense_DJI_G
# GSM430342.CEL	GSM430342.CEL	antisense_DJI_G
# GSM430343.CEL	GSM430343.CEL	s_control_m
# GSM430344.CEL	GSM430344.CEL	s_control_m
# GSM430345.CEL	GSM430345.CEL	s_control_h
# GSM430346.CEL	GSM430346.CEL	s_control_h

####################################################### End "phenodata.txt" file creation ################################################################################## 

############################## Begin Normalisation Block ########################################################################################################

celfiles <- read.affy(covdesc="phenodata.txt", path=inputDirectory) # creates affybatch object
expression <- expresso(celfiles, normalize.method='quantiles', bgcorrect.method='rma', pmcorrect.method='pmonly',summary.method='medianpolish')

name="MAplot"
pdf(name)
MAplot(celfiles,pairs=TRUE,plot.method="smoothScatter", cex=.7)
dev.off()

name="MAplot2"
pdf(name)
MAplot(celfiles,pairs=FALSE,plot.method="smoothScatter")
dev.off()

exprs <- exprs(expression) # extracts logfold change (expression) data
conditions <- pData(expression)$Target
#arrayQualityMetrics(expressionset = expression, outdir = "GSE4788_QAnorm", force = FALSE, do.logtransform = TRUE, intgroup = colnames(pData(expression)))

name = "boxplotnorm.pdf"
pdf(name)
BiocGenerics::boxplot(exprs,col='red',names=pData(expression)$File_Name)
dev.off()


allunique_Conditions <- unique(conditions)
otherconditions <- allunique_Conditions[allunique_Conditions !=baseline]
reversetwo <- rev(otherconditions)
allotherconditions <- c(otherconditions, reversetwo)
nocontrols <- paste(otherconditions[1],"-", otherconditions[2])

############################## End Normalisation Block ###########################################################################################################

################ Selection of Differentially Expressed Genes ################


## Limma DEG calculations are done here:
gc()
f <- factor(conditions, levels=allunique_Conditions)
design <- model.matrix(~0+f)
colnames(design) <- allunique_Conditions
fit <- lmFit(expression, design)
contrast <- paste(allotherconditions,"-", c(baseline,baseline))
contrast.matrix <- makeContrasts(contrast[1],contrast[2], nocontrols,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

##############################################################################################################

#link gene symbols to affymatrix
goi_symbols <- read.table("parkinson_goi.txt")
#goi_symbols <- apply(goi_symbols, 2, function(x) toupper(x)) # use to change lowercase to uppercase  
goi_symbols <- data.frame(goi_symbols)
genesofinterest <- goi_symbols$V1
genesofinterest <- as.character(genesofinterest)
mygenelist <- mget(genesofinterest, revmap(mgu74aSYMBOL), ifnotfound=NA) #use for human datasets only! 
mygenedf <- do.call("rbind", lapply(mygenelist,data.frame)) # creates dataframe from gene and symbol list
which(is.na(mygenelist)) # Gives names of genes of interest  not included in affymatrix
u <- mygenedf$X..i..
mygenedf_nona <- mygenedf[!is.na(u),] # removes NA's from dataframe of genes
mygenedf_nona <- as.character(mygenedf_nona)

# link probeid's to genesymbols
goi_linkedprobes <- select(mgu74a.db, keys= mygenedf_nona, columns = "SYMBOL", keytype = "PROBEID") #for human dataset only!
affys <- goi_linkedprobes$PROBEID


# loop through limmma coef's to create 3 limma tables for the 3 comparisons
tops <- list()
symbol_linkedtogenes <- list()


for (i in 1:3){
  tops[[i]] <- topTable(fit2, number= Inf,coef=i)[affys,] #adjust="none", p.value=ttestPVal, lfc=log2(minFC))
  symbol_linkedtogenes[[i]] <- merge(tops[[i]], goi_linkedprobes, by.x=0, by.y="PROBEID")
  col_1st_goi <- grep("SYMBOL", names(symbol_linkedtogenes[[i]]))
  symbol_linkedtogenes[[i]] <- symbol_linkedtogenes[[i]][,c(col_1st_goi, (1:ncol(symbol_linkedtogenes[[i]]))[-col_1st_goi])]
  colnames(symbol_linkedtogenes[[i]])[2] <- "ProbesetIds"
  entrezs <- base::unlist(AnnotationDbi::as.list(mgu74aENTREZID[affys]))
  ddd <- ldply(entrezs)
  ddd <- ddd[order(ddd$.id),]
  symbol_linkedtogenes[[i]]$ENTREZID<-ddd$V1 # Adds ENTREZID column
}

# check that forloop is giving right result

check1 <- topTable(fit2, number= Inf,coef=3)[affys,] #adjust="none", p.value=ttestPVal, lfc=log2(minFC))
linkedtogenes1 <- merge(check1, goi_linkedprobes, by.x=0, by.y="PROBEID")
col_1st_goi1 <- grep("SYMBOL", names(linkedtogenes1))
linkedtogenes1 <- linkedtogenes1[,c(col_1st_goi1, (1:ncol(linkedtogenes1))[-col_1st_goi1])]
colnames(linkedtogenes1)[2] <- "ProbesetIds"
entrezs <- unlist(AnnotationDbi::as.list(mgu74aENTREZID[affys]))
ddd <- ldply(entrezs)
ddd <- ddd[order(ddd$.id),]
linkedtogenes1$ENTREZID<-ddd$V1 # Adds ENTREZID column

all.equal(symbol_linkedtogenes[[3]], linkedtogenes1)


# Final versions of Limma tables
tops1 <- symbol_linkedtogenes[[1]]
tops2 <- symbol_linkedtogenes[[2]]
tops3 <- symbol_linkedtogenes[[3]]
#tops4 <- symbol_linkedtogenes[[4]]
#tops5 <- symbol_linkedtogenes[[5]]

write.xlsx(tops1, file="limma_table_GSE4788_MMEvcontrol.xlsx")
write.xlsx(tops1, file="limma_table_GSE4788_MMFvcontrol.xlsx")
write.xlsx(tops1, file="limma_table_GSE4788_MMEvMMF.xlsx")

# create expression tables for heatmap

exprs=exprs(expression)[affys,]
affyid=rownames(exprs)
symbol2probe=mgu74aSYMBOL[affyid]
annots=toTable(symbol2probe)
str(annots)
#check to see if dataframes are same
df1 <- sort(annots$probe_id, decreasing=TRUE) 
df2 <- sort(goi_linkedprobes$PROBEID, decreasing=TRUE)
all.equal(df1,df2)
#back to main workflow
exprs=exprs[annots$probe_id,]

#if multiple probe sets map to a gene, select the one with maximal IQR
iqrs=apply(exprs, 1, IQR)
sel.rn=tapply(1:nrow(annots), annots$symbol, function(x){
  x[which.max(iqrs[x])]
})
exprs.symb=exprs[sel.rn,] #expression data with single gene names
rownames(exprs.symb)=names(sel.rn)

mexprs.symb <- exprs
rownames(mexprs.symb)=annots$symbol
mexprs.symb <- mexprs.symb[order(rownames(mexprs.symb)),] #expression data where there are more than one of the smae gene

################################################################# End Selection of Differentially Expressed Genes ############

sum(topTable(fit2, number=Inf, coef=1)$adj.P.Val < 0.05)
sigtable1 <- topTable(fit2, number=Inf, coef=1, p.value = 0.05)
sigaffyid <- rownames(sigtable1)

exprs_sig <- exprs[sigaffyid,]






################################################################# Begin Heatmap Block ##########################################



# heatmap for single genes
pdfOutputFile <- paste(outputDirectory, "/", runName, "_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Rowv=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=1, cexCol=.55)
dev.off()

# heatmap for multiple genes
pdfOutputFile <- paste(outputDirectory, "/", runName2, "_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(mexprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Rowv=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.7, cexCol=.55)
dev.off()

# heatmap for single genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName , "genecluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=1, cexCol=.55)
dev.off()

# heatmap for multiple genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2 , "_mgenecluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(mexprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.7, cexCol=.55)
dev.off()

# heatmap for total genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2 , "_tgenecluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.01, cexCol=.01)
dev.off()

# heatmap for significant genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2 , "_siggenecluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs_sig, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.01, cexCol=.01)
dev.off()

# heatmap for total genes with no gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2, "_tnoclust_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Rowv=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.01, cexCol=.01)
dev.off()

# heatmap for significant genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2 , "_siggene_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="control") "#FF0000" else {if (Target=="MME") "#0000FF" else  "#2dbdfc"} }
patientcolors <- unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs_sig, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Rowv=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.01, cexCol=.01)
dev.off()
################################### End Heatmap Block ##################################################################

################################### Begin Single gene Limma tables ###########################################################

sannots <- annots[sel.rn,]

# Final versions of Limma tables
stops1 <- tops1[rownames(sannots),] 
stops2 <- tops2[rownames(sannots),]
stops3 <- tops3[rownames(sannots),]

write.xlsx(stops1, file="slimma_table_GSE4788_MMEvcontrol.xlsx")
write.xlsx(stops2, file="slimma_table_GSE4788_MMFvcontrol.xlsx")
write.xlsx(stops3, file="slimma_table_GSE4788_MMEvMMF.xlsx")

################################## End Single gene Limma tables #########################################################

x <- org.Mm.egACCNUM
mapped_genes <- mappedkeys(x)
xx <- AnnotationDbi::as.list(x[mapped_genes])
geneUniverse <- (unique(names(xx)))


# Gene Ontology Categories for condition  that were shown to be relatively Higher (more expressed) 
# for (1) "antisense_DJI_B than s_control_m" 
# for (2) "antisense_DJI_G than s_control_h" 
# for (3) "antisense_DJI_B than s_control_m" 
# for (4) "antisense_DJI_G than s_control_h"
# for (5) "antisense_DJI_B than antisense_DJI_G"

sig_goi <- list()
upreg_goi <- list()
upgenesLinkedtoEntrezIds <- list()
GOstats_upgenes <- list()
paramsBPH <- list()
paramsCCH <- list()
paramsMFH <- list()
BP_GOdataH <- list()
BP_SummH <- list()
CC_GOdataH <- list()
CC_SummH <- list()
MF_GOdataH <- list()
MF_SummH <- list()

for(i in 1:3){
  
  sig_goi[[i]]<-subset(symbol_linkedtogenes[[i]],symbol_linkedtogenes[[i]][,7]<0.05)
  upreg_goi[[i]] <- subset(sig_goi[[i]], sig_goi[[i]][,3]>0)
  up_gostat <- upreg_goi[[i]][1]
  upgenesCHR <- up_gostat$SYMBOL
  upgenesLinkedtoEntrezIds[[i]] <- select(mgu74a.db, keys= upgenesCHR, "ENTREZID", "SYMBOL")
  GOstats_upgenes[[i]] <- upgenesLinkedtoEntrezIds[[i]][,2]
  
  paramsBPH[[i]] <- new("GOHyperGParams", geneIds=GOstats_upgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsCCH[[i]] <- new("GOHyperGParams", geneIds=GOstats_upgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Mm.eg.db", ontology="CC", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsMFH[[i]] <- new("GOHyperGParams", geneIds=GOstats_upgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Mm.eg.db", ontology="MF", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  BP_GOdataH[[i]] <- hyperGTest(paramsBPH[[i]]) # Biological Process.
  BP_SummH[[i]] <- summary(BP_GOdataH[[i]])
  
  CC_GOdataH[[i]] <- hyperGTest(paramsCCH[[i]]) # Cellular Component.
  CC_SummH[[i]] <- summary(CC_GOdataH[[i]])
  
  MF_GOdataH[[i]] <- hyperGTest(paramsMFH[[i]]) # Molecular Function.
  MF_SummH[[i]] <- summary(MF_GOdataH[[i]])
  
}

# Gene Ontology Categories for condition  that were shown to be relatively Lower (less expressed)

downreg_goi <- list()
downgenesLinkedtoEntrezIds <- list()
GOstats_downgenes <- list()
paramsBPL <- list()
paramsCCL <- list()
paramsMFL <- list()
BP_GOdataL <- list()
BP_SummL <- list()
CC_GOdataL <- list()
CC_SummL <- list()
MF_GOdataL <- list()
MF_SummL <- list()

for(i in 1:3){
  
  
  downreg_goi[[i]] <- subset(sig_goi[[i]], sig_goi[[i]][,3]<0)
  down_gostat <- downreg_goi[[i]][1]
  downgenesCHR <- down_gostat$SYMBOL
  downgenesLinkedtoEntrezIds[[i]] <- select(mgu74a.db, keys= downgenesCHR, "ENTREZID", "SYMBOL")
  GOstats_downgenes[[i]] <- downgenesLinkedtoEntrezIds[[i]][,2]
  
  paramsBPL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsCCL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Mm.eg.db", ontology="CC", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsMFL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Mm.eg.db", ontology="MF", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  BP_GOdataL[[i]] <- hyperGTest(paramsBPL[[i]]) # Biological Process.
  BP_SummL[[i]] <- summary(BP_GOdataL[[i]])
  
  CC_GOdataL[[i]] <- hyperGTest(paramsCCL[[i]]) # Cellular Component.
  CC_SummL[[i]] <- summary(CC_GOdataL[[i]])
  
  MF_GOdataL[[i]] <- hyperGTest(paramsMFL[[i]]) # Molecular Function.
  MF_SummL[[i]] <- summary(MF_GOdataL[[i]])
  
}


## Write out the results
source_https('https://raw.github.com/bobthecat/codebox/master/write.GOhyper.r')

bph1GO <- write.GOhyper(BP_GOdataH[[1]], filename="hBP_GSE4788_MMEvcontrol.xls")
bph2GO <- write.GOhyper(BP_GOdataH[[2]], filename="hBP_GSE4788_MMFvcontrol.xls")
#bph3GO <- write.GOhyper(BP_GOdataH[[3]], filename="hBP_GSE4788_MMEvMMF.xls")

cch1GO <- write.GOhyper(CC_GOdataH[[1]], filename="hCC_GSE4788_MMEvcontrol.xls")
cch2GO <- write.GOhyper(CC_GOdataH[[2]], filename="hCC_GSE4788_MMFvcontrol.xls")
#cch3GO <- write.GOhyper(CC_GOdataH[[3]], filename="hCC_GSE4788_MMEvMMF.xls")

mfh1GO <- write.GOhyper(MF_GOdataH[[1]], filename="hMF_GSE4788_MMEvcontrol.xls")
mfh2GO <- write.GOhyper(MF_GOdataH[[2]], filename="hMF_GSE4788_MMFvcontrol.xls")
#mfh3GO <- write.GOhyper(MF_GOdataH[[3]], filename="hMF_GSE4788_MMEvMMF.xls")

bpl1GO <- write.GOhyper(BP_GOdataL[[1]], filename="lBP_GSE4788_MMEvcontrol.xls")
bpl2GO <- write.GOhyper(BP_GOdataL[[2]], filename="lBP_GSE4788_MMFvcontrol.xls")
#bpl3GO <- write.GOhyper(BP_GOdataL[[3]], filename="lBP_GSE4788_MMEvMMF.xls")


ccl1GO <- write.GOhyper(CC_GOdataL[[1]], filename="lCC_GSE4788_MMEvcontrol.xls")
ccl2GO <- write.GOhyper(CC_GOdataL[[2]], filename="lCC_GSE4788_MMFvcontrol.xls")
#ccl3GO <- write.GOhyper(CC_GOdataL[[3]], filename="lCC_GSE4788_MMEvMMF.xls")


mfl1GO <- write.GOhyper(MF_GOdataL[[1]], filename="lMF_GSE4788_MMEvcontrol.xls")
mfl2GO <- write.GOhyper(MF_GOdataL[[2]], filename="lMF_GSE4788_MMFvcontrol.xls")
#mfl3GO <- write.GOhyper(MF_GOdataL[[3]], filename="lMF_GSE4788_MMEvMMF.xls")

###################################### Colour Map #############################################

# color scale
c1 <- rainbow(32,v=seq(0.5,1,length=32),s=seq(1,0.3,length=32),start=4/6,end=4.0001/6);
pie(rep(1,32),col=c1);   # show off these colors
c2 <- rainbow(32,v=seq(0.5,1,length=32),s=seq(1,0.3,length=32),start=1/6,end=1.0001/6);
pie(rep(1,32),col=c2);   # show off these colors
c3 <- c(c1,rev(c2));  # rev reverses the list
#rm(c1,c2)

######################################## Correlation Analysis ###################################
 
# Here the degree of correlation between different gene ontologies is analysed based on the 
# the differentailly expressed genes that they have in common.

# The following pieces of code take the objects created from the GOstats analysis as input. It then creates a matrix of 
# correlation values for gene ontologies and uses it to create a correlogram showing the percenatge of genes shared 
# between gene ontologies in the form a clustered heatmap. It also creates a table of gene ontologies showing the over
# represented genes in each ontology. Rows are ordered according to the clustering heatmap of correlation values. 
 
#upregulated pathways

bph1square <- makebeta1_pathwaygenetable(bph1GO)
hr <- hclust(cor.dist(bph1square[[1]]))
hc <- hclust(cor.dist(bph1square[[1]]))
pdf(file='hBP_correlmap_GSE4788_controlvMME.pdf', height=13, width=13)
heatmap.2(bph1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.1, cexCol=0.1)
dev.off()
bph1clustered <- bph1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bph1clustered)
clustered <- bph1square[[2]][match(clst, bph1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="hBP_genesinpathways_GSE4788_controlvMME.xlsx")

cch1square <- makebeta1_pathwaygenetable(cch1GO)
hr <- hclust(cor.dist(cch1square[[1]]))
hc <- hclust(cor.dist(cch1square[[1]]))
pdf(file='hCC_correlmap_GSE4788_controlvMME.pdf', height=13, width=13)
heatmap.2(cch1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
cch1clustered <- cch1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(cch1clustered)
clustered <- cch1square[[2]][match(clst, cch1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="hCC_genesinpathways_GSE4788_controlvMME.xlsx")

mfh1square <- makebeta1_pathwaygenetable(mfh1GO)
hr <- hclust(cor.dist(mfh1square[[1]]))
hc <- hclust(cor.dist(mfh1square[[1]]))
pdf(file='hMF_correlmap_GSE4788_controlvMME.pdf', height=13, width=13)
heatmap.2(mfh1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
mfh1clustered <- mfh1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(mfh1clustered)
clustered <- mfh1square[[2]][match(clst, mfh1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="hMF_genesinpathways_GSE4788_controlvMME.xlsx")

#################################################################
######## down reg

bpl1square <- makebeta1down_pathwaygenetable(bpl1GO)
hr <- hclust(cor.dist(bpl1square[[1]]))
hc <- hclust(cor.dist(bpl1square[[1]]))
pdf(file='lBP_correlmap_GSE4788_controlvMME.pdf', height=13, width=13)
heatmap.2(bpl1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.1, cexCol=0.1)
dev.off()
bpl1clustered <- bpl1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bpl1clustered)
clustered <- bpl1square[[2]][match(clst, bpl1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lBP_genesinpathways_GSE4788_controlvMME.xlsx")


#This one is wrong
bpl2square <- makebeta2down_pathwaygenetable(bpl2GO)
hr <- hclust(cor.dist(bpl2square[[1]]))
hc <- hclust(cor.dist(bpl2square[[1]]))
pdf(file='lBP_correlmap_GSE4788_controlvMMF.pdf', height=13, width=13)
heatmap.2(bpl2square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
bpl2clustered <- bpl2square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bpl2clustered)
clustered <- bpl2square[[2]][match(clst, bpl2square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lBP_genesinpathways_GSE4788_controlvMMF.xlsx")

ccl1square <- makebeta1down_pathwaygenetable(ccl1GO)
hr <- hclust(cor.dist(ccl1square[[1]]))
hc <- hclust(cor.dist(ccl1square[[1]]))
pdf(file='lCC_correlmap_GSE4788_controlvMME.pdf', height=13, width=13)
heatmap.2(ccl1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
ccl1clustered <- ccl1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(ccl1clustered)
clustered <- ccl1square[[2]][match(clst, ccl1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lCC_genesinpathways_GSE4788_controlvMME.xlsx")

mfl1square <- makebeta1down_pathwaygenetable(mfl1GO)
hr <- hclust(cor.dist(mfl1square[[1]]))
hc <- hclust(cor.dist(mfl1square[[1]]))
pdf(file='lMF_correlmap_GSE4788_controlvMME.pdf', height=13, width=13)
heatmap.2(mfl1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
mfl1clustered <- mfl1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(mfl1clustered)
clustered <- mfl1square[[2]][match(clst, mfl1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lMF_genesinpathways_GSE4788_controlvMME.xlsx")


#this one is wrong
mfl2square <- makebeta2down_pathwaygenetable(mfl2GO)
hr <- hclust(cor.dist(mfl2square[[1]]))
hc <- hclust(cor.dist(mfl2square[[1]]))
pdf(file='lMF_correlmap_GSE4788_controlvMMF.pdf', height=13, width=13)
heatmap.2(mfl1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
mfl2clustered <- mfl2square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(mfl2clustered)
clustered <- mfl2square[[2]][match(clst, mfl2square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lMF_genesinpathways_GSE4788_controlvMMF.xlsx")
