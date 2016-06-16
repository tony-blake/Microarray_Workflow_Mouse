# Microarray_Workflow_Mouse

Workflow to analyse data from the mgu74a platform

In this repository there are 5 files that comprise a workflow for analysing the GSE17204 dataset.

"Microarray_workflow_mouse.R" - downloads the data from GEO and performs LIMMA, GOstats, and correlation of over rperesented gene ontologies analyses.

"make_pathwaygenetable.R" - function callled in correlation analysis section of "Microarray_workflow_mouse.R". It creates a matrix of correlation values for gene ontologies and uses it to create a correlogram showing the percenatge of genes shared between gene ontologies in the form a clustered heatmap. It also creates a table of gene ontologies showing the over represented genes in each ontology. Rows are ordered according to the clustering heatmap of correlation values.

"make_betamatrix2d.R" - function called in the "make_pathwaygenetable.R" file. The function creates the matrix of values that represent the percentage of genes shared between gene ontologies.

"source_https.R" - function called in "workflowGSE17204.R" that can download and run other functions from a different users Github account

"phenodata.txt" - meta data file. This needs to placed in the same directory as the .CEL files once they have been downlaoded.

It should be noted also that the main file "MIcroarray_workflow_mouse.R" needs to be run in sections, each section being a code block. The reason for doing so is primarily to catch bugs but also because the metadata text file for the LIMMA needs to be created (directions for doing this are in the comments in "Microarray_workflow_mouse.R") After the .CEL files have been downlaoded they are placed in the data directory created from running the codeblock that downloaded the GEO data. It is in this directory that the metadata textfile needs to be placed.By doing this the next block of code will create an affybatch object for LIMMA.
