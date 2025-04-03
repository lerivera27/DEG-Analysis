# DEG-Analysis-Lab
# Task: Finding Differentially Expressed Genes (DEGs) specifically analyzing BLCA to identify genes with significant expression changes between normal and tumorous tissue. 

## Objective: Use the DESeq2 R package within the R-Studio environment to identify DEGs between normal GTEX and tumor TCGA RNA-Seq profiles in a gene expression matrix (GEM). 


## Upon running the DESeq2 the following steps must be conducted: 
1. Start by setting up the DESeq2 package in your Linux environment. 
2. Download RNA-seq data for both normal and tumor for your desired cancer subset (we will use BLCA-bladder tissue samples). 
    a. Then, process and combine them into a single Gene Expression Matrix (GEM) for analysis. 
3. Set up a comparison matrix that will clearly define normal vs. tumor. 


# DESeq2 Installation in R: 
``` 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

## Set the working directory: 
``` 
setwd("/scratch/leiarar/deg-analysis-lab") #information to be used as an example! 

```

# After installation in R Studio: 
4. Use DESeq2 to analyze the data and identify genes that show significant differences in expression between the two conditions. 
5. After processing both conditions, merge and preprocess the matrices to explore similarities and differences across cancers. 

# Extract Data and Run DESeq2 
## Extract a sub-matrix of counts for each group: 
``` 
subgem <- function(gem, anot, group) {
  datalist = list()
  subanot = subset(anot, Comparison == group)
  for (id in subanot$Sample) {
    ind = which(colnames(gem) == id)
    genes = gem[0]
    exp = gem[,ind]
    datalist[[id]] <- exp
  }
  subcounts = cbind(genes, datalist)
  return(subcounts)
}
```

## Extract a subset of the sample annotation matrix for each group: 
``` 
subanot <- function(anot, group) {
  datalist = list()
  subanot = subset(anot, Comparison == group)
  return(subanot)
}

```
# Running the Script: 
To execute the DEG analysis, follow these steps: 
1. Ensure all required files are prepared: 
* 'bladder-gtex-tcga-clean.csv' which is the cleaned processed GEM file. 
* 'bladder-gtex-tcga.comparison.csv' which is the comparison of the sample groups normal versus tumorous. 
2. Open R-Studio and navigate to the working directory containing these files. 
3. Load DESeq2. 
4. Run the script. 
5. The output file 'bladder-gtex-tcga-degs.csv' will be generated containing all the differentially expressed resutls of the two conditions: normal vs. tumor for BLCA. 


## Run DESeq2: 

``` 
run_deseq <- function(counts, annotation) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = annotation,
                                design = ~ Group)
  dds <- dds[rowSums(counts(dds)) >= 50,]
  dds <- DESeq(dds)
  norm = fpm(dds)
  
  conditionA = which(annotation[2] == "BLCA_GTEX")
  conditionB = which(annotation[2] == "BLCA_TCGA")
  norm = subset(norm, select=c(conditionA, conditionB))
  
  res <- results(dds, contrast=c("Group", "BLCA_GTEX", "BLCA_TCGA"))
  res <- cbind(res, norm)
  return(res)
}
``` 
## Print Results to a separate output file: 
``` 
main <- function(countfile, anotfile, outfile) {
  counts = read.delim(countfile, sep=',', header=TRUE, row.names='Hugo_Symbol')
  samples = read.delim(anotfile, sep=',', row.names = NULL, check.names=FALSE)
  groups = unique(samples$Comparison)
  
  for (t in groups) {
    subcounts = subgem(counts, samples, t)
    subannotation = subanot(samples, t)
    results = run_deseq(subcounts, subannotation)
    
    f_results = subset(results, padj < 0.05)
    o_results = f_results[order(f_results$padj),]
    write.csv(o_results, outfile, row.names = TRUE)
  }
  return(results)
} 
``` 

## Run the Analysis using the following: 
``` 
main('bladder-gtex-tcga-clean.csv', 'bladder-gtex-tcga.comparison.csv', 'bladder-gtex-tcga-degs.csv')
``` 

## License - this project is licensed under the MIT License.

### Contact information: 
Leiara Rivera 
* leiarar@clemson.edu 
* lrivera@email.sc.edu 
