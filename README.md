# DEG-Analysis-Lab
# Task: Finding Differentially Expressed Genes (DEGs), specifically analyzing BRCA to identify genes with significant expression changes between normal and tumorous tissue. 

## Objective: Use the DESeq2 R package within the R-Studio environment to identify DEGs between normal GTEX and tumor TCGA RNA-Seq profiles in a gene expression matrix (GEM). 


## Upon running the DESeq2, the following steps must be conducted: 
1. Start by setting up the DESeq2 package in your Linux environment. 
2. Download RNA-seq data for both normal and tumor for your desired cancer subset (we will use BRCA-breast tissue samples). 
    a. Then, process and combine them into a single Gene Expression Matrix (GEM) for analysis. 
3. Set up a comparison matrix to clearly define normal vs. tumor. 


## Step 1: DESeq2 Installation in R: 
``` 

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2") 

```

## Step 2: Set the working directory: 
``` 

setwd("/scratch/leiarar/BRCA-deg-analysis") #information to be used as an example! 


```


# After installation in R Studio: 
4. Use DESeq2 to analyze the data and identify genes that show significant differences in expression between the two conditions. 
5. After processing both conditions, merge and preprocess the matrices to explore similarities and differences across cancers. 

# Running the Script: 
To execute the DEG analysis, follow these steps: 
1. Ensure all required files are prepared: 
* 'breast-gtex-tcga-clean.csv' is the cleaned processed GEM file. 
* 'breast-gtex-tcga.comparison.csv' is the comparison of the sample groups normal versus tumorous. 
2. Open R-Studio and navigate to the working directory containing these files. 
3. Load DESeq2. 
4. Run the script. 
5. The output file 'breast-gtex-tcga-degs.csv' will be generated containing all the differentially expressed results of the two conditions: normal vs. tumor for BRCA. 


## Step 3: Extract Data and Run DESeq2 
### Extract a sub-matrix of counts for each group: 
```
subgem <- function(gem, anot, group) {
  group_samples <- subset(anot, Comparison == group)$Sample
  group_samples <- intersect(colnames(gem), group_samples)
  
  if (length(group_samples) == 0) {
    stop("No matching samples in count matrix for group: ", group)
  }

  subcounts <- gem[, group_samples, drop = FALSE]
  return(subcounts)
}
```

### Extract a subset of the sample annotation matrix for each group:
```
subanot <- function(anot, group) {
  return(subset(anot, Comparison == group))
}
```
### Run DESeq2 and return DEG stats & normalized expression: 
``` 
run_deseq <- function(counts, annotation) {
  annotation <- annotation[!duplicated(annotation$Sample), ]
  rownames(annotation) <- annotation$Sample

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = annotation,
                                design = ~ Group)
  
  dds <- dds[rowSums(counts(dds)) >= 50, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("Group", "BRCA_GTEX", "BRCA_TCGA"))
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)
```

  ### Get normalized counts (FPM):
 ```
norm <- as.data.frame(fpm(dds))
  norm$GeneID <- rownames(norm)
  ```

  ### Merge columns and specify their specific headers: 
```
combined <- merge(res_df, norm, by = "GeneID")
  combined <- combined[, c("GeneID", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                           setdiff(colnames(combined), c("GeneID", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")))]
  return(combined)
}
```
### Main defines the input functions 'gene_counts_file', 'sample_annotation_file', 'output_file': 
```
main <- function(gene_counts_file, sample_annotation_file, output_file) {         #three input arguments are gene_counts_file, sample_annotation_file, and deseq_results_file
  message("Reading count matrix...")
  counts <- read.csv(gene_counts_file, row.names = "Hugo_Symbol", check.names = FALSE)
  colnames(counts) <- trimws(colnames(counts))    #helps to remove any whitespace from column names and helps to prevent mismatches 

  message("Reading sample annotation...")    
  samples <- read.csv(sample_annotation_file, check.names = FALSE)
  samples$Sample <- trimws(samples$Sample)

  matching_samples <- intersect(colnames(counts), samples$Sample)   #helps to ensure that sample IDs are present in both the count matrix and the annotation 
  if (length(matching_samples) == 0) stop("No overlapping samples found.")

  groups <- unique(samples$Comparison)
  results_list <- list()

  for (group in groups) {
    message("Processing group: ", group)
    subcounts <- subgem(counts, samples, group)
    subannot <- subanot(samples, group)
    subannot <- subannot[subannot$Sample %in% colnames(subcounts), ]

    message("Running DESeq2 for group: ", group)   #prints number of samples and genes for the group before processing DESeq2
    res <- run_deseq(subcounts, subannot)    #will return the results for each gene

    res_filtered <- subset(res, padj < 0.05)
    res_sorted <- res_filtered[order(res_filtered$padj), ]

    results_list[[group]] <- res_sorted
  }
```
  ### Combine all group results for processing of the data: 
```
if (length(results_list) > 0) {
    all_results <- do.call(rbind, results_list)
    message("Writing results to file: ", output_file)
    write.table(all_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {       #save results in a '.tsv' file 
    message("No significant DEGs found.")
  }

  return(invisible(results_list))
}
```
## Step 4:  Run the analysis (output_file = 'BRCA-deg-analysis.tsv'): 
```
main("breast-gtex-tcga-clean.csv", "breast-gtex-tcga.comparison.csv", "BRCA-deg-analysis.tsv")

```

## Step 5: Analyze results

## License - this project is licensed under the MIT License.

### Contact information: 
Leiara Rivera 
* leiarar@clemson.edu 
* lrivera@email.sc.edu 
