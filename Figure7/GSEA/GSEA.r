library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(svglite)
library(org.Hs.eg.db)

# Prompt user for the x value (cluster number)
x <- as.integer(readline(prompt = "Enter the cluster number (e.g., 8): "))

# Prompt user for the input file paths
input_file_24h <- readline(prompt = "Enter the path to the 24h log2fc file (e.g., '~/desktop/DEG/cluster8_24log2fc'): ")
df_24h <- read.delim(input_file_24h, header = FALSE, sep = "\t")
original_gene_list_24h <- df_24h$V2
names(original_gene_list_24h) <- df_24h$V1
gene_list_24h = sort(original_gene_list_24h, decreasing = TRUE)

# Perform GO enrichment analysis for 24h
gse_24h <- gseGO(geneList=gene_list_24h, 
                 ont ="BP", 
                 keyType = "SYMBOL", 
                 minGSSize = 3, 
                 maxGSSize = 800,
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "none")

geneRatio_24h <- sapply(strsplit(gse_24h$core_enrichment, "/"), function(x) length(x)) / gse_24h$setSize
selected_columns_24h <- gse_24h[, c("Description", "pvalue", "NES")]
selected_columns_24h$GeneRatio <- geneRatio_24h
selected_columns_24h$Sign <- ifelse(selected_columns_24h$NES > 0, "Positive", "Negative")

# Prompt user for output file for 24h
output_file_24h <- readline(prompt = "Enter the path to save the 24h output (e.g., '~/desktop/cluster8_24h.xlsx'): ")
write.csv(selected_columns_24h, file = output_file_24h, row.names = FALSE)

###########################

# Prompt user for the input file paths for 3h log2fc
input_file_3h <- readline(prompt = "Enter the path to the 3h log2fc file (e.g., '~/desktop/DEG/cluster8_3log2fc'): ")
df_3h <- read.delim(input_file_3h, header = FALSE, sep = "\t")
original_gene_list_3h <- df_3h$V2
names(original_gene_list_3h) <- df_3h$V1
gene_list_3h = sort(original_gene_list_3h, decreasing = TRUE)

# Perform GO enrichment analysis for 3h
gse_3h <- gseGO(geneList=gene_list_3h, 
                ont ="BP", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800,
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "none")

geneRatio_3h <- sapply(strsplit(gse_3h$core_enrichment, "/"), function(x) length(x)) / gse_3h$setSize
selected_columns_3h <- gse_3h[, c("Description", "pvalue", "NES")]
selected_columns_3h$GeneRatio <- geneRatio_3h
selected_columns_3h$Sign <- ifelse(selected_columns_3h$NES > 0, "Positive", "Negative")

# Prompt user for output file for 3h
output_file_3h <- readline(prompt = "Enter the path to save the 3h output (e.g., '~/desktop/cluster8_3h.xlsx'): ")
write.csv(selected_columns_3h, file = output_file_3h, row.names = FALSE)

##########################################

# Prompt user for the input file paths for 1h log2fc
input_file_1h <- readline(prompt = "Enter the path to the 1h log2fc file (e.g., '~/desktop/DEG/cluster8_1log2fc'): ")
df_1h <- read.delim(input_file_1h, header = FALSE, sep = "\t")
original_gene_list_1h <- df_1h$V2
names(original_gene_list_1h) <- df_1h$V1
gene_list_1h = sort(original_gene_list_1h, decreasing = TRUE)

# Perform GO enrichment analysis for 1h
gse_1h <- gseGO(geneList=gene_list_1h, 
                ont ="BP", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800,
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "none")

geneRatio_1h <- sapply(strsplit(gse_1h$core_enrichment, "/"), function(x) length(x)) / gse_1h$setSize
selected_columns_1h <- gse_1h[, c("Description", "pvalue", "NES")]
selected_columns_1h$GeneRatio <- geneRatio_1h
selected_columns_1h$Sign <- ifelse(selected_columns_1h$NES > 0, "Positive", "Negative")

# Prompt user for output file for 1h
output_file_1h <- readline(prompt = "Enter the path to save the 1h output (e.g., '~/desktop/cluster8_1h.xlsx'): ")
write.csv(selected_columns_1h, file = output_file_1h, row.names = FALSE)
