library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(pheatmap)

setwd('/Users/juanjovel/jj/data_analysis/RNA_training_careem/data')
data_file = 'control_vs_treatment_10up.tsv'
metadata_file = 'control_vs_treatment_meta.tsv'

data = read.table(data_file, sep = '\t', header = T, row.names = 1)
metadata = read.table(metadata_file, sep = '\t', header = T, row.names = 1)

data_matrix = as.matrix(round(data, digits = 0))
class(data_matrix)

# Create DESeq2 object to contain data and metadata
dss = DESeqDataSetFromMatrix(countData = data_matrix, colData = metadata,
                             design =~ group )

# Calculate size factor (library size) for normalization
dss = estimateSizeFactors(dss)

# Apply regularized logarithmic transformation
rld = rlog(dss, fitType = "local")

# # Apply variance stabilizing transformation
vsd = vst(dss, fitType = "local")

# Calculate Euclidean distances from transposed
# normalized data (either rld or vsd)
distances = dist(t(assay(rld)))

# Convert dist object to a matrix (required by pheatmap)
dist_matrix = as.matrix(distances)

# OPTIONAL, create a vector with formatted names to be
# deployed on heatmap
new_labels <- c(rep("Treatment", 5), rep("Control", 5))

# Pass vectors to appear as labels on columns
# and rows of heatmap
colnames(dist_matrix) = metadata$label
rownames(dist_matrix) = new_labels

# Create a custom color palette (255 + 1 is the number of breaks
# in the color scale)
#colors = colorRampPalette(c("red","white","blue"))(255)

# Use a predefined color palette from the package RColorBrewer
# 9 indicates the intensity of the colors (ranges from 3 to 11)
colors = colorRampPalette(brewer.pal(9, "Reds"))(255)

# Create hierarchichal clustering heatmap
pheatmap(dist_matrix, 
         clustering_distance_cols = distances,
         clustering_distance_rows = distances,
         col=colors,
         cex=1.1)

# Let's conduct Principal Component Analyses (PCA)

# Calculate the first two Principal Components (Eigen's vectors)
data = plotPCA(rld, intgroup = c("group","label"), returnData = T)

# Calculate percentage of variance for each principal component
percentVar = round(100 * attr(data, "percentVar"))

# OPTIONAL: Create a formatted vector for colors
Group = c(rep("Control", 5), rep("Treatment", 5))

# Create first layer of plot which will include the x and y values as well
# as the factor used for coloring (x and y labels are defined too)
PCA_EucDist = ggplot(data, aes(x=data$PC1, y=data$PC2, 
                     color = Group)) +
                       xlab(paste("PC1 :", percentVar[1], "% variance")) +
                        ylab(paste("PC2 :", percentVar[2], "% variance"))
  
# Add the points (dots) layer, and define shape and size of dots
PCA_EucDist + geom_point(shape = 20, size = 7) +
  # Add labels to each sample in the plot
  geom_text(aes(label=metadata$label), hjust=0.5, vjust=2, color="black")


##### DIFFERENTIAL EXPRESSION ANALYSIS ####
dss = DESeq(dss)
resultsNames(dss)

de_res = results(dss)

sign_res = subset(de_res, padj < 0.05)
head(sign_res)


sign_data = cbind(Transcript=rownames(sign_res), sign_res)
head(sign_data)
# Export all DE expression results + normalized data
all_results_file = gsub('.tsv', '_allRes.tsv', data_file)
all_results = merge(as.data.frame(de_res), as.data.frame(counts(dss, normalized=T)),
                    by='row.names', sort=F)

# Export results to a file
write.table(all_results, all_results_file, sep = '\t', quote = F, row.names = F)


##### Added on 230706 #####

# Adding annotations to the DE results from Biomart

library(biomaRt)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
transcripts = sign_data$Transcript
transcripts_clean = gsub("\\.\\d+$", "", transcripts)
attributes = c("ensembl_transcript_id", 
                "hgnc_symbol", 
                "wikigene_description", 
                "name_1006", 
                "definition_1006")

records <- getBM(attributes = attributes, filters = "ensembl_transcript_id", 
                values = transcripts_clean, mart = mart)

first_annotations_df <- data.frame()

# Loop through each transcript ID
for (transcript_id in transcripts_clean) {
  # Filter the records for the current transcript ID
  filtered_records <- records[records$ensembl_transcript_id == transcript_id, ]
  
  # Get the first annotation record for the current transcript ID
  first_annotation_df     <- as.data.frame(filtered_records[1, ])
  #first_annotation_df  <- as.data.frame(first_annotation)
  first_annotations_df <- rbind(first_annotations_df, first_annotation_df)
  
}

# Merge unique records with your original table
full_table <- cbind(sign_data, first_annotations_df)

# Export DE results
results_file = gsub('.tsv', '_q0.05_withAnnotations.tsv', data_file)

write.table(full_table, results_file, sep = '\t', quote = F, row.names = F)


library(EnhancedVolcano)

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  full_table$log2FoldChange < -1, 'forestgreen',
  ifelse(full_table$log2FoldChange > 1, 'firebrick2',
         'dodgerblue'))
names(keyvals)[keyvals == 'firebrick2'] <- 'Upregulated'
names(keyvals)[keyvals == 'dodgerblue'] <- '-1 < log2FC< 1'
names(keyvals)[keyvals == 'forestgreen'] <- 'Downregulated'


EnhancedVolcano(full_table,
                lab = full_table$hgnc_symbol,
                selectLab = c('IRF4','CUX1','SH3BP2'),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 4,
                colCustom = keyvals,
                colAlpha = 1,
                shape = 19,
                legendLabels=c('Not sig.','Log2FC','p-value',
                               'p-value & Log2FC'),
                legendPosition = 'top',
                legendLabSize = 13,
                legendIconSize = 3.0
                )


