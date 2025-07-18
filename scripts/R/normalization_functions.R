# Comprehensive Microbiome Data Normalization Function
# Author: Juan Jovel (Use them at your own risk)
# Description: Normalizes microbiome abundance data using multiple methods

normalize_microbiome_data <- function(df, 
                                    methods = c("relative", "rarefied", "clr", "tss", "css", "tmm", "deseq2_vst"),
                                    rarefaction_depth = NULL,
                                    seed = 123,
                                    return_list = TRUE,
                                    verbose = TRUE) {
  
  # Load required libraries
  required_packages <- c("vegan", "compositions", "DESeq2", "edgeR", "metagenomeSeq")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (verbose) {
        cat("Installing missing package:", pkg, "\n")
      }
      
      # Handle Bioconductor packages
      if (pkg %in% c("DESeq2", "edgeR", "metagenomeSeq")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, quiet = TRUE)
      } else {
        install.packages(pkg, quiet = TRUE)
      }
      
      library(pkg, character.only = TRUE)
    }
  }
  
  # Input validation
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame")
  }
  
  if (ncol(df) < 2) {
    stop("Data frame must have at least 2 columns (taxa names + at least 1 sample)")
  }
  
  # Extract taxa names and count matrix
  taxa_names <- df[, 1]
  count_matrix <- as.matrix(df[, -1])
  rownames(count_matrix) <- taxa_names
  
  # Validate count data
  if (any(count_matrix < 0, na.rm = TRUE)) {
    stop("Count matrix contains negative values")
  }
  
  if (any(!is.finite(count_matrix), na.rm = TRUE)) {
    stop("Count matrix contains non-finite values")
  }
  
  # Remove taxa with zero counts across all samples
  non_zero_taxa <- rowSums(count_matrix, na.rm = TRUE) > 0
  if (sum(non_zero_taxa) < nrow(count_matrix)) {
    if (verbose) {
      cat("Removing", sum(!non_zero_taxa), "taxa with zero counts across all samples\n")
    }
    count_matrix <- count_matrix[non_zero_taxa, ]
    taxa_names <- taxa_names[non_zero_taxa]
  }
  
  # Initialize results list
  normalized_data <- list()
  
  if (verbose) {
    cat("Normalizing microbiome data using methods:", paste(methods, collapse = ", "), "\n")
    cat("Original data dimensions:", nrow(count_matrix), "taxa x", ncol(count_matrix), "samples\n")
  }
  
  # 1. Relative abundance (proportional) normalization
  if ("relative" %in% methods) {
    if (verbose) cat("Computing relative abundance normalization...\n")
    
    rel_abundance <- sweep(count_matrix, 2, colSums(count_matrix, na.rm = TRUE), "/")
    rel_abundance[is.na(rel_abundance)] <- 0
    
    rel_df <- data.frame(Taxa = rownames(rel_abundance), rel_abundance, stringsAsFactors = FALSE)
    normalized_data[["relative_abundance"]] <- rel_df
  }
  
  # 2. Rarefaction (subsampling) normalization
  if ("rarefied" %in% methods) {
    if (verbose) cat("Computing rarefaction normalization...\n")
    
    # Determine rarefaction depth
    if (is.null(rarefaction_depth)) {
      rarefaction_depth <- min(colSums(count_matrix, na.rm = TRUE))
      if (verbose) cat("Using minimum sample depth for rarefaction:", rarefaction_depth, "\n")
    }
    
    # Check if rarefaction depth is feasible
    min_depth <- min(colSums(count_matrix, na.rm = TRUE))
    if (rarefaction_depth > min_depth) {
      warning("Rarefaction depth (", rarefaction_depth, ") exceeds minimum sample depth (", 
              min_depth, "). Using minimum depth.")
      rarefaction_depth <- min_depth
    }
    
    set.seed(seed)
    # Transpose for vegan::rrarefy (expects samples as rows)
    rarefied_t <- vegan::rrarefy(t(count_matrix), sample = rarefaction_depth)
    rarefied_matrix <- t(rarefied_t)
    
    rarefied_df <- data.frame(Taxa = rownames(rarefied_matrix), rarefied_matrix, stringsAsFactors = FALSE)
    normalized_data[["rarefied"]] <- rarefied_df
  }
  
  # 3. Centered Log-Ratio (CLR) transformation
  if ("clr" %in% methods) {
    if (verbose) cat("Computing CLR transformation...\n")
    
    # Add pseudocount to avoid log(0)
    pseudocount <- 1
    clr_matrix <- count_matrix + pseudocount
    
    # Apply CLR transformation
    clr_transformed <- compositions::clr(t(clr_matrix))
    clr_matrix <- t(clr_transformed)
    
    clr_df <- data.frame(Taxa = rownames(clr_matrix), clr_matrix, stringsAsFactors = FALSE)
    normalized_data[["clr"]] <- clr_df
  }
  
  # 4. Total Sum Scaling (TSS) - same as relative abundance but with scaling factor
  if ("tss" %in% methods) {
    if (verbose) cat("Computing TSS normalization...\n")
    
    # Scale to 1 million reads (common convention)
    scaling_factor <- 1e6
    tss_matrix <- sweep(count_matrix, 2, colSums(count_matrix, na.rm = TRUE), "/") * scaling_factor
    tss_matrix[is.na(tss_matrix)] <- 0
    
    tss_df <- data.frame(Taxa = rownames(tss_matrix), tss_matrix, stringsAsFactors = FALSE)
    normalized_data[["tss"]] <- tss_df
  }
  
  # 5. Cumulative Sum Scaling (CSS) from metagenomeSeq
  if ("css" %in% methods) {
    if (verbose) cat("Computing CSS normalization...\n")
    
    tryCatch({
      # Create MRexperiment object
      phenoData <- AnnotatedDataFrame(data.frame(row.names = colnames(count_matrix)))
      mr_obj <- newMRexperiment(count_matrix, phenoData = phenoData)
      
      # Calculate normalization factors
      p <- cumNormStatFast(mr_obj)
      mr_obj <- cumNorm(mr_obj, p = p)
      
      # Extract normalized counts
      css_matrix <- MRcounts(mr_obj, norm = TRUE, log = FALSE)
      
      css_df <- data.frame(Taxa = rownames(css_matrix), css_matrix, stringsAsFactors = FALSE)
      normalized_data[["css"]] <- css_df
    }, error = function(e) {
      if (verbose) cat("CSS normalization failed:", e$message, "\n")
    })
  }
  
  # 6. TMM (Trimmed Mean of M-values) normalization from edgeR
  if ("tmm" %in% methods) {
    if (verbose) cat("Computing TMM normalization...\n")
    
    tryCatch({
      # Create DGEList object
      dge <- DGEList(counts = count_matrix)
      
      # Calculate normalization factors
      dge <- calcNormFactors(dge, method = "TMM")
      
      # Get normalized counts (CPM with normalization factors)
      tmm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
      
      tmm_df <- data.frame(Taxa = rownames(tmm_matrix), tmm_matrix, stringsAsFactors = FALSE)
      normalized_data[["tmm"]] <- tmm_df
    }, error = function(e) {
      if (verbose) cat("TMM normalization failed:", e$message, "\n")
    })
  }
  
  # 7. DESeq2 Variance Stabilizing Transformation
  if ("deseq2_vst" %in% methods) {
    if (verbose) cat("Computing DESeq2 VST...\n")
    
    tryCatch({
      # Create sample metadata
      coldata <- data.frame(condition = rep("A", ncol(count_matrix)), 
                           row.names = colnames(count_matrix))
      
      # Create DESeq2 object
      dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                   colData = coldata,
                                   design = ~ 1)  # No design for normalization only
      
      # Apply variance stabilizing transformation
      vst_matrix <- assay(vst(dds, blind = TRUE))
      
      vst_df <- data.frame(Taxa = rownames(vst_matrix), vst_matrix, stringsAsFactors = FALSE)
      normalized_data[["deseq2_vst"]] <- vst_df
    }, error = function(e) {
      if (verbose) cat("DESeq2 VST failed:", e$message, "\n")
    })
  }
  
  # Print summary statistics
  if (verbose) {
    cat("\n=== Normalization Summary ===\n")
    for (method_name in names(normalized_data)) {
      norm_matrix <- as.matrix(normalized_data[[method_name]][, -1])
      cat(sprintf("%-15s: Range [%.3f, %.3f], Mean = %.3f\n", 
                  method_name, 
                  min(norm_matrix, na.rm = TRUE), 
                  max(norm_matrix, na.rm = TRUE),
                  mean(norm_matrix, na.rm = TRUE)))
    }
  }
  
  # Return results
  if (return_list) {
    return(normalized_data)
  } else {
    # If only one method requested, return single data frame
    if (length(normalized_data) == 1) {
      return(normalized_data[[1]])
    } else {
      return(normalized_data)
    }
  }
}

# Helper function to create example data
create_example_microbiome_data <- function(n_taxa = 50, n_samples = 20, seed = 123) {
  set.seed(seed)
  
  # Generate taxa names
  taxa_names <- paste0("Taxa_", sprintf("%03d", 1:n_taxa))
  
  # Generate sample names
  sample_names <- paste0("Sample_", sprintf("%02d", 1:n_samples))
  
  # Generate count data with realistic microbiome distribution
  # Using negative binomial to simulate overdispersed count data
  count_data <- matrix(0, nrow = n_taxa, ncol = n_samples)
  
  for (i in 1:n_taxa) {
    for (j in 1:n_samples) {
      # Simulate realistic microbiome counts with high variation
      mean_count <- rpois(1, lambda = runif(1, 10, 1000))
      count_data[i, j] <- rnbinom(1, size = 2, mu = mean_count)
    }
  }
  
  # Add some zeros (common in microbiome data)
  zero_indices <- sample(1:(n_taxa * n_samples), size = n_taxa * n_samples * 0.3)
  count_data[zero_indices] <- 0
  
  # Create data frame
  df <- data.frame(Taxa = taxa_names, count_data, stringsAsFactors = FALSE)
  colnames(df) <- c("Taxa", sample_names)
  
  return(df)
}

# Example usage and demonstration
if (FALSE) {  # Set to TRUE to run examples
  
  # Create example data
  example_data <- create_example_microbiome_data(n_taxa = 30, n_samples = 10)
  print("Example data (first 5 rows):")
  print(example_data[1:5, 1:6])
  
  # Normalize using all methods
  normalized_results <- normalize_microbiome_data(
    example_data,
    methods = c("relative", "rarefied", "clr", "tss", "tmm"),
    verbose = TRUE
  )
  
  # Access specific normalized data
  relative_abundance <- normalized_results$relative_abundance
  rarefied_data <- normalized_results$rarefied
  
  # Normalize using single method
  clr_only <- normalize_microbiome_data(
    example_data,
    methods = "clr",
    return_list = FALSE
  )
  
  print("CLR normalized data (first 5 rows):")
  print(clr_only[1:5, 1:6])
}
