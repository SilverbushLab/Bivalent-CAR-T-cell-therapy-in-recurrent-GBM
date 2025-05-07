library(tidyr)
library(purrr)

"%!in%" <- Negate("%in%")

cd4cd8_car_present = function(
  obj
){
  cd8cd4 = FetchData(obj, c("CD8A", "CD8B", "CD4"), slot = "counts")
  ct.cd8cd4 = cd8cd4 %>% mutate(
    CD4CD8_BY_EXPRS = case_when(
      CD4 > 0 & CD8A == 0 & CD8B == 0 ~ "CD4+CD8-",
      CD4 == 0 & (CD8A > 0 | CD8B > 0) ~ "CD4-CD8+",
      CD4 == 0 & CD8A == 0 & CD8B == 0 ~ "CD4-CD8-",
      CD4 > 0 & (CD8A > 0 | CD8B > 0) ~ "CD4+CD8+",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD4CD8_BY_EXPRS)
  rownames(ct.cd8cd4) = rownames(cd8cd4)
  obj = AddMetaData(obj, ct.cd8cd4)

  cd3 = FetchData(obj, c("CD3D", "CD3E", "CD3G"), slot = "counts")
  cd3 = cd3 %>% mutate(
    CD3_BY_EXPRS = case_when(
      CD3D > 0 | CD3E > 0 | CD3G > 0 ~ "CD3",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD3_BY_EXPRS)
  obj = AddMetaData(obj, cd3)

  if("CAR" %in% rownames(GetAssayData(obj, slot = c("data"), assay = "SCT"))){
    car.ftr = FetchData(obj, c("CAR"), slot = "data", assay = "SCT")
    obj$CAR_BY_EXPRS = as.factor(car.ftr$`CAR` > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }

  obj
}


list_to_df <- function(named_list) {
  max_length <- max(sapply(named_list, length))
  padded_list <- lapply(named_list, function(x) c(x, rep('', max_length - length(x))))
  data.frame(padded_list, check.names = FALSE)
}
                        

#' Generate and Export Top Markers for Clusters
#'
#' This function filters markers based on adjusted p-value, selects top markers for each cluster,
#' and exports the results to a CSV file.
#'
#' @param oupMarker A dataframe containing marker information with columns 'cluster', 'gene', 'p_val_adj', and 'pct.diff'.
#' @param p_val_threshold Numeric. The adjusted p-value threshold for filtering markers. Default is 0.05.
#' @param top_n Integer. The number of top markers to select for each cluster. Default is 30.
#' @param output_file Character. The file path for the output CSV. Default is "results/top_cluster_markers.csv".
#'
#' @return A dataframe of top markers for each cluster in wide format.
#'
#' @importFrom dplyr filter group_by slice_max mutate ungroup across
#' @importFrom tidyr pivot_wider
#' @importFrom purrr map_chr
#' @export
generate_top_cluster_markers <- function(oupMarker, 
                                         p_val_threshold = 0.05, order_by = 'avg_log2FC', asc_desc=1,
                                         top_n = 30, print_table = TRUE,
                                         output_file = "results/top_cluster_markers.csv") {
  
  all_markers <- oupMarker %>% filter(p_val_adj < p_val_threshold)
  
  df_wide_top <- all_markers %>%
    dplyr::select(cluster, gene, !!sym(order_by)) %>%
    group_by(cluster) %>%
    slice_max(order_by = asc_desc * !!sym(order_by), n = top_n, with_ties = FALSE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  # Optionally print cluster distribution
  if (print_table) {
      print(table(df_wide_top$cluster))
      print(as.data.frame(df_wide_top[c('cluster','gene')] %>% 
        pivot_wider(names_from = cluster, values_from = gene, values_fn = ~ paste(.x, collapse = ", "))))
  }
  
  df_wide_top <- df_wide_top[, c('cluster', 'gene', 'rank')] %>%
    pivot_wider(id_cols = rank, 
                names_from = cluster, 
                values_from = gene, 
                values_fn = list) %>%
    mutate(across(-rank, ~map_chr(., ~paste(.x, collapse = ", "))))
  
  df_wide_top <- df_wide_top[-c(1)]
  
  # Export to CSV
  write.csv(df_wide_top, file = output_file, row.names = FALSE)
  
  return(df_wide_top)
}



#' Perform NMF Analysis on Seurat Object and Visualize Results
#'
#' This function performs Non-negative Matrix Factorization (NMF) analysis on a Seurat object,
#' generates metaprograms, visualizes results, and performs GSEA on top programs.
#'
#' @param seurat_obj A Seurat object to analyze
#' @param split_by Character. The column name to split the object by. Default is "sample_id".
#' @param nmf_k Numeric vector. Range of k values for NMF. Default is 5:15.
#' @param nmf_nfeatures Integer. Number of features to use in NMF. Default is 5000.
#' @param meta_nprograms Integer. Number of metaprograms to generate. Default is 10.
#' @param meta_max_genes Integer. Maximum number of genes per metaprogram. Default is 50.
#' @param reduction String. Name of the reduction to use for visualization. Default is "harmony_UMAP".
#' @param output_file Character. File path for saving NMF programs. Default is "results/nmf_programs.csv".
#'
#' @return A list containing the updated Seurat object and NMF results
#'
#' @import Seurat
#' @import SeuratWrappers
#' @import ggplot2
#' @importFrom SeuratObject AddModuleScore_UCell
#'
#' @export
perform_nmf_analysis <- function(seurat_obj,
                                 split_by = "sample_id",
                                 nmf_k = 5:15,
                                 nmf_nfeatures = 5000,
                                 meta_nprograms = 10,
                                 meta_max_genes = 50,
                                 reduction = "harmony_UMAP",
                                 output_file = "results/nmf_programs.csv") {
  
  # Split object
  seu.list <- SplitObject(seurat_obj, split.by = split_by)
  
  # Perform multiNMF
  geneNMF.programs <- multiNMF(seu.list, slot="data", k=nmf_k, L1=c(0,0), 
                               do_centering=TRUE, nfeatures = nmf_nfeatures)
  
  # Generate metaprograms
  geneNMF.metaprograms <- getMetaPrograms(
    geneNMF.programs, nprograms=meta_nprograms, max.genes=meta_max_genes, 
    hclust.method="ward.D2", min.confidence=0.1
  )
  
  # Print metaprogram metrics
  print(geneNMF.metaprograms[['metaprograms.metrics']])
  
  # Save metaprograms to CSV
  mps <- list_to_df(geneNMF.metaprograms$metaprograms.genes)
  write.csv(mps, file = output_file, row.names=FALSE)
  
  # Perform GSEA on top programs
  top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
    runGSEA(program, universe=rownames(seurat_obj), category = "C5", subcategory='BP')
  })
  
  # Add module scores
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features = mp.genes, ncores=4, name = "")
  
  # Visualize feature plots
  options(repr.plot.width = 12, repr.plot.height = 12, repr.plot.res=150)
  print(FeaturePlot(seurat_obj, features = names(mp.genes), reduction = reduction, 
                    ncol=3, min.cutoff='q25', max.cutoff='q75') &
          theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank()) & 
          xlab('') & ylab(''))
  
  # Return results
  return(list(seurat_obj = seurat_obj, 
              nmf_results = geneNMF.metaprograms, 
              gsea_results = top_p))
}



# Function modified from Gregory M Chen, Written By Alice Wang

FindDoublets <- function(seurat.rna, PCs = 1:50, exp_rate = 0.02, sct = FALSE, multiome=FALSE, v5=FALSE){
  # sct--do SCTransform or not
  DefaultAssay(seurat.rna) <- "RNA"
  
  ## pK identification
  sweep.res.list <- paramSweep(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_infer <- 0.09
  
  if (multiome == TRUE)  {
    pK_infer <- bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric), ]$pK
    pK_infer <- as.numeric(levels(pK_infer))[pK_infer]
  }
  
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)          
  exp_rate <- (dim(seurat.rna)[2] / 2000) * exp_rate
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = pk_infer, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = pk_infer, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  return(seurat.rna)
}

add_totals_to_plot <- function(p, column_name) {
  # Calculate totals
  totals <- aggregate(Freq ~ get(column_name), data = p$data, FUN = sum)
  names(totals)[1] <- column_name

  # Create new plot with totals
  p_with_totals <- p +
    geom_text(data = totals, 
              aes(x = .data[[column_name]], y = 1.01, label = Freq), 
              size = 4,
              inherit.aes = FALSE,
              position = position_dodge(width = 0),
              vjust = 0) +
    scale_y_continuous(labels = scales::percent_format(), 
                       limits = c(0, 1.05),
                       expand = c(0, 0)) +
    theme( 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p_with_totals)
}



analyze_gene_expression <- function(seurat_object, gene, cluster_threshold = 20, expression_threshold = 0) {
  # Set default assay
  DefaultAssay(seurat_object) <- 'RNA'
  
  # Initialize Gene_Expression column
  seurat_object$Gene_Expression <- paste(gene, '-')
  
  # Mark cells as gene+ if expression is > 0
  seurat_object$Gene_Expression[GetAssayData(object = seurat_object, layer='data')[gene,] > expression_threshold] <- paste(gene, '+')
  
  # Create boolean column for gene positive cells
  seurat_object$gene_positive <- seurat_object$Gene_Expression == paste(gene, '+')
  
  # Calculate percentage of gene+ cells per cluster
  cluster_percentages <- aggregate(gene_positive ~ seurat_clusters, 
                                   data = seurat_object@meta.data, 
                                   FUN = function(x) mean(x) * 100)
  
  # Assign cluster gene status based on threshold
  cluster_percentages$cluster_gene_status <- ifelse(cluster_percentages$gene_positive > cluster_threshold, 
                                                    paste(gene, '+'), paste(gene, '-'))
  
  # Set row names for easier indexing
  rownames(cluster_percentages) <- cluster_percentages$seurat_clusters
  
  # Add cluster percentages to Seurat object metadata
  seurat_object <- AddMetaData(seurat_object, 
                               cluster_percentages[as.character(seurat_object@meta.data[['seurat_clusters']]),])
  
  # Create plot
  plot <- ggplot(cluster_percentages, aes(x = as.factor(seurat_clusters), y = gene_positive, fill = cluster_gene_status)) +
    geom_bar(stat = "identity") + 
    ylim(0,100) +
    geom_hline(yintercept = cluster_threshold, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(x = "Seurat Cluster", y = paste("Percentage of", gene, "+ Cells"), fill = "Cluster Status") +
    ggtitle(paste(gene)) + 
    mytheme  # Note: ensure mytheme is defined in your environment
  
  # Return a list containing the updated Seurat object and the plot
  return(list(seurat_object = seurat_object, plot = plot))
}



#' Preprocess a Seurat object for single-cell RNA-seq analysis
#'
#' This function performs a series of preprocessing steps on a Seurat object,
#' including normalization, variable feature selection, cell cycle scoring,
#' and dimensionality reduction.
#'
#' @param obj A Seurat object to be processed.
#' @param g2m_genes A vector of G2M phase genes for cell cycle scoring.
#' @param s_genes A vector of S phase genes for cell cycle scoring.
#' @param nfeatures The number of variable features to select (default: 4000).
#'
#' @return A list containing two elements:
#'   \item{obj}{The processed Seurat object}
#'   \item{var_minus_cycling}{Variable features with cell cycle genes removed}
#'
#' @details
#' The function performs the following steps:
#' 2. Normalizes data
#' 3. Finds variable features
#' 4. Performs TCR gene processing
#' 5. Scores cell cycle phases
#' 6. Removes cell cycle effects from variable features
#' 7. Scales data
#' 8. Runs PCA
#'
#' @note This function requires the Seurat and SingleCellExperiment packages to be installed.
#'
#' @examples
#' \dontrun{
#' result <- preprocess_seurat_object(seurat_obj, g2m_genes, s_genes)
#' processed_obj <- result$obj
#' var_features <- result$var_minus_cycling
#' }
#'
#' @export
                               
preprocess_seurat_object <- function(obj, g2m_genes, s_genes,  nfeatures = 4000, phase_threshold=1) {

  # Normalize data
  obj <- NormalizeData(obj)
  
  # Find variable features
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  
  # Create a copy for TCR gene processing
  objv5 <- obj
  objv5[["RNA"]] <- as(object = objv5[["RNA"]], Class = "Assay")
  objv5 <- quietTCRgenes(objv5, assay = 'RNA')

  per.batch.var.feats <- VariableFeatures(objv5)
  objv5 <- FindVariableFeatures(objv5, nfeatures = nfeatures)
  objv5 <- quietTCRgenes(objv5, assay = 'RNA')
  joint.var.feats <- VariableFeatures(objv5)
  var.feats <- union(joint.var.feats, per.batch.var.feats)
  
  # Perform cell cycle scoring
  objv5 <- CellCycleScoring(objv5, g2m.features = g2m_genes, s.features = s_genes, verbose = FALSE)
  
  # Convert to SingleCellExperiment object
  obj.sce <- as.SingleCellExperiment(objv5)
  
  # Get variance explained by cell cycle phase
  diff <- getVarianceExplained(obj.sce, "Phase")
  discard <- diff > phase_threshold
  
  # Filter out cell cycle genes from variable features
  var_minus_cycling <- var.feats[var.feats %!in% rownames(objv5)[which(discard)]]
  
  # Clean up memory
  rm(objv5)

  VariableFeatures(obj) <- var_minus_cycling
  
  # Scale data
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  
  # Run PCA
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  
  # Return processed object and filtered variable features
  return(obj)
}



#' Create Cell Type Labels Based on Gene Expression
#'
#' This function processes a Seurat object to create cell type labels based on gene expression,
#' marking cells as positive or negative for a specific gene.
#'
#' @param seurat_obj A Seurat object containing the single-cell data.
#' @param gene Character string specifying the gene to analyze.
#' @param assay Character string specifying which assay to use. Default is 'RNA'.
#' @param slot Character string specifying which slot to use ('counts', 'data', or 'scale.data').
#'             Default is 'data'.
#' @param threshold Numeric value specifying the expression threshold for positive cells.
#'                 Default is 0.
#'
#' @return A list containing:
#'        \itemize{
#'          \item counts: A table of marker type counts by sample
#'          \item target_ct: The name of the positive cell type (e.g., "CD8B+")
#'          \item seurat_obj: The modified Seurat object with new marker_type column
#'        }
#'
#' @details The function adds a 'marker_type' column to the Seurat object's metadata,
#' labeling cells as either positive or negative for the specified gene based on
#' the expression threshold.
#'
#' @examples
#' # Process Seurat object for CD8B expression
#' results <- create_marker_types(carneg.tcells, "CD8B")
#'
#' # Access results
#' counts <- results$counts
#' target_ct <- results$target_ct
#' modified_seurat <- results$seurat_obj
#'
#' @import Seurat
#'
#' @export
create_marker_types <- function(seurat_obj, gene, slot = "data", threshold = 0) {
    
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (!is.character(gene) || length(gene) != 1) {
    stop("gene must be a single character string")
  }
  
  if (!gene %in% rownames(seurat_obj)) {
    stop(sprintf("gene '%s' not found in the Seurat object", gene))
  }
  
  if (!slot %in% c("counts", "data", "scale.data")) {
    stop("slot must be one of: 'counts', 'data', or 'scale.data'")
  }
  
  # Get expression data
  data <- GetAssayData(object = seurat_obj, layer = slot)
  
  # Add to Seurat object
  seurat_obj$marker_type <- paste0(gene, '-')
  seurat_obj$marker_type[data[gene,] > threshold] <- paste0(gene, '+')
  
  # Create counts table
  counts <- table(seurat_obj$marker_type, seurat_obj$sample_id)
  
  # Define target cell type
  target_ct <- paste0(gene, '+')
  
  # Return results
  return(list(counts = counts, target_ct = target_ct, seurat_obj = seurat_obj))
}
                                   


#' Analyze and Visualize Cell Proportions Between Timepoints
#'
#' This function analyzes the proportions of specific cell types between two timepoints (D0 and D7),
#' performs statistical tests, and creates a visualization of the changes.
#'
#' @param counts A matrix or table object containing cell counts, with cell types as rows
#'               and sample IDs as columns.
#' @param target_ct Character string specifying the target cell type to analyze (must be
#'                 a row name in the counts matrix).
#' @param timepoint_pairs A data frame with three columns:
#'        \itemize{
#'          \item Patient: Character vector of patient IDs
#'          \item D0: Character vector of sample IDs for Day 0
#'          \item D7: Character vector of sample IDs for Day 7
#'        }
#' @param p_value_y Numeric value specifying the y-position of the p-value annotation.
#'                  Default is 48.
#' @param y_max Numeric value specifying the maximum value for y-axis.
#'             Default is 50.
#'
#' @return A list containing:
#'        \itemize{
#'          \item plot: A ggplot2 object showing the proportion changes
#'          \item t_test: Results of paired t-test on D7/D0 ratios
#'          \item geometric_mean: Geometric mean of the D7/D0 ratios
#'          \item data: Data frame containing the analyzed data
#'        }
#'
#' @details The function calculates the proportion of target cells relative to total cells
#' for each timepoint, performs a two-sided t-test on the D7/D0 ratios (null hypothesis: ratio = 1),
#' and creates a line plot showing the changes in proportions for each patient.
#'
#' @examples
#' # Create example timepoint pairs
#' timepoint_pairs <- data.frame(
#'   Patient = c("P3", "P6", "P2", "P5"),
#'   D0 = c("P3D0", "P6D0", "P2D0", "P5D0_run2"),
#'   D7 = c("P3D7", "P6D7", "P2D7", "P5D7")
#' )
#'
#' # Run analysis with custom y-axis settings
#' # results <- analyze_cell_proportions(counts, target_ct, timepoint_pairs,
#' #                                    p_value_y = 30, y_max = 35)
#'
#' @import ggplot2
#' @importFrom stats t.test
#'
#' @export
analyze_cell_proportions <- function(counts, target_ct, timepoint_pairs, 
                                   p_value_y = 48, y_max = 50) {
  # Input validation
  if (!is.character(target_ct) || length(target_ct) != 1) {
    stop("target_ct must be a single character string")
  }
  
  if (!target_ct %in% rownames(counts)) {
    stop("target_ct must be a valid row name in counts matrix")
  }
  
  if (!all(c("Patient", "D0", "D7") %in% colnames(timepoint_pairs))) {
    stop("timepoint_pairs must have columns: Patient, D0, and D7")
  }
  
  if (!all(c(timepoint_pairs$D0, timepoint_pairs$D7) %in% colnames(counts))) {
    stop("Some sample IDs in timepoint_pairs are not found in counts matrix")
  }
  
  if (!is.numeric(p_value_y)) {
    stop("p_value_y must be numeric")
  }
  
  if (!is.numeric(y_max)) {
    stop("y_max must be numeric")
  }
  
  if (p_value_y >= y_max) {
    warning("p_value_y should be less than y_max for optimal visualization")
  }
  
  # Calculate proportions
  target_celltypes <- counts[target_ct, ]
  total_t_cells <- colSums(counts)
  proportions <- target_celltypes / total_t_cells
  
  # Create data frame for analysis
  data <- data.frame(
    Patient = timepoint_pairs$Patient,
    D0 = proportions[timepoint_pairs$D0],
    D7 = proportions[timepoint_pairs$D7]
  )
  
  # Calculate ratios
  data$ratio <- data$D7 / data$D0
  
  # Perform statistical tests
  result <- t.test(data$ratio, mu = 1, alternative = "two.sided")
  geom_mean <- exp(mean(log(data$ratio)))
  
  # Create plot data
  plot_data <- data.frame(
    Patient = rep(data$Patient, times = 2),
    TimePoint = rep(c("D0", "D7"), each = nrow(data)),
    Proportion = c(data$D0, data$D7)*100
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = TimePoint, y = Proportion, group = Patient, color = Patient)) +
    geom_line() +
    geom_point(size = 2) +
    labs(title = paste0(target_ct, " cells"),
         y = "Percent of CAR- T and NK Cells",
         x = "") +
    annotate("text", x = 0.8, y = p_value_y,
             label = sprintf('italic(p) == %.3f', result$p.value),
             hjust = 0, vjust = 1, parse=TRUE) + 
    scale_color_brewer(palette = "Set1") + 
    theme_classic() +  
    ylim(0, y_max) + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Return results as a list
  return(list(
    plot = p,
    t_test = result,
    geometric_mean = geom_mean,
    data = data
  ))
}



#' Create Correlation Analysis Plot with Statistics
#'
#' This function creates a scatter plot with correlation statistics for two variables,
#' including a trend line and statistical annotations.
#'
#' @param data A data frame containing the required columns.
#' @param x_var Character string specifying the x-axis variable name.
#' @param y_var Character string specifying the y-axis variable name.
#' @param label_var Character string specifying the column to use for point labels.
#' @param x_label Character string for x-axis label. 
#' @param y_label Character string for y-axis label.
#' @param padding Numeric value for axis padding (fraction of range). Default is 0.2.
#' @param stats_x,stats_y Numeric values for positioning correlation statistics.
#'                       If NULL, will be automatically positioned.
#'
#' @return A list containing:
#'        \itemize{
#'          \item plot: The ggplot object
#'          \item correlation: Pearson correlation coefficient
#'          \item p_value: P-value from correlation test
#'          \item lm_model: Linear model object
#'        }
#'
#' @import ggplot2
#' @import ggsci
#'
#' @examples
#' data <- data.frame(
#'   Patient = c("P1", "P2", "P3"),
#'   difference = c(-5, 0, 5),
#'   change_in_tumor = c(10, 20, 30)
#' )
#' 
#' results <- create_correlation_plot(
#'   data,
#'   x_var = "difference",
#'   y_var = "change_in_tumor",
#'   label_var = "Patient",
#'   x_label = "Change in Treg Proportion\n(% of CAR- T and NK cells, D0 vs D7)",
#'   y_label = "Reduction in Tumor Size (%)"
#' )
#'
#' @export
create_correlation_plot <- function(data, x_var, y_var, label_var,
                                  x_label, y_label,
                                  padding = 0.2,
                                  stats_x = NULL, stats_y = NULL) {
  # Input validation
  required_cols <- c(x_var, y_var, label_var)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Calculate axis limits with padding
  x_range <- range(data[[x_var]], na.rm = TRUE)
  y_range <- range(data[[y_var]], na.rm = TRUE)
  
  x_padding <- padding * diff(x_range)
  y_padding <- padding * diff(y_range)
  
  x_limits <- c(x_range[1] - x_padding, x_range[2] + x_padding)
  y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)
  
  # Calculate correlation statistics
  cor_test <- cor.test(data[[x_var]], data[[y_var]], 
                       method = "pearson")
  correlation <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # Fit linear model
  lm_model <- lm(reformulate(x_var, y_var), data = data)
  new_x <- data.frame(x = x_limits)
  names(new_x) <- x_var
  predicted_y <- predict(lm_model, newdata = new_x)
  
  # Create the plot
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(size = 2, color = pal_npg()(2)[[2]]) +
    geom_text(aes_string(label = label_var), 
              vjust = -1, hjust = 0.5) +
    geom_line(data = data.frame(x = x_limits, y = predicted_y),
              aes(x = x, y = y), 
              color = pal_npg()(2)[[1]], 
              linetype = "dashed") +
    labs(x = x_label,
         y = y_label) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(5, 5, 5, 5),
      aspect.ratio = 1
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)
  
  # Add correlation statistics
  if (is.null(stats_x)) stats_x <- x_limits[1] + x_padding/2
  if (is.null(stats_y)) stats_y <- y_limits[1] + y_padding/2
  
  p <- p + annotate(
    "text", 
    x = stats_x, 
    y = stats_y,
    label = sprintf(
      'atop(italic(r) == %.3f * phantom("."), italic(p) == %.3f * phantom("...."))',
      correlation, 
      p_value
    ),
    hjust = 0, 
    vjust = 1, 
    parse = TRUE
  )
  
  # Return results
  return(list(
    plot = p,
    correlation = correlation,
    p_value = p_value,
    lm_model = lm_model
  ))
}


                                   
do_VolcanoPlot <- function(sample,
                           de_genes,
                           genes.up = NULL, genes.down = NULL,
                           pval_cutoff = 0.05,
                           FC_cutoff = 2,
                           pt.size = 2,
                           border.size = 1.5,
                           border.color = "black",
                           font.size = 14,
                           font.type = "sans",
                           plot.title = NULL,
                           plot.subtitle = NULL,
                           plot.caption = NULL,
                           plot_lines = TRUE,
                           line_color = "grey75",
                           line_size = 0.5,
                           add_gene_tags = TRUE,
                           add_tag_side = "both",
                           order_tags_by = "both",
                           n_genes = 5,
                           use_labels = FALSE,
                           colors.use = "steelblue",
                           min.segment.length=0.14, 
                           nudge_x=0.5, 
                           force=0.1,
                           plot.title.face = "plain",
                           plot.subtitle.face = "plain",
                           plot.caption.face = "italic",
                           axis.title.face = "plain",
                           axis.text.face = "plain",
                           legend.title.face = "plain",
                           legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))


  # Check logical parameters.
  logical_list <- list("add_gene_tags" = add_gene_tags,
                       "plot_lines" = plot_lines,
                       "use_labels" = use_labels)
  # Check numeric parameters.
  numeric_list <- list("pval_cutoff" = pval_cutoff,
                       "FC_cutoff" = FC_cutoff,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "line_size" = line_size,
                       "n_genes" = n_genes)
  # Check character parameters.
  character_list <- list("border.color" = border.color,
                         "font.type" = font.type,
                         "line_color" = line_color,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "add_tag_side" = add_tag_side,
                         "order_tags_by" = order_tags_by,
                         "colors.use" = colors.use,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
 

  `%>%` <- magrittr::`%>%`
  colors <- do_ColorPalette(colors.use, tetradic = TRUE)
  names(colors) <- c("A", "C", "B", "D")

  if (!("gene" %in% colnames(de_genes))){
    data <- de_genes %>%
            tibble::rownames_to_column(var = "gene")
  } else {
    data <- de_genes
  }


  data <- data %>%
          tibble::as_tibble() %>%
          dplyr::select(c("p_val_adj", "avg_log2FC", "gene")) %>%
          dplyr::mutate("p_val_adj" = replace(.data$p_val_adj, .data$p_val_adj == 0, .Machine$double.xmin)) %>%
          dplyr::mutate(log_p = -log10(.data$p_val_adj)) %>%
          dplyr::select(-"p_val_adj")

  pval_cutoff <- -log10(pval_cutoff)
  data$color <- NA
  data$color[data$avg_log2FC > FC_cutoff & data$log_p > pval_cutoff] <- "A"
  data$color[data$avg_log2FC < -FC_cutoff & data$log_p > pval_cutoff] <- "B"
  data$color[(data$log_p < pval_cutoff) | (abs(data$avg_log2FC) < FC_cutoff)] <- "C"

  max_value <- max(abs(c(min(data$avg_log2FC), max(data$avg_log2FC)))) + 0.6
  x_lims <- c(-max_value, max_value)

  # Shuffle the data.
  data <- data[sample(rownames(data), nrow(data)), ]
  p <- data %>%
       ggplot2::ggplot(mapping = ggplot2::aes(x = .data$avg_log2FC,
                                              y = .data$log_p)) +
       ggplot2::geom_point(size = pt.size * border.size, alpha = 0.5,
                           color = border.color) +
       ggplot2::geom_point(mapping = ggplot2::aes(color = .data$color),  alpha = 0.9,
                           size = pt.size) +
       ggplot2::labs(title = plot.title,
                     subtitle = plot.subtitle,
                     caption = plot.caption) +
       ggplot2::scale_color_manual(values = c(pal_nejm()(5)[1], pal_nejm()(5)[2],'grey')) +
       ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4),
                                                     title.position = "top",
                                                     title.hjust = 0.5)) +
       ggplot2::xlim(x_lims) +
       ggplot2::xlab(expression(paste(log["2"], "(fold change)"))) +
       ggplot2::ylab(expression(paste("-", log["10"], "(P value adjusted)")))

  if (isTRUE(plot_lines)){
    p <- p +
         ggplot2::geom_hline(yintercept = pval_cutoff,
                             color = line_color,
                             linewidth = line_size,
                             linetype = "dashed") +
         ggplot2::geom_vline(xintercept = FC_cutoff,
                             color = line_color,
                             linewidth = line_size,
                             linetype = "dashed") +
         ggplot2::geom_vline(xintercept = -FC_cutoff,
                             color = line_color,
                             linewidth = line_size,
                             linetype = "dashed")
  }

  if (isTRUE(add_gene_tags)){
    if (order_tags_by == "custom"){
       data.up <- data[data$gene %in% genes.up,] %>%  as.data.frame()
       data.down <- data[data$gene %in% genes.down,] %>%  as.data.frame()
    } else if (order_tags_by == "both"){
      data.up <- data %>%
                 dplyr::arrange(dplyr::desc(.data$log_p),
                                dplyr::desc(.data$avg_log2FC)) %>%
                 as.data.frame() %>%
                 utils::head(n_genes)

      data.down <- data %>%
                   dplyr::arrange(dplyr::desc(.data$log_p),
                                  .data$avg_log2FC) %>%
                   as.data.frame() %>%
                   utils::head(n_genes)
    } else if (order_tags_by == "pvalue"){
      data.up <- data %>%
                 dplyr::filter(.data$avg_log2FC > 0) %>%
                 dplyr::arrange(dplyr::desc(.data$log_p),
                                dplyr::desc(.data$avg_log2FC)) %>%
                 as.data.frame() %>%
                 utils::head(n_genes)

      data.down <- data %>%
                   dplyr::filter(.data$avg_log2FC < 0) %>%
                   dplyr::arrange(dplyr::desc(.data$log_p)) %>%
                   as.data.frame() %>%
                   utils::head(n_genes)
    } else if (order_tags_by == "logfc"){
      data.up <- data %>%
                 dplyr::filter(.data$log_p > pval_cutoff) %>%
                 dplyr::arrange(dplyr::desc(.data$avg_log2FC)) %>%
                 as.data.frame() %>%
                 utils::head(n_genes)

      data.down <- data %>%
                   dplyr::filter(.data$log_p > pval_cutoff) %>%
                   dplyr::arrange(.data$avg_log2FC) %>%
                   as.data.frame() %>%
                   utils::head(n_genes)

      data.else <- data %>%
                   dplyr::filter(.data$avg_log2FC < FC_cutoff) %>%
                   dplyr::filter(.data$avg_log2FC > -FC_cutoff) %>%
                   as.data.frame()
      data.else$gene <- ' '
    }

    if (add_tag_side == "both") {
      data.label <- dplyr::bind_rows(data.up, data.down)

    } else if (add_tag_side == "positive") {
      data.label <- data.up

    } else if (add_tag_side == "negative") {
      data.label <- data.down

    }

    if (base::isFALSE(use_labels)){
      p <- p +
           ggrepel::geom_text_repel(data = data.up,
                                    mapping = ggplot2::aes(label = .data$gene),
                                    max.overlaps = 1000, min.segment.length=min.segment.length, nudge_x=nudge_x, force=force,
                                    color = "black",
                                    fontface = "italic") +
         ggrepel::geom_text_repel(data = data.down,
                                    mapping = ggplot2::aes(label = .data$gene),
                                    max.overlaps = 1000, min.segment.length=min.segment.length, nudge_x=-nudge_x, force=force,
                                    color = "black",
                                    fontface = "italic")
    } else if (isTRUE(use_labels)){
      p <- p +
           ggrepel::geom_label_repel(data = data.label,
                                     mapping = ggplot2::aes(label = .data$gene),
                                     max.overlaps = 1000,
                                     color = "black",
                                     fontface = "italic")
    }

  }
  p <- p +
       ggplot2::theme_classic(base_size = font.size) +
       ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                      plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                      plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                      plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                      legend.text = ggplot2::element_text(face = legend.text.face),
                      legend.title = ggplot2::element_text(face = legend.title.face),
                      panel.grid = ggplot2::element_blank(),
                      legend.position = "none",
                      legend.justification = "center",
                      axis.title.x = ggplot2::element_text(face = axis.title.face, color = "black"),
                      axis.title.y = ggplot2::element_text(face = axis.title.face, angle = 90, color = "black"),
                      axis.text = ggplot2::element_text(face = axis.text.face, color = "black"),
                      # axis.line = ggplot2::element_line(color = "black"),
                      # axis.ticks = ggplot2::element_line(color = "black"),
                      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                      panel.border = element_rect(fill=NA, linewidth=1),
                      axis.line = element_line(colour = 'black', linewidth = 0),
                      legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  return(p)
}



#' Filter Genes Based on Ensembl Biotype
#'
#' This function retrieves gene biotype information from Ensembl and filters
#' genes based on specified criteria.
#'
#' @param genes Character vector of gene symbols to query.
#' @param biotype Character vector of accepted biotypes. Default is "protein_coding".
#' @param ensembl_dataset Character string specifying the Ensembl dataset.
#'                       Default is "hsapiens_gene_ensembl".
#' @param mart_name Character string specifying the Mart name. Default is "ensembl".
#' @param attributes Character vector of Ensembl attributes to retrieve.
#'                  Default is c("hgnc_symbol", "gene_biotype").
#'
#' @return A list containing:
#'        \itemize{
#'          \item keep: Logical vector indicating which genes pass the filter
#'          \item filtered_genes: Character vector of genes that pass the filter
#'          \item gene_info: Data frame with complete gene information from Ensembl
#'        }
#'
#' @details The function queries Ensembl using biomaRt to get gene biotype
#' information and filters genes based on specified biotypes. It maintains
#' the original gene order for the output.
#'
#' @examples
#' # Filter for protein-coding genes
#' results <- filter_genes_by_biotype(rownames(tcells))
#'
#' # Filter for multiple biotypes
#' results <- filter_genes_by_biotype(
#'   rownames(tcells),
#'   biotype = c("protein_coding", "lincRNA")
#' )
#'
#' @import biomaRt
#'
#' @export
filter_genes_by_biotype <- function(genes,
                                  biotype = "protein_coding",
                                  ensembl_dataset = "hsapiens_gene_ensembl",
                                  mart_name = "ensembl",
                                  attributes = c("hgnc_symbol", "gene_biotype")) {
  
  # Input validation
  if (!is.character(genes)) {
    stop("genes must be a character vector")
  }
  
  if (!is.character(biotype)) {
    stop("biotype must be a character vector")
  }
  
  # Remove any NA or empty values from genes
  genes <- genes[!is.na(genes) & genes != ""]
  
  # Connect to Ensembl
  tryCatch({
    ensembl <- useMart(mart_name, dataset = ensembl_dataset)
  }, error = function(e) {
    stop("Failed to connect to Ensembl: ", e$message)
  })
  
  # Query Ensembl
  gene_info <- tryCatch({
    getBM(attributes = attributes,
          filters = "hgnc_symbol",
          values = genes,
          mart = ensembl)
  }, error = function(e) {
    stop("Failed to retrieve gene information from Ensembl: ", e$message)
  })
  
  # Match to original gene order
  gene_info <- gene_info[match(genes, gene_info$hgnc_symbol), ]
  rownames(gene_info) <- NULL  # Clean up row names
  
  # Create filter
  keep <- (!is.na(gene_info$hgnc_symbol) & 
            gene_info$gene_biotype %in% biotype)
  
  # Get filtered genes
  filtered_genes <- genes[keep]
  
  # Return results
  return(list(
    keep = keep,
    filtered_genes = filtered_genes,
    gene_info = gene_info
  ))
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}





# source: https://github.com/naomihabiblab/BEYOND_DLPFC/blob/main/1.%20Library%20preprocessing/utils/run.doublet.finder.R
#' Search for doublets using DoubletFinder
#'
#' Run DoubletFinder on given object. Store DoubletFinder's `proportion of Artificial Nearest Neighbors` as a meta.data
#' column with specified name. If requested, the simulated doublet parent names as well as parent identity distribution
#' is stored under \code{Tool(object = object, "RunDoubletFinder")}
#'
#' @param object Seurat object to run DoubletFinder over
#' @param pN proportion of simulated doublets from object after adding these doublets. 
#' @param pK proportion of dataset size to be used as number of neighbors. If NULL DoubletFinder's \code{paramSweep_v3}
#' is executed. If `k` is the desired number of neighbors then pass `k/(ncol(object)/(1-.25))`
#' @param pcs A vector PCs to use for the run of DoubletFinder. If NULL uses the `dims` logged under the `FindNeighbors`
#' command of the specified `assay` and `reduction`
#'
RunDoubletFinder <- function(object, assay = DefaultAssay(object = object), pN = .25, pK = NULL, pcs = NULL, sct = F, 
                             score.name="doublet.score", reduction="PCA", batch.size=5000, sim.idents = NULL,
                             compute.parent.ident.distribution=T, ...) {
  require(DoubletFinder); require(data.table)
  
  if(is.null(pcs)) {
    pcs <- Command(object, paste("FindNeighbors", assay, reduction, sep = "."))$dims
  }
  if (is.null(pK)) {
    sweep.res.list <- DoubletFinder::paramSweep_v3(object, PCs = pcs, sct = sct)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <-DoubletFinder::find.pK(sweep.stats)
    pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), 
    ]$pK))
  }
  
  df <- doubletFinder_v3(object, PCs = pcs, pN = pN, pK = pK, nExp = 0, idents = sim.idents,
                         reuse.pANN = F, sct = sct, batch.size = batch.size, 
                         get.neighbor.doublets = compute.parent.ident.distribution)
  object@meta.data[,score.name] <- df@meta.data[, grep("pANN", colnames(df@meta.data))]
  
  if(compute.parent.ident.distribution) {
    message("Computing Doublet Neighbors Parents Ident Distribution...")
    tool <- df@tools$doubletFinder_v3
    
    tool$neighbor.doublets <- tool$neighbor.doublets[sapply(tool$neighbor.doublets, function(l) length(l) != 0)]
    n <- names(tool$neighbor.doublets)
    tool$neighbor.doublets <- as.data.frame(rbindlist(lapply(tool$neighbor.doublets, function(l) data.table(t(l))), fill=T))
    rownames(tool$neighbor.doublets) <- n
    
    tool$parent.ident.distribution <- ComputeDoubletNeighborParentsDistribution(object, tool, score.column = score.name, ...)
    Tool(object) <- tool
  }
  
  object <- Seurat::LogSeuratCommand(object)
  return(object)
}




                                                             #' Perform Differential Expression Analysis Using Muscat
#'
#' This function conducts differential expression analysis on a Seurat object
#' using the muscat package methodology. It processes single-cell RNA sequencing
#' data and identifies differentially expressed genes between groups within clusters.
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA sequencing data
#' @param cluster_col Character string specifying the column name in meta.data containing cluster information
#' @param group_col Character string specifying the column name in meta.data containing group information
#' @param sample_col Character string specifying the column name in meta.data containing sample information
#' @param method Character string specifying the DE method to use (e.g., "DESeq2", "edgeR", "limma")
#'
#' @return A data frame containing differential expression results
#'
#' @examples
#' \dontrun{
#' results <- muscat.de.genes(
#'   seurat_obj = my_seurat_object,
#'   cluster_col = "seurat_clusters",
#'   group_col = "condition",
#'   sample_col = "sample_id",
#'   method = "DESeq2"
#' )}
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom SingleCellExperiment as.SingleCellExperiment
#' @importFrom muscat prepSCE aggregateData
#'
muscat.de.genes <- function(seurat_obj, cluster_col, group_col, sample_col, method, min_cells = 10, filter = 'both') {
    # Input validation
    if (!inherits(seurat_obj, "Seurat")) {
        stop("Input must be a Seurat object")
    }
    
    required_cols <- c(cluster_col, group_col, sample_col)
    missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns in meta.data: %s", 
                    paste(missing_cols, collapse = ", ")))
    }
    
    # Set default assay to RNA
    DefaultAssay(seurat_obj) <- "RNA"
    
    # Convert to SingleCellExperiment
    tryCatch({
        sce <- as.SingleCellExperiment(seurat_obj)
        sce <- prepSCE(sce, 
                      kid = cluster_col, 
                      gid = group_col, 
                      sid = sample_col, 
                      drop = TRUE)
    }, error = function(e) {
        stop("Error in converting to SingleCellExperiment: ", e$message)
    })
    
    # Aggregate data
    pb <- aggregateData(sce, 
                       assay = "counts", 
                       fun = "sum", 
                       by = c("cluster_id", "sample_id"))
    
    return(muscat_pbDS(pb, sce, method, min_cells, filter))
}

#' Process Pseudo-bulk Differential Expression Results
#'
#' This helper function processes the results from pseudo-bulk differential
#' expression analysis, filtering for significant genes and calculating
#' summary statistics.
#'
#' @param pb Pseudo-bulk data object from aggregateData
#' @param sce SingleCellExperiment object
#' @param method Character string specifying the DE method
#'
#' @return A data frame containing filtered and processed differential expression results
#'
#' @importFrom dplyr arrange filter bind_rows
#' @importFrom muscat pbDS
#'
muscat_pbDS <- function(pb, sce, method, min_cells = 10, filter = 'both') {
    # Perform differential expression analysis
    res <- pbDS(pb, 
                method = method,
                filter = filter,
                min_cells = min_cells,
                verbose = TRUE)
    
    # Process significant results
    significant <- lapply(res$table[[1]], function(u) {
        arrange(filter(u, p_adj.loc <= 0.05), p_adj.loc)
    })
    
    # Calculate summary statistics
    n_de <- vapply(significant, nrow, numeric(1))
    p_de <- format(n_de / nrow(sce) * 100, digits = 3)
    
    # Print summary
    summary_df <- data.frame(
        "Num DE.genes" = n_de,
        "% DE.genes" = p_de,
        check.names = FALSE
    )
    print(summary_df)
    
    # Return processed results without column 10 (typically metadata)
    return(bind_rows(res$table[[1]])[-c(10)])
}



no.labs <- labs(x=NULL, y=NULL, color=NULL, fill=NULL)
no.axes <- theme(axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_blank())
theme_embedding <- theme_classic() + no.axes