suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
    library(scCustomize) 
    library(SCP)
    library(ggplot2)
    library(dplyr)
    library(knitr)
    library(readr) 
    library(ggsci)
    library(scater)
    library(DoubletFinder)
    library(Trex)    
    library(SCpubr)
    library(biomaRt)
    library(data.table)
    library(genekitr)
    library(Azimuth)
    library(UCell)
    library(harmony)

    library(SummarizedExperiment)
    library(limma)
    library(muscat)
    library(purrr)
    library(scater)
    library(ggpubr)
    library(fgsea)
    
    
})


options(future.globals.maxSize = 128*1024**3)
plan(strategy = "multicore", workers = 1)
plan()

set.seed(123)


panel.path = 'Figures/'



"%!in%" <- Negate("%in%")

load("../../data/cycle.rda") # tirosh cycling genes doi/10.1126/science.aad0501


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



# Adapted fro SCpubr https://enblacar.github.io/SCpubr-book/                                 
do_VolcanoPlot <- function(sample, de_genes, genes.up = NULL, genes.down = NULL, pval_cutoff = 0.05,
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
                           legend.text.face = "plain") {
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

  max_value <- max(abs(c(min(data$avg_log2FC), max(data$avg_log2FC)))) +0.6
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


preprocess_seurat_object <- function(obj, g2m_genes, s_genes,  nfeatures = 4000, phase_threshold=1) {

  obj <- NormalizeData(obj) 
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  
  objv5 <- obj
  objv5[["RNA"]] <- as(object = objv5[["RNA"]], Class = "Assay")
  objv5 <- quietTCRgenes(objv5, assay = 'RNA')

  per.batch.var.feats <- VariableFeatures(objv5)
  objv5 <- FindVariableFeatures(objv5, nfeatures = nfeatures)
  objv5 <- quietTCRgenes(objv5, assay = 'RNA')
  joint.var.feats <- VariableFeatures(objv5)
  var.feats <- union(joint.var.feats, per.batch.var.feats)
  
  objv5 <- CellCycleScoring(objv5, g2m.features = g2m_genes, s.features = s_genes, verbose = FALSE)
  
  obj.sce <- as.SingleCellExperiment(objv5)
  
  diff <- getVarianceExplained(obj.sce, "Phase")
  discard <- diff > phase_threshold
  
  var_minus_cycling <- var.feats[var.feats %!in% rownames(objv5)[which(discard)]]
  
  rm(objv5)

  VariableFeatures(obj) <- var_minus_cycling
  
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  
  return(obj)
}

cell_type_levels = c(
  "CD4+ CAR T",
  "CD4+ T",
  "CD8+ CAR T",
  "CD8+ T",
  "Cycling CD8+ T",
    
  "Mixed CAR+/CAR- Tfh",
  "Treg",
  "NK",
  "MAIT",
    
  "Monocyte",  
  "Macrophage",
  "cDC1",
  "cDC2",
  "mregDC",
  "pDC",
    
  "B",
  "Plasma"
)


cell_type_cols <- c(
  "CD4+ T"               = "#A6CEE3",
  "Macrophage"           = "#4992C2",
  "CD4+ CAR T"           = "#569EA4",
  "Plasma"               = "#AADB84",
  "Treg"                 = "#52AF43",
  "cDC2"                 = "#8BAD49",
  "cDC1"                 = "#FDB761",
  "Mixed CAR+/CAR- Tfh"  = "#FE8B14",
  "mregDC"               = "#FD8C4C",
  "Monocyte"             = "#F68181",
  "B"                    = "#E5292B",
  "pDC"                  = "#70449D",
  "Cycling CD8+ T"       = "#C7B699",
  "CD8+ CAR T"           = "#D46F84",
  "CD8+ T"               = "#B294C7",
  "NK"                   = "#E6CB75",
  "MAIT"                 = "#B15928"
)


