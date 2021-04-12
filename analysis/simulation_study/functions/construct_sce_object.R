library(SingleCellExperiment)
library(scater)
library(batchelor)

#' Constrói um objeto SCE incluindo dados corrigidos, PCA, UMAP e TSNE.

construct_sce <- function(simul_obj, get_redDim = TRUE, get_cluster = FALSE,
                          verbose = TRUE, vst = FALSE, corrected = TRUE) {
   
  colnames(simul_obj$counts) <- paste0("Cell_", 1:ncol(simul_obj$counts))
  rownames(simul_obj$counts) <- paste0("Gene_", 1:nrow(simul_obj$counts))
  sce <- SingleCellExperiment(assays = list(counts = simul_obj$counts))
  colData(sce)$group <- simul_obj$design$group
  colData(sce)$individuo <- simul_obj$design$individuo
  colData(sce)$cell_type <- simul_obj$design$cell_type
  colData(sce)$comb <- simul_obj$design$comb
  rowData(sce)$is_de <- simul_obj$is_de

  # Filtra células com contagem nula
  zero_cell_count <- colSums(assay(sce, "counts")) == 0
  if (sum(zero_cell_count) > 0) {
    sce <- sce[, !zero_cell_count]
  }
  # Filtra genes com contagem nula
  zero_gene_count <- rowSums(assay(sce, "counts")) == 0
  if (sum(zero_gene_count) > 0) {
    sce <- sce[!zero_gene_count, ]
  }
  # Normalizaçao usando scran
  dummy_clust <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, cluster = dummy_clust)
  sce <- scuttle::logNormCounts(sce)
  
  # Correção nos valores da contagem do efeito de individuo
  if (corrected) {
    sce_corrected <- rescaleBatches(sce, batch = sce$individuo)
    assay(sce, "corrected") <- assay(sce_corrected, "corrected")
  }
  
  # Normalização via variance stabilizing transformation for UMI 
  if (vst) {
    y <- sctransform::vst(umi = counts(sce), min_cells = 0, verbosity = 0)$y
    assay(sce, "vstresiduals") <- y
  }

  # Gene e células com contagem zero pós correção
  zero_gene_corrected <- rowSums(assay(sce, "corrected")) == 0
  zero_cell_corrected <- colSums(assay(sce, "corrected")) == 0
  
  if (get_redDim) {
    sce <- runPCA(sce)
    sce <- runUMAP(sce)
    # Redução da Dimensionalidade nos valores corrigidos
    sce <- runPCA(sce, exprs_values = "corrected", name = "PCA_corrected")
    sce <- runUMAP(sce, dimred = "PCA_corrected", name = "UMAP_corrected", 
                   exprs_values = "corrected")
    # Redução da Dimensionalidade nos valores do resíduo
    sce <- runPCA(sce, exprs_values = "vstresiduals", name = "PCA_vst")
    sce <- runUMAP(sce, dimred = "PCA_vst", name = "UMAP_vst", 
                   exprs_values = "vst")
    
  }
  
  if (get_cluster) {
    clust_cell_type <- kmeans(reducedDim(sce, "PCA"), iter.max = 50, 
                              centers = length(unique(sce$cell_type)))
    sce$clust_cell_type <- factor(clust_cell_type$cluster)
  }

  if (verbose) {
    cat("Total de célula", ncol(sce), "\n")
    cat("Total de células com contagem nula", sum(zero_cell_count), "\n")
    cat("Total de genes com contagem nula", sum(zero_gene_count), "\n")
    cat("Total de células com contagem nula pós correção", 
        sum(zero_cell_corrected), "\n")
    cat("Total de genes com contagem nula pós correção", 
        sum(zero_gene_corrected), "\n")
  }
  
  return(sce)
}
