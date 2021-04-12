rm(list = ls())

wd <- "../latex/figuras/04_capitulo4/02_data_analysis"
setwd("./covid/")

# Pkgs --------------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## Read the data
suppressMessages(library(HDF5Array))
## Feature selection
suppressMessages(library(scran)) 
## Batch correction
suppressMessages(library(batchelor))

# Colors to plot
col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]


# Read the merged sce -----------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
sce

# Compute PCA, tSNE and UMAP ----------------------------------------------

sce <- runPCA(
  sce, 
  ncomponents = 100,
  exprs_values = "logcounts",
  BSPARAM = BiocSingular::RandomParam()
)

set.seed(666)
sce <- runTSNE(sce, dimred = "PCA")
sce <- runUMAP(sce, dimred = "PCA")
reducedDimNames(sce) <- paste0(reducedDimNames(sce), "_uncorrected")

# Remove batch effect -----------------------------------------------------

sce_corrected <- rescaleBatches(sce, batch = sce$sample)
assay(sce, "logcounts_corrected") <- assay(sce_corrected, "corrected")

set.seed(666)
ini <- proc.time()
sce <- runPCA(
  x = sce, 
  exprs_values = "logcounts_corrected",
  ncomponents = 100,
  name = "PCA_corrected"
)
sce <- runTSNE(
  x = sce, 
  dimred = "PCA_corrected",
  exprs_values = "logcounts_corrected",
  name = "TSNE_corrected",
  external_neighbors = TRUE, 
  BSPARAM = BiocSingular::RandomParam()
)
sce <- runUMAP(
  x = sce, 
  dimred = "PCA_corrected",
  exprs_values = "logcounts_corrected",
  name = "UMAP_corrected"
)
fim <- proc.time() - ini
fim[3]/60



# Salvando o objeto -------------------------------------------------------

quickResaveHDF5SummarizedExperiment(sce)


# Visualizing -------------------------------------------------------------

x11(); plotReducedDim(sce, dimred = "PCA_uncorrected", colour_by = "sample")
x11(); plotReducedDim(sce, dimred = "PCA_corrected", colour_by = "sample")

x11(); plotReducedDim(sce, dimred = "TSNE_uncorrected", colour_by = "sample")
x11(); plotReducedDim(sce, dimred = "TSNE_corrected", colour_by = "sample")

p <- plotReducedDim(sce, dimred = "UMAP_uncorrected") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave(filename = "umap_vizul.png", path = wd,
       plot = p, device = "png", width = 6, height = 4)
x11(); plotReducedDim(sce, dimred = "UMAP_corrected", colour_by = "sample")

x11(); plotReducedDim(sce, dimred = "UMAP_uncorrected", colour_by = "group")
x11(); plotReducedDim(sce, dimred = "UMAP_corrected", colour_by = "group")


# Efeito da correção na expressão de alguns genes -------------------------

plot_gene_expression <- function(sce, gene_name, colour_by, save=TRUE) {
  
  id_gene <- which(rownames(sce) == gene_name)  
  
  tb <- tibble::tibble(
    counts = counts(sce)[id_gene, ],
    log_counts = logcounts(sce)[id_gene, ],
    log_corrected = assay(sce, "logcounts_corrected")[id_gene, ],
    colour_by = sce[[colour_by]]
  )
  
  if (length(unique(tb$colour_by)) == 3) {
    col_H <- col_H[3]
    col_M <- col_M[3]
    col_S <- col_S[5]
  }
  
  p1 <- ggplot(tb, aes(x = counts, fill = colour_by)) +
    geom_density(alpha = 0.6) +
    labs(x = "Contagem observada", y = "Densidade", fill = "") +
    scale_fill_manual(values = c(col_H, col_M, col_S))
  
  p2 <- ggplot(tb, aes(x = log_counts, fill = colour_by)) +
    geom_density(alpha = 0.6) +
    labs(x = expression(Log[2]~"Contagem normalizada"), y = "Densidade") +
    scale_fill_manual(values = c(col_H, col_M, col_S))
  
  p3 <- ggplot(tb, aes(x = log_corrected, fill = colour_by)) +
    geom_density(alpha = 0.6) +
    labs(x = expression(Log[2]~"Contagem normalizada e corrigida"), y = "Densidade") +
    scale_fill_manual(values = c(col_H, col_M, col_S))
  
  legend_p <- get_legend(p1)
  prow <- plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    align = 'hv',
    hjust = -1,
    nrow = 3
  )
  p_legended <- plot_grid(prow, legend_p, rel_widths = c(3, .4))
  my_title <- ggdraw() +
    draw_label(paste0("Expressão do gene ", gene_name), x = 0.05, hjust = 0)
  p_out <- plot_grid(my_title, p_legended, ncol = 1, rel_heights = c(0.1, 1))
  
  if (save) {
    wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/cap5"
    ff <- file.path(wd, paste0("expr_", colour_by, "_", gene_name, ".pdf"))
    ggsave(filename = ff, plot = p_out, device = "pdf", width = 11, height = 8)
  } 
  else {
    return(p_out)
  }
}

plot_gene_expression(sce, gene_name = "S100A6", colour_by = "group")
plot_gene_expression(sce, gene_name = "HLA-DRA", colour_by = "group")

plot_gene_expression(sce, gene_name = "S100A6", colour_by = "sample")
plot_gene_expression(sce, gene_name = "HLA-DRA", colour_by = "sample")
plot_gene_expression(sce, gene_name = "APOE", colour_by = "sample")
plot_gene_expression(sce, gene_name = "C1QA", colour_by = "sample")
plot_gene_expression(sce, gene_name = "FTL", colour_by = "sample")
plot_gene_expression(sce, gene_name = "LYZ", colour_by = "sample")



# Generating beautiful plots ----------------------------------------------

plotDim <- function(sce, dimred, colour_by, save = TRUE) {
  
  dimreds <- paste0(dimred, c("_uncorrected", "_corrected"))
  
  tb_DimRed <- reducedDim(sce, dimreds[1])[, 1:2] %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(colour_by = factor(sce[[colour_by]]),
                  fonte = "Não corrigido") %>% 
    dplyr::bind_rows(
      reducedDim(sce, dimreds[2])[, 1:2] %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(colour_by = factor(sce[[colour_by]]),
                      fonte = "Corrigido")
    ) %>% 
    dplyr::mutate(
      fonte = forcats::fct_relevel(fonte, "Não corrigido")
    )
  
  names(tb_DimRed) <- c("x", "y", "colour_by", "fonte")
  
  col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
  col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
  col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]
  
  if (length(unique(tb_DimRed$colour_by)) == 3) {
    col_H <- col_H[length(col_H)]
    col_M <- col_M[length(col_M)]
    col_S <- col_S[length(col_S)]
  }
  
  p_dim <- ggplot(tb_DimRed, aes(x=x, y=y, col=colour_by)) +
    facet_wrap(~fonte, scales = "free") + 
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = c(col_H, col_M, col_S)) +
    labs(x = paste0(dimred, " 1"), y = paste0(dimred, " 2"), col = "") +
    theme(legend.position = "top") 
  
  if (save) {
    ff <- paste0(dimred, ".png")
    ff <- file.path(wd, ff)
    ggsave(filename = ff, plot = p_dim, device = "png", width = 14, height = 10)
  } else {
    return(p_dim)
  }
}


## Plots
plotDim(sce, dimred = "PCA", colour_by = "sample")
plotDim(sce, dimred = "TSNE", colour_by = "sample")
plotDim(sce, dimred = "UMAP", colour_by = "sample")

plotDim(sce, dimred = "PCA", colour_by = "group")
plotDim(sce, dimred = "TSNE", colour_by = "group")
plotDim(sce, dimred = "UMAP", colour_by = "group")

