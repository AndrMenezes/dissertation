rm(list = ls())

wd <- "../latex/figuras/04_capitulo4/02_data_analysis"
setwd("./covid/")
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

# Packages ----------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## Read the data
suppressMessages(library(HDF5Array))

## For some convenience functions
suppressMessages(library(SingleCellExperiment))

## For cluster annotation 
suppressMessages(library(clustifyr))
suppressMessages(library(SingleR))

# Colors to plot
col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]

# Read the merged sce -----------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
sce


# Anotação ----------------------------------------------------------------

## Carregando matriz de referência
ref_data <- celldex::HumanPrimaryCellAtlasData()
ref_data
assay(ref_data)[1:10, 1:10]

## Identificando genes presentes nas duas bases
genes_sce <- rownames(sce)
genes_ref <- rownames(ref_data)
interseccao <- intersect(genes_sce, genes_ref)
apenas_sce  <- genes_sce[-which(genes_sce %in% interseccao)]
## Total de genes que serão utilizados
length(interseccao)

## Genes perdidos
length(apenas_sce)


## Função para anotar os clusters

get_annotation <- function(sce, ref_data, genes, cluster, method = "clustifyr") {
  
  # Total de células em cada cluster
  tb_ncells <- table(cluster) %>% 
    as.data.frame() %>% 
    dplyr::rename(ncells=Freq)
  
  # Tipo da célula para cada cluster
  df <- as.data.frame(colData(ref_data))
  tb_type <- dplyr::as_tibble(df) %>% 
    dplyr::mutate(type = rownames(df))
  
  # Utilizando clustifyr
  if (method == "clustifyr") {
    X    <- as.matrix(assay(sce, "logcounts_corrected"))
    umap <- reducedDim(sce, "UMAP_corrected")
    meta <- data.frame(cluster = cluster, umap = umap)
    ref_mat <- assay(ref_data)
    colnames(ref_mat) <- ref_data$label.fine
    
    pred_clustifyr <- clustify(
      input = X,
      metadata = meta,
      cluster_col = "cluster", 
      ref_mat = assay(ref_data), 
      query_genes = genes
    )
    
    tb_clustifyr <- cor_to_call(
      cor_mat = pred_clustifyr,
      cluster_col = "cluster" 
    ) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(tb_type, by = "type") %>% 
      dplyr::left_join(tb_ncells, by = "cluster")
    #%>% 
    #dplyr::rename_with(~paste0(.x, "_clustifyr"), -dplyr::starts_with("cluster"))
    
    return(tb_clustifyr)
    
  }
  
  
  ## Utilizando o SingleR
  if(method == "SingleR"){
    pred_SingleR <- SingleR(
      test = sce,
      ref = ref_data,
      method = "cluster",
      clusters = cluster,
      labels = colnames(ref_data),
      assay.type.test = "logcounts_corrected"
    )
    tb_SingleR <- dplyr::tibble(
      cluster = rownames(pred_SingleR),
      type = pred_SingleR$labels,
      r = pred_SingleR$tuning.scores$first
    ) %>% 
      dplyr::left_join(tb_type, by = "type") %>% 
      dplyr::rename_with(~paste0(.x, "_SingleR"), -dplyr::starts_with("cluster")) %>% 
      dplyr::arrange(cluster)
    
    return(tb_SingleR)
  }
  
}
ks <- paste0("kmeans_", c(6, 9, 11, 13, 15, 17, 20))
system.time(
  lt_annotation <- lapply(ks, function(z) {
    cat(z, "\n")
    get_annotation(sce, ref_data, genes = interseccao, cluster = sce[[z]])
  })
)
names(lt_annotation) <- c(6, 9, 11, 13, 15, 17, 20)
tb_annotation <- dplyr::bind_rows(lt_annotation, .id = "kmeans") %>% 
  dplyr::mutate(kmeans = forcats::fct_relevel(kmeans, "6", "9"))
saveRDS(tb_annotation, "./data/tb_annotation_clustifyr.rds")

# Anotação principal
tb_annotation_main <- tb_annotation %>% 
  dplyr::group_by(kmeans, label.main) %>% 
  dplyr::summarise(
    ncluster = dplyr::n(),
    ncells = sum(ncells), .groups = "drop"
  ) 
fill_breaks <- round(seq(min(tb_annotation_main$ncells), max(tb_annotation_main$ncells), l = 4))
fill_labels <- prettyNum(fill_breaks, big.mark = ".", decimal.mark = ",")
ggplot(tb_annotation_main, aes(x=kmeans, y=label.main, fill=ncells)) +
  geom_tile(col = "black", alpha = 0.85) +
  geom_text(aes(label = ncluster), size = 6) +
  labs(x = "Número de grupos (k)", y = "Anotação principal do grupo", fill = "Total de células") +
  colorspace::scale_fill_continuous_sequential(palette="Blues 2", breaks=fill_breaks, labels=fill_labels) +
  cowplot::theme_cowplot() + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1)) +
  theme(legend.position = "top")

ff <- file.path(wd, "anotacao_principal_kmeans.pdf")
ggsave(ff, device = "pdf", width = 9, height = 6)


# Anotação específica

tb_annotation_fine <- tb_annotation %>% 
  dplyr::group_by(kmeans, label.fine) %>% 
  dplyr::summarise(
    ncluster = dplyr::n(),
    ncells = sum(ncells), .groups = "drop"
  ) 
fill_breaks <- round(seq(min(tb_annotation_fine$ncells), max(tb_annotation_fine$ncells), l = 4))
fill_labels <- prettyNum(fill_breaks, big.mark = ".", decimal.mark = ",")
ggplot(tb_annotation_fine, aes(x=kmeans, y=label.fine, fill=ncells)) +
  geom_tile(col = "black", alpha = 0.85) +
  geom_text(aes(label = ncluster), size = 6) +
  labs(x = "Número de grupos (k)", y = "Anotação específica do grupo", fill = "Total de células") +
  colorspace::scale_fill_continuous_sequential(palette="Blues 2", breaks=fill_breaks, labels=fill_labels) +
  cowplot::theme_cowplot() + 
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1)) +
  theme(legend.position = "top")

ff <- file.path(wd, "anotacao_especifica_kmeans.pdf")
ggsave(ff, device = "pdf", width = 9, height = 6)


## Tabela com as anotações do grupo 13
options(OutDec = ",")
mat <- as.data.frame(colData(sce)[, c("kmeans_13", "r", "label.main", "label.fine")]) %>% 
  dplyr::as_tibble() %>% 
  dplyr::group_by(kmeans_13, r, label.main, label.fine) %>%
  dplyr::count(name = "ncells") %>%
  dplyr::mutate(ncells = prettyNum(ncells, big.mark = ".", decimal.mark = ','),
                r = FF(r)) %>% 
  dplyr::select(kmeans_13, ncells, r, label.main, label.fine) %>% 
  as.matrix()

printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")

ff <- file.path(wd, "tab_annotation_kmeans.tex")
file.remove(ff)
sink(ff, append = T)
cat('\\begin{table}[H]', "\n")
cat("\\caption{Anotação dos tipo de células predominantes para cada grupo.}", "\n")
cat("\\onehalfspacing \\centering \\begin{tabular}{ccccc} \\toprule 
    Grupo & Total de células &$\\rho$ & Principal & Específico \\\\ \\midrule  ")
invisible(apply(mat, 1, printmrow))
cat('\\bottomrule \\end{tabular}', "\n")
cat('\\end{table}', "\n")
sink()


# Incluí anotação do cluster e salva sce -----------------------------------
annotation_clust <- tb_annotation[tb_annotation$kmeans == 13,]
colData(sce) <- merge(colData(sce), annotation_clust[, -c(1, 3, 8)], 
                      by.x = "kmeans_13", by.y = "cluster")
quickResaveHDF5SummarizedExperiment(sce)         
