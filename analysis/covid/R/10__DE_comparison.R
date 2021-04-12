rm(list = ls())

setwd("./covid/")
wd <- "D:/dissertation/latex/figuras/04_capitulo4/02_data_analysis"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,
                                            format="f"))}


# Packages ----------------------------------------------------------------

## Read the data
suppressMessages(library(HDF5Array))

## For utils
suppressMessages(library(SingleCellExperiment))

## For aggregate
suppressMessages(library(scuttle))

## For latex tables
library(xtable)

## For %>% 
library(magrittr)

## For plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## For complex heatmap
library(ComplexHeatmap)

## For save the heatmap
library(GetoptLong)

## Colors to plot
col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]



# Data --------------------------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
table(sce$label.main, sce$group)
sce <- sce[, sce$label.main != "Neutrophils"]
# Total de amostras por indivíduo e cluster
table(sce$sample, sce$label.main)
# Exclui indivíduo HC2 no cluster Pre-B_cell_CD34-
# remove <- !(sce$label.main == "Pre-B_cell_CD34-" & sce$sample == "HC2")
# sce <- sce[, remove]

# Results of DE analysis --------------------------------------------------

df_mm <- readRDS("./data/DE_mlm_logcounts.rds")
cols <- Reduce(intersect, lapply(df_mm, names))
df_mm <- do.call("rbind", lapply(df_mm, function(x) x[, cols] )) %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::mutate(fdr = p.adjust(p.value, method = "BH")) %>% 
  dplyr::ungroup()

df_edgeR <- readRDS("./data/DE_results_edgeR.rds") %>%
  dplyr::filter(!is.na(PValue)) %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::mutate(fdr = p.adjust(PValue, "BH"))

head(df_mm)
head(df_edgeR)

df_mm %>% 
  dplyr::filter(gene == "FABP4")
df_edgeR %>% 
  dplyr::filter(gene == "FABP4")



# Tabela com gene mais significativo --------------------------------------

options(OutDec = ",")
ntop <- 1

## Modelo linear misto
top_genes_mm <- df_mm %>% 
  dplyr::group_by(label.main, contrast) %>% 
  dplyr::slice_min(order_by = p.value, n = ntop)

tab_mm <- top_genes_mm[, c(9, 2, 1, 3, 6, 10)]
print(xtable(tab_mm, digits = c(rep(1, 4), 4, -2, -2)), include.rownames=FALSE)

## Modelo pseudo-bulk edgeR
top_genes_e <- df_edgeR %>% 
  dplyr::group_by(label.main, contrast) %>% 
  dplyr::slice_min(order_by = PValue, n = ntop)

## Tabela
tab_e <- top_genes_e[, c(7, 1, 6, 2, 5, 8)]
print(xtable(tab_e, digits = c(rep(1, 4), 4, -2, -2)), include.rownames=FALSE)

options(OutDec = ".")


# Volcano plot ------------------------------------------------------------

## Modelo linear misto
tb_volcano <- df_mm %>% 
  dplyr::group_by(contrast, label.main) %>% 
  dplyr::mutate(
    label = ifelse(fdr <= quantile(fdr, probs = 0.01), gene, ""),
    contrast  = gsub(" - ", " versus ", contrast),
    label = ifelse(label.main == "Epithelial_cells", "", label)
  ) %>% 
  dplyr::ungroup()
tb_volcano %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::count(label != "")

ggplot(tb_volcano, aes(x = estimate, y = -log10(fdr), label = label,
                       col = label.main)) +
  facet_wrap(~contrast, scales = "free") +
  geom_point(alpha = 0.4, size = 3) +
  # geom_text(aes(y = -log10(fdr) + 0.3),  col = "black") +
  ggrepel::geom_text_repel(col = "black", max.time = 60, max.overlaps = 100) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = bquote(~Log[2]~ 'fold change'), y = bquote("-" ~Log[10]~ "FDR"), 
       col = "") +
  theme(legend.position = "top")
ggsave(file.path(wd, "DE_volcano_mixed_model.pdf"), device = "pdf",
       width = 14, height = 8)

## Modelo edgeR

tb_volcano <- df_edgeR %>%
  dplyr::group_by(contrast, label.main) %>% 
  dplyr::mutate(
    label = ifelse(fdr <= quantile(fdr, probs = 0.01), gene, ""),
    contrast  = gsub(" - ", " versus ", contrast),
    label = ifelse(label.main %in% c("Epithelial_cells", "Macrophage"), "", 
                   label)
  )
tb_volcano %>% 
  dplyr::group_by(label.main, contrast) %>% 
  dplyr::count(m = label != "") %>% 
  dplyr::filter(m)

ggplot(tb_volcano, aes(x = logFC, y = -log10(fdr), label = label,
                       col = label.main)) +
  facet_wrap(~contrast, scales = "free_x") +
  geom_point(alpha = 0.4, size = 3) +
  # geom_text(aes(y = -log10(fdr) + 0.1, label = label),  col = "black") +
  ggrepel::geom_text_repel(col = "black", max.time = 60) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = bquote(~Log[2]~ 'fold change'), y = bquote("-" ~Log[10]~ "FDR"), 
       col = "") +
  theme(legend.position = "top") +
  panel_border()
ggsave(file.path(wd, "DE_volcano_pseudo_bulk.pdf"), device = "pdf",
       width = 14, height = 8)



# Comparando a significancia de cada método -------------------------------

counter <- 1L
out_tot <- out_prop <- genes_mm <- genes_edgeR <- list()

for (label in sort(unique(df_mm$label.main))) {
  
  for (con in sort(unique(df_mm$contrast))) {
    
    chosen_mm <- (df_mm$contrast == con) & (df_mm$label.main == label)
    chosen_edgeR <- (df_edgeR$contrast == con) & (df_edgeR$label.main == label)
    cur_mm <- df_mm[chosen_mm, ]
    cur_edgeR <- df_edgeR[chosen_edgeR, ]
    
    cur_mm$rank <- rank(cur_mm$fdr, ties.method = "first")
    cur_edgeR$rank <- rank(cur_edgeR$fdr, ties.method = "first")
    cur_mm <- cur_mm[order(cur_mm$rank), ]
    cur_edgeR <- cur_edgeR[order(cur_edgeR$rank), ]
      
    ## Total number of markers genes detected for group g at a p-value of 5%
    sig_mm <- cur_mm[cur_mm$fdr <= 0.05, ]
    sig_edgeR <- cur_edgeR[cur_edgeR$fdr <= 0.05, ]
    out_tot[[counter]] <- data.frame(label = label, contrast = con, 
                                     MLM = nrow(sig_mm), edgeR = nrow(sig_edgeR),
                                     down_mm = sum(sig_mm$estimate < 0),
                                     down_edgeR = sum(sig_edgeR$logFC < 0),
                                     up_edgeR = sum(sig_edgeR$logFC > 0),
                                     up_mm = sum(sig_mm$estimate > 0)) 
    # Seleciona os 5 genes com menor fdr
    genes_mm[[counter]] <- data.frame(cell_type = label, contrast = con,
                                      gene = cur_mm$gene, 
                                      pvalue = cur_mm$p.value, 
                                      rank = rank(cur_mm$p.value,
                                                  ties.method = "first"))
    genes_edgeR[[counter]] <- data.frame(cell_type = label, contrast = con,
                                         gene = cur_edgeR$gene,
                                         pvalue = cur_edgeR$PValue, 
                                         rank = rank(cur_edgeR$PValue,
                                                     ties.method = "first"))
    ## Proportion of the top-ranking markers genes shared between methods.
    ## Top genes were defined as those with the smallest 20, 50, or 100 p-values 
    ## in each comparison
    df <- data.frame()
    for (top in c(20, 50, 100)) {
      top_mm <- cur_mm[cur_mm$rank <= top, ]$gene
      top_edgeR <- cur_edgeR[cur_edgeR$rank <= top, ]$gene
      aux <- data.frame(label = label, contrast = con, 
                        top = top,
                        shared = length(intersect(top_mm, top_edgeR)) / top)
      df <- rbind(df, aux)
    }
    out_prop[[counter]] <- df
    
    counter <- counter + 1L
  }
}
df_tot <- do.call("rbind", out_tot)
df_prop <- do.call("rbind", out_prop)
genes_mm <- do.call("rbind", genes_mm)
genes_edgeR <- do.call("rbind", genes_edgeR)

head(genes_edgeR)
head(genes_mm)



# Ranqueando os genes conforme o valor-p ----------------------------------

min_rank_mm <- genes_mm %>% 
  dplyr::group_by(cell_type, gene) %>% 
  dplyr::summarise(min_rank = min(rank), .groups = 'drop') %>% 
  dplyr::arrange(cell_type, min_rank)
min_rank_edgeR <- genes_edgeR %>% 
  dplyr::group_by(cell_type, gene) %>% 
  dplyr::summarise(min_rank = min(rank), .groups = 'drop') %>% 
  dplyr::arrange(cell_type, min_rank)
  

# Gráfico comparando o total de genes detectado com fdr <= 5% -------------

tb_long <- df_tot[, 1:4] %>% 
  tidyr::pivot_longer(cols = -c(label, contrast))
tb_long %>% 
  dplyr::filter(label != "Epithelial_cells") %>% 
  ggplot(aes(x = contrast, y = value, fill = name)) +
  facet_wrap(~label) +
  geom_col(position = position_dodge(), alpha = 0.8, col = "black") +
  colorspace::scale_fill_discrete_qualitative() +
  scale_y_continuous(limits = c(0,max(df_tot$MLM)+5), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Comparação", y = "Total de genes significativos ", fill = "") +
  theme(legend.position = "top")

ggsave(filename = file.path(wd, "total_genes_DE.pdf"), device = "pdf", 
       width = 10, height = 8)


# Gráfico com total de genes compatilhados --------------------------------

ggplot(df_prop, aes(x = contrast, y = shared, fill = factor(top))) +
  facet_grid(~label) +
  geom_col(position = position_dodge(), alpha = 0.8, col = "black") +
  coord_flip() +
  colorspace::scale_fill_discrete_qualitative() +
  scale_y_continuous(limits = c(0, 0.4), expand = expansion(mult = c(0, 0.08)),
                     breaks = scales::pretty_breaks(6)) +
  labs(x = "Grupo", y = "Proporção de genes DE compartilhados entre os métodos", 
       fill = "") +
  theme(legend.position = "top") +
  panel_border()
ggsave(filename = file.path(wd, "prop_genes_DE.pdf"), device = "pdf", 
       width = 12, height = 8)


# Visualização no nível dos indivíduos ------------------------------------

visul_sample_level <- function(sce, genes_res, n_top, save = TRUE, method) {
  
  ## Informação dos top n genes mais DE
  df_top_n <- genes_res %>% 
    dplyr::group_by(cell_type) %>% 
    dplyr::slice_head(n = n_top)
  
  ## Agregando as contagens em média para cada amostras dentro de cada grupo 
  ## de célula
  top_genes <- df_top_n$gene
  sce_agg <- aggregateAcrossCells(
    x = sce, 
    ids = colData(sce)[, c("sample", "label.main")],
    statistics = "mean",
    subset.row = top_genes,
    use.assay.type = "logcounts")
  colData(sce_agg) <- colData(sce_agg)[, c("sample", "group", "label.main",
                                           "ncells")]
  sce_agg$sample <- factor(sce_agg$sample)
  sce_agg$group <- factor(sce_agg$group)
  sce_agg$label.main <- factor(sce_agg$label.main)
  
  ## Reorganizando a matriz de expressão e os metadados
  mat <- assay(sce_agg)
  lt_mats <- lt_meta <- list()
  counter <- 1L
  for (k in sort(unique(df_top_n$cell_type))) {
    
    chosen_rows <- df_top_n[df_top_n$cell_type == k, ]$gene
    chosen_cols <- sce_agg$label.main == k
    lt_mats[[counter]] <- mat[chosen_rows, chosen_cols]
    lt_meta[[counter]] <- data.frame(sample = sce_agg$sample[chosen_cols],
                                     group = sce_agg$group[chosen_cols],
                                     cell_type = k) 
    counter <- counter + 1L
  }
  df_metadata <- do.call("rbind", lt_meta)
  df_metadata$cell_type <- factor(df_metadata$cell_type)
  df_metadata$group <- factor(df_metadata$group)
  mat_expr <- do.call("rbind", lt_mats)
  colnames(mat_expr) <- levels(df_metadata$sample)
  
  ## Inclui valor NA para o indivíduo HC2 no tipo de célula Pre-B_cell_CD34-,
  ## pois há somente uma observação para agregar nesse nível
  HC2_B_cell <- (df_metadata$cell_type == "Pre-B_cell_CD34-" &
                   df_metadata$sample == "HC2")
  mat_expr[19:24, "HC2"] <- NA_real_
   
  ## Organizando as informações para gerar o gráfico
  lgd_aes <- list(labels_gp = gpar(fontsize = 6),
                  title_gp = gpar(fontface = "bold", fontsize = 8))
  
  ## Anotação das linas
  cols_rows <- setNames(c("#882E72", "#FF7F00", "#E7298A", "#A6761D"),
                        levels(df_metadata$cell_type))
  row_split <- rep(levels(df_metadata$cell_type), each = 6)
  row_anno <- rowAnnotation(
    df = data.frame(`População` = row_split),
    col = list(`População` = cols_rows),
    gp = gpar(col = "white"),
    show_annotation_name = FALSE,
    annotation_legend_param = lgd_aes)
  
  ## Anotação das colunas
  cols_columns <- setNames(c(col_H[2], col_M[2], col_S[3]), 
                           levels(df_metadata$group))
  col_split <- rep(c("HC", "M", "S"), c(3, 3, 6))
  col_anno <- columnAnnotation(
    df = data.frame(Grupo = col_split),
    col = list(Grupo = cols_columns),
    gp = gpar(col = "white"),
    show_annotation_name = FALSE,
    annotation_legend_param = lgd_aes)
  
  ## Normalização das expressões médias
  mat_norm <- t(apply(mat_expr, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) 
  }))
  # mat_norm
  # rowMeans(mat_norm)
  # rowSds(mat_norm)
  
  ## Enfim plotando...
  hp <- Heatmap(
    matrix = mat_norm, 
    col = rev(RColorBrewer::brewer.pal(11, "RdBu")), 
    name = "Expressão\nnormalizada",
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    left_annotation = row_anno, 
    top_annotation = col_anno, 
    column_names_rot = TRUE,
    row_split = row_split, 
    column_split = col_split,
    heatmap_legend_param = lgd_aes,
    row_names_gp = gpar(fontsize = 6), 
    column_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontface = "bold", fontsize = 10),
    use_raster = TRUE, 
    raster_device = "CairoPNG")
  if (save) {
    pdf(file.path(wd, paste0("heatmap_expr_mean_", method,".pdf")), width = 10,
        height = 8)
    draw(hp)
    dev.off()
  }
  return(list(df_top_n = df_top_n, mat_expr = mat_expr))  
}

lt_tops_mm <- visul_sample_level(sce, genes_res = min_rank_mm, n_top = 6, 
                                 method = "MM", save = FALSE)

lt_tops_edgeR <- visul_sample_level(sce, genes_res = min_rank_edgeR, n_top = 6, 
                                    method = "edgeR", save = FALSE)



# Visualização no nível das células ---------------------------------------

visul_cell_level <- function(sce, df_top_n, save = TRUE, method) {
  
  ## Matriz com as expressões gênicas
  mat_expr_sc <- logcounts(sce)
  
  ## Organizando os dados no formato long para plotar
  lt_tbs <- list()
  counter <- 1L
  for (k in sort(unique(df_top_n$cell_type))) {
    chosen_cols <- sce$label.main == k
    chosen_rows <- df_top_n[df_top_n$cell_type == k, ]$gene
    expr_data <- mat_expr_sc[chosen_rows, chosen_cols] %>% 
      t() %>% 
      dplyr::as_tibble() %>% 
      dplyr::mutate(label = k)
    col_data <- colData(sce)[chosen_cols, c("sample", "group")] %>% 
      dplyr::as_tibble()
    lt_tbs[[counter]] <- dplyr::bind_cols(expr_data, col_data) %>% 
      tidyr::pivot_longer(cols = -c(label, sample, group), names_to = "gene",
                          values_to = "expr")
    counter <- counter + 1L
  }
  tb_long <- do.call("rbind", lt_tbs)
  
  ## Para cada tipo de célula plota a distribuição dos genes mais sig.
  for (k in sort(unique(df_top_n$cell_type))) {
    p_export <- tb_long %>% 
      dplyr::filter(label == k) %>% 
      ggplot(aes(x = sample, y = expr, col = group)) +
      facet_wrap(. ~ gene, scales="free_y") +
      ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.4, groupOnX = TRUE, 
                                   bandwidth = 1) +
      scale_color_manual(values = c(col_H[2], col_M[2], col_S[3])) +
      labs(x = "", y = expression(Log[2]~"expressão normalizada"), col = "") +
      panel_border() +
      theme_cowplot() +
      theme(legend.position = "top")
    ff <- paste0("expr_", k, "_", method, ".png")
    ggsave(filename = ff, path = wd, plot = p_export, width = 14, 
           height = 8)
  }
  return(tb_long)
}
tb_long_mm <- visul_cell_level(sce, df_top_n = lt_tops_mm$df_top_n, 
                               method = "MM")
tb_long_edgeR <- visul_cell_level(sce, df_top_n = lt_tops_edgeR$df_top_n, 
                                  method = "edgeR")


# Visulalização das amostras bulk -----------------------------------------

visul_bulk_level <- function(sce, df_top_n, method, save = TRUE) {
  
  
  for (k in sort(unique(df_top_n$cell_type))) {
    chosen_row <- df_top_n[df_top_n$cell_type == k, ]$gene
    chosen_col <- sce$label.main == k
  
    sce_pseudo <- aggregateAcrossCells(
      x = sce, ids = colData(sce)[, c("sample", "kmeans_13")],
      subset.col = chosen_col)
    colData(sce_pseudo) <- colData(sce_pseudo)[, c("sample", "group", 
                                                   "label.main", "ncells")]
    sce_pseudo <- logNormCounts(sce_pseudo)
    scater::plotExpression(sce_pseudo, features = chosen_row, x = "group")
    y <- DGEList(counts(sce_pseudo), samples = as.data.frame(cur_data))
    
  }
  
  colData(sce_pseudo)
  table(sce_pseudo$label.main, sce_pseudo$group)
  
  
  
}

# End ---------------------------------------------------------------------
sessionInfo()