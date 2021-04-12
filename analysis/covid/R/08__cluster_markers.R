rm(list = ls())

setwd("./covid/")
wd <- "../latex/figuras/04_capitulo4/02_data_analysis"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")

# Packages ----------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## Read the data
suppressMessages(library(HDF5Array))

## For some convenience functions
suppressMessages(library(SingleCellExperiment))

## For marker gene detection
suppressMessages(library(scran))

## For heatmap
library(pheatmap)

# Read the merged sce -----------------------------------------------------
sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
sce


# Get markers using Welch t-test ------------------------------------------
markers_welch <- findMarkers(
  sce,
  groups = factor(sce$kmeans_13),
  block = factor(sce$sample),
  test.type = "t", 
  pval.type = "any",
  direction = "up",
  assay.type = "logcounts_corrected"
)
head(markers_welch[[1]], 6)[,1:3]
head(markers_welch[[2]], 6)[,1:3]


# Get markers using Wilcoxon test -----------------------------------------
markers_wilcoxon <- findMarkers(
  sce,
  groups = factor(sce$kmeans_13),
  block = factor(sce$sample),
  test.type = "wilcox", 
  pval.type = "any",
  direction = "up",
  assay.type = "logcounts_corrected"
)
head(markers_wilcoxon[[1]], 6)[,1:3]
head(markers_wilcoxon[[2]], 6)[,1:3]


# Get markers using linear model ------------------------------------------
design <- model.matrix(~sce$sample)
design <- design[,-1,drop=FALSE]
head(design)

markers_lm <- findMarkers(
  sce,
  groups = factor(sce$kmeans_13),
  design = design,
  pval.type = "any",
  direction = "up",
  assay.type = "logcounts_corrected"
)

head(markers_lm[[1]], 6)[,1:4]
head(markers_lm[[2]], 6)[,1:5]
head(markers_lm[[3]], 6)[,1:5]


# Comparing methods -------------------------------------------------------
out_tot <- list()
out_prop <- list()
for (g in seq_along(markers_lm)) {
  
  markers_g_t <- markers_welch[[g]]
  markers_g_w <- markers_wilcoxon[[g]]
  markers_g_lm <- markers_lm[[g]]
  
  ## Total number of markers genes detected for group g at a FDR of 5%
  sig_t <- nrow(markers_g_t[markers_g_t$FDR <= 0.05, ])
  sig_w <- nrow(markers_g_w[markers_g_w$FDR <= 0.05, ])
  sig_lm <- nrow(markers_g_lm[markers_g_lm$FDR <= 0.05, ])
  
  df_tot <- data.frame(group = g, t = sig_t, w = sig_w, lm = sig_lm) 
  out_tot[[g]] <- df_tot
  
  ## Proportion of the top-ranking markers genes shared between methods.
  ## Top genes were defined as those with the smallest 20, 50, or 100 p-values 
  ## in each comparison
  
  markers_g_t$rank <- rank(markers_g_t$FDR, ties.method = "first")
  markers_g_w$rank <- rank(markers_g_w$FDR, ties.method = "first")
  markers_g_lm$rank <- rank(markers_g_lm$FDR, ties.method = "first")
  # 
  # markers_g_t <- markers_g_t[order(rownames(markers_g_t)), ]
  # markers_g_w <- markers_g_w[order(rownames(markers_g_w)), ]
  # markers_g_lm <- markers_g_lm[order(rownames(markers_g_lm)), ]
  # 
  # df_rank <- data.frame(gene = rownames(markers_g_t), 
  #                       rank_t = markers_g_t$rank, 
  #                       rank_w = markers_g_w$rank,
  #                       rank_lm = markers_g_lm$rank)
  # df_rank$min_rank <- ifelse(df_rank$rank_t < df_rank$rank_w, df_rank$rank_t, 
  #                            df_rank$rank_w)
  # df_rank$dif_rank <- df_rank$rank_t - df_rank$rank_w
  # plot(df_rank$min_rank, df_rank$dif_rank)
  
  df_prop <- data.frame()
  for (top in c(20, 50, 100)) {
    top_t <- rownames(markers_g_t[markers_g_t$rank <= top, ])
    top_w <- rownames(markers_g_w[markers_g_w$rank <= top, ])
    top_lm <- rownames(markers_g_lm[markers_g_lm$rank <= top, ])
    
    r <- ifelse(top == 20, paste0("\\multirow{3}{*}{", g, "}"), "")
    aux <- data.frame(group = g, group_lab = r, top = top,
                      t_w = length(intersect(top_t, top_w)) / top,
                      t_lm = length(intersect(top_t, top_lm)) / top,
                      w_lm = length(intersect(top_w, top_lm)) / top)
    df_prop <- rbind(df_prop, aux)
  }
  out_prop[[g]] <- df_prop
}
out_tot <- do.call("rbind", out_tot)
out_prop <- do.call("rbind", out_prop)

## Gráficos
tb_long <- out_tot %>% 
  tidyr::pivot_longer(cols = -group) %>% 
  dplyr::mutate(
    label = dplyr::case_when(
      name == "t" ~ "Welch",
      name == "w" ~ "Wilcoxon",
      name == "lm" ~ "Modelo linear"), 
    group = factor(group), label = factor(label)
  )
ggplot(tb_long, aes(x = group, y = value, fill = label)) +
  geom_col(position = position_dodge(), alpha = 0.8, col = "black") +
  colorspace::scale_fill_discrete_qualitative() +
  scale_y_continuous(limits = c(0,1800), expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::pretty_breaks(8)) +
  labs(x = "Grupo", y = "Total de genes marcadores", fill = "") +
  theme(legend.position = "top")

ggsave(filename = file.path(wd, "total_markers.pdf"), device = "pdf", 
       width = 10, height = 8)

tb_long <- out_prop %>% 
  dplyr::select(-group_lab) %>% 
  tidyr::pivot_longer(cols = -c(group, top)) %>% 
  dplyr::mutate(
    label = dplyr::case_when(
      name == "t_w" ~ "Welch e Wilcoxon",
      name == "t_lm" ~ "Welch e Modelo linear",
      name == "w_lm" ~ "Wilcoxon e Modelo linear"), 
    group = factor(group), 
    top = forcats::fct_relevel(factor(paste0("Top ", top)), "Top 20", "Top 50")
  )
ggplot(tb_long, aes(x = group, y = value, fill = label)) +
  facet_grid(~top) +
  geom_col(position = position_dodge(), alpha = 0.8, col = "black") +
  coord_flip() +
  colorspace::scale_fill_discrete_qualitative() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::pretty_breaks(8)) +
  labs(x = "Grupo", y = "Proporção de genes marcadores compartilhados entre os métodos", 
       fill = "") +
  theme(legend.position = "top")
ggsave(filename = file.path(wd, "prop_markers.pdf"), device = "pdf", 
       width = 12, height = 8)


## Tabelas
out_tot[, -1] <- apply(out_tot[, -1], 2, prettyNum, decimal.mark = ",", 
                       big.mark = '.')
out_prop[, 3:5] <- apply(out_prop[, 3:5], 2, prettyNum, decimal.mark = ",",
                         nsmall = 2) 

ff <- file.path(wd, "total_markers.tex")
file.remove(ff)
sink(ff, append = T)
cat('\\begin{table}[H]', "\n")
cat("\\caption{Total de genes marcadores considerando FDR de 5\\% em cada agrupamento conforme método.}", "\n")
cat("\\onehalfspacing \\centering \\begin{tabular}{cccc} \\toprule 
    Grupo & Welch & Wilcoxon & Modelo linear \\\\ \\midrule  ")
invisible(apply(out_tot, 1, printmrow))
cat('\\bottomrule \\end{tabular}', "\n")
cat('\\end{table}', "\n")
sink()

ff <- file.path(wd, "prop_markers.tex")
file.remove(ff)
sink(ff, append = T)
cat('\\begin{table}[H]', "\n")
cat("\\caption{Proporção dos top genes marcadores compartilhados entre métodos. Os top genes foram definidos como aqueles com os menores 20, 50 ou 100 valores p das comparações múltiplas.}", "\n")
cat("\\onehalfspacing \\centering \\begin{tabular}{cccc} \\toprule 
    Grupo & Top & Welch e Wilcoxon & Welch e Modelo linear & Wilcoxon e Modelo linear \\\\ \\midrule  ")
invisible(apply(out_prop, 1, printmrow))
cat('\\bottomrule \\end{tabular}', "\n")
cat('\\end{table}', "\n")
sink()




# Tables with the gene marker and logFC for all groups --------------------

options(OutDec = ",")
groups <- seq_along(markers_lm)
for (g in groups) {
  tab_g <- markers_lm[[g]]
  tab_g <- tab_g[rank(tab_g$FDR, ties.method = "first") <= 8, -c(1:4)]
  logFCs <- getMarkerEffects(tab_g)
  
  p_heat <- pheatmap(logFCs, cluster_cols = TRUE, cluster_rows = TRUE,
                     display_numbers = TRUE, number_color = "black", 
                     number_format = "%.3f", fontsize_number = 14, 
                     main = paste0("Grupo ", g),
                     angle_col = 0, 
                     fontsize = 20,
                     fontsize_row = 16,
                     fontsize_col = 16,
                     breaks = seq(-2, 2, length=101))
  pfull <- ggdraw() +
    draw_plot(gridExtra::arrangeGrob(p_heat[[4]]), x = 0, y = 0)
  
  pdf(file.path(wd, paste0("heatmap_logFC_group_", g, ".pdf")), width = 12, 
      height = 8)
  print(pfull)
  dev.off()
  
  tab_g <- cbind(gene = rownames(tab_g), tab_g)
  tab_g[, -1] <- apply(tab_g[, -1], 2, FF)
  tab_g <- t(as.matrix(tab_g))
  col <- tab_g[1, ]
  col <- c("Comparação", col)
  col <- c(paste0(col[-6], sep = " &"), col[6])
  tab_g <- tab_g[-1, ]
  
  row <- paste0(g, " vs. ", groups[-g])
  tab_g <- cbind(row, tab_g)
  
  ff <- file.path(wd, paste0("markers_group_", g,".tex"))
  file.remove(ff)
  sink(ff, append = T)
  cat('\\begin{table}[H]', "\n")
  cat("\\caption{Log fold change dos top oito genes ranqueados pelo valor-p das comparações múltiplas entre grupos --- 
      Agrupamento", g,".}", "\n")
  cat("\\onehalfspacing \n \\scalefont{0.9} \n \\centering \n \\begin{tabular}{crrrrrrrr} \\toprule \n")
  cat(col, "\\\\ \\midrule \n")
  invisible(apply(tab_g, 1, printmrow))
  cat('\\bottomrule \\end{tabular}', "\n")
  cat('\\end{table}', "\n")
  sink()
  cat(g, "\n")  
}
options(OutDec = ".")




# Heatmap with the mean expression level of markers genes -----------------
plot_markers <- function(markers, top = 5, method) {
  
  my_breaks <- seq(-2, 2, length=101)
  prefix <- "logFC"
  if (method == "wilcoxon") {
    my_breaks <- seq(0, 1, length.out = 101)
    prefix = "AUC"
  } 
  
  
  for (g in seq_along(markers)) {

    tab_g <- markers[[g]][1:top, ]
    # tab_g <- tab_g[tab_g$Top <= top, ]
    logFCs <- getMarkerEffects(tab_g, prefix = prefix)

    # cat("Total de marcadores no grupo", g, nrow(tab_g), "\n")
    
    p_heat <- pheatmap(logFCs, cluster_cols = TRUE, cluster_rows = TRUE,
                       display_numbers = TRUE, number_color = "black", 
                       number_format = "%.3f", fontsize_number = 14, 
                       main = paste0("Grupo ", g),
                       angle_col = 0, 
                       fontsize = 20,
                       fontsize_row = 16,
                       fontsize_col = 16,
                       breaks = my_breaks)
    pfull <- ggdraw() +
      draw_plot(gridExtra::arrangeGrob(p_heat[[4]]), x = 0, y = 0)
    
    pdf(file.path(wd, paste0("heatmap_effect_group_", g, "_", method, ".pdf")), 
        width = 12, height = 8)
    print(pfull)
    dev.off()
  }

  # Heatmap with the mean expression level of markers genes
  
  top_markers <- lapply(markers, function(z) rownames(z[1:top, ]))
  gene_ids <- unique(unlist(top_markers))

  sce_pseudo <- scater::sumCountsAcrossCells(
    x = sce[gene_ids,],
    ids = DataFrame(label=sce$kmeans_13),
    average = TRUE,
    exprs_values = "logcounts_corrected"
  )
  
  mat <- assay(sce_pseudo)
  colnames(mat) <- sce_pseudo$label
  
  p_heat <- pheatmap::pheatmap(
    mat = mat,
    scale = "row",
    cluster_rows = TRUE,
    angle_col = 0, 
    fontsize = 20,
    fontsize_row = 16,
    fontsize_col = 16,
    color = rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )
  
  pfull <- cowplot::ggdraw() +
    cowplot::draw_plot(gridExtra::arrangeGrob(p_heat[[4]]), x = 0, y = 0)
  
  pdf(paste0(wd, "/heatmap_mean_markers_", method, ".pdf"), width = 16, 
      height = 10)
  print(pfull)
  dev.off()

    
} 
plot_markers(markers = markers_lm, top = 5, method = "lm")
plot_markers(markers = markers_welch, top = 5, method = "welch")
plot_markers(markers = markers_wilcoxon, top = 5, method = "wilcoxon")
