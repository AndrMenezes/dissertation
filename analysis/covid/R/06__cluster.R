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

## For clusters methods
suppressMessages(library(bluster))

## Compute Gap statistics
library(cluster)

## Heatmaps
library(pheatmap)

# Colors to plot
col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]

# Read the merged sce -----------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
sce

mat <- reducedDim(sce, "PCA_corrected")
dim(mat)


# Compute graph cluster ---------------------------------------------------

# system.time(
#  clust_walktrap <- clusterRows(mat[, 1:50], NNGraphParam(cluster.fun = "walktrap"))
# )
# saveRDS(clust_walktrap, "./data/clust_walktrap")
# clust_louvain <- clusterRows(mat[, 1:50], NNGraphParam(cluster.fun = "louvain"))

# Compute k-means clustering for k=1...25 ---------------------------------

k_seq <- 1:25
set.seed(69)
time_start <- proc.time()
km_list <- lapply(k_seq, function(k) {
  cat(k, "\n")
  kmeans(mat, iter.max = 50, centers = k)
})
time_end <- proc.time()
time <- time_end - time_start
time[3] / 60
# saveRDS(km_res, file="./data/06_cluster_kmeans.rds")


## Plot WSS
size <- sapply(km_list, "[[", "size")
names(size) <- k_seq
wss <- sapply(km_list, "[[", "tot.withinss")
best <- which.min(wss)
df_wss <- tibble::tibble(wss=wss, k=k_seq)
ggplot(df_wss, aes(x=k, y=wss)) +
  geom_line(linetype = "dashed") + 
  geom_point(shape = 21, fill = "grey", col = "black", size = 3, alpha = 0.6) +
  scale_x_continuous(breaks = seq(1, 25, by = 2), limits = c(1, 25)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(y = "WSS", x = "Número de grupos (k)")
ff <- file.path(wd, "wss_kmeans.pdf")
ggsave(ff, device = "pdf", width = 8, height = 6)




# Incluí cluster e salva sce ----------------------------------------------

sce$kmeans_6 <- factor(km_list[[6]]$cluster)
sce$kmeans_9 <- factor(km_list[[9]]$cluster)
sce$kmeans_11 <- factor(km_list[[11]]$cluster)
sce$kmeans_13 <- factor(km_list[[13]]$cluster)
sce$kmeans_15 <- factor(km_list[[15]]$cluster)
sce$kmeans_17 <- factor(km_list[[17]]$cluster)
sce$kmeans_20 <-  factor(km_list[[20]]$cluster)

table(sce$kmeans_6)
table(sce$kmeans_11)
table(sce$kmeans_13)
table(sce$kmeans_15)
table(sce$kmeans_20)

quickResaveHDF5SummarizedExperiment(sce)

# Evaluating cluster stability --------------------------------------------

set.seed(66)
ari <- bootstrapStability(mat, clusters = sce$kmeans_13, 
                          mode = "index", 
                          BLUSPARAM = KmeansParam(13, iter.max = 50),
                          iterations = 50)
ari
# 0.7484193

set.seed(66)
ratio <- bootstrapStability(mat, clusters = sce$kmeans_13,
                            mode = "ratio", 
                            iterations = 50,
                            BLUSPARAM = KmeansParam(13, iter.max = 50))


seq_hcl <- colorspace::sequential_hcl(100, palette = "BluYl")

mat_numbers <- ratio
mat_numbers <- FF(mat_numbers, 3)
mat_numbers[lower.tri(mat_numbers)] <- ""

p_heat <- pheatmap(ratio, cluster_cols = FALSE, cluster_rows = FALSE,
         #display_numbers = mat_numbers, number_color = "black", 
         #number_format = "%.3f", fontsize_number = 14,
         color=seq_hcl, breaks = seq(-1, 1, length=101))
pfull <- ggdraw() +
  draw_plot(gridExtra::arrangeGrob(p_heat[[4]]), x = 0, y = 0)

ff <- file.path(wd, "heatmap_ari.pdf")
pdf(ff, width = 15, height = 8)
print(pfull)
dev.off()

set.seed(66)
prob <- bootstrapStability(mat, clusters = sce$kmeans_13, 
                           iterations = 50,
                           BLUSPARAM = KmeansParam(13, iter.max = 50),
                           compare = scran::coassignProb)

mat_numbers <- prob
mat_numbers <- FF(mat_numbers, 3)
mat_numbers[lower.tri(mat_numbers)] <- ""

p_heat <- pheatmap(prob, cluster_cols = FALSE, cluster_rows = FALSE,
                   display_numbers = mat_numbers, number_color = "black", 
                   number_format = "%.3f", fontsize_number = 14,
                   color = rev(seq_hcl))
pfull <- ggdraw() +
  draw_plot(gridExtra::arrangeGrob(p_heat[[4]]), x = 0, y = 0)

ff <- file.path(wd, "heatmap_coassign.pdf")
pdf(ff, width = 15, height = 8)
print(pfull)
dev.off()

# Computing the silhouette width ------------------------------------------

tb_sil <- approxSilhouette(mat, sce$kmeans_13) %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(cl_re = factor(paste0("Grupo ", sce$kmeans_13)),
                cl_re = forcats::fct_relevel(cl_re, paste0("Grupo ", 1:13)))
tb_sil %>% 
  dplyr::group_by(cl_re) %>% 
  dplyr::summarise(silhueta_media = mean(width))

ggplot(tb_sil, aes(x = factor(cl_re), y = width )) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red") +
  labs(x = "Agrupamentos", y = "Silhueta de cada observação")

ff <- file.path(wd, "silhueta.pdf")
ggsave(ff, device = "pdf", width = 8, height = 6)

# Relationship between the clusters obtained ------------------------------

cent_tree <- hclust(dist(km_list[[13]]$centers), "ward.D2")
cent_tree <- as.dendrogram(cent_tree)
plot(cent_tree, type = "rectangle")


# Root-mean-squared deviation ---------------------------------------------

## represents the spread of cells within each cluster. A cluster is more 
## likely to have a low RMSD if it has no internal structure and is 
## separated from other clusters

ncells <- tabulate(kmeans_clust$cluster)
tab <- data.frame(wcss=kmeans_clust$withinss, ncells=ncells)
tab$rms <- sqrt(tab$wcss / tab$ncells)
tab


# Total de células por indivíduo em cada cluster --------------------------

tb_n_cluster <- colData(sce) %>% 
  dplyr::as_tibble() %>% 
  dplyr::select(sample, kmeans_13) %>% 
  dplyr::rename(cluster = kmeans_13) %>% 
  dplyr::group_by(sample, cluster) %>% 
  dplyr::count() %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::mutate(pct = n/sum(n)) %>% 
  dplyr::ungroup() %>%  
  dplyr::mutate(
    cl_re = gsub("kmeans_", "",cluster),
    cl_re = forcats::fct_relevel(cl_re, function(x) as.character(sort(as.integer(x))))
  )

tb_n_cluster %>% 
  ggplot(aes(x=cl_re, y=100*pct, fill=factor(sample))) +
  geom_col(col="white") +
  scale_fill_manual(values = c(col_H, col_M, col_S)) +
  cowplot::theme_half_open(12)+
  labs(x = "Grupo", y="% de células", fill="") +
  scale_y_continuous(
    limits=c(0,100),
    expand=expansion(mult = c(0, 0.05))
  ) +
  theme(text = element_text(size = 14), 
        axis.text = element_text(size = 14))

ff <- file.path(wd, "pct_cell_cluster.pdf")
ggsave(filename = ff, device = "pdf", width = 12, height = 8)


# Tabela cruzada grupo versus individuo -----------------------------------
tab_freq <- table(sce$kmeans_13, sce$group) %>% 
  as.data.frame() %>% 
  dplyr::rename(grupo = Var1, individuo = Var2) %>% 
  dplyr::group_by(grupo) %>% 
  dplyr::mutate(prop = 100 * Freq / sum(Freq),
                prop = prettyNum(round(prop, 2), decimal.mark = ",", nsmall = 2),
                Freq = prettyNum(Freq, big.mark = ".", decimal.mark = ","),
                lab = paste0(Freq, " (", prop, "\\%)")) %>% 
  dplyr::arrange(grupo) %>%
  dplyr::select(-c(Freq, prop)) %>% 
  tidyr::pivot_wider(names_from = c(individuo), values_from = lab)
library(xtable)
print.xtable(xtable(tab_freq), include.rownames = FALSE, 
             sanitize.text.function = force)



# UMAP conforme os clusters -----------------------------------------------

tb_DimRed <- reducedDim(sce, "UMAP_corrected")[, 1:2] %>% 
  dplyr::as_tibble() %>% 
  dplyr::rename(x=V1, y=V2) %>% 
  dplyr::mutate(
    cl_re = forcats::fct_relevel(sce$kmeans_13, function(x) as.character(sort(as.integer(x))))
  )

p_dim <- ggplot(tb_DimRed, aes(x=x, y=y, col=cl_re)) +
  geom_point(size = 2) +
  labs(x = "UMAP 1", y = "UMAP 2", col = "") +
  # scale_color_manual(values = clusterExperiment::bigPalette) +
  theme_cowplot() +
  background_grid() +
  theme(legend.position = "top")

x11(); p_dim

sce <- schex::make_hexbin(sce, nbins = 100,
                          dimension_reduction = "UMAP_corrected", 
                          use_dims = c(1,2))

p_hex <- schex::plot_hexbin_meta(sce, 
                                 col = "kmeans_13",
                                 action = "majority",
                                 colors = clusterExperiment::bigPalette)
label_df <- schex::make_hexbin_label(sce, col = "kmeans_13")
p_hex <- p_hex +
  ggtitle("") +
  labs(x= "UMAP 1", y = "UMAP 2", fill = "") +
  theme_cowplot() +
  background_grid() +
  theme(legend.position = "top") +
  ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label), 
                            colour="black",  label.size = NA, fill = NA, size = 5)
x11(); p_hex

ff <- file.path(wd, "phex_cluster.pdf")
ggsave(filename = ff, plot = p_hex, device = "pdf", width = 9, height = 7)


