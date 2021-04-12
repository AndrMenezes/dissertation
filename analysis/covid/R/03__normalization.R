rm(list = ls())

setwd("./covid/")
wd <- "../latex/figuras/04_capitulo4/02_data_analysis"

# Pkgs --------------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## Normalization by deconvolution
suppressMessages(library(scran))
## Apply logNormcounts in sce object
suppressMessages(library(scuttle))
# Read the data
suppressMessages(library(HDF5Array))

# Colors to plot
col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]

# Read the datas ----------------------------------------------------------

samples <- c(paste0("HC", 1:3), paste0("M", 1:3), paste0("S", 1:6))
dir_to_read <- paste0("./data/", samples)
l_sce <- lapply(dir_to_read, function(d) loadHDF5SummarizedExperiment(dir = d))

# Normalization -----------------------------------------------------------

set.seed(66)
l_sce <- lapply(seq_along(l_sce), function(z) {
  sce <- l_sce[[z]]
  
  ini <- proc.time()
  clust <- quickCluster(sce)
  sce <- computeSumFactors(sce, cluster = clust, min.mean = 0.1)
  sce <- logNormCounts(sce)
  fim <- proc.time() - ini
  
  cat("Sample", sce$sample[1], "finished. Time elapsed", fim[3]/60, "minutes.\n")
  return(sce)
})




# librarySizeFactor versus sizeFactor -------------------------------------

tb_lib_size <- lapply(seq_along(l_sce), function(z){
  tb <- tibble::tibble(
    sf_deconvolution = sizeFactors(l_sce[[z]]),
    sf_library = librarySizeFactors(l_sce[[z]]),
    sample = l_sce[[z]]$sample[1]
  )
  return(tb)
})
tb_lib_size <- do.call("rbind", tb_lib_size)
ggplot(tb_lib_size, aes(x = log(sf_library), y = log(sf_deconvolution))) +
  facet_wrap(~sample) +
  geom_point(fill = "grey", alpha = 0.5, shape = 21, size = 1.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Fator de escala tamanho da biblioteca", y = "Fator de escala deconvolution") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6))
ff <- file.path(wd, "size_factor.png")
ggsave(filename = ff, device = "png", width = 14, height = 8, dpi = 200)

tb_lib_size %>% 
  dplyr::filter(sample == "S6") %>% 
  ggplot(aes(x = log(sf_library), y = log(sf_deconvolution))) +
  geom_point(fill = "grey", alpha = 0.5, shape = 21, size = 2.0) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Fator de escala tamanho da biblioteca", y = "Fator de escala deconvolution") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6))
ff <- file.path(wd, "one_size_factor.png")
ggsave(filename = ff, device = "png", width = 8, height = 6, dpi = 200)


# Distribution of one normalized cell -------------------------------------
cells_sample <- lapply(l_sce, function(z) colnames(z))
set.seed(69)
cells_chosen <- sapply(cells_sample, function(z) sample(z, size = 1))

sce_cells <- do.call("cbind", lapply(l_sce, function(z) z[, colnames(z) %in% cells_chosen]))
tb_cells <- logcounts(sce_cells) %>% 
  as.matrix() %>% 
  dplyr::as_tibble() %>% 
  tidyr::pivot_longer(cols = dplyr::everything(), 
                      names_to = "Barcode",
                      values_to = "logNormCounts") %>% 
  dplyr::left_join(
    dplyr::as_tibble(colData(sce_cells)[, c("Barcode", "sample", "group")]),
    by = "Barcode"
  ) %>% 
  dplyr::mutate(
    group = forcats::fct_relevel(group, "HC", "M", "S"),
  )


ggplot(tb_cells, aes(x=logNormCounts, fill=factor(group))) +
  facet_wrap(~Barcode, scales = "free") +
  geom_density(alpha=0.6) +
  labs(x = "Logaritmo das expressão gênica normalizada",
       y = "Densidade",
       fill = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  theme(legend.position = "top")
ff <- file.path(wd, "log_expr_norm_cells.pdf")
ggsave(filename = ff, device = "pdf", width = 14, height = 8)


# Distribution of normalized genes ----------------------------------------
genes <- c("HLA-DRA", "S100A6", "FTL", "C1QA")
sce_chosen <- do.call("cbind", lapply(l_sce, function(x) x[genes, ]))

tb <- logcounts(sce_chosen) %>% 
  as.matrix() %>% 
  t() %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(sample = sce_chosen$sample,
                group = sce_chosen$group)

tb_long <- tb %>% 
  tidyr::pivot_longer(cols = -c(sample, group), 
                      names_to = "gene", 
                      values_to = "logNormCounts")

tb_long %>% 
  dplyr::filter(gene != "APOE") %>% 
  ggplot(aes(x = logNormCounts, fill = factor(group))) +
  facet_grid(gene~sample, scales = "free_y") +
  geom_density(alpha = 0.9) +
  labs(x = expression(log[2]~"Contagem normalizada"),
       y = "Densidade",
       fill = "") +
  scale_fill_manual(values = c(col_H[3], col_M[3], col_S[6])) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  theme(legend.position = "top")


ff <- file.path(wd, "log_expr_norm_genes.pdf")
ggsave(filename = ff, device = "pdf", width = 14, height = 8)


# Saving ------------------------------------------------------------------
for (i in seq_along(l_sce)) {
  quickResaveHDF5SummarizedExperiment(l_sce[[i]])
  cat(i, "\n")
}

