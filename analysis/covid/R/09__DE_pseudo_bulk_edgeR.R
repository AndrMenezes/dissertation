rm(list = ls())

setwd("./covid/")
wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/04_capitulo4/02_data_analysis"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,
                                            format="f"))}

# Packages ----------------------------------------------------------------

## Read the data
suppressMessages(library(HDF5Array))

## For utilities
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))

## For plot
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## For export tables
library(xtable)

## For fit the models
source("./R/edgeR_DE.R")


# Data --------------------------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
table(sce$label.main, sce$group)
sce <- sce[, sce$label.main != "Neutrophils"]

# Gera pseudo bulk amostras -----------------------------------------------

sce_pseudo <- aggregateAcrossCells(
  x = sce, ids = colData(sce)[, c("sample", "kmeans_13")])
colData(sce_pseudo) <- colData(sce_pseudo)[, c("sample", "group", "kmeans_13", 
                                               "label.main", "ncells")]
colData(sce_pseudo)
table(sce_pseudo$label.main, sce_pseudo$group)

tab_freq <- table(sce_pseudo$label.main, sce_pseudo$group) %>% 
  as.data.frame() %>% 
  dplyr::rename(grupo = Var1, individuo = Var2) %>% 
  dplyr::group_by(grupo) %>% 
  dplyr::mutate(prop = 100 * Freq / sum(Freq),
                prop = prettyNum(round(prop, 2), decimal.mark = ",", nsmall = 2),
                Freq = prettyNum(Freq, big.mark = ".", decimal.mark = ","),
                lab = paste0(Freq, " (", prop, "\\%)")) %>% 
  dplyr::arrange(grupo) %>%
  dplyr::select(-c(Freq, prop)) %>% 
  tidyr::pivot_wider(names_from = c(individuo), values_from = lab) %>% 
  dplyr::ungroup()
library(xtable)
print.xtable(xtable(tab_freq), include.rownames = FALSE, 
             sanitize.text.function = force)



# Realiza a análise DE utilizando a abordagem do pacote edgeR -------------
## Contrastes para realizar as comparações múltiplas entre HC, M e S
my_contrast <- matrix(c(1, -1, 0, 
                        1, 0, -1, 
                        0, 1, -1), ncol = 3)

res <- DE_edgeR(x = sce_pseudo, label = sce_pseudo$label.main, 
                design = ~0 + group, contrast = my_contrast)

dplyr::glimpse(res)



# Gráficos com as dispersões estimadas ------------------------------------

l_fit <- lapply(res, "[[", 2)
x <- l_fit[[1]]
plot(x$AveLogCPM, sqrt(x$trended.dispersion))

gplots <- lapply(l_fit, function(y) {
  tb <- dplyr::tibble(aveLogCPM = y$AveLogCPM,
                      disp_gene = y$tagwise.dispersion,
                      disp_trend = y$trended.dispersion)

  ggplot(tb, aes(x = aveLogCPM, y = sqrt(disp_gene))) +
    geom_point(aes(color = "Específica"), size = 1.0) +
    geom_line(aes(y = sqrt(disp_trend), color = "Tendência"), 
              size = 1.4) +
    geom_hline(aes(color = "Comum",
                   yintercept = sqrt(y$common.dispersion)), size = 1.4) +
    scale_color_manual(values = c("red", "black", "dodgerblue"), 
                       name = "Dispersão") +
    theme(legend.position = "top") +
    labs(x = "Média do log CPM", y = "Coeficiente de variação biológico") +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    ggtitle(y$samples$label.main[1])

})
grid_plots <- gridExtra::grid.arrange(gplots[[1]], gplots[[2]], gplots[[3]],
                                      gplots[[4]])
ggsave(filename = "edgeR_dispersion_estimate.png", plot = grid_plots,
       path = wd, width = 12, height = 8)
ggsave(filename = "edgeR_dispersion_one_estimate.png", plot = gplots[[3]] + ggtitle(""),
       path = wd, width = 12, height = 8)



# Resultados dos testes ---------------------------------------------------

df_res <- lapply(res, function(x) as.data.frame(x$tabs)) %>% 
  dplyr::bind_rows(.id = "label.main") %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(
    contrast  = gsub("1\\*", "", contrast),
    contrast  = gsub("-", "- ", contrast)
  )
head(df_res)

# saveRDS(df_res, file = "./data/DE_results_edgeR.rds")
df_res <- readRDS(file = "./data/DE_results_edgeR.rds")

# Volcano plot ------------------------------------------------------------

tb_volcano <- df_res %>%
  dplyr::filter(!is.na(PValue)) %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::mutate(
    fdr = p.adjust(PValue, "BH"),
    label = ifelse(fdr <= 0.001, gene, ""),
    contrast  = gsub("1\\*", "", contrast),
    contrast  = gsub("-", "versus ", contrast)
  )

ggplot(tb_volcano, aes(x = logFC, y = -log10(fdr), label = label,
                       col = label.main)) +
  facet_wrap(~contrast, scales = "free_x") +
  geom_point(alpha = 0.4, size = 3) +
  # geom_text(aes(y = -log10(fdr) + 0.1, label = label),  col = "black") +
  ggrepel::geom_text_repel(col = "black") +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = bquote(~Log[2]~ 'fold change'), y = bquote("-" ~Log[10]~ "FDR"), 
       col = "") +
  theme(legend.position = "top") +
  panel_border()
ggsave(file.path(wd, "DE_volcano_pseudo_bulk.pdf"), device = "pdf",
       width = 14, height = 8)


# Analisando expressão de alguns genes ------------------------------------

## Genes que foram mais significativos por tipo de célula e comparação

ntop <- 1
top_genes <- df_res %>% 
  dplyr::filter(!is.na(PValue)) %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::mutate(fdr = p.adjust(PValue, "BH"),
                contrast  = gsub("1\\*", "", contrast),
                contrast  = gsub("-", "- ", contrast)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(label.main, contrast) %>% 
  dplyr::slice_min(order_by = PValue, n = ntop)

## Tabela
tab <- top_genes[, c(7, 1, 6, 2, 5, 8)]

options(OutDec = ",")
print(xtable(tab, digits = c(rep(1, 4), 4, -2, -2)), include.rownames=FALSE)




## Gráfico da expressão gênica
lt <- list()
for (i in 1:nrow(top_genes)) {
  
  con <- as.character(top_genes$contrast[i])
  label <- as.character(top_genes$label.main[i])
  gene <- top_genes$gene[i]
  group <- strsplit(con, split = " - ")[[1]]
  chosen <- (sce$group %in% group) & (sce$label.main == label)
  sce_chosen <- sce[, chosen]
  log_counts <- assay(sce_chosen, "logcounts_corrected")[gene, ]
  df_aux <- data.frame(gene = gene, contrast = con, log_counts = log_counts, 
                       label.main = sce_chosen$label.main, 
                       group = sce_chosen$group)
  lt[[i]] <- df_aux
}
df_to_plot <- do.call("rbind", lt)
head(df_to_plot)


x11()
ggplot(df_to_plot, aes(x = log_counts, fill = group)) +
  facet_wrap(contrast ~ label.main, scales = "free") +
  geom_density(alpha = 0.8) +
  scale_fill_manual(values = c(col_H[3], col_M[3], col_S[6]))
ggsave(file.path(wd, "DE_genes_density.pdf"), device = "pdf",
       width = 14, height = 8)
