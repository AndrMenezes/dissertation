rm(list = ls())

setwd("./covid/")
wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/cap5"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,
                                            format="f"))}

# Packages ----------------------------------------------------------------

## Read the data
suppressMessages(library(HDF5Array))

suppressMessages(library(SingleCellExperiment))

## For linear mixed model
library(lme4)

## For multiple comparison tests
library(emmeans)

emm_options(lmerTest.limit = 1e10)

library(magrittr)

## For plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

# Colors to plot
col_H <- colorspace::sequential_hcl(4, palette = "Blues 3")[-4]
col_M <- colorspace::sequential_hcl(4, palette = "Greens 3")[-4]
col_S <- colorspace::sequential_hcl(7, palette = "Reds 3")[-7]

# Data --------------------------------------------------------------------

sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")
table(sce$label.main, sce$group)
sce <- sce[, sce$label.main != "Neutrophils"]


# Tabela com frequencia das observações -----------------------------------

tab_freq <- table(sce$label.main, sce$group) %>% 
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



# Seleciona os genes que tem menos que 90% de contagem igual a 0 ----------

pct_zero_counts <- DelayedMatrixStats::rowMeans2(assay(sce, "counts") == 0)
df_info <- data.frame(gene = rownames(sce), pct_zero_counts=pct_zero_counts)
df_info <- df_info[order(df_info$pct_zero), ]
gkeep <- df_info[df_info$pct_zero < 0.9, ]$gene
length(gkeep)
100 * length(gkeep) / nrow(df_info)

sce <- sce[gkeep, ]



# Modelo linear misto com os pesos do voom --------------------------------

## Dataframe com as variáveis do modelo
df <- as.data.frame(colData(sce)[, c("group", "label.main", "sample")]) 
df$group <- factor(df$group)
df$label.main <- factor(df$label.main)
df$sample <- factor(df$sample)
df$logcounts <- 0

## Estimação dos pesos para cada gene
voom_fit <- limma::voom(counts = counts(sce), 
                        design = model.matrix(~group * label.main, df))
wi <- voom_fit$weights
dim(wi)
hist(wi[1, ])
hist(wi[10, ])

## Matriz de expressão gênica
expr_mat <- logcounts(sce)
dim(expr_mat)

## Estimação do modelo misto para cada gene
ini <- proc.time()
out_fits <- lapply(seq_len(nrow(sce)), function(g) {
  df$logcounts <- expr_mat[g, ]
  tryCatch(
    fit <- lmer(formula = logcounts ~ 0 + group * label.main + (1 | sample),
                data = df, weights = wi[g ,]), 
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    emm <- pairs(emmeans(fit, specs = "group", by = "label.main", 
                         lmer.df = "satterthwaite"), 
                 adjust = "none")
    df_out <- as.data.frame(summary(emm))
    vc <- as.data.frame(VarCorr(fit))
    df_out$sigma2_ranef <- vc$vcov[1]
    df_out$sigma2 <- vc$vcov[2]
    df_out$gene <- rownames(sce)[g]
    cat(g, rownames(sce)[g], "OK! \n")
  } 
  else {
    df_out <- NULL
    cat(g, rownames(sce)[g], "Fail! \n")
  }
  return(df_out)
})
fim <- proc.time() - ini
fim[3] / 60 / 60
90/60
# 2h e 58 minutos

sum(sapply(out_fits, is.null))
saveRDS(out_fits, file = "./data/DE_mlm_logcounts.rds")
# out_fits <- readRDS(file = "./data/DE_mlm_logcounts.rds")

cols <- Reduce(intersect, lapply(out_fits, names))
df_res <- do.call("rbind", lapply(out_fits, function(x) x[, cols] ))
df_res <- df_res %>% 
  dplyr::group_by(gene) %>% 
  dplyr::mutate(fdr = p.adjust(p.value, method = "BH"),
                icc = sigma2_ranef / (sigma2_ranef + sigma2)) %>% 
  dplyr::arrange(-icc) %>% 
  dplyr::ungroup()
saveRDS(df_res, file = "./data/DE_mixed_model.rds")

# df_res <- readRDS(file = "./data/DE_mixed_model.rds")


# Número de ajustes que o efeito aleatório estimado foi zero
df_res %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>% 
  dplyr::count(zero_sigma2_ranef = sigma2_ranef == 0)
df_res %>% 
  dplyr::filter(sigma2_ranef == 0) %>% 
  dplyr::distinct(gene)


# Distribuição da variância estimada do efeito aleatório ------------------

range(df_res$icc)
range(df_res$sigma2_ranef)

df_res %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>% 
  dplyr::filter(sigma2_ranef != 0) %>% 
  ggplot(aes(x = icc)) +
  geom_density(fill = "grey78")
df_res %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>% 
  dplyr::filter(sigma2_ranef != 0) %>% 
  ggplot(aes(x = sigma2_ranef)) +
  geom_density(fill = "grey78")


# Volcano plot ------------------------------------------------------------

tb_volcano <- df_res %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::mutate(
    label = ifelse(fdr <= quantile(fdr, probs = 0.01), gene, ""),
    contrast  = gsub(" - ", " versus ", contrast)
  )
tb_volcano %>% 
  dplyr::group_by(label.main) %>% 
  dplyr::count(label != "")

ggplot(tb_volcano, aes(x = estimate, y = -log10(fdr), label = label,
                       col = label.main)) +
  facet_wrap(~contrast, scales = "free_y") +
  geom_point(alpha = 0.4, size = 3) +
  # geom_text(aes(y = -log10(fdr) + 0.3),  col = "black") +
  ggrepel::geom_text_repel(col = "black", max.time = 60) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = bquote(~Log[2]~ 'fold change'), y = bquote("-" ~Log[10]~ "FDR"), 
       col = "") +
  theme(legend.position = "top")
ggsave(file.path(wd, "DE_volcano_mixed_model.pdf"), device = "pdf",
       width = 14, height = 8)


# Analisando expressão de alguns genes ------------------------------------

## Genes que foram mais significativos por tipo de célula e comparação

ntop <- 1
top_genes <- df_res %>% 
  dplyr::group_by(label.main, contrast) %>% 
  dplyr::slice_min(order_by = fdr, n = ntop)

## Tabela
tab <- top_genes[, c(9, 2, 1, 3, 6, 10)]

options(OutDec = ",")
print(xtable(tab, digits = c(rep(1, 4), 4, -2, -2)), include.rownames=FALSE)


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


x11()
ggplot(df_to_plot, aes(x = group, y = log_counts, col = group)) +
  facet_wrap(contrast ~ label.main, scales = "free") +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, width = 0.4, groupOnX = TRUE, 
                               bandwidth = 1) +
  scale_color_manual(values = c(col_H[3], col_M[3], col_S[6]))
ggsave(file.path(wd, "DE_genes_point.pdf"), device = "pdf",
       width = 14, height = 8)

x11()
ggplot(df_to_plot, aes(y = group, x = log_counts, fill = group)) +
  facet_wrap(contrast ~ label.main, scales = "free") +
  ggridges::geom_density_ridges(alpha = 0.9) +
  scale_fill_manual(values = c(col_H[3], col_M[3], col_S[6]))
ggsave(file.path(wd, "DE_genes_ridges.pdf"), device = "pdf",
       width = 14, height = 8)



#######################
scater::plotExpression(sce[, sce$label.main == "Monocyte"], 
                       features = "CSTA", x = "group", 
                       exprs_values = "logcounts")


sce_filt <- sce[, sce$label.main == "Monocyte"]
tb <- assay(sce_filt, "logcounts_corrected")[c("CSTA", "MT-CO3"), ] %>% 
  t() %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(group = sce_filt$group, sample = sce_filt$sample)

ggplot(tb, aes(x = CSTA, fill = group)) +
  geom_density(alpha = 0.6)

ggplot(tb, aes(x = CSTA, y = sample, fill = group)) +
  ggridges::geom_density_ridges(alpha = 0.7)

ggplot(tb, aes(x = `MT-CO3`, col = group)) +
  stat_ecdf()
