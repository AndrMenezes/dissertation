rm(list = ls())

wd <- "../latex/figuras/04_capitulo4/02_data_analysis"
setwd("./covid/")

# Pkgs --------------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## Feature selection
suppressMessages(library(scran))
# Read the data
suppressMessages(library(HDF5Array))

# Read the list of sces ---------------------------------------------------

samples <- c(paste0("HC", 1:3), paste0("M", 1:3), paste0("S", 1:6))
dir_to_read <- paste0("./data/", samples)
l_sce <- lapply(dir_to_read, function(d) loadHDF5SummarizedExperiment(dir = d))

# Feature selection -------------------------------------------------------

system.time(
  l_dec <- lapply(l_sce, modelGeneVar)
)

### Visualizing the fit
l_fit <- lapply(seq_along(l_dec), function(i) {
  
  meta_fit <- metadata(l_dec[[i]])
  out <- dplyr::tibble(
    gene = rownames(l_dec[[i]]),
    mean = meta_fit$mean,
    var = meta_fit$var,
    fit = meta_fit$trend(mean),
    sig_fdr = ifelse(l_dec[[i]]$FDR < 0.01, "Alta variância", "Baixa variância"),
    sig_bio = ifelse(l_dec[[i]]$bio > 0, "Alta variância", "Baixa variância"),
    sample = l_sce[[i]]$sample[1]
  )
  out
})

tb_fit <- do.call("rbind", l_fit)

tb_fit %>% 
  # dplyr::filter(sample == "HC1") %>% 
  dplyr::group_by(sample) %>% 
  dplyr::count(sig_fdr, sig_bio) %>% 
  dplyr::filter(sig_fdr == "Baixa variância", sig_bio == "Baixa variância")


p_fit <- tb_fit %>% 
  # dplyr::filter(sample == "HC1") %>% 
  ggplot(aes(x = mean, y = var)) +
  facet_wrap(~sample, scales = "free_x", nrow = 4) +
  geom_point(size = 2.2, alpha = 0.6, fill = "grey", shape = 21) +
  geom_line(aes(y = fit), size = 1.4, col = "blue") +
  labs(x = "Média do log(expressão normalizada)", 
       y = "Variância do log(expressão normalizada)") +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  theme(legend.position = "top")

ff <- file.path(wd, "feature_selection_samples.png")
ggsave(filename = ff, plot = p_fit, device = "png", width = 15, height = 10)

one_fit <- tb_fit %>% 
  dplyr::filter(sample == "S2") %>% 
  ggplot(aes(x = mean, y = var)) +
  geom_point(size = 2.2, alpha = 0.6, fill = "grey", shape = 21) +
  geom_line(aes(y = fit), size = 1.4, col = "blue") +
  labs(x = "Média do log(expressão normalizada)", 
       y = "Variância do log(expressão normalizada)") +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(4)) +
  theme(legend.position = "top")

ff <- file.path(wd, "one_feature_selection.png")
ggsave(filename = ff, plot = one_fit, device = "png", width = 12, height = 7)

## Combining
hvgs <- combineVar(l_dec)
hvgs <- hvgs[order(hvgs$bio, decreasing=TRUE),] 
hvgs
xtable(as.data.frame(hvgs[1:10, c(1:4, 6)]), digits=c(1, rep(4, 4), -2))

## Creating tables
one_dec <- l_dec[[1]]
one_dec <- one_dec[order(one_dec$bio, decreasing=TRUE),] 
xtable(as.data.frame(one_dec[1:10, c(1:4, 6)]), digits=c(1, rep(4, 4), -2))


# Select features and concatenate sces objects ----------------------------

chosen <- which(hvgs$FDR <= 0.05)
length(chosen)
length(chosen) / nrow(hvgs)

sce <- do.call("cbind", lapply(l_sce, function(x) x[chosen,]))
sce


# Saving ------------------------------------------------------------------
counts(sce) <- Matrix(as.matrix(counts(sce)), sparse=TRUE)
logcounts(sce) <- Matrix(as.matrix(logcounts(sce)), sparse=TRUE)
saveHDF5SummarizedExperiment(sce, dir = "./data/all_samples", verbose = TRUE)

