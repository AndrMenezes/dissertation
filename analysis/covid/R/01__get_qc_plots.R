rm(list = ls())

# Pkgs --------------------------------------------------------------------

library(magrittr)
## Read and QC
suppressMessages(library(DropletUtils))
suppressMessages(library(scater))
## To find the mitocondrial genes
suppressMessages(library(EnsDb.Hsapiens.v86))

# Metafile and arqs -------------------------------------------------------

setwd("./covid/")
meta <- read.delim("./data/meta.txt")
ff <- list.files(path = "./data", pattern = ".h5")

# Read the data and perform standard QC -----------------------------------

getQC <- function(ff, sce=FALSE) {
  
  tmp <- read10xCounts(paste0("./data/",ff), col.names = TRUE)
  
  rownames(tmp) <- uniquifyFeatureNames(
    rowData(tmp)$ID, 
    rowData(tmp)$Symbol
  )
  
  ## Exclude genes with zero counts for all cells
  no_zero <- rowSums(counts(tmp)) > 0
  tmp <- tmp[which(no_zero), ]

  location <- mapIds(
    EnsDb.Hsapiens.v86, 
    keys = rowData(tmp)$ID,
    column = "SEQNAME", 
    keytype = "GENEID"
  )
  
  aux <- stringr::str_remove(
    string = ff, 
    pattern = ".h5"
  )
  
  id <- meta %>% 
    dplyr::filter(sample == aux) 
    
  
  stats <- perCellQCMetrics(
    x = tmp,
    subsets = list(mito = which(location == "MT"))
  )
  

  thresholds <- rbind(
    attr(isOutlier(stats$sum, log=TRUE, type="lower"), "thresholds"),
    attr(isOutlier(stats$detected, log=TRUE, type="lower"), "thresholds"),
    attr(isOutlier(stats$subsets_mito_percent, type="higher"), "thresholds")
  )
  
  qc <- quickPerCellQC(
    df = stats,
    n_features = "detected",
    lib_size = "sum",
    percent_subsets = "subsets_mito_percent"
  )
  
  out <- cbind(
    stats, qc,
    sample = id$sample_new, 
    thresholds_libsize = thresholds[1,1], 
    thresholds_feature = thresholds[2,1],
    thresholds_mito = thresholds[3,2]
  ) %>% 
    dplyr::as_tibble()

  out
}


# Plot QC metrics ---------------------------------------------------------

plotQC <- function(tb, save=FALSE) {

  df <- tb %>% 
    dplyr::rename(my_discard = discard) %>% 
    dplyr::mutate(
      id = dplyr::case_when(
        sample %in% paste0("HC", 1:3) ~ "HC",
        sample %in% paste0("M", 1:3) ~ "M",
        sample %in% paste0("S", 1:6) ~ "S"
      ),
      their_discard = ifelse(sum < 1000 | detected < 200 | detected > 6000 | subsets_mito_percent > 10, TRUE, FALSE),
      discard = dplyr::case_when(
        my_discard==TRUE & their_discard==TRUE ~"both",
        my_discard==TRUE & their_discard==FALSE ~ "my",  
        my_discard==FALSE & their_discard==TRUE ~ "their",
        my_discard==FALSE & their_discard==FALSE ~ "none"  
      )
    ) %>% 
    dplyr::select(
      sample, sum, detected, subsets_mito_percent, 
      thresholds_libsize, thresholds_feature, thresholds_mito,
      my_discard, their_discard, discard
   )
  
  
  tab <- tb %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise(
      thresholds_libsize = max(thresholds_libsize),
      thresholds_feature = max(thresholds_feature),
      thresholds_mito = max(thresholds_mito),
      size = sum(low_lib_size), 
      features = sum(low_n_features), 
      mito = sum(high_subsets_mito_percent), 
      discard = sum(discard),
      total = dplyr::n(), 
      pct = discard/total*100
    )
  
  tab_compar <- df %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise(
      my_thresholds_libsize = max(thresholds_libsize),
      their_thresholds_libsize = 1000,
      
      my_thresholds_feature = max(thresholds_feature),
      their_thresholds_feature = "<200 e >6000",
      
      my_thresholds_mito = max(thresholds_mito),
      their_thresholds_mito = 10,
      
      my_discard = sum(my_discard),
      their_discard = sum(their_discard),
      
      total = dplyr::n(),
      my_pct = my_discard/total*100,
      their_pct = their_discard/total*100
    )

  aux <- theme(
    legend.position = "top", 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) 
  
  df <- df %>% 
    dplyr::mutate(their_discard = ifelse(their_discard, "Sim", "Não"))
  
  ## UMIs
  p_umis <- ggplot(df, aes(x = sample, y = sum)) +
    facet_wrap(~sample, scales = "free") +
    ggbeeswarm::geom_quasirandom(aes(col = their_discard), alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1) +
    labs(x = "", y = "Tamanho da biblioteca", col = "Baixa qualidade?") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    colorspace::scale_color_discrete_diverging(palette = "Berlin") +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    aux
  ## Genes
  p_genes <- ggplot(df, aes(x = sample, y = detected)) +
    facet_wrap(~sample, scales = "free") +
    ggbeeswarm::geom_quasirandom(aes(col = their_discard), alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1) +
    labs(x = "", y = "Genes detectados", col = "Baixa qualidade?") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    colorspace::scale_color_discrete_diverging(palette = "Berlin") +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    aux
  ## Mito
  p_mito <- ggplot(df, aes(x = sample, y = subsets_mito_percent )) +
    facet_wrap(~sample, scales = "free") +
    geom_violin(colour = "grey69", size = 1, alpha = 0.2, scale = "width", width = 0.8) +
    ggbeeswarm::geom_quasirandom(aes(col = their_discard), alpha = 0.6, width=0.4, groupOnX=TRUE, bandwidth=1) +
    labs(x = "", y = "% de genes mitocondriais", col = "Baixa qualidade?") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    colorspace::scale_color_discrete_diverging(palette = "Berlin") +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    aux
  ## Mito vs umis
  p_mito_umis <- ggplot(df, aes(x = log(sum), y = subsets_mito_percent, col = their_discard)) +
    geom_point(alpha = 0.6) +
    facet_wrap(~sample, scales = "free") +
    labs(x = "Tamanho da biblioteca (log)", y = "% genes mitocondriais", col = "Baixa qualidade?") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    colorspace::scale_color_discrete_diverging(palette = "Berlin") +
    cowplot::theme_cowplot() +
    cowplot::background_grid() +
    theme(legend.position = "top")
  
  if (save) {
    wd <- "../latex/figuras/04_capitulo4/02_data_analysis"
    ggsave(paste0(wd, "umis.png"), plot = p_umis, device = "png", width = 16, height = 12)
    ggsave(paste0(wd, "genes.png"), plot = p_genes, width = 16, height = 12)
    ggsave(paste0(wd, "mito.png"), plot = p_mito, width = 16, height = 12)
    ggsave(paste0(wd, "mito_umis.png"), plot = p_mito_umis, width = 16, height = 12)
    saveRDS(tab, file = "./data/qc_tab.rds")
    saveRDS(tab_compar, file = "./data/qc_tab_compar.rds")
  } 
  else {
    list(
      p_umis = p_umis,
      p_genes = p_genes,
      p_mito = p_mito,
      tab = tab
    )
  }
}

# Run and save ------------------------------------------------------------

tb <- purrr::map_df(ff, getQC)
plotQC(tb, save = TRUE)

tab <- readRDS(file = "./data/qc_tab.rds")
tab <- tb %>% 
  dplyr::select(sample, sum, detected, subsets_mito_percent) %>% 
  dplyr::mutate(discard_lib = sum < 1000,
                discard_gene = detected < 200 | detected > 6000,
                discard_mito = subsets_mito_percent > 10,
                discard = discard_lib | discard_gene | discard_mito) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(
    size = sum(discard_lib), 
    features = sum(discard_gene), 
    mito = sum(discard_mito), 
    discard = sum(discard),
    total = dplyr::n(), 
    pct = discard/total*100
  ) %>% 
  dplyr::select(-total)

xt <- xtable(tab)
print.xtable(xt, include.rownames = FALSE, sanitize.text.function = force)
