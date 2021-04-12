
# Pkgs --------------------------------------------------------------------

library(magrittr)
## Read and QC
suppressMessages(library(DropletUtils))
suppressMessages(library(scater))
## To find the mitocondrial genes
suppressMessages(library(EnsDb.Hsapiens.v86))


setwd("./covid/")

# Metafile and arqs -------------------------------------------------------

meta <- read.delim("./data/meta.txt")
ff <- list.files(path = "./data", pattern = ".h5")

# To read the data --------------------------------------------------------

getDescr <- function(ff) {
  tmp <- read10xCounts(paste0("./data/",ff), col.names = TRUE)
  
  rownames(tmp) <- uniquifyFeatureNames(
    rowData(tmp)$ID, 
    rowData(tmp)$Symbol
  )
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
    dplyr::filter(sample == aux) %>% 
    dplyr::select(sample_new) %>% 
    dplyr::pull() %>% 
    as.character()
  
  stats <- perCellQCMetrics(
    x = tmp,
    subsets = list(mito = which(location == "MT"))
  ) %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(sum, detected, subsets_mito_sum, subsets_mito_detected) %>%
    dplyr::mutate(sample = id)
  
  gen0 <- sum(rowSums(counts(tmp)) == 0)
  mt_sz <- sum(stats$subsets_mito_sum)
  descr <- dplyr::tibble(
    id = id,
    cells = ncol(tmp),
    umi_size = sum(stats$sum),
    genes = nrow(tmp),
    zero_genes = paste0(
      prettyNum(gen0, big.mark = ".", decimal.mark = ","), " (",
      prettyNum(round(100 * gen0/nrow(tmp), 1), big.mark = ".",
                decimal.mark = ",", nsmall = 1), "%)" 
    ),
    mito = sum(location == "MT", na.rm = TRUE),
    mito_size = paste0(
      prettyNum(mt_sz, big.mark = ".", decimal.mark = ","), " (",
      prettyNum(round(100 * mt_sz/umi_size, 1), big.mark = ".",
                decimal.mark = ",", nsmall = 1), "%)" 
    )
  ) %>% 
    dplyr::mutate_if(is.integer, prettyNum, big.mark = ".", decimal.mark = ",")
  
  list(
    descr = descr,
    stats = stats
  )
}


# To plot the data --------------------------------------------------------

plotDescr <- function(lt) {
  tmp <- do.call("rbind", purrr::map(lt, "stats")) %>% 
    dplyr::arrange(sample) %>% 
    dplyr::mutate(
      id = dplyr::case_when(
        sample %in% paste0("HC", 1:3) ~ "HC",
        sample %in% paste0("M", 1:3) ~ "M",
        sample %in% paste0("S", 1:6) ~ "S"
      )
    )
  
  p_umis <- ggplot(tmp, aes(x = sample, y = log1p(sum), fill = id)) +
    geom_violin(show.legend = FALSE) +
    cowplot::theme_cowplot() +
    ggtitle("Contagem dos UMIs")
  
  p_genes <- ggplot(tmp, aes(x = sample, y = log1p(detected), fill = id)) +
    geom_violin(show.legend = FALSE) +
    cowplot::theme_cowplot() +
    ggtitle("Genes detectados")
  
  p_mito <- ggplot(tmp, aes(x = sample, y = log1p(subsets_mito_sum), fill = id)) +
    geom_violin(show.legend = FALSE) +
    cowplot::theme_cowplot() +
    ggtitle("% de UMIs em genes da mitocondria")
  
  list(
    p_umis = p_umis,
    p_genes = p_genes,
    p_mito = p_mito
  )
}


# Run and save ------------------------------------------------------------

if (!any(list.files(path = "./data/") == "list_samples.Rdata")) {
  lt <- purrr::map(ff, getDescr)
  saveRDS(lt, file = "./data/list_samples.rds")
} else {
  lt <- readRDS(file = "./data/list_samples.rds")
}


if (TRUE) {
  tab <- do.call("rbind", purrr::map(lt, "descr")) %>% 
    dplyr::arrange(id)
  saveRDS(tab, file = "./data/tab_descr.rds")
  
  library(xtable) 
  xt <- xtable(tab)
  print.xtable(xt, include.rownames = F, sanitize.text.function = force)
  
}




