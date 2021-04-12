rm(list = ls())

wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/04_capitulo4/01_simulation_study"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,
                                            format="f"))}
# Simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", 
       encoding = "utf8")

# Construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", 
       encoding = "utf8")


# Packages ----------------------------------------------------------------

## For utils
suppressMessages(library(SingleCellExperiment))

## For PCA and UMAP
suppressMessages(library(scater))

## For %>% 
library(magrittr)

## For plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

## For tex tables
library(xtable)

# Cenários da simulação para erro do tipo I -------------------------------

n_samples <- c(3, 6)
n_group <- c(2, 3)
var_effect <- c(0.1, 1e6)
df_scenario <- expand.grid(n_group = n_group, n_samples = n_samples, 
                           var_effect = var_effect)
df_scenario$scenario <- 1:nrow(df_scenario)
head(df_scenario)


# Leitura dos dados, calcula PCA e UMAP para cada cenário -----------------

l_sces <- lapply(1:nrow(df_scenario), function(i) {
  x <- readRDS(paste0("./simulation_study/data/type_I_error/sce_scenario_", i, ".rds"))
  x <- runPCA(x)
  x <- runUMAP(x)
  cat(i, "\n")
  return(x)
})


# Tabela com informações dos cenários -------------------------------------
df_scenario$ncells <- prettyNum(sapply(l_sces, ncol), big.mark = ".", 
                                decimal.mark = ",")
df_scenario$var_effect <- ifelse(df_scenario$var_effect == 0.1, "Sim", "Não")
print(xtable(df_scenario[, c(4, 3, 1, 2, 5)],  digits = 0), 
      include.rownames = FALSE)


# Gráficos UMAP  ----------------------------------------------------------

l_umap <- lapply(l_sces, function(x) {
  tb <- reducedDim(x, "UMAP") %>% 
    dplyr::as_tibble()
  names(tb) <- c("UMAP_1", "UMAP_2")
  tb$group <- x$group
  tb$individuo <- x$individuo
  tb$cell_type <- x$cell_type
  return(tb)
})
names(l_umap) <- 1:length(l_umap)
tb_umap <- dplyr::bind_rows(l_umap, .id = "scenario")
tb_umap$scenario <- paste0("Cenário ", tb_umap$scenario)

p_umaps <- ggplot(tb_umap, aes(x = UMAP_1, y = UMAP_2, col = individuo, 
                               shape = cell_type)) +
  facet_wrap(~scenario, nrow = 2) +
  geom_point(alpha = 0.6, size = 3.0) +
  labs(x = "UMAP 1", y = "UMAP 2", col = "Indivíduo", shape = "Tipo de célula") +
  theme(legend.position = "top")
x11(); p_umaps
ggsave(filename = "umaps_type_I.png", plot = p_umaps, path = wd, width = 10,
       height = 8)


