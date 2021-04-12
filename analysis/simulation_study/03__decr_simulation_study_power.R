rm(list = ls())

wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/cap5/simulation_study"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,
                                            format="f"))}
# Simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", 
       encoding = "utf8")

# Construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", 
       encoding = "utf8")

#  Dados com os parâmetros estimados
parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")

# Cenários para avaliar o poder -------------------------------------------

n_samples <- 3
n_group <- 2
var_effect <- c(0.1, 1e6)
pct_de <- c(0.05, 0.1, 0.2, 0.5)
fc <- c(2, 4)
n_genes <- 2000

df_scenario <- expand.grid(n_genes = n_genes, n_group = n_group,
                           n_samples = n_samples, pct_de = pct_de, fc = fc,
                           var_effect = var_effect)
df_scenario$scenario <- 1:nrow(df_scenario)


# Geração dos dados -------------------------------------------------------

set.seed(6669)
l_sces_power <- lapply(1:nrow(df_scenario), function(i) {
  sim_data <- simulate_sce(estimated_parameters = parms_3_samples, 
                           n_genes = df_scenario$n_genes[i], 
                           n_cell_type = 1, range_n_cells = c(200, 500), 
                           n_group = df_scenario$n_group[i], 
                           var_effect = df_scenario$var_effect[i], 
                           pct_de = df_scenario$pct_de[i],
                           fc = df_scenario$fc[i])
  x <- construct_sce(simul_obj = sim_data, get_redDim = FALSE, 
                     get_cluster = FALSE, verbose = FALSE)
  x <- runPCA(x)
  x <- runUMAP(x)
  saveRDS(x, paste0("./simulation_study/data/power/sce_scenario_", i, ".rds"))
  cat(i, "\n")
  return(x)
})


# Tabela com informações dos cenários do poder ----------------------------
df_scenario$var_effect <- ifelse(df_scenario$var_effect == 0.5, "Sim", "Não")
df_scenario$ncells <- prettyNum(sapply(l_sces_power, ncol), big.mark = ".",
                                decimal.mark = ",")
df_scenario$pct_de <- paste0(prettyNum(10 * df_scenario$pct_de, big.mark = ".",
                                       decimal.mark = ",", nsmall = 1),
                             "%")
print(xtable(df_scenario[, c(7, 6, 4, 5, 8)],  digits = c(rep(0, 4), 0, 0)), 
      include.rownames = FALSE)

# Gráficos UMAP  ----------------------------------------------------------

l_umap <- lapply(l_sces_power, function(x) {
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
tb_umap$scenario <- forcats::fct_relevel(
  paste0("Cenário ", tb_umap$scenario),
  paste0("Cenário ", 1:nrow(df_scenario))
)
                                         

p_umaps <- ggplot(tb_umap, aes(x = UMAP_1, y = UMAP_2, col = individuo, 
                               shape = group)) +
  facet_wrap(~scenario, nrow = 4, scales = "free") +
  geom_point(alpha = 0.6, size = 3.0) +
  labs(x = "UMAP 1", y = "UMAP 2", col = "Indivíduo", shape = "Grupo") +
  theme(legend.position = "top")
x11(); p_umaps
ggsave(filename = "umaps_power.pdf", plot = p_umaps, path = wd, width = 10,
       height = 8)


