rm(list = ls())

# Função para simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", encoding = "utf8")

# Função para construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", encoding = "utf8")

# Parâmetros estimados utilizando informação de 3 indivíduos

parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")
str(parms_3_samples)

# Dados 1: sem diferença de grupos ----------------------------------------

set.seed(6969)
simul_not_de <- simulate_sce(estimated_parameters = parms_3_samples, 
                             n_genes = 1000, n_group = 2, var_effect = 1e6,
                             n_cell_type = 2, range_n_cells = c(200, 500))
sce_not_de <- construct_sce(simul_obj = simul_not_de)


gridExtra::grid.arrange(
  plotReducedDim(sce_not_de, dimred = "UMAP", colour_by = "group") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_not_de, dimred = "UMAP_corrected", colour_by = "group") +
    ggtitle("Com correção")
)

gridExtra::grid.arrange(
  plotReducedDim(sce_not_de, dimred = "UMAP", colour_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_not_de, dimred = "UMAP_corrected", colour_by = "individuo") +
    ggtitle("Com correção")
)



# Dados 2: com alguns genes diferentes ------------------------------------

set.seed(6969)
simul_de <- simulate_sce(estimated_parameters = parms_3_samples, 
                         n_genes = 1000, n_group = 2,
                         pct_de = 0.01, var_effect = 0.07,
                         fc = 3, range_n_cells = c(200, 500))
sce_de <- construct_sce(simul_obj = simul_de)
ncol(sce_de)

gridExtra::grid.arrange(
  plotReducedDim(sce_de, dimred = "UMAP", colour_by = "group") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_de, dimred = "UMAP_corrected", colour_by = "group") +
    ggtitle("Com correção")
)

gridExtra::grid.arrange(
  plotReducedDim(sce_de, dimred = "UMAP", colour_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_de, dimred = "UMAP_corrected", colour_by = "individuo") +
    ggtitle("Com correção")
)



# Dados 3: 3 grupos sem genes diferentes ----------------------------------

set.seed(6969)
simul_not_de_3 <- simulate_sce(estimated_parameters = parms_3_samples, 
                               n_genes = 1000, n_group = 3, 
                               range_n_cells = c(200, 500))
sce_not_de_3 <- construct_sce(simul_obj = simul_not_de_3)
ncol(sce_not_de_3)

gridExtra::grid.arrange(
  plotReducedDim(sce_not_de_3, dimred = "UMAP", colour_by = "group") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_not_de_3, dimred = "UMAP_corrected", colour_by = "group") +
    ggtitle("Com correção")
)

gridExtra::grid.arrange(
  plotReducedDim(sce_not_de_3, dimred = "UMAP", colour_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_not_de_3, dimred = "UMAP_corrected", colour_by = "individuo") +
    ggtitle("Com correção")
)




# Dados 4: 3 grupos com alguns genes diferentes ---------------------------

set.seed(6969)
simul_de_3 <- simulate_sce(estimated_parameters = parms_3_samples, 
                           n_genes = 1000, n_group = 3,
                           n_de = 400, fc = 3,
                           range_n_cells = c(200, 500))
sce_de_3 <- construct_sce(simul_obj = simul_de_3)


gridExtra::grid.arrange(
  plotReducedDim(sce_de_3, dimred = "UMAP", colour_by = "group") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_de_3, dimred = "UMAP_corrected", colour_by = "group") +
    ggtitle("Com correção")
)

gridExtra::grid.arrange(
  plotReducedDim(sce_de_3, dimred = "UMAP", colour_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_de_3, dimred = "UMAP_corrected", colour_by = "individuo") +
    ggtitle("Com correção")
)





# Dados 5: 4 grupos sem genes diferentes ----------------------------------

set.seed(6969)
simul_not_de_4 <- simulate_sce(estimated_parameters = parms_3_samples, 
                               n_genes = 1000, n_group = 4,
                               range_n_cells = c(200, 500))
sce_not_de_4 <- construct_sce(simul_obj = simul_not_de_4)


gridExtra::grid.arrange(
  plotReducedDim(sce_not_de_4, dimred = "UMAP", colour_by = "group") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_not_de_4, dimred = "UMAP_corrected", colour_by = "group") +
    ggtitle("Com correção")
)

gridExtra::grid.arrange(
  plotReducedDim(sce_not_de_4, dimred = "UMAP", colour_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_not_de_4, dimred = "UMAP_corrected", colour_by = "individuo") +
    ggtitle("Com correção")
)




# Dados 6: 4 grupos com metade dos genes diferentes -----------------------

set.seed(6969)
simul_de_4 <- simulate_sce(estimated_parameters = parms_3_samples, 
                           n_genes = 1000, n_group = 4, 
                           n_de = 400, fc = 3,
                           range_n_cells = c(200, 500))
sce_de_4 <- construct_sce(simul_obj = simul_de_4)


gridExtra::grid.arrange(
  plotReducedDim(sce_de_4, dimred = "UMAP", colour_by = "group",
                 shape_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_de_4, dimred = "UMAP_corrected", colour_by = "group") +
    ggtitle("Com correção") 
)

gridExtra::grid.arrange(
  plotReducedDim(sce_de_4, dimred = "UMAP", colour_by = "individuo") +
    ggtitle("Sem correção"),
  plotReducedDim(sce_de_4, dimred = "UMAP_corrected", colour_by = "individuo") +
    ggtitle("Com correção")
)


# Salvando os dados simulados ---------------------------------------------

path <- "./simulation_study/data"
# 2 grupos
saveRDS(sce_not_de, file.path(path, "test_samples_3_groups_2_nde.rds"))
saveRDS(sce_de, file.path(path, "test_samples_3_groups_2_de"))
# 3 grupos
saveRDS(sce_not_de_3, file.path(path, "test_samples_3_groups_3_nde.rds"))
saveRDS(sce_de_3, file.path(path, "test_samples_3_groups_3_de.rds"))
# 4 grupos
saveRDS(sce_not_de_4, file.path(path, "test_samples_3_groups_4_nde.rds"))
saveRDS(sce_de_4, file.path(path, "test_samples_3_groups_4_de.rds"))


