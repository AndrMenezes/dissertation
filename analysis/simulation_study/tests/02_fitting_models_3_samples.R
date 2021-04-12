rm(list = ls())

# Função para simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", encoding = "utf8")

# Função para construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", encoding = "utf8")

# Script com funções para ajuste dos modelos
source("./simulation_study/functions/de_models_current.R", encoding = "utf8")

# Parâmetros estimados utilizando informação de 3 indivíduos

parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")
str(parms_3_samples)


# Gerando dados

set.seed(66)
simul_nde <- simulate_sce(estimated_parameters = parms_3_samples,
                          n_genes = 2000, n_group = 3, var_effect = 1e6,
                          range_n_cells = c(200, 500))
sce <- construct_sce(simul_obj = simul_nde, get_redDim = FALSE)


# edgeR single-cell
system.time(
  edgeR_sc <- fit_models(sce, model = "edgeR_sc")
)
mean(edgeR_sc <= 0.05)

# edgeR pseudo bulk
system.time(
  edgeR_bulk <- fit_models(sce, model = "edgeR_bulk")
)
mean(edgeR_bulk$pvalue <= 0.05)

# Mixed linear model in log-corrected single-cell
system.time(
  mlm_sc <- fit_models(sce[1:100, ], model = "lmm")
)
mean(mlm_sc$pvalue <= 0.05)

