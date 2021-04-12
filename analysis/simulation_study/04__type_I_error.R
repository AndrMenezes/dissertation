rm(list = ls())
library(foreach)

# Simulação para estimar o erro do tipo I.

# Simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", 
       encoding = "utf8")

# Construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", 
       encoding = "utf8")

# Modelos que serão avaliados
source("./simulation_study/functions/de_models_current.R", encoding = "utf8")

# Leitura dos parâmetros estimados dos dados
parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")
parms_6_samples <- readRDS("./simulation_study/data/parms_6_samples.rds")

# Arquivo para exportar resultados
log_file <- "res_type_I_error.txt"

# Remove arquivos para não sobreescrever resultados
setwd("./simulation_study/data")
rm_logs <- dir(pattern = "log_")
rm_sces <- dir(pattern = "sce_")
junks <- c(log_file, rm_sces, rm_logs)
file.remove(junks)

# Iterações monte carlo
M <- 5
# Modelos
MODELS <- c("edgeR_bulk", "edgeR_sc", "lmm")

# Gera o arquivo com os resultados
save_fun <- function(label, scenario, x, mc, log_file) {
  alphas <- c(0.05, 0.01, 0.001)
  discard <- is.na(x)
  total <- sum(!discard)
  for (alpha in alphas) {
    # Protect against NAs
    output <- data.frame(label, scenario, mc, alpha, 
                         sum(x <= alpha & !discard) / total)
    write.table(output, file = log_file, sep = "\t", append = TRUE, 
                col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
}

# Cenários
n_samples <- c(3, 6)
n_group <- c(2, 3)
var_effect <- c(0.1, 1e6)
n_genes <- 2000

df_scenario <- expand.grid(n_genes = n_genes, n_group = n_group,
                           n_samples = n_samples, var_effect = var_effect)
df_scenario$id <- 1:nrow(df_scenario)
head(df_scenario)


# Cria cluster
# cl <- parallel::makeCluster(nrow(df_scenario))
cl <- parallel::makeForkCluster(nrow(df_scenario))
doParallel::registerDoParallel(cl)

# Run!
ini <- proc.time()
set.seed(6669)
foreach(i = 1:nrow(df_scenario)) %dopar% {
# for (i in 1:nrow(df_scenario)) {
  
  # Inicializa os parâmetros
  n_ge <- df_scenario$n_genes[i]
  n_gr <- df_scenario$n_group[i]
  v_ef <- df_scenario$var_effect[i]
  id    <- i
  
  if (df_scenario$n_samples[i] == 3) {
    est_parms <- parms_3_samples
  } else {
    est_parms <- parms_6_samples
  }
  
  # Repete a simulação M vezes
  for(j in 1:M) {
    # Simula dados
    sim_data <- simulate_sce(estimated_parameters = est_parms, n_genes = n_ge, 
                             range_n_cells = c(200, 500), n_group = n_gr, 
                             var_effect = v_ef, n_cell_type = 2)
    
    # Obtém os valores corrigidos utilizando batchelor::rescaleBatches
    sce <- construct_sce(simul_obj = sim_data, get_redDim = FALSE, 
                         get_cluster = FALSE, verbose = FALSE)

    if (j == 1) {
      saveRDS(sce, file = paste0("./simulation_study/data/sce_scenario_", 
                                     id, ".rds"))
    }
    
    # Ajusta os modelos e salva resultados
    for (m in MODELS) {
      write.table(data.frame(id, j, m), 
                  file = paste0("./simulation_study/data/log_scenario_", 
                                id, ".txt"), 
                  sep = "\t", append = TRUE, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE)
      # cat("Cenário:", id, "| Iteração:", j, "| Método:", m, "\n")
      fit <- fit_models(sce, model = m)
      save_fun(label = m, scenario = id, mc = j, x = fit$pvalue,
               log_file = log_file)
    }
  }
}
fim <- proc.time() - ini
parallel::stopCluster(cl)
