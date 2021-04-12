rm(list = ls())
library(foreach)
library(magrittr)

# Simulação para estimar o erro do tipo I.

# Simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", 
       encoding = "utf8")

# Construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", 
       encoding = "utf8")

# Modelos que serão avaliados
source("./simulation_study/functions/de_models_current.R", encoding = "utf8")

# Aréa sob a curva utilizando a regra do trapézio
trapz <- function(x, y) {
  m <- length(x)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2*m
  p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
  p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]
  return(0.5*(p1-p2))
}

# Calcula curva ROC
compute_roc <- function(x, is_de) {
  obj_rocr <- ROCR::prediction(predictions = x, labels = is_de)
  roc  <- ROCR::performance(obj_rocr, measure = "tpr", x.measure = "fpr")
  esp  <- 1 - roc@x.values[[1]]
  sen  <- 1 - roc@y.values[[1]]
  return(data.frame(esp = esp, sens = sen))
}

# Calcula o FDR com valor-p corrigido por BH
compute_fdr <- function(x, is_de, threshold = 0.05) {
  sig <- p.adjust(x, method = "BH") <= threshold
  if (any(sig, na.rm=TRUE)) { 
    return(sum(sig & !is_de, na.rm = TRUE) / sum(sig, na.rm = TRUE))
  } else {
    return(0)
  }
}

# Gera o arquivo com os resultados
save_fun <- function(x, label, scenario, mc, log_file) {
  discard <- is.na(x$pvalue)
  pvalue <- x$pvalue[!discard]
  is_de <- x$is_de[!discard]
  roc_out <- compute_roc(pvalue, is_de)
  output <- data.frame(label = label, scenario = scenario, mc = mc,
                       fdr = compute_fdr(pvalue, is_de),
                       auc = trapz(1 - roc_out$esp, roc_out$sens))
  write.table(output, file = log_file, sep = "\t", append = TRUE, 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}


# Leitura dos parâmetros estimados dos dados
parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")
parms_6_samples <- readRDS("./simulation_study/data/parms_6_samples.rds")

# Arquivos para exportar resultados
log_file <- "res_power.txt"

# Remove arquivos para não sobreescrever resultados
setwd("./simulation_study/data")
rm_logs <- dir(pattern = "log_power_")
rm_sces <- dir(pattern = "sce_power_")
junks <- c(log_file, rm_sces, rm_logs)

if (any(file.exists(junks))) {
  file.remove(junks[which(file.exists(junks))])  
}

# Iterações monte carlo
M <- 4
# Modelos
MODELS <- c("edgeR_bulk", "edgeR_sc", "lmm")

# Cenários
n_samples <- 3
n_group <- 2
var_effect <- c(0.5, 1e6)
pct_de <- c(0.05, 0.1, 0.2)
fc <- c(2, 4)
n_genes <- 2000

df_scenario <- expand.grid(n_genes = n_genes, n_group = n_group,
                           n_samples = n_samples, var_effect = var_effect,
                           pct_de = pct_de, fc = fc)
df_scenario$id <- 1:nrow(df_scenario)
head(df_scenario)


# Cria cluster
# cl <- parallel::makeCluster(nrow(df_scenario))
cl <- parallel::makeForkCluster(nrow(df_scenario))
doParallel::registerDoParallel(cl)

# Run!
ini <- proc.time()
set.seed(6669)
# foreach(i = 1:nrow(df_scenario)) %dopar% {
for (i in 1:nrow(df_scenario)) {
  
  # Inicializa os parâmetros
  n_ge <- df_scenario$n_genes[i]
  n_gr <- df_scenario$n_group[i]
  v_ef <- df_scenario$var_effect[i]
  p_de <- df_scenario$pct_de[i]
  fc <- df_scenario$fc[i]
  
  if (df_scenario$n_samples[i] == 3) {
    est_parms <- parms_3_samples
  } else {
    est_parms <- parms_6_samples
  }
  
  # Repete a simulação M vezes
  lt_roc <- list()
  counter <- 1L
  for(j in 1:M) {
    # Simula dados
    sim_data <- simulate_sce(estimated_parameters = est_parms, n_genes = n_ge, 
                             range_n_cells = c(200, 500), n_group = n_gr, 
                             var_effect = v_ef, n_cell_type = 2, 
                             pct_de = p_de, fc = fc)
    
    # Obt�m objeto sce
    sce <- construct_sce(simul_obj = sim_data, get_redDim = FALSE, 
                         get_cluster = FALSE, verbose = FALSE,
                         corrected = FALSE)

    # if (j == 1) {
    #   saveRDS(sce, file = paste0("sce_power_", i, ".rds"))
    # }
    
    # Ajusta os modelos e salva resultados
    for (m in MODELS) {
      write.table(data.frame(i, j, m), 
                  file = paste0("log_power_", i, ".txt"), 
                  sep = "\t", append = TRUE, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE)
      # cat("Cenário:", i, "| Iteração:", j, "| Método:", m, "\n")
      fit <- fit_models(sce, model = m)
      save_fun(x = fit, label = m, scenario = i, mc = j, log_file = log_file)
      df_roc <- compute_roc(fit$pvalue, fit$is_de)
      df_roc$method <- m
      df_roc$id <- 1:nrow(df_roc)
      lt_roc[[counter]] <- df_roc
      counter <- counter + 1L
    }
  }
  
  res_roc <- do.call("rbind", lt_roc) %>% 
    dplyr::group_by(method, id) %>% 
    dplyr::summarise(esp = mean(esp, na.rm = TRUE), 
                     sens = mean(sens, na.rm = TRUE))
  
  saveRDS(res_roc, file = paste0("roc_res_", i, ".rds"))
}
fim <- proc.time() - ini
parallel::stopCluster(cl)
