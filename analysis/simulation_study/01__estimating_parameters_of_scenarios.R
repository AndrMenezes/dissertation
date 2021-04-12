
# Ler os dados
library(HDF5Array)

# Função para estimar os parâmetros
source("./simulation_study/functions/parameters_estimation.R")

# Leitura do objeto sce
sce <- loadHDF5SummarizedExperiment(dir = "./data/all_samples")

# Lista com os indívudos considerados
samples <- list(c("HC1", "M1", "S1"), 
                c("HC1", "M1", "S1", "HC2", "M2", "S2"))

for (j in seq_along(samples)) {
  
  ff <- paste0("parms_", length(samples[[j]]), "_samples.rds")
  
  cat("   Gerando parâmetros considerando", 
      length(samples[[j]]), " individuos:  \n \n \n")
  
  sce_chosen <- sce[, sce$sample %in% samples[[j]]]
  estimate_parameters(counts = counts(sce_chosen), sample_id = sce_chosen$sample,
                      sf = sce_chosen$sizeFactor, file_name = ff)
  
  cat("\n \n \n ")
}

