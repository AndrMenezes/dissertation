library(scran)
library(edgeR)
library(parallel)
library(lme4)

# Estima parâmetros dos dados para utilizar na simulação.
# Utiliza os dados do experimento scRNA-seq de covid.

estimate_parameters <- function(counts, sample_id, sf = NULL, export = TRUE, 
                                path_name = "./simulation_study/data", file_name) {
  
  # Cria um objeto do pacote edgeR ------------------------------------------
  
  # Este objeto já contém o tamanho da biblioteca que serão utilizadas na simulação
  
  time_ini <- proc.time()
  
  cat("Criando objeto edgeR \n \n")
  
  edgeR_obj <- DGEList(counts)
  
  # Estima fator de escala usando deconvolution 

  if (is.null(sf)) {
    cat("Estimando fator de escala via deconvolution \n")
    ini <- proc.time()
    cluster_id <- quickCluster(edgeR_obj$counts,
                               BPPARAM = BiocParallel::MulticoreParam(6))
    sf <- calculateSumFactors(edgeR_obj$counts, clusters = cluster,
                              BPPARAM = BiocParallel::MulticoreParam(6))
    fim <- proc.time() - ini
    cat("Estimação do fator de escala concluída em ", fim[3]/60, " minutos \n \n")
  }  
  edgeR_obj$samples$norm.factors <- sf / edgeR_obj$samples$lib.size
  edgeR_obj$samples$norm.factors <- edgeR_obj$samples$norm.factors / 
    exp(mean(log(edgeR_obj$samples$norm.factors))) 
  
  
  # Os fatores de escala entram nos modelos como offset para estimação
  # dos parâmetros da BN de cada gene
  
  
  # Estima para cada gene o parâmetro do dispersão da BN 
  
  cat("Estimando parâmetros de dispersão \n")
  
  ini <- proc.time()
  design_disp <- model.matrix(~sample_id)
  edgeR_obj <- estimateDisp(edgeR_obj, design_disp, prior.df = 0, trend = "none")
  fim <- proc.time() - ini
  
  cat("Estimação da dispersão concluída em ", fim[3]/60, " minutos \n \n")
  
  # Estima para cada gene a média da BN 
  
  cat("Estimando parâmetros de média \n")
  
  ini <- proc.time()
  centered_off <- getOffset(edgeR_obj)
  centered_off <- centered_off - mean(centered_off)
  logmeans <- mglmOneGroup(edgeR_obj$counts, offset = centered_off, 
                           dispersion = edgeR_obj$tagwise.dispersion)
  fim <- proc.time() - ini
  
  cat("Estimação da média concluída em ", fim[3]/60, " minutos \n \n")
  
  # Estima a variancia do efeito aleatório de individuo 
  
  cat("Estimando variância do efeito aleatório \n")
  
  ini <- proc.time()
  sigma2_collected <- mclapply(seq_len(nrow(edgeR_obj)), function(i) {
    output <- NA_real_
    try({ 
      df <- data.frame(Counts = as.integer(edgeR_obj$counts[i, ]), 
                       sample = sample_id)
      fit <- glmer(formula = Counts ~ offset(log(sf)) + (1 | sample), 
                   data = df, family = negative.binomial(
                     1 / edgeR_obj$tagwise.dispersion[i]))
      output <- unlist(VarCorr(fit))
    })
    return(output)
  }, mc.cores = 10)
  fim <- proc.time() - ini
  
  cat("Estimação da variância concluída em ", fim[3]/60, " minutos \n \n")
  
  # Lista com parâmetros estimados
  
  estimated_parameters <- list(
    cells_library_size = edgeR_obj$samples$lib.size,
    cells_size_factor = sf,
    genes_mean = exp(logmeans),
    genes_dispersion = edgeR_obj$tagwise.dispersion,
    sigma2_ranef = unlist(sigma2_collected),
    sample_id = sample_id
  )
  
  time_fim <- proc.time() - time_ini
  
  cat("Tempo total percorrido ", time_fim[3]/60, " minutos \n\n")
  
  # Exporta resultados
  
  if (export) {
    saveRDS(estimated_parameters, file = file.path(path_name, file_name))
  } else {
    return(estimated_parameters)
  }
}
