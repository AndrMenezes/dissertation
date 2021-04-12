rm(list = ls())

# Simulação para estimar o erro do tipo I.

# Simular dados de scRNAseq
source("./simulation_study/functions/simulating_scRNA_counts.R", 
       encoding = "utf8")

# Construir um objeto SCE
source("./simulation_study/functions/construct_sce_object.R", 
       encoding = "utf8")

# Modelos que serão avaliados
source("./simulation_study/functions/de_models_power.R", encoding = "utf8")

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


set.seed(123)
parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")

sim_data <- simulate_sce(estimated_parameters = parms_3_samples, 
                         n_genes = 2000, range_n_cells = c(200, 500), 
                         n_group = 2, var_effect = 0.1, n_cell_type = 1, 
                         pct_de = 0.2, fc = 2)
sce <- construct_sce(simul_obj = sim_data, get_redDim = FALSE, 
                     get_cluster = FALSE, verbose = FALSE, vst = FALSE)


# genes <- rownames(sce)[which(rowData(sce)$is_de)]
# plotExpression(sce, features = genes[1:10], x = "group")


# Abordagens edgeR usando pseudo-bulk -------------------------------------
edgeR_LRT <- DE_edgeR(sce, type = "bulk", eb = FALSE)
compute_fdr(x = edgeR_LRT$pvalue, is_de = edgeR_LRT$is_de)
table(p.adjust(edgeR_LRT$pvalue, "BH") <= 0.05, edgeR_LRT$is_de)

edgeR_QL <- DE_edgeR(sce, type = "bulk", eb = TRUE)
compute_fdr(x = edgeR_QL$pvalue, is_de = edgeR_QL$is_de)
table(p.adjust(edgeR_QL$pvalue, "BH") <= 0.05, edgeR_QL$is_de)


# Abordagem regressão quantílica ------------------------------------------

system.time(qr <- quant_reg(sce))
compute_fdr(x = qr$pvalue, is_de = qr$is_de)
table(p.adjust(qr$pvalue, "BH") <= 0.05, qr$is_de)


# Abordagem modelo linear misto -------------------------------------------

## Utilizando as contagens log-normalizadas pelo método deconvolution
lmm_lc  <- lmm(sce, name_assay = "logcounts", weight = TRUE)
compute_fdr(x = lmm_lc$pvalue, is_de = lmm_lc$is_de, threshold = 0.01)
hist(lmm_lc $sigma2_ranef, breaks = 50)
mean(lmm_lc $sigma2_ranef == 0)
icc_lmm_lc <- lmm_lc$sigma2_ranef / (lmm_lc $sigma2_ranef + lmm_lc $sigma2)
range(icc_lmm_lc)

## Utilizando as contagens log-normalizadas pelo método deconvolution e corrigidas
lmm_lcc <- lmm(sce, name_assay = "logcounts")
hist(lmm_lcc$sigma2_ranef, breaks = 50)
mean(lmm_lcc$sigma2_ranef == 0)
icc_lmm_lcc <- lmm_lcc$sigma2_ranef / (lmm_lcc$sigma2_ranef + lmm_lcc$sigma2)
range(icc_lmm_lcc)
compute_fdr(x = lmm_lcc$pvalue, is_de = lmm_lcc$is_de)

out5 <- lmm(sce, name_assay = "vst")
hist(out5$sigma2_ranef, breaks = 50)
mean(out5$sigma2_ranef == 0)
icc5 <- out5$sigma2_ranef / (out5$sigma2_ranef + out5$sigma2)
hist(icc5)
compute_fdr(x = out5$pvalue, is_de = out5$is_de, threshold = 0.01)



########################
# Tentativa de modelagem usando a distribuição gamma


