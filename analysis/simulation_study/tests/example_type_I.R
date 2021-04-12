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

set.seed(123)
parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")

sim_data <- simulate_sce(estimated_parameters = parms_3_samples, 
                         n_genes = 2000, range_n_cells = c(200, 500), 
                         n_group = 2, var_effect = 0.1, n_cell_type = 1)
sce <- construct_sce(simul_obj = sim_data, get_redDim = FALSE, 
                     get_cluster = FALSE, verbose = FALSE, vst = FALSE)


# Abordagens edgeR usando pseudo-bulk -------------------------------------
edgeR_LRT <- DE_edgeR(sce, type = "bulk", eb = FALSE)
table(edgeR_LRT$pvalue <= 0.05, edgeR_LRT$is_de)
mean(edgeR_LRT$pvalue <= 0.05)
mean(edgeR_LRT$pvalue <= 0.01)

edgeR_QL <- DE_edgeR(sce, type = "bulk", eb = TRUE)
mean(edgeR_QL$pvalue <= 0.05)
mean(edgeR_QL$pvalue <= 0.01)


# Abordagem regressão quantilica ------------------------------------------

system.time(qr <- lqr(sce))
mean(qr <= 0.05) # 0.055
mean(qr <= 0.01) # 0.0295


# Abordagem modelo linear misto -------------------------------------------

## Utilizando as contagens log-normalizadas pelo método deconvolution
lmm_lc  <- lmm(sce, name_assay = "logcounts", weight = TRUE)
mean(lmm_lc$pvalue <= 0.05)
mean(lmm_lc$pvalue <= 0.01)
hist(lmm_lc $sigma2_ranef, breaks = 50)
mean(lmm_lc $sigma2_ranef == 0)
icc_lmm_lc <- lmm_lc$sigma2_ranef / (lmm_lc $sigma2_ranef + lmm_lc $sigma2)
range(icc_lmm_lc)

## Utilizando as contagens log-normalizadas pelo método deconvolution e corrigidas
lmm_lcc <- lmm(sce, name_assay = "corrected")
mean(lmm_lcc$pvalue <= 0.05)
mean(lmm_lcc$pvalue <= 0.01)
hist(lmm_lcc$sigma2_ranef, breaks = 50)
mean(lmm_lcc$sigma2_ranef == 0)
icc_lmm_lcc <- lmm_lcc$sigma2_ranef / (lmm_lcc$sigma2_ranef + lmm_lcc$sigma2)
range(icc_lmm_lcc)
