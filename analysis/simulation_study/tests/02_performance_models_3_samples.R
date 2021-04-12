
# Simulação 1: 2 grupos com e var_effect = 10 (sigma2 = 2)
df_1 <- readRDS("./simulation_study/data/test_1_simul_samples_3_nde.rds")
df_1 %>% 
  dplyr::group_by(method, comparison) %>% 
  dplyr::summarise(
    n_fails = sum(is.na(pvalue)),
    alpha_0.05 = mean(pvalue <= 0.05, na.rm=TRUE),
    alpha_0.025 = mean(pvalue <= 0.025, na.rm=TRUE),
    alpha_0.01 = mean(pvalue <= 0.01, na.rm=TRUE)
  )

# Simulação 2: 2 grupos com e var_effect = 1 (sigma2 = 20)
df_2 <- readRDS("./simulation_study/data/test_2_simul_samples_3_nde.rds")
df_2 %>% 
  dplyr::group_by(method, comparison) %>% 
  dplyr::summarise(
    n_fails = sum(is.na(pvalue)),
    alpha_0.05 = mean(pvalue <= 0.05, na.rm=TRUE),
    alpha_0.025 = mean(pvalue <= 0.025, na.rm=TRUE),
    alpha_0.01 = mean(pvalue <= 0.01, na.rm=TRUE)
  )


# Simulação 3: 3 grupos com e var_effect = 10 (sigma2 = 2)
df_3 <- readRDS("./simulation_study/data/test_3_simul_samples_3_nde.rds")
df_3 %>% 
  dplyr::group_by(method) %>% 
  dplyr::summarise(
    n_fails = sum(is.na(pvalue)),
    alpha_0.05 = mean(pvalue <= 0.05, na.rm=TRUE),
    alpha_0.025 = mean(pvalue <= 0.025, na.rm=TRUE),
    alpha_0.01 = mean(pvalue <= 0.01, na.rm=TRUE)
  )
