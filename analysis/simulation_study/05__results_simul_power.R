rm(list = ls())
wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/cap5/simulation_study"

# Pacotes -----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())


# Funções úteis -----------------------------------------------------------
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")


# Cenários da simulação ---------------------------------------------------

n_samples <- 3
n_group <- 2
var_effect <- c(0.5, 1e6)
pct_de <- c(0.05, 0.1, 0.2, 0.5)
fc <- c(2, 4)
n_genes <- 2000

df_scenario <- expand.grid(n_genes = n_genes, n_group = n_group,
                           n_samples = n_samples, var_effect = var_effect,
                           pct_de = pct_de, fc = fc)
df_scenario$scenario <- 1:nrow(df_scenario)
head(df_scenario)


# Resultados com FDR e AUC ------------------------------------------------

res_power <- read.table("./simulation_study/data/power/res_power_1cell.txt", 
                        sep = "\t")
names(res_power) <- c("method", "scenario", "mc", "fdr", "auc")
res_power <- res_power[res_power$method != "lmm", ]
head(res_power)


res_power_mm <- read.table("./simulation_study/data/power/res_power_mm.txt", 
                        sep = "\t")
names(res_power_mm) <- c("method", "scenario", "mc", "fdr", "auc", 
                         "sigma2_ranef", "sigma2")
head(res_power_mm)

res_power <- rbind(res_power, res_power_mm[, -c(6, 7)])


# Resumindo as iterações Monte Carlo
df_power <- res_power %>% 
  group_by(method, scenario) %>% 
  summarise(fdr = mean(fdr), auc = mean(auc), .groups = "drop") %>% 
  arrange(fdr)

# Incluindo parâmetros do cenário
df_power <- df_power %>% 
  left_join(df_scenario, by = "scenario") %>% 
  arrange(scenario)

df_power %>% 
  filter(method == "lmm")
df_power %>% 
  filter(method == "edgeR_bulk")
df_power %>% 
  filter(method == "edgeR_sc")


# Gráfico resumindo o fdr estimado para alpha = 0.05
df_power %>%
  mutate(var_effect = ifelse(var_effect == 0.5, "Efeito de indivíduo", 
                             "Sem efeito de indivíduo"),
         pct_de = paste0(pct_de * 10, "% de genes DE"),
         fc = paste0("fold change de ", fc),
         method = gsub("_", " ", method),
         method = ifelse(method == "lmm", "MLM", method)) %>% 
  ggplot(aes(x = method, y = fdr, col = factor(var_effect))) +
  facet_grid(fc ~ pct_de) +
  geom_point(size = 3, shape = 15, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.05, col = "dodgerblue") +
  scale_y_continuous(breaks = c(0, 0.05, 0.15, 0.25, 0.35, 0.42)) +
  colorspace::scale_color_discrete_qualitative() +
  panel_border() +
  labs(x = "Modelo", y = "FDR estimado", col = "") +
  theme(legend.position = "top")
ggsave(filename = file.path(wd, "simul_fdr.pdf"), width = 12, height = 7)

# Gráfico resumindo a AUC estimada 
df_power %>%
  mutate(var_effect = ifelse(var_effect == 0.5, "Efeito de indivíduo", 
                             "Sem efeito de indivíduo"),
         pct_de = paste0(pct_de * 10, "% de genes DE"),
         fc = paste0("fold change de ", fc),
         method = gsub("_", " ", method),
         method = ifelse(method == "lmm", "MLM", method)) %>% 
  ggplot(aes(x = method, y = auc, col = factor(var_effect))) +
  facet_grid(fc ~ pct_de) +
  geom_point(size = 3, shape = 15, position = position_dodge(width = 0.4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  colorspace::scale_color_discrete_qualitative() +
  panel_border() +
  labs(x = "Modelo", y = "AUC estimada", col = "") +
  theme(legend.position = "top")
ggsave(filename = file.path(wd, "simul_auc.pdf"), width = 12, height = 7)


# Estimativa da variância do efeito aleatório
res_power_mm %>% 
  group_by(scenario) %>% 
  summarise(sigma2_ranef = mean(sigma2_ranef), sigma2 = mean(sigma2)) %>% 
  left_join(df_scenario, by = "scenario") %>% 
  arrange(scenario)


# Dados da curva ROC ------------------------------------------------------

wd_res <- "./simulation_study/data/power"

lt_rocs <- lapply(1:nrow(df_scenario), function(i) {
  roc_res <- readRDS(file.path(wd_res, paste0("roc_res_", i, ".rds")))
  roc_res <- roc_res[roc_res$method != "lmm", ]
  roc_res_mm <- readRDS(file.path(wd_res, paste0("roc_res_mm_", i, ".rds")))
  roc_res <- rbind(roc_res, roc_res_mm)
  roc_res$scenario <- i
  return(roc_res)
})
df_roc <- do.call("rbind", lt_rocs) %>% 
  left_join(df_scenario, by = "scenario") %>% 
  arrange(scenario)


df_roc %>%
  filter(fc == 2) %>% 
  mutate(var_effect = ifelse(var_effect == 0.5, "Efeito de indivíduo", 
                             "Sem efeito de indivíduo"),
         pct_de = paste0(pct_de * 10, "% de genes DE"),
         method = gsub("_", " ", method),
         method = ifelse(method == "lmm", "MLM", method)) %>% 
  ggplot(aes(x = esp, y = sens, col = factor(method))) +
  facet_grid(var_effect ~ pct_de) +
  geom_line(size = 0.8) +
  scale_y_continuous(breaks = scales::pretty_breaks(6),
                     expand = expansion(mult = c(0, 0.1))) +
  colorspace::scale_color_discrete_qualitative() +
  panel_border() +
  labs(x = "1 - Especificidade", y = "Sensibilidade", 
       col = "") +
  theme(legend.position = "top")
ggsave(filename = file.path(wd, "simul_roc.pdf"), width = 12, height = 7)


