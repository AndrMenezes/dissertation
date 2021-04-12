rm(list = ls())

library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())
wd <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/cap5/simulation_study"
FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}
printmrow <- function(x) cat(cat(x,sep=" & "),"\\\\ \n")

# Cenários da simulação
n_samples <- c(3, 6)
n_group <- c(2, 3)
var_effect <- c(0.1, 1e6)
df_scenario <- expand.grid(n_group = n_group, n_samples = n_samples, 
                           var_effect = var_effect)
df_scenario$scenario <- 1:nrow(df_scenario)
head(df_scenario)

# Resultados da simulação
res_type_I <- read.table("./simulation_study/data/type_I_error/res_type_I_error.txt", 
                         sep = "\t")
names(res_type_I) <- c("method", "scenario", "mc", "alpha", "type_I")
head(res_type_I)

# Resumindo as iterações Monte Carlo
df_error <- res_type_I %>% 
  arrange(method, scenario, mc, alpha) %>% 
  group_by(method, scenario, alpha) %>% 
  summarise(mean_error_rate = mean(type_I),
            sd_error_rate = sd(type_I), .groups = "drop")



# Incluindo parâmetros do cenário
df_error <- df_error %>% 
  left_join(df_scenario, by = "scenario")


df_error %>% 
  filter(alpha == 0.05) %>% 
  filter(method == "lmm")
df_error %>% 
  filter(alpha == 0.05) %>% 
  filter(method == "edgeR_bulk")

# Tabela com todos os resultados
options(OutDec = ",")
tab <- df_error %>% 
  filter(alpha != 0.001) %>% 
  mutate(prt = FF(mean_error_rate),
         var_effect = ifelse(var_effect == 0.1, "Sim", "Não")) %>% 
  select(var_effect, n_samples, n_group, alpha, method, prt) %>% 
  tidyr::pivot_wider(names_from = method, values_from = prt)
tab$var_effect <- c(rep(c("\\multirow{8}{*}{Sim}", ""), c(1, 7)),
                    rep(c("\\multirow{8}{*}{Não}", ""), c(1, 7)))
tab$n_samples <- rep(c(rep(c("\\multirow{4}{*}{3}", ""), c(1, 3)),
                       rep(c("\\multirow{4}{*}{6}", ""), c(1, 3))), 2)
tab$n_group <- rep(c("\\multirow{2}{*}{2}", "", "\\multirow{2}{*}{3}", ""), 4)

ff <- file.path(wd, "simul_error_rate.tex")
file.remove(ff)
sink(ff, append = T)
cat('\\begin{table}[H]', "\n")
cat("\\caption{.}", "\n")
cat("\\onehalfspacing 
     \\centering 
     \\begin{tabular}{ccccccc} \\toprule 
     Efeito de indivíduo & Indivíduos & Grupos & $\\alpha$ & edgeR \\emph{bulk} 
     & edgeR sc & MLM
     \\\\ \\midrule  ")
invisible(apply(tab, 1, printmrow))
cat('\\bottomrule \\end{tabular}', "\n")
cat('\\end{table}', "\n")
sink()

options(OutDec = ".")

# Gráfico resumindo o erro do tipo I estimado para alpha = 0.01
df_error %>%
  filter(alpha == 0.05) %>% 
  mutate(var_effect = ifelse(var_effect == 0.1, "Efeito de indivíduo", 
                             "Sem efeito de indivíduo"),
         n_group = paste0(n_group, " grupos"),
         n_samples = paste0(n_samples, " indivíduos"),
         method = gsub("_", " ", method),
         method = ifelse(method == "lmm", "MLM", method)) %>% 
  ggplot(aes(x = method, y = mean_error_rate, col = factor(var_effect))) +
  facet_grid(n_group~n_samples, scales = "free_y") +
  geom_point(size = 3, shape = 15, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0.05, size = 0.7, col = "dodgerblue") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  colorspace::scale_color_discrete_qualitative() +
  panel_border() +
  labs(x = "Modelo", y = "Erro do tipo I estimado", col = "") +
  theme(legend.position = "top")
ggsave(filename = file.path(wd, "simul_error_rate.pdf"),
       width = 10, height = 7)




