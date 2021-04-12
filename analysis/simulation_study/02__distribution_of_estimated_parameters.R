rm(list = ls())
library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid())

# Análise dos parâmetros estimados nos dois cenários

# path_fig <- "./simulation_study/figures/"
path_fig <- "C:/Users/User/Dropbox/Unicamp/dissertacao/notas/figuras/cap5/simulation_study"

# Cenário com 3 indivíduos ------------------------------------------------

parms_3_samples <- readRDS("./simulation_study/data/parms_3_samples.rds")
str(parms_3_samples)
table(parms_3_samples$sample_id)

# Tamanho da biblioteca

tb_lib <- dplyr::tibble(lib_size = parms_3_samples$cells_library_size,
                        size_factor = parms_3_samples$cells_size_factor,
                        sample_id = parms_3_samples$sample_id) %>% 
  dplyr::mutate(lib_size = lib_size / mean(lib_size))

ggplot(tb_lib, aes(x=log2(lib_size))) +
  geom_density()
ggplot(tb_lib, aes(x=log2(lib_size), fill = sample_id)) +
  geom_density(alpha = 0.6) +
  labs(x = expression(Log[2]~"Tamanho da biblioteca"), y = "Densidade", fill = "")
ggsave(filename = "library_size_3_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Fator de escala
tb_lib %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(mean = mean(lib_size))

tapply(tb_lib$sample_id, tb_lib$lib_size, mean)
mean(tb_lib$size_factor)

ggplot(tb_lib, aes(x=log2(size_factor))) +
  geom_density()
ggplot(tb_lib, aes(x=log2(size_factor), fill = sample_id)) +
  geom_density(alpha = 0.6) +
  labs(x = expression(Log[2]~"Fator de escala"), y = "Densidade", fill = "")
ggsave(filename = "size_factor_3_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Distribuição dos parâmetros dos genes

tb_genes <- dplyr::tibble(average_BN = parms_3_samples$genes_mean,
                          dispersion_BN = parms_3_samples$genes_dispersion,
                          variance_ranef = parms_3_samples$sigma2_ranef)

# Médias estimadas para cada gene usando a BN
ggplot(tb_genes, aes(x=log2(average_BN))) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  labs(x = expression(Log[2]~"Contagem média estimada"), 
       y = "Número de genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
ggsave(filename = "average_count_3_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Parâmetro de disepersão estimados para cada gene usando a BN
ggplot(tb_genes, aes(x=log2(dispersion_BN))) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  labs(x = expression(Log[2]~"Dispersão da binomial negativa"),
       y = "Número de genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
ggsave(filename = "dispersion_3_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Dispersão versus médias estimadas para cada gene usando a BN
tb_genes %>% 
  dplyr::filter(dispersion_BN > 0.001) %>% 
  ggplot(aes(x = log2(average_BN), y = log2(dispersion_BN))) +
  stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = 1),
                 contour = FALSE, show.legend = FALSE) + 
  geom_point(size = 0.5, alpha = 0.4, color = "grey20") +
  stat_density2d(geom = "tile", 
                 aes(fill = ..density..^0.25,
                     alpha = ifelse(..density..^0.25<0.4, 0, 1)), 
                 contour = FALSE, show.legend = FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
  labs(x = expression(Log[2]~"Contagem média estimada"), 
       y = expression(Log[2]~"Dispersão da binomial negativa")) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  theme_cowplot()
ggsave(filename = "smooth_3_samples.pdf", 
       path = path_fig, width = 8, height = 6)


# Variância do efeito aleatório
tb_genes %>% 
#  dplyr::filter(variance_ranef != 0) %>% 
  ggplot(aes(x=variance_ranef)) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  geom_vline(xintercept = mean(tb_genes$variance_ranef), col = "red") +
  labs(x = expression(sigma^2~" estimado"),
       y = "Número de genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
ggsave(filename = "sigma2_3_samples.pdf", path = path_fig, width = 8, 
       height = 6)

# Cenário com 6 indivíduos ------------------------------------------------

parms_6_samples <- readRDS("./simulation_study/data/parms_6_samples.rds")
str(parms_6_samples)
table(parms_6_samples$sample_id)

# Tamanho da biblioteca

tb_lib <- dplyr::tibble(lib_size = parms_6_samples$cells_library_size,
                        size_factor = parms_6_samples$cells_size_factor,
                        sample_id = parms_6_samples$sample_id) %>% 
  dplyr::mutate(lib_size = lib_size / mean(lib_size))

ggplot(tb_lib, aes(x=log2(lib_size))) +
  geom_density()
ggplot(tb_lib, aes(x=log2(lib_size), fill = sample_id)) +
  geom_density(alpha = 0.6) +
  labs(x = expression(Log[2]~"Tamanho da biblioteca"), y = "Densidade", fill = "")
ggsave(filename = "library_size_6_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Fator de escala
ggplot(tb_lib, aes(x=log2(size_factor))) +
  geom_density()
ggplot(tb_lib, aes(x=log2(size_factor), fill = sample_id)) +
  geom_density(alpha = 0.6) +
  labs(x = expression(Log[2]~"Fator de escala"), y = "Densidade", fill = "")
ggsave(filename = "size_factor_6_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Distribuição dos parâmetros dos genes

tb_genes <- dplyr::tibble(average_BN = parms_6_samples$genes_mean,
                          dispersion_BN = parms_6_samples$genes_dispersion,
                          variance_ranef = parms_6_samples$sigma2_ranef)

# Médias estimadas para cada gene usando a BN
ggplot(tb_genes, aes(x=log2(average_BN))) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  labs(x = expression(Log[2]~"Contagem média estimada"), 
       y = "Número de genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
ggsave(filename = "average_count_6_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Parâmetro de disepersão estimados para cada gene usando a BN
ggplot(tb_genes, aes(x=log2(dispersion_BN))) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  labs(x = expression(Log[2]~"Dispersão da binomial negativa"),
       y = "Número de genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
ggsave(filename = "dispersion_6_samples.pdf", 
       path = path_fig, width = 8, height = 6)

# Dispersão versus médias estimadas para cada gene usando a BN
tb_genes %>% 
  dplyr::filter(dispersion_BN > 0.001) %>% 
  ggplot(aes(x = log2(average_BN), y = log2(dispersion_BN))) +
  stat_density2d(geom = "tile", aes(fill = ..density..^0.25, alpha = 1),
                 contour = FALSE, show.legend = FALSE) + 
  geom_point(size = 0.5, alpha = 0.4, color = "grey20") +
  stat_density2d(geom = "tile", 
                 aes(fill = ..density..^0.25,
                     alpha = ifelse(..density..^0.25<0.4, 0, 1)), 
                 contour = FALSE, show.legend = FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) +
  labs(x = expression(Log[2]~"Contagem média estimada"), 
       y = expression(Log[2]~"Dispersão da binomial negativa")) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  theme_cowplot()
ggsave(filename = "smooth_6_samples.pdf", 
       path = path_fig, width = 8, height = 6)


# Variância do efeito aleatório
tb_genes %>% 
#  dplyr::filter(variance_ranef != 0) %>% 
  ggplot(aes(x=variance_ranef)) +
  geom_histogram(bins = 30, color = "black", fill = "grey") +
  geom_vline(xintercept = mean(tb_genes$variance_ranef), col = "red") +
  labs(x = expression(sigma^2~" estimado"),
       y = "Número de genes") +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
ggsave(filename = "sigma2_6_samples.pdf", path = path_fig, width = 8,
       height = 6)
mean(parms_3_samples$sigma2_ranef)
mean(parms_6_samples$sigma2_ranef)
