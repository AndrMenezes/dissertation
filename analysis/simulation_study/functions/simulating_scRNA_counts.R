
# Simula dados de scRNA-seq com estrutura hierarquica do indivíduo
# Código é baseado no artigo de Lun e Marioni (2016). 
# Ver https://github.com/MarioniLab/PlateEffects2016
# Os autores simularam uma estrutura de dependecia de plates 
# (lugar onde a célula foi sequenciada)


# estimated_parameters <- readRDS("./simulation_study/data/parms_3_samples.rds")
# range_n_cells <- c(200, 500)
# n_group <- 2
# n_cell_type <- 2
# pct_de <- 0.01
# n_genes <- 2000

simulate_sce <- function(estimated_parameters, n_genes, 
                         range_n_cells = c(200, 500),
                         n_group = 2, var_effect = 1, 
                         pct_de = 0, fc = 2, n_cell_type = 2) {
  
  # Inicializando os parâmetros estimados
  
  # Média e parâmetro de dispersão dos genes
  genes_means <- estimated_parameters$genes_mean
  genes_dispersion <- estimated_parameters$genes_dispersion
  
  # Tamanho da biblioteca das células
  cells_library_size <- estimated_parameters$cells_library_size
  cells_library_size <- cells_library_size / mean(cells_library_size)
  
  # Vetor indicando a origem (do individuo) da célula 
  individuo_origin <- factor(estimated_parameters$sample_id)
  
  # Variância estimada do efeito aleatório. Será utilizado para simular o efeito
  # do indivíduo no gene.
  sigma2 <- mean(estimated_parameters$sigma2_ranef)
  sigma2 <- sigma2 / var_effect
  
  
  # Número de indivíduos e as condições experimentais (grupos por individuos)
  n_individuo <- nlevels(individuo_origin)
  conditions <- rep(LETTERS[1:n_group], each = n_individuo)
  names(conditions) <- rep(levels(individuo_origin), n_group)
  
  n_de <- n_genes * pct_de
  chosen_de <- sample(n_genes, n_de, replace = FALSE)

  # Repete a geração dos dados para simular um efeito de estado/tipo de célula
  lt_counts <- list()
  lt_design <- list()
  for (k in seq_len(n_cell_type)) {
    
    # Gera aleatoriamente a quantidade de célula por grupo e individuo
    n_cells_by_group_individuo <- round(
      runif(length(conditions), range_n_cells[1], range_n_cells[2]))
    names(n_cells_by_group_individuo) <- rep(levels(individuo_origin), n_group)
    
    # Seleciona aleatoriamente as células dos individuos
    cells_lib_size <- sapply(seq_along(n_cells_by_group_individuo), function(i){
      id <- individuo_origin == names(n_cells_by_group_individuo)[i]
      sample(seq_along(cells_library_size[id]), n_cells_by_group_individuo[i], 
             replace = TRUE)
    })
    cells_lib_size   <- unlist(cells_lib_size)
    cells_grouping   <- factor(rep(conditions, n_cells_by_group_individuo))
    cells_individuos <- factor(rep(rep(levels(individuo_origin), n_group), 
                                   n_cells_by_group_individuo))
    
    # Cria a matriz de experimento com informações da simulação
    design <- data.frame(group = cells_grouping, individuo = cells_individuos)
    design$comb <- factor(paste0(design$individuo, "_", design$group))
    order_design <- factor(paste0(names(conditions), "_", conditions), 
                           levels = paste0(names(conditions), "_", conditions) ,
                           ordered = TRUE)
    design$comb <- factor(design$comb, levels = order_design, ordered = TRUE)
    
    # Escolhendo aleatoriamente os parâmetros para os genes
    chosen_genes <- sample(length(genes_means), n_genes, replace=TRUE)
    curent_mean <- genes_means[chosen_genes]
    curent_disp <- genes_dispersion[chosen_genes]
    
    # Simulando efeito do individuo
    genes_individuo_effect <- curent_mean * matrix(
      data = exp(rnorm(n_genes * n_individuo, mean = -sigma2/2, sd = sqrt(sigma2))), 
      nrow = n_genes, ncol = n_individuo
    )
    genes_individuo_effect <- do.call(
      cbind, replicate(n_group, genes_individuo_effect, simplify=FALSE))
    
    # Testando a funçãoe
    # colnames(genes_individuo_effect) <- paste0(names(conditions), "_", conditions) 
    # head(genes_individuo_effect)
    # genes_group_effect <- do.call(cbind, replicate(n_group-1, rnorm(n_genes, sd = 1) * genes_individuo_effect, 
    #                                               simplify=FALSE))
    # genes_effect <- cbind(genes_individuo_effect, genes_group_effect)
    # colnames(genes_effect) <- paste0(names(conditions), "_", conditions)
    
    # Adicionando genes DE
    if (n_de > 0) {
      ingroup <- unique(conditions)[1:(n_group - 1)]
      is_up <- rep(c(TRUE, FALSE), length.out = n_de)
      # fc <- sqrt(fc)
      for (j in 1:length(ingroup)) {
        fc <- j * fc
        m <- which(conditions %in% ingroup[j])
        genes_individuo_effect[chosen_de[is_up],m] <- genes_individuo_effect[chosen_de[is_up],m] * fc
        genes_individuo_effect[chosen_de[is_up],-m] <- genes_individuo_effect[chosen_de[is_up],-m] / fc
        genes_individuo_effect[chosen_de[!is_up],m] <- genes_individuo_effect[chosen_de[!is_up],m] / fc
        genes_individuo_effect[chosen_de[!is_up],-m] <- genes_individuo_effect[chosen_de[!is_up],-m] * fc
      }
    } else {
      chosen_de <- integer(0)
    }
    is_de <- logical(n_genes)
    is_de[chosen_de] <- TRUE
    
    
    # Expandindo e adicionando o efeito do tamanho da biblioteca
    expanded <- as.integer(design$comb)
    cell_means <- genes_individuo_effect[, expanded]
    cell_means <- t(t(cell_means) * cells_lib_size)
    
    # Simulando as contagens para cada célula e gene
    n_cells <- ncol(cell_means)
    counts <- matrix(rnbinom(n_cells * n_genes, mu = cell_means, size = 1 / curent_disp),
                     ncol = n_cells, nrow = n_genes)
    design <- dplyr::arrange(design, comb)
    design$cell_type <- factor(rev(LETTERS)[k])
    
    # Guarda resultados
    lt_counts[[k]] <- counts
    lt_design[[k]] <- design
  }
  
  design <- do.call("rbind", lt_design)
  counts <- do.call("cbind", lt_counts)

  lt_out <- list(counts = counts, design = design, is_de = is_de)
  return(lt_out)
}


