library(edgeR)

#' @description Realizar a análise DE utilizando edgeR
#'
#' @param x objeto SingleCellExperiment contendo o coldata e as expressões
#' @param label vetor especificando as populações das células
#' @param design formula especificando o modelo
#' @param contrast contrastes especificando as comparações
#' @param filter_gene logico indicando se realiza a filtragem dos genes segundo
#' a abordagem da função `edgeR::filterByExpr`
#' @param assay.type indica o nome do assay com as contagens

DE_edgeR <- function(x, label, design, contrast = NULL, filter_gene = TRUE, 
                     eb = FALSE, assay.type = 1) {
  
  expr_mat <- assay(x, i = assay.type)
  col_data <- colData(x)

  counter <- 1L  
  lt_labels <- list()
  for (l in sort(unique(label))) {
    chosen <- l == label
    cur_expr <- expr_mat[, chosen, drop = FALSE]
    cur_data <- col_data[chosen,, drop = FALSE]
    y <- DGEList(cur_expr, samples = as.data.frame(cur_data))
    cur_design <- model.matrix(design, data = cur_data)
    
    lt_labels[[counter]] <- pairwise_edgeR(y, cur_design, contrast, 
                                           filter_gene, eb)
    counter <- counter + 1L
  }
  names(lt_labels) <- unique(label)
  
  return(lt_labels)
  
}

#' @description Realizar as comparações conforme o contraste especificando 
#' modelo BN do pacote `edgeR`
#'
#' @param y objeto da classe DGEList contendo o coldata e as expressões gênicas
#' @param cur_design formula especificando o modelo
#' @param contrast contrastes especificando as comparações
#' @param filter_gene logico indicando se realiza a filtragem dos genes segundo
#' a abordagem da função `edgeR::filterByExpr`
#'

pairwise_edgeR <- function(y, cur_design, contrast, filter_gene = TRUE, 
                           eb = FALSE) {
  ngenes <- nrow(y)
  genes_names <- rownames(y)
  
  # Se contraste nulo cria contrastes para comparaçao múltipla
  if (is.null(contrast)) {
    lev <- factor(1:ncol(cur_design))
    lfm <- diag(length(lev))
    cbn <- utils::combn(seq_along(lev), 2)
    contrast <- t(lfm[cbn[1, ], ] - lfm[cbn[2, ], ])
  }
  
  # Filtra genes que tem contagens suficientes para a análise estatística.
  if (filter_gene) {
    gkeep <- filterByExpr(y, design = cur_design)
    y <- y[gkeep, ]
  } else {
    gkeep <- c()
  }
  # Normaliza e 
  y <- calcNormFactors(y)
  y <- scaleOffset(y, getOffset(y))
  
  # Estima parâmetro de dispersão e o modelo via QL
  y <- estimateDisp(y, cur_design)
  if (eb) {
    fit <- glmQLFit(y, cur_design, robust = TRUE)
  } else {
    fit <- glmFit(y, cur_design, dispersion = y$tagwise.dispersion)
  }
  
  # Realiza as comparações múltiplas
  lt_tabs <- list()
  for (j in 1:ncol(contrast)) {
    if (eb) {
      res <- glmQLFTest(fit, contrast = contrast[, j])
    } else {
      res <- glmLRT(fit, contrast = contrast[, j])
    }
    tab <- topTags(res, n=Inf, sort.by="none")
    expander <- match(seq_len(ngenes), which(gkeep))
    df <- DataFrame(tab$table[expander,,drop=FALSE])[, -5]
    df$contrast <- gsub("group", "", tab$comparison)
    rownames(df) <- genes_names
    lt_tabs[[j]] <- df
  }
  # Exporta data.frame
  tabs <- do.call(rbind, lt_tabs)
  tabs$gene <- rownames(tabs)
  out <- list(tabs = tabs, fit = y)
  return(out)
}




