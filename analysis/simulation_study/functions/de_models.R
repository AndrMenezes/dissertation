library(lme4)
library(emmeans)
library(edgeR)
library(quantreg)

emm_options(lmerTest.limit = 1e8)

# Funções para ajustar modelos para análise DE

# edgeR -------------------------------------------------------------------

DE_edgeR <- function(x, type, assay.type = 1, eb = FALSE) {
  
  if (type == "bulk") {
    x <- scater::aggregateAcrossCells(
      x = x, ids = colData(x)[, c("individuo", "group", "cell_type")])
    colData(x) <- colData(x)[, -c(1:5)]
  }
  
  label <- x$cell_type
  df <- ~0 + group
  
  n_group <- length(unique(x$group))
  if (n_group == 2) {
    contrast <- matrix(c(1, -1), ncol = 1)
  } else {
    lev <- factor(1:n_group)
    lfm <- diag(length(lev))
    cbn <- utils::combn(seq_along(lev), 2)
    contrast <- t(lfm[cbn[1, ], ] - lfm[cbn[2, ], ])
  }
  
  expr_mat <- assay(x, i = assay.type)
  col_data <- colData(x)
  
  counter <- 1L  
  lt_labels <- list()
  for (l in sort(unique(label))) {
    chosen <- l == label
    cur_expr <- expr_mat[, chosen, drop = FALSE]
    cur_data <- col_data[chosen,, drop = FALSE]
    y <- DGEList(cur_expr, samples = as.data.frame(cur_data), 
                 genes = rowData(x))
    cur_df <- model.matrix(df, data = cur_data)
    lt_labels[[counter]] <- pairwise_edgeR(y, cur_df, contrast, eb)
    
    # cat(counter, '\n')
    counter <- counter + 1L
  }
  names(lt_labels) <- unique(label)
  pvalue <- do.call("c",lapply(lt_labels, "[[", "pvalue"))
  is_de <- do.call("c",lapply(lt_labels, "[[", "is_de"))
  output <- list(pvalue = pvalue, is_de = is_de)
  return(output)
}

pairwise_edgeR <- function(y, cur_df, contrast, eb) {
  ngenes <- nrow(y)
  genes_names <- rownames(y)
  
  # Normaliza e escala offset
  y <- calcNormFactors(y)
  y <- scaleOffset(y, getOffset(y))
  # Estima parâmetro de dispersão e o modelo via QL
  y <- estimateDisp(y, cur_df)

  if (eb) {
    # Estima o modelo NB via Quasi-verossimilhança com empircal bayes
    fit <- glmQLFit(y, cur_df, robust = TRUE)
    fun_test <- get("glmQLFTest")
  } else {
    # Estima o modelo NB via verossimilhança
    fit <- glmFit(y, cur_df, dispersion = y$tagwise.dispersion)
    fun_test <- get("glmLRT")
  }
  
  
  # Realiza as comparações múltiplas
  lt_pvalues <- list()
  lt_genes <- list()
  for (j in 1:ncol(contrast)) {
    res <- fun_test(fit, contrast = contrast[, j])
    lt_pvalues[[j]] <- res$table$PValue
    lt_genes[[j]] <- res$genes$is_de
  }
  # Exporta data.frame
  output <- list(pvalue = unname(do.call("c", lt_pvalues)), 
                 is_de = unname(do.call("c", lt_genes)))

  return(output)
}

# Modelo misto ------------------------------------------------------------
lmm <- function(x, name_assay = "logcounts", weight = TRUE) {
  
  expr_matrix <- assay(x, name_assay)
  df <- data.frame(group = x$group, individuo = x$individuo,
                       cell_type = x$cell_type, logcounts = 0)
  if (weight) {
    wi <- voom(counts(x), design = model.matrix(~group, df))$weights
  } else {
    wi <- matrix(rep(1, ncol(x) * nrow(x)), ncol = ncol(x), nrow = nrow(x))
  }
  
  lt_results <- lapply(1:nrow(expr_matrix), function(g) {
    df$logcounts <- expr_matrix[g, ]
    suppressMessages(tryCatch(
      fit <- lmer(formula = logcounts ~ group + (1 | individuo), data = df,
                  weights = wi[g, ]),
      error = function(e) NULL))
    if (!is.null(fit)) {
      emm <- pairs(emmeans(fit, specs = "group", lmer.df = "asymptotic"), 
                   adjust = "none")
      vc <- as.data.frame(VarCorr(fit))
      output <- data.frame(pvalue = summary(emm)$p.value)
      output$sigma2_ranef <- vc$vcov[1]
      output$sigma2 <- vc$vcov[2]
      output$gene <- rep(rownames(x)[g], nrow(output))
      output$is_de <- rowData(x)$is_de[g]
      cat(g, "OK! \n")
    }
    else {
      output <- NULL
      cat(g, "Error! \n")
    }
    return(output)
  })
  
  out <- do.call("rbind", lt_results)
  
  
  return(out)
}

glmm <- function(x, name_assay = "counts", family = poisson(link = "log")) {
  
  expr_matrix <- assay(x, name_assay)
  df <- data.frame(group = x$group, individuo = x$individuo,
                   cell_type = x$cell_type, sf = x$sizeFactor,
                   counts = 0)
  
  ff <- counts ~ group + (1 | individuo) + offset(log(sf))

  lt_results <- lapply(1:nrow(expr_matrix), function(g) {
    
    df$counts <- expr_matrix[g, ]
    suppressMessages(tryCatch(fit <- glmer(formula = ff, data = df, 
                                           family = family), 
             error = function(e) NULL))
    if (!is.null(fit)) {
      emm <- pairs(emmeans(fit, specs = "group", lmer.df = "satterthwaite"), 
                   adjust = "none")
      vc <- as.data.frame(VarCorr(fit))
      output <- data.frame(pvalue = summary(emm)$p.value)
      output$sigma2_ranef <- vc$vcov[1]
      output$gene <- rep(rownames(x)[g], nrow(output))
      output$is_de <- rowData(x)$is_de[g]
      cat(g, "OK! \n")
    }
    else {
      output <- NULL
      cat(g, "Error! \n")
    }
    return(output)
  })
  
  out <- do.call("rbind", lt_results)
  
  
  return(out)
}


# Fit all models ----------------------------------------------------------

fit_models <- function(sce, model) {
  
  out <- switch (model,
    "lmm" = lmm(x = sce, name_assay = "logcounts"),
    "edgeR_sc" = DE_edgeR(x = sce, type = "sc"),
    "edgeR_bulk" = DE_edgeR(x = sce, type = "bulk")
  )
  
  return(out)  
}
