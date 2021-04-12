rm(list = ls())

setwd("./covid/")

# Pkgs --------------------------------------------------------------------

library(magrittr)
## Read .h5 objects
suppressMessages(library(DropletUtils))
## QC
suppressMessages(library(scater))
## Save h5 in disk
suppressMessages(library(HDF5Array))
## To find the mitocondrial genes
suppressMessages(library(EnsDb.Hsapiens.v86))

# Metafile and arqs -------------------------------------------------------
meta <- read.delim("./data/meta.txt")
ff <- list.files(path = "./data", pattern = ".h5")[-13]

# Apply quality control ---------------------------------------------------

get_qc <- function(sce) {
  
  location <- mapIds(
    EnsDb.Hsapiens.v86, 
    keys = rowData(sce)$ID,
    column = "SEQNAME", 
    keytype = "GENEID"
  )
  
  stats <- perCellQCMetrics(
    x = sce,
    subsets = list(mito = which(location == "MT"))
  )
  
  ## Apply filter according to the Authors
  discard <- stats$sum < 1000 | 
    stats$detected < 200 | stats$detected > 6000 | 
    stats$subsets_mito_percent > 10 
  
  sce <- sce[, !discard]
  
  sce
}


# Apply QC and return data ------------------------------------------------

get_data <- function(onefile) {
  
  tmp <- read10xCounts(paste0("./data/", onefile), col.names = TRUE)
  
  # ?quickResaveHDF5SummarizedExperiment(tmp)
  
  # saveHDF5SummarizedExperiment(tmp, dir = "h5", prefix = "")
  
  rownames(tmp) <- uniquifyFeatureNames(
    rowData(tmp)$ID, 
    rowData(tmp)$Symbol
  )
  metadata(tmp) <- list()
  sce <- get_qc(tmp)

  ## Including colData information
  filt <- gsub(pattern = ".h5", replacement = "", x = onefile)
  id <- dplyr::filter(meta, sample == filt)
  
  sce$sample  <- id$sample_new
  sce$disease <- id$disease
  sce$group   <- id$group
  
  sce
}



# Get data in a list ------------------------------------------------------

l_sce <- lapply(ff, get_data)
names(l_sce) <- sapply(l_sce, function (x) x$sample[1])

# Extract and filter genes that are zero for all samples ------------------
zero_genes <- lapply(l_sce, {function(x)
  names(which(rowSums(counts(x)) == 0))
})
zero_genes <- unique(unlist(zero_genes))
length(zero_genes)

l_sce <- lapply(l_sce, function(x) {
  x[-which(rownames(x) %in% zero_genes), ]
})
# Total de genes em cada paciente (amostra)
sapply(l_sce, nrow)
# Total de células (observações) em cada paciente (amostra)
sapply(l_sce, ncol)


# Save data using .h5 format ----------------------------------------------

for (i in seq_along(l_sce)) {
  dir_to_save <- paste0("./data/", l_sce[[i]]$sample[1])
  saveHDF5SummarizedExperiment (l_sce[[i]], dir = dir_to_save)
  cat(i, "\n")
}

# saveRDS(l_sce, file = "./data/01_list_sces.rds")
