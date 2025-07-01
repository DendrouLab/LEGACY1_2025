options(bitmapType='cairo-png')

library(dplyr)
library(SingleR)
library(SingleCellExperiment)
library(Seurat)
library(qs)
library(dittoSeq)
library(BiocParallel)

# Directories
legacy_dat_dir <- "sce.qs"
ellebedy_dat_dir <- "ellebedy.qs"
out_dir <- ""

# Import and convert to sce
dat_tcell <- qread(ellebedy_dat_dir)

dat_tcell$ellebedy_annotations <- dat_tcell$Tfh_type_fine
dat_tcell$ellebedy_annotations <- as.character(dat_tcell$ellebedy_annotations)
dat_tcell$ellebedy_annotations[is.na(dat_tcell$ellebedy_annotations)] <- dat_tcell$ellebedy_seurat_clusters[is.na(dat_tcell$ellebedy_annotations)]
dat_tcell$ellebedy_annotations <- as.factor(dat_tcell$ellebedy_annotations)

message("Loading legacy data")
legacy_dat <- qread(legacy_dat_dir)

# By single-cell (ellebedy = reference)
message("Starting ref mapping")
message("Wilcox, lgocounts, clusters = NONE")
pred.ellebedy1 <- SingleR(test = legacy_dat, ref = dat_tcell, de.method = "wilcox", assay.type.test="logcounts", assay.type.ref="logcounts", labels = dat_tcell$ellebedy_annotations,BPPARAM=MulticoreParam(8))
 
message("Done ref mapping. Saving results")
qsave(pred.ellebedy1, "LEGACYdata_refEllebedySingleCells.qs")

# By single-cell (legacy = reference)
pred.ellebedy2 <- SingleR(test = dat_tcell, ref = legacy_dat, de.method = "wilcox", assay.type.test="logcounts", assay.type.ref="logcounts", labels = legacy_dat$celltype,BPPARAM=MulticoreParam(8))
 
message("Done ref mapping. Saving results")
qsave(pred.ellebedy2, "Ellebedydata_refLEGACYSingleCells.qs")

# Plotting LEGACY object with Ellebedy ref mapped annotations
legacy_sce$barcode <- rownames(colData(legacy_sce))
res_singlecell <- data.frame("cellbarcodes" = row.names(pred.singlecell), "labels" = pred.singlecell$labels) 
legacy_sce$singleR_singlecellRes <- unlist(lapply(legacy_sce$barcode, function(x) res_singlecell$labels[match(x, res_singlecell$cellbarcodes)]))

DA_pl <- dittoBarPlot(legacy_sce, "singleR_singlecellRes", scale = "percent", group.by = "visitid_vaxtype",  ylab = "Prop of Cells")
DA_pl
ggsave(filename=file.path(out_dir, paste0("StackedBarPlot_perVisitIDVaxType_Ellebedyref_LEGACYdata.pdf")), width = 5, height = 5)


# Plotting Ellebedy object with LEGACY ref mapped annotations
res_singlecell <- data.frame("cellbarcodes" = row.names(pred.singlecell), "labels" = pred.singlecell$labels) 

DA_p2 <- dittoBarPlot(dat_tcell, "singleR_singlecellRes", scale = "percent", group.by = "day", split.by = "donor", ylab = "Prop of Cells")
DA_p2
ggsave(filename=file.path(out_dir, paste0("StackedBarPlot_perSample_LEGACYref_Ellebedydata.pdf")), width = 7, height = 5)
```