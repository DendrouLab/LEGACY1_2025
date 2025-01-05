# scRepertoire

# Libraries
library('scRepertoire') #R 4.3.3, version 2.0.0
#BiocManager::install("scRepertoire", version = "3.18")
library(ggplot2)
library(qs)
library(dplyr)
library(randomcoloR)
library(SingleCellExperiment)

source("scRepertoire_utils.R")

# Functions
loadContigs2 <- function (input, format = "10X") {
  if (inherits(x = input, what = "character")) {
    format.list <- list(WAT3R = "barcode_results.csv", `10X` = "filtered_contig_annotations.csv", 
                        AIRR = "airr_rearrangement.tsv", Immcantation = "_data.tsv", 
                        MiXCR = "clones.tsv", JSON = ".json", TRUST4 = "barcode_report.tsv", 
                        BD = "Contigs_AIRR.tsv", Omniscope = c("_OSB.csv", 
                                                               "_OST.csv"))
    file.pattern <- format.list[[format]]
    #input <- "C:/Users/jacquelines/OneDrive - Nexus365/Science/Projects/LEGACY/Analysis_LEGACY/005_year1multi_analysis/immcantation/Immcantation/"
    #file.pattern <- "_data.tsv"
    contig.files <- list.files(input, paste0(file.pattern, 
                                             collapse = "|"), recursive = TRUE, full.names = TRUE)
    if (format %in% c("10X", "WAT3R", "Omniscope")) {
      df <- lapply(contig.files, read.csv)
    }
    else if (format %in% c("json")) {
      df <- lapply(contig.files, function(x) {
        tmp <- as.data.frame(fromJSON(x))
      })
    }
    else {
      df <- lapply(contig.files, read.delim)
    }
  }
  else if (inherits(x = input, what = "list") | inherits(x = input, 
                                                         what = "data.frame")) {
    df <- .checkList(input)
  }
  loadFunc <- switch(format, `10X` = .parse10x, AIRR = .parseAIRR, 
                     JSON = .parseJSON, MiXCR = .parseMiXCR, TRUST4 = .parseTRUST4, 
                     BD = .parseBD, WAT3R = .parseWAT3R, Omniscope = .parseOmniscope, 
                     Immcantation = .parseImmcantation, stop("Invalid format provided"))
  df <- loadFunc(df)
  return(df)
}

.parseImmcantation<- function(df) {
  for (i in seq_along(df)) {
    df[[i]][df[[i]] == ""] <- NA
    df[[i]] <- as.data.frame(df[[i]])
    df[[i]] <- df[[i]][,c("sequence_id", "locus", "consensus_count",  "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa"),]
    colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
    df[[i]]$barcode <- stringr::str_split(df[[i]][,"barcode"], "_", simplify = TRUE)[,1]
  }
  return(df)
}

.theCall <- function(df, x, check.df = TRUE) {
  x <- .convertClonecall(x)
  if(check.df) {
    if(inherits(df, "list") & !any(colnames(df[[1]]) %in% x)) {
      stop("Check the clonal variable (cloneCall) being used in the function, it does not appear in the data provided.")
    } else if (inherits(df, "data.frame") & !any(colnames(df) %in% x)) {
      stop("Check the clonal variable (cloneCall) being used in the function, it does not appear in the data provided.")
    }
  }
  return(x)
}

# helper for .theCall
.convertClonecall <- function(x) {
  
  clonecall_dictionary <- hash::hash(
    "gene" = "CTgene",
    "genes" = "CTgene",
    "ctgene" = "CTgene",
    "ctstrict" = "CTstrict",
    "nt" = "CTnt",
    "nucleotide" = "CTnt",
    "nucleotides" = "CTnt",
    "ctnt" = "CTnt",
    "aa" = "CTaa",
    "amino" = "CTaa",
    "ctaa" = "CTaa",
    "gene+nt" = "CTstrict",
    "strict" = "CTstrict",
    "ctstrict" = "CTstrict"
  )
  
  x <- tolower(x)
  
  if (!is.null(clonecall_dictionary[[x]])) {
    return(clonecall_dictionary[[x]])
  }
  else {
    warning("A custom variable ", x, " will be used to call clones")
    return(x)
  }
}

# Import TCR contigs
in_dir <- "Immcantation/"
contig_list <- loadContigs2(input = in_dir, format = "Immcantation")

sample_names <- list.files(in_dir, paste0("_data.tsv", collapse = "|"), recursive = TRUE, full.names = TRUE)
sample_names <- gsub("", "", sample_names)
sample_names <- gsub("_tr_db-pass_cleaned_gexonly_data.tsv", "", sample_names)

combined.TCR <- combineTCR(contig_list, 
                           samples = sample_names,
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

## Add metadata
metadat <- data.frame(sample_names)
metadat$donor_id <- gsub("-.*", "", metadat$sample_names)
metadat$sample_type <- gsub(".*-", "", metadat$sample_names)
metadat$visit_id <- paste0("V", substring(metadat$sample_type, 1, 1))
metadat$sample_side <- substring(metadat$sample_type, nchar(metadat$sample_type), nchar(metadat$sample_type))

clinical_metadat <- read.csv("LEGACYmetadata.csv")
columns_keep <- c("donor_id", "days_V4afterV3", "Vaccine_arm")
clinical_metadat <- clinical_metadat[columns_keep]
clinical_metadat$Vaccine_arm <- gsub("H", "", clinical_metadat$Vaccine_arm)

metadat <- merge(metadat, clinical_metadat, by = "donor_id", all.x = TRUE)

metadat <- transform(metadat, vaccine_type = ifelse(Vaccine_arm==sample_side, "Ipsilateral", "Contralateral"))
metadat$visitid_vaxtype <- paste0(metadat$visit_id, "_", metadat$vaccine_type)

metadat <- metadat[match(sample_names, metadat$sample_names),]

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "visit_id", 
                            variables = metadat$visit_id)

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "vaccine_type", 
                            variables = metadat$vaccine_type)

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "visitid_vaxtype", 
                            variables = metadat$visitid_vaxtype)

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "days_V4afterV3", 
                            variables = metadat$days_V4afterV3)

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "donor_id", 
                            variables = metadat$donor_id)

cloneabundance_dat <- clonalAbundance(combined.TCR, 
                                      cloneCall = "strict", 
                                      scale = FALSE, 
                                      exportTable = TRUE)

cloneabundance_summary <- cloneabundance_dat %>% group_by(values) %>% summarise(unique_clonotype_count = n_distinct(CTstrict))
cloneabundance_summary$visitid <- gsub(".*-", "", cloneabundance_summary$values)


filtercutoff <- 250
toremove <- subset(cloneabundance_summary, unique_clonotype_count < filtercutoff, select = values)

## Remove samples that have too few clonotypes
combined.TCR <- combined.TCR[names(combined.TCR) %in% toremove$values == FALSE]

#### Importantly, the major requirement for the attachment is matching contig cell barcodes and barcodes in the row names of the meta data of the Seurat or Single-Cell Experiment object. If these do not match, the attachment will fail. Based on ease, we suggest making changes to the single-cell object barcodes

# Adding in SCE object
sce_filename <- file.path("sce_object.qs")

sce <- qread(sce_filename)

## Change barcode to match scRepertoire (ie. LEG1002-20R_ACATACGCATGTAGTC-1)
barcodes_original <- colnames(sce)

sampleid_barcode <- gsub(".*-1-", "", barcodes_original)
barcode_barcode <- gsub("-LEG.*", "", barcodes_original)

barcodes_new <- paste0(sampleid_barcode, "_", barcode_barcode)

colnames(sce) <- barcodes_new

## Add metadata
clinical_metadat <- read.csv("LEGACY_metadata.csv")
columns_keep <- c("donor_id", "days_V4afterV3", "Vaccine_arm")
clinical_metadat <- clinical_metadat[columns_keep]
clinical_metadat$Vaccine_arm <- paste0(clinical_metadat$Vaccine_arm, " FNA")

sce$vaccineside <- unlist(lapply(sce$donor_id, function(x) clinical_metadat$Vaccine_arm[match(x, clinical_metadat$donor_id)]))

sce$vaccine_type <- ifelse(sce$vaccineside==sce$samplesite, "Ipsilateral", "Contralateral")
sce$visitid_vaxtype <- paste0(sce$visit_id, "_", sce$vaccine_type)

## Based off of 250 unique clone cutoff
remove_samples <- c("LEG1028-20L", "LEG1017-20L", "LEG1006-20L", "LEG1006-20R", "LEG1019-20L")

sce <- subset(sce, , (!(sample_id %in% remove_samples)))
sce$sample_id <- droplevels(sce$sample_id)

## Adding CD4 and CD8
medium_resolution <- read.csv("LEGACY_annotation.csv") ## created from SCE object
sce$medium_celltype <- unlist(lapply(sce$celltype, function(x) medium_resolution$mediumbroad_resolution[match(x, medium_resolution$annotation)]))


## Subset for only CD4 and CD8 T cells (take out MAIT cells, NK cells )
sce <- subset(sce, , ((medium_celltype %in% c("CD4", "CD8"))))

sce$celltype <- droplevels(sce$celltype)

## Add clone type based on barcode
tcrlocus_barcode <- read.csv("QC_alphabetapairings.csv")

tcrlocus_barcode$new_barcode <- paste0(tcrlocus_barcode$sample_id, "_", tcrlocus_barcode$barcode)

sce$new_barcode <- colnames(sce)

sce$summarytra_trb <- unlist(lapply(sce$new_barcode, function(x) tcrlocus_barcode$summarytra_trb[match(x, tcrlocus_barcode$new_barcode)]))


# Add TCR info into SCE
sce <- combineExpression(combined.TCR, 
                         sce, 
                         cloneCall="strict",
                         chain = "both",
                         group.by = "donor_id", 
                         proportion = TRUE, 
                         addLabel = TRUE)

# CD4 expansion with only alpha beta
manuscript_outdir <- file.path("")


test <- sapply(combined.TCR, nrow)
sum_test <- sum(test)

cd4_sce <- subset(sce, , ((medium_celltype %in% c("CD4"))))

cd4_sce$celltype <- droplevels(cd4_sce$celltype)
cd4_sce$cluster <- cd4_sce$celltype

# Subset for only alpha-beta pairings
cd4_sce <- subset(cd4_sce, , ((summarytra_trb %in% c("tra1_trb1"))))

## Calculating clone_id
meta <- data.frame(colData(cd4_sce))
rownames(meta) <- cd4_sce@colData@rownames
clu <- which(colnames(meta) == "label") # as set by colLabels()
colnames(meta)[clu] <- "ident"

df <- meta

cloneCall <- "strict" 
cloneCall <- .theCall(df, cloneCall) #converts it to "CTstrict
barcodes <- rownames(df)
colnames(df)[ncol(df)] <- "majorCluster"

group.by <- "donor_id"
group.levels <- unique(df[, group.by])

df2 <- df %>% group_by(df[, group.by], df[, cloneCall]) %>% dplyr::mutate(n = n()) %>% as.data.frame()
rownames(df2) <- barcodes
remove.pos <- which(df2[, cloneCall] %in% c("", NA))
df2 <- df2[-remove.pos, ]
df2[, "clone.status"] <- ifelse(df2[, "n"] > 1, "Yes", "No")


# Statistics by donor and visitid
expan_summary <- df2 %>% group_by(donor_id, visit_id, majorCluster, clone.status) %>% summarise(expan_n = n()) %>% mutate(expan_freq = expan_n / sum(expan_n))

expan_summary %>% filter(clone.status=="Yes") %>% 
  ggplot(aes(x = majorCluster, y = expan_freq, fill = visitid_vaxtype, col = visitid_vaxtype)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +  
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(subtitle = "expansion > 1")

## Paired anova + posthoc (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/)
library(tidyverse)
library(ggpubr)
library(rstatix)

expan_summary <- df2 %>% group_by(donor_id, majorCluster, clone.status) %>% summarise(expan_n = n()) %>% mutate(expan_freq = expan_n / sum(expan_n))

expan_summary_stats <- expan_summary %>% filter(clone.status=="Yes")

ggboxplot(expan_summary_stats, x = "majorCluster", y = "expan_freq", add = "jitter",
          color = "majorCluster", palette = "jco")+
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova", label.y = 0.65)+      # Add global p-value
  stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = ".all.") +
  geom_hline(yintercept = median(expan_summary_stats$expan_freq), linetype = 2) +
  scale_color_manual(values=COLOURS$celltype)
ggsave(filename=file.path(manuscript_outdir, "CD4expanFreq_ANOVAwithWilcoxTestagainstMedian_pvalues_alphabetaonly.pdf"), width = 10, height = 6)

ggboxplot(expan_summary_stats, x = "majorCluster", y = "expan_freq", add = "jitter",
          color = "majorCluster")+
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova", label.y = 0.65)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = TRUE) +
  geom_hline(yintercept = median(expan_summary_stats$expan_freq), linetype = 2) +
  scale_color_manual(values=COLOURS$celltype)

ggsave(filename=file.path(manuscript_outdir, "CD4expanFreq_ANOVAwithWilcoxTestagainstMedian_pstars_alphabetaonly.pdf"), width = 6, height = 6)

write.csv(expan_summary_stats, file = file.path(manuscript_outdir, "CD4_expansion_freq_alphabetaonly.csv"))

# CD8 expansion with only alpha beta
test <- sapply(combined.TCR, nrow)
sum_test <- sum(test)

cd8_sce <- subset(sce, , ((medium_celltype %in% c("CD8"))))

cd8_sce$celltype <- droplevels(cd8_sce$celltype)
cd8_sce$cluster <- cd8_sce$celltype

# Subset for only alpha-beta pairings
cd8_sce <- subset(cd8_sce, , ((summarytra_trb %in% c("tra1_trb1"))))

## Calculating clone_id
meta <- data.frame(colData(cd8_sce))
rownames(meta) <- cd8_sce@colData@rownames
clu <- which(colnames(meta) == "label") # as set by colLabels()
colnames(meta)[clu] <- "ident"

df <- meta

cloneCall <- "strict" 
cloneCall <- .theCall(df, cloneCall) #converts it to "CTstrict
barcodes <- rownames(df)
colnames(df)[ncol(df)] <- "majorCluster"

group.by <- "donor_id"
group.levels <- unique(df[, group.by])

df2 <- df %>% group_by(df[, group.by], df[, cloneCall]) %>% dplyr::mutate(n = n()) %>% as.data.frame()
rownames(df2) <- barcodes
remove.pos <- which(df2[, cloneCall] %in% c("", NA))
df2 <- df2[-remove.pos, ]
df2[, "clone.status"] <- ifelse(df2[, "n"] > 1, "Yes", "No")


## Statistics by donor and visitid
expan_summary <- df2 %>% group_by(donor_id, visitid_vaxtype, majorCluster, clone.status) %>% summarise(expan_n = n()) %>% mutate(expan_freq = expan_n / sum(expan_n))

expan_summary %>% filter(clone.status=="Yes") %>% 
  ggplot(aes(x = majorCluster, y = expan_freq, fill = visitid_vaxtype, col = visitid_vaxtype)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +  
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(subtitle = "expansion > 1")

## Paired anova + posthoc (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/)
library(tidyverse)
library(ggpubr)
library(rstatix)

expan_summary <- df2 %>% group_by(donor_id, majorCluster, clone.status) %>% summarise(expan_n = n()) %>% mutate(expan_freq = expan_n / sum(expan_n))

expan_summary_stats <- expan_summary %>% filter(clone.status=="Yes")

ggboxplot(expan_summary_stats, x = "majorCluster", y = "expan_freq", add = "jitter",
          color = "majorCluster", palette = "jco")+
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova", label.y = 0.65)+      # Add global p-value
  stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = ".all.") +
  geom_hline(yintercept = median(expan_summary_stats$expan_freq), linetype = 2) +
  scale_color_manual(values=COLOURS$celltype)
ggsave(filename=file.path(manuscript_outdir, "CD8expanFreq_ANOVAwithWilcoxTestagainstMedian_pvalues_alphabetaonly.pdf"), width = 10, height = 6)

ggboxplot(expan_summary_stats, x = "majorCluster", y = "expan_freq", add = "jitter",
          color = "majorCluster")+
  rotate_x_text(angle = 45)+
  stat_compare_means(method = "anova", label.y = 0.65)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = TRUE) +
  geom_hline(yintercept = median(expan_summary_stats$expan_freq), linetype = 2) +
  scale_color_manual(values=COLOURS$celltype)

ggsave(filename=file.path(manuscript_outdir, "CD8expanFreq_ANOVAwithWilcoxTestagainstMedian_pstars_alphabetaonly.pdf"), width = 6, height = 6)

# Normalised Shannon Diversity
dat_test <- data.frame(sce$new_barcode, sce$sample_id, sce$celltype, sce$summarytra_trb, sce$visit_id, sce$visitid_vaxtype, sce$donor_id)

sce_counts <- dat_test %>% group_by("donor_id" = sce.donor_id, "celltype" = sce.celltype) %>% summarise(n = n())
sce_counts$celltype_donorid <- paste0(sce_counts$celltype, "__", sce_counts$donor_id)

sce$celltype_donorid <- paste0(sce$celltype, "__", sce$donor_id)

sce_clean <- subset(sce, , ((summarytra_trb %in% c("tra1_trb1"))))

diversity_dat2 <- clonalDiversity(sce_clean, group.by = "celltype_donorid", x.axis = "donor_id", cloneCall = "strict", exportTable = TRUE, skip.boots = T)

diversity_dat2$celltype <- gsub("__.*", "", diversity_dat2$celltype_donorid)
diversity_dat2$donor_id <- gsub(".*__", "", diversity_dat2$celltype_donorid)

diversity_dat2 <- merge(diversity_dat2, sce_counts[,c("celltype_donorid", "n")], by = "celltype_donorid", all.x = T)

diversity_dat2 <- merge(diversity_dat2, sample_metadat, by = "donor_id", all.x = T)

diversity_dat2 <- transform(diversity_dat2, u20_cells = ifelse(n<20, "U20_cells", "O20_cells"))
diversity_dat2 <- transform(diversity_dat2, u50_cells = ifelse(n<50, "U50_cells", "O50_cells"))

diversity_dat2 %>% filter(!(celltype %in% c("CD8_GZMB", "CD8_GZMK", "CD8_Naive"))) %>% ggplot(aes(x = as.factor(celltype), y = norm.entropy)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), aes(shape = u20_cells)) +
  #facet_wrap(~celltype, scales = "free_y") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette = "RdYlBu")

diversity_dat2 %>% filter(!(celltype %in% c("CD8_GZMB", "CD8_GZMK", "CD8_Naive"))) %>% ggplot(aes(x = as.factor(celltype), y = norm.entropy)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  #facet_wrap(~celltype, scales = "free_y") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette = "RdYlBu")

diversity_dat2 %>% filter((celltype %in% c("CD8_GZMB", "CD8_GZMK", "CD8_Naive"))) %>% ggplot(aes(x = as.factor(celltype), y = norm.entropy)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  #facet_wrap(~celltype, scales = "free_y") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_brewer(palette = "RdYlBu")

library(ggpubr)
diversity_dat2 %>% filter(!(celltype %in% c("CD8_GZMB", "CD8_GZMK", "CD8_Naive"))) %>%
  ggboxplot(x = "celltype", y = "norm.entropy",
            color = "celltype", palette = "simpsons")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_hline(yintercept = 1) +
  stat_compare_means(method = "anova", label.y = 1.01)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") 