# Figure 1 with ultrasound LN diameter, PCA

library(ggpubr)
library(rstatix)
library(factoextra)
library("FactoMineR")
library("factoextra")
library("missMDA") #imputeMFA
library(dplyr)
library(purrr)
library(viridis)

out_dir <- file.path(".")

# Import data
metadat <- read.csv("LN_USmeasurements_withArea.csv")

columns_remove2 <- c("X", "X.1", "V2_Ipsi_EllipseArea_mm2", "V2_WidthIpsiLN_mm", "V2_Contra_EllipseArea_mm2", "V2_WidthContraLN_mm", "V4_Contra_EllipseArea_mm2", "V4_WidthIpsiLN_mm", "V4_Ipsi_EllipseArea_mm2", "V4_WidthContraLN_mm")
metadat <- metadat[,!(names(metadat) %in% columns_remove2)]

# Adding in sampling total counts
totalcount_sample <- read.csv("sampleprocessing_LNtotalcounts.csv")
totalcount_sample$vaxtype <- ifelse(totalcount_sample$sampletype == totalcount_sample$vaccine_side, "Ipsi", "Contra")

totalcount_sample <- unique(totalcount_sample[c("donor_id", "visit_id", "sampletype", "totalnumbercells", "livedead", "vaxtype")])
totalcount_sample$total_count <- totalcount_sample$totalnumbercells*(totalcount_sample$livedead/100)
totalcount_sample$visit_id <- gsub("V2", "PreVax", totalcount_sample$visit_id)
totalcount_sample$visit_id <- gsub("V4", "PostVax", totalcount_sample$visit_id)
totalcount_sample$count_visitid_vaxtype <- paste0("totalcount_", totalcount_sample$visit_id, totalcount_sample$vaxtype)
totalcount_sample <- totalcount_sample[c("donor_id", "total_count", "count_visitid_vaxtype")]
totalcount_sample <- tidyr::pivot_wider(totalcount_sample, names_from = count_visitid_vaxtype, values_from = total_count)

# Merging data
pca_dat <- merge(metadat, totalcount_sample, by = "donor_id")

rownames(pca_dat) <- pca_dat$donor_id
pca_dat <- pca_dat[,!names(pca_dat) %in% c("donor_id")]

# Impute missing data
pca_dat <- as.data.frame(pca_dat)
res.impute2 <- imputePCA(pca_dat, ncp = length(pca_dat)/2) 

# PCA
res.pca2 <- PCA(res.impute2$completeObs)

# Visualise PCA results

plt1 <- fviz_contrib(res.pca2, choice = "var", axes = 1, top = 10)
ggsave(plt1, filename = file.path(out_dir, "Contribution_PCA1.pdf"), width = 8, height = 4)

plt2 <- fviz_contrib(res.pca2, choice = "var", axes = 2, top = 10)
ggsave(plt2, filename = file.path(out_dir, "Contribution_PCA2.pdf"), width = 8, height = 4)

variables_plot <- colnames(pca_dat)
for (variableplotting in variables_plot){
  #variableplotting <- "BMI"
  fviz_pca_ind(res.pca2, col.ind = "black") + geom_point(aes(color=res.pca2$call$X$BMI), size=3) + scale_color_viridis() + labs(col = variableplotting)
  ggsave(filename = file.path(out_dir, paste0("PCA_byDonor_", variableplotting,".pdf")), width = 10, height = 8)
}

## LN Diameter vs days post vax
COLOURS <- list()

COLOURS[['visitid_vaxtype']] <- c("V2_Contra" = "#A6CEE3", 
                                  "V2_Ipsi" = "#FDBF6F",
                                  "V4_Contra" = "#1F78B4", 
                                  "V4_Ipsi" = "#FF7F00")

COLOURS[['visitid_vaxtype2']] <- c("PreVaxContra" = "#A6CEE3", 
                                   "PreVaxIpsi" = "#FDBF6F",
                                   "PostVaxContra" = "#1F78B4", 
                                   "PostVaxIpsi" = "#FF7F00")

dat <- read.csv(file.path(in_dir, "LymphNodeSize_Metadata.csv"))

columns_remove <- c("X")

dat <- dat[,!(names(dat) %in% columns_remove)]

width_dat <- dat[,names(dat) %in% c("donor_id", "days_V4afterV3", "V2_LengthIpsiLN_mm", "V2_LengthContraLN_mm", "V4_LengthIpsiLN_mm", "V4_LengthContraLN_mm")]

width_dat <- width_dat %>% pivot_longer(!c(donor_id, days_V4afterV3), names_to = "sampletype", values_to = "diameter_mm")
width_dat$sampletype <- sub("LN_mm", "", width_dat$sampletype)
width_dat$sampletype <- sub("Length", "", width_dat$sampletype)
width_dat$timepoint <- sub("_.*", "", width_dat$sampletype)
width_dat$sampleside <- sub(".*_", "", width_dat$sampletype)

compare_means(diameter_mm ~ sampletype, data = width_dat)
my_comparisons <- list(c("V2_Ipsi", "V2_Contra"), c("V2_Ipsi", "V4_Ipsi"), c("V2_Contra", "V4_Contra"), c("V4_Ipsi", "V4_Contra"))

width_dat$sampletype <- factor(width_dat$sampletype, levels = c("V2_Contra", "V2_Ipsi", "V4_Contra", "V4_Ipsi"))
ggboxplot(width_dat, x = "sampletype", y = "diameter_mm", color = "sampletype", add = "jitter") +
  scale_color_manual(values = COLOURS$visitid_vaxtype) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 30)                   # Add global p-value
ggsave(filename=file.path(out_dir, "UltrasoundMeasurement_Diameter.pdf"), width = 5, height = 6)


## Days of post-vax
width_dat %>% filter(sampletype == "V4_Ipsi") %>% ggplot(aes(x = days_V4afterV3, y = diameter_mm)) + geom_jitter(width = 0.05) +
  geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
  stat_cor(label.x = 3, label.y = 25, size = 4) +
  theme_pubr() +
  labs(subtitle = "V4_Ipsi") +
  ylim(5, 25)
ggsave(filename=file.path(out_dir, "V4IpsiDiameter_daysV4afterV3.pdf"), width = 6, height = 5)


width_dat %>% filter(sampletype == "V4_Contra") %>% ggplot(aes(x = days_V4afterV3, y = diameter_mm)) + geom_jitter(width = 0.05) +
  geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
  stat_cor(label.x = 3, label.y = 25, size = 4) +
  theme_pubr() +
  labs(subtitle = "V4_Contra") +
  ylim(5, 25)
ggsave(filename=file.path(out_dir, "V4ContraDiameter_daysV4afterV3.pdf"), width = 6, height = 5)
