# Figure 1 with ultrasound area & STERR, PCA

library(ggpubr)
library(rstatix)
library(factoextra)
library("FactoMineR")
library("factoextra")
library("missMDA") #imputeMFA
library(dplyr)
library(purrr)
library(viridis)

out_dir <- "."

# Import data
metadat <- read.csv("LN_USmeasurements_withArea.csv")

metadat <- transform(metadat, V2_Ipsi_EllipseArea_mm2 = pi*(V2_LengthIpsiLN_mm/2)*(V2_WidthIpsiLN_mm/2))
metadat <- transform(metadat, V2_Contra_EllipseArea_mm2 = pi*(V2_LengthContraLN_mm/2)*(V2_WidthContraLN_mm/2))
metadat <- transform(metadat, V4_Ipsi_EllipseArea_mm2 = pi*(V4_LengthIpsiLN_mm/2)*(V4_WidthIpsiLN_mm/2))
metadat <- transform(metadat, V4_Contra_EllipseArea_mm2 = pi*(V4_LengthContraLN_mm/2)*(V4_WidthContraLN_mm/2))

columns_remove2 <- c("V2_LengthIpsiLN_mm", "V2_WidthIpsiLN_mm", "V2_LengthContraLN_mm", "V2_WidthContraLN_mm", "V4_LengthIpsiLN_mm", "V4_WidthIpsiLN_mm", "V4_LengthContraLN_mm", "V4_WidthContraLN_mm")
metadat <- metadat[,!(names(metadat) %in% columns_remove2)]

# Calculate ellipse area uncertaintity
ellipse_stdev <- function(dat, width_name, length_name){
  width_dat = dat[,width_name]
  length_dat = dat[,length_name]
  
  num_samples = length(width_dat)
  
  width_mean = mean(width_dat)
  length_mean = mean(length_dat)
  
  width_sd = sd(width_dat)
  length_sd = sd(length_dat)
  
  ellipse_stdev = (pi/4)*((width_mean*length_mean))*sqrt((width_sd^2/width_mean^2)+(length_sd^2/length_mean^2))
  return(ellipse_stdev)
}

summary_dat <- metadat %>% group_by(sampletype, timepoint, sampleside) %>% summarise(mean_area = mean(ellipseArea_mm2), 
                                                                                      n_number = n())
PreVax_Ipsi <- ellipse_stdev(metadat, "V2_WidthIpsiLN_mm", "V2_LengthIpsiLN_mm")
PreVax_Contra <- ellipse_stdev(metadat, "V2_WidthContraLN_mm", "V2_LengthContraLN_mm")
PostVax_Ipsi <- ellipse_stdev(metadat, "V4_WidthIpsiLN_mm", "V4_LengthIpsiLN_mm")
PostVax_Contra <- ellipse_stdev(metadat, "V4_WidthContraLN_mm", "V4_LengthContraLN_mm")

stdev_ellipsearea <- c(PreVax_Contra, PreVax_Ipsi, PostVax_Contra, PostVax_Ipsi)
names(stdev_ellipsearea) <- c("PreVax_Contra", "PreVax_Ipsi", "PostVax_Contra", "PostVax_Ipsi")

summary_dat$stdev_area <- stdev_ellipsearea
summary_dat$sterror_area <- summary_dat$stdev_area/sqrt(summary_dat$n_number)

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

# Comparing processing cell numbers with number days post-vax
totalcount_sample_save <- tidyr::pivot_wider(totalcount_sample, names_from = count_visitid_vaxtype, values_from = total_count)

us_dat <- merge(metadat, totalcount_sample_save, by = "donor_id")

Q_V4ipsi <- quantile(us_dat$totalcount_PostVaxIpsi, probs=c(.25, .75), na.rm = FALSE)
iqr_V4ipsi <- IQR(us_dat$totalcount_PostVaxIpsi)
no_outliers_V4ipsi <- subset(us_dat, us_dat$totalcount_PostVaxIpsi > (Q_V4ipsi[1] - 1.5*iqr_V4ipsi) & us_dat$totalcount_PostVaxIpsi < (Q_V4ipsi[2]+1.5*iqr_V4ipsi))

Q_V4contra <- quantile(us_dat$totalcount_PostVaxContra, probs=c(.25, .75), na.rm = FALSE)
iqr_V4contra <- IQR(us_dat$totalcount_PostVaxContra)
no_outliers_V4contra <- subset(us_dat, us_dat$totalcount_PostVaxContra > (Q_V4contra[1] - 1.5*iqr_V4contra) & us_dat$totalcount_PostVaxContra < (Q_V4contra[2]+1.5*iqr_V4contra))


ggplot(no_outliers_V4ipsi,aes(days_V4afterV3, totalcount_PostVaxIpsi)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
  stat_cor(label.x = 3, label.y = 1E7, size = 4) +
  stat_regline_equation(label.x = 3, label.y = 1.2E7, size = 4) +
  theme_pubr()
ggsave(filename = file.path(dir, "V4TotalCellCountIpsi_daysV4afterV3_noOutliers.pdf"))

ggplot(no_outliers_V4contra,aes(days_V4afterV3, totalcount_PostVaxContra)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
  stat_cor(label.x = 3, label.y = 2500000, size = 4) +
  stat_regline_equation(label.x = 3, label.y = 2000000, size = 4) +
  theme_pubr()
ggsave(filename = file.path(dir, "V4TotalCellCountContra_daysV4afterV3_noOutliers.pdf"))