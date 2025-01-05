# 
# Based on https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html
#### singularity shell -B /well/legacy/:/well/legacy/  /well/legacy/users/pps914/packages/immcantation/immcantation_suite-4.4.0.sif
#### Rscript run_immcantation_cleanup.R --infile whatever.csv --outfile whatever.csv

##TO DO
### MAY NEED TO REMOVE disagreements between the C-region primers and the reference alignment
#ParseDb.py select -d db.tsv -f v_call j_call c_call -u "IGH" \
#    --logic all --regex --outname heavy
#ParseDb.py select -d db.tsv -f v_call j_call c_call -u "IG[LK]" \
#    --logic all --regex --outname light

library(alakazam)
library(data.table)
library(dowser)
library(dplyr)
library(ggplot2)
library(scoper)
library(shazam)
library(optparse)

option_list <- list(
    make_option(c("--vdj_type"), default="none",
                help="ig or tr"),
    make_option(c("--samplename"), default="none",
                help="LEG1018-20L_ig"),
    make_option(c("--outdir"), default="none",
                help=""),
    make_option(c("--infile"), default="none",
                help="Input file name such as 02_MakeDb/LEG1018-20L_ig_db-pass.tsv"),
  	make_option(c("--outfile"), default="none",help = "Output file name such as 02_MakeDb/LEG1018-20L_ig_db-pass_cleaned.tsv"),
    make_option(c("--csvfilename"), default="none",
                help="02_MakeDb/LEG1018-20L_ig_db-pass_cleaned.csv")
    )

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")

print(opt)

#data <- readChangeoDb("/well/legacy/users/pps914/output/005_year1multi/005_year1multi_immarch/02_MakeDb/LEG1002-4FL_tr_db-pass.tsv")
data <- readChangeoDb(opt$infile)
message(cat(paste("There are", nrow(data), "rows in the data.\n")))
numbers_pre <- nrow(data)
data <- data %>% filter(productive == TRUE)

message(cat(paste("There are", nrow(data), "rows in the data, keeping only productive.\n")))
numbers_production <- nrow(data)

gexbarcodeids <- read.csv("/well/legacy/users/pps914/output/005_year1multi/005_year1multi_subclustering/final_annotation/final005_year1multi_finalannotations_withredcap_01122023/Finalannotations19102023_withmetadata01122023.csv")

if (opt$vdj_type == "ig"){
  # remove cells with multiple heavy chain
  multi_heavy <- table(filter(data, locus == "IGH")$cell_id)
  multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]

  data <- filter(data, !cell_id %in% multi_heavy_cells)
  message(cat(paste("There are", nrow(data), "rows in the data after filtering out cells with multiple heavy chains.\n")))
  numbers_nomultiheavy <- nrow(data)

  # split cells by heavy and light chains
  heavy_cells <- filter(data, locus == "IGH")$cell_id
  light_cells <- filter(data, locus == "IGK" | locus == "IGL")$cell_id
  no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]

  data <- filter(data, !cell_id %in% no_heavy_cells)
  message(cat(paste("There are", nrow(data), "rows in the data after filtering out cells without heavy chains.")))
  numbers_haveheavy <- nrow(data)
  
  # Export version 1 (not matched with GEX)
  writeChangeoDb(data, file = opt$outfile)
  
  # Filter for GEX matching only
  #opt$samplename <- "LEG1002-4FL"
  gexbarcodeids <- subset(gexbarcodeids, sample_id == opt$samplename)
  gexbarcodeids$cell_id <- gsub(paste0("-", opt$samplename), "", gexbarcodeids$index)
  
  withgex <- table(filter(data, cell_id %in% gexbarcodeids$cell_id)$cell_id)
  withgex_cells <- names(withgex)[withgex > 1]
  
  data <- filter(data, cell_id %in% withgex_cells)
  numbers_havegex <- nrow(data)
  message(cat(paste("There are", nrow(data), "rows in the data after filtering out cells without GEX. \n")))
  
  # Export
  message(cat(paste("Exporting: ", opt$samplename)))
  numbers <- data.frame("Pre-filtering" = numbers_pre, "Productive_only" = numbers_production, "No_multiHeavyChains" = numbers_nomultiheavy, "Have_HeavyChain" = numbers_haveheavy, "Have_GEX" = numbers_havegex)
  write.csv(numbers, file = opt$csvfilename, row.names=FALSE)
  filename2 <- gsub(".tsv", "_gexonly.tsv", opt$outfile)
  writeChangeoDb(data, file = filename2)
  
} else if (opt$vdj_type == "tr") {
  
  # Filtering
  ## Remove TRD and TRG
  data <- data %>% filter(locus %in% c("TRA", "TRB"))
  message(cat(paste("There are", nrow(data), "cells. Removed locus = TRD or TRG")))
  numbers_noTRDG <- nrow(data)

  ## Remove any cell with 3+ TRA
  multi_tra <- table(filter(data, locus == "TRA")$cell_id)
  multi_tra_cells <- names(multi_tra)[multi_tra > 2]

  data <- filter(data, !cell_id %in% multi_tra_cells)
  message(cat(paste("There are", nrow(data), "rows in the data after filtering out cells with 3+ TRA chains.\n")))
  numbers_no3TRA <- nrow(data)

  ## Remove any cell with 3+ TRB
  multi_trb <- table(filter(data, locus == "TRB")$cell_id)
  multi_trb_cells <- names(multi_trb)[multi_trb > 2]

  data <- filter(data, !cell_id %in% multi_trb_cells)
  message(cat(paste("There are", nrow(data), "rows in the data after filtering out cells with 3+ TRB chains.\n")))
  numbers_no3TRB <- nrow(data)

  # Convert to ImmunArch-friendly, AIRR format
  airr_dat <- data
  airr_dat$cdr1_aa <- translateDNA(airr_dat$cdr1, trim = TRUE)
  airr_dat$cdr2_aa <- translateDNA(airr_dat$cdr2, trim = TRUE)
  airr_dat$cdr3_aa <- translateDNA(airr_dat$cdr3, trim = TRUE)
  airr_dat$cdr3_aa <- translateDNA(airr_dat$cdr3, trim = TRUE)
  airr_dat$fwr1_aa <- translateDNA(airr_dat$fwr1, trim = TRUE)
  airr_dat$fwr2_aa <- translateDNA(airr_dat$fwr2, trim = TRUE)
  airr_dat$fwr3_aa <- translateDNA(airr_dat$fwr3, trim = TRUE)
  airr_dat$fwr4_aa <- translateDNA(airr_dat$fwr4, trim = TRUE)

  airr_dat$duplicate_count <- airr_dat$umi_count #Immunarch expects column called duplicate_count, not umi_count
 
  # Export without GEX cleanup
  writeChangeoDb(data, file = opt$outfile)
  airr::write_rearrangement(airr_dat, file = paste0(opt$outdir, "/immunarch_airr/", opt$samplename, ".tsv"))
  
  # Filter for GEX matching only
  #opt$samplename <- "LEG1002-4FL"
  gexbarcodeids <- subset(gexbarcodeids, sample_id == opt$samplename)
  gexbarcodeids$cell_id <- gsub(paste0("-", opt$samplename), "", gexbarcodeids$index)
  
  withgex <- table(filter(data, cell_id %in% gexbarcodeids$cell_id)$cell_id)
  withgex_cells <- names(withgex)[withgex > 1]
  
  data <- filter(data, cell_id %in% withgex_cells)
  numbers_havegex <- nrow(data)
  message(cat(paste("There are", nrow(data), "rows in the data after filtering out cells without GEX. \n")))
  
  # Adding in CDR3_aa column for scRepertoire
  data$cdr3_aa <- translateDNA(data$cdr3, trim = TRUE)

  # Exporting with GEX cleanup
  message(cat(paste("Exporting: ", opt$samplename)))
  numbers <- data.frame("Pre-filtering" = numbers_pre, "Productive_only" = numbers_production, "No_TRD_TRG" = numbers_noTRDG, "Lessthan3_TRA" = numbers_no3TRA, "Lessthan3_TRB" = numbers_no3TRB, "Have_GEX" = numbers_havegex)
  write.csv(numbers, file = opt$csvfilename, row.names=FALSE)
  filename2 <- gsub(".tsv", "_gexonly_data.tsv", opt$outfile)
  writeChangeoDb(data, file = filename2)
  airr::write_rearrangement(airr_dat, file = paste0(opt$outdir, "/immunarch_airr_withGEXonly/", opt$samplename, ".tsv"))

} else {
  message(cat(paste("This sample was skipped: ", opt$samplename)))
}

message("Done")

