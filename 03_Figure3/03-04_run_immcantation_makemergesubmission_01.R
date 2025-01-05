### Updated May 9, 2023
### Generate mergedonor_submissionfile.tsv and merge_files.csv

# Libraries
library(dplyr)
library(optparse)


option_list <- list(
    make_option(c("--infile"), default="/well/legacy/users/pps914/output/005_year1multi/005_year1multi_immarch/01032024_immcantationsubmissionfile.tsv",
                help="Input file name SOMETHING.tsv"),
  	make_option(c("--outdir"), default="/well/legacy/users/pps914/output/005_year1multi/005_year1multi_immarch/04_MergebyDonor",
                help = "Output directory"), 
    make_option(c("--filelocations"), default="/well/legacy/users/pps914/output/005_year1multi/005_year1multi_immarch/03_CleanedDb/",
                help = "Where the DB cleaned up files are"),
    make_option(c("--fileextension"), default="_db-pass_cleaned.tsv",
                help = "What the DB cleaned up files are called")
    )

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")

print(opt)

# Import
dat <- read.table("/well/legacy/users/pps914/output/005_year1multi/005_year1multi_immarch/01032024_immcantationsubmissionfile.tsv", sep = '\t', header = TRUE)
#dat$file_name <- paste0(dat$sample_id, "_", dat$vdj_type)

# Generate .csv per sample, and a submission .tsv file per run
donor_ids <- unique(dat$donor_id)

out_submissiontsv <- {}

for (i in 1:length(donor_ids)){
  dat_temp <- subset(dat, donor_id == donor_ids[i])
  message(cat(print(paste(donor_ids[i], " | All records | Rows:",nrow(dat_temp)))))
  
  message(cat(paste("There are", sum(dat_temp$vdj_type == "ig"), "ig samples in", donor_ids[i])))
  message(cat(paste("There are", sum(dat_temp$vdj_type == "tr"), "tr samples in", donor_ids[i])))

  # Library Output
  out_librarycsv <- data.frame("donor_id" = dat_temp$donor_id, "sample_id" = dat_temp$sample_id, "vdj_type" = dat_temp$vdj_type, "data_dir" = paste0(opt$filelocations, dat_temp$sample_id, "_", dat_temp$vdj_type, opt$fileextension))
  write.table(out_librarycsv, file.path(opt$outdir,"groupedbydonor", paste0(donor_ids[i], "_db-pass_cleaned_filenames.tsv")), row.names=FALSE, quote=FALSE, sep = '\t')
    
  # Submission Output
  submission_temp <- data.frame("donor_id" = donor_ids[i], "filenames_csv" = file.path(opt$outdir,"groupedbydonor", paste0(donor_ids[i], "_db-pass_cleaned_filenames.tsv")), "end" = "end")
  out_submissiontsv <- rbind(out_submissiontsv, submission_temp)
}

write.table(out_submissiontsv, file.path(opt$outdir, "04_MergeByDonor_submissionfile.tsv"), quote=FALSE, sep='\t', row.names = FALSE)