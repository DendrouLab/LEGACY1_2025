### Updated May 9, 2023
### Merge .tsv AIRR files by donor
### singularity shell -B /well/legacy/:/well/legacy/  /well/legacy/users/pps914/packages/immcantation/immcantation_suite-4.4.0.sif
### Rscript 04_02_mergeperdonor.R --donorfilenames /well/legacy/users/pps914/scripts/immcantation_bydonor/04_MergebyDonor/groupedbydonor/LEG1018_db-pass_cleaned_filenames.tsv \


# Libraries
library(dplyr)
library(optparse)
library(alakazam)


option_list <- list(
    make_option(c("--donorfilenames"), default="",
                help="Input file with all the filenames DONORID_db-pass_cleaned_filenames.tsv"),
  	make_option(c("--outdir"), default="",
                help = "Output directory"), 
    make_option(c("--savename"), default="_db-pass_cleaned_combinedbydonor.tsv",
                help = "The suffix of the created .tsv file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")

print(opt)

# Import
read_file <- read.table(opt$donorfilenames, sep = '\t', header = TRUE)
donor_id <- unique(read_file$donor_id)

# IG merge
ig_files <- which(read_file$vdj_type == "ig")

merge_file <- {}
for (i in ig_files){
    temp <- readChangeoDb(read_file$data_dir[i])
    temp$sample_id <- read_file$sample_id[i]
    temp$cellbarcode <- paste0(gsub("_contig_[1-20]*", "", temp$sequence_id), "-", temp$sample_id) #different contigs refer to different locus (after cleanup)

    merge_file <- rbind(merge_file, temp)
}

message(capture.output((table(as.factor(merge_file$sample_id)))))
message(capture.output((table(as.factor(merge_file$locus)))))

writeChangeoDb(merge_file, file = file.path(opt$outdir, paste0(donor_id, opt$savename)))

