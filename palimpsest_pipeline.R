# INSTRUCTIONS
# 1. This pipeline is designed for mutation signature analysis using Palimpsest and works best with data from cBioPortal.
# 2. Ensure that mutation and clinical data are based on the **hg19** reference genome, as hg38 may cause errors in annotation steps.
# 3. Prepare mutation data to include a `PATIENT_ID` column that matches with clinical data.
# 4. Clinical data should include `PATIENT_ID` and `SUBTYPE` columns, and all column names should be in uppercase.
# 5. Update all placeholder paths and variables (gene_name, subtype, file paths, etc.) before running the script.
# 6. Run each section sequentially to ensure proper setup and avoid runtime errors.

# Libraries
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(Palimpsest)

# CHANGE THESE VARIABLES
gene_name <- "^name_of_the_target_gene"  # Gene of interest (do not remove "^")
subtype <- "cancer_subtype as mentioned in the clinical file"  # Specific cancer subtype
mut_location <- "write/the/mutation/file/path/here"  # Path to mutation data file
clinical_location <- "write/the/clinical/file/path/here"  # Path to clinical data file
num_of_denovo_sings <- 5  # Number of de novo signatures to extract (modify as needed)
parent_subfile_name <- "subfile_name_here"
parent_file_name <- "/write/the/parent/file/path/here/"

print_statement <- "Total unique patients"  # Optional; no need to change

# RESULT DIRECTORY SETUP
# Creates a parent and a gene-specific result directory if not already present
file_name <- paste0("/", gene_name, "_", subtype, "_Palimpsest")
result_parent_dir <- paste0(parent_file_name, parent_subfile_name) if (!file.exists(result_parent_dir)) dir.create(result_parent_dir)
result_dir <- paste0(result_parent_dir, file_name) if (!file.exists(result_dir)) dir.create(result_dir)

# FILE TO LOG DETAILS
details_file <- paste0(result_dir, "/", gene_name, "_", subtype, "details.txt") file.create(details_file)

#-LOAD MUTATION DATA----
mutation_df <- read.delim2(mut_location, header = TRUE, sep = ",", stringsAsFactors = FALSE)

#-FIXING THE MUTATION DATASET - Clean the mutation file as needed ----
# Ensure mutation dataset has a `PATIENT_ID` column aligned with clinical data.
mut_clean_df <- mutation_df  # Save cleaned mutation data here

# usual fixes
# mut_clean_df$PATIENT_ID <- gsub(".{3}$", "", mut_clean_df$Tumor_Sample_Barcode)
# mut_clean_df$PATIENT_ID <- gsub("-01", "", mut_clean_df$Tumor_Sample_Barcode)

# log a summary
write(paste(print_statement, "in mutation dataframe are:", length(unique(mut_clean_df$PATIENT_ID))), file = details_file, append = TRUE)

#-Clinical Data & further filtration----
# Do the necessary data cleaning as needed. column names should be in UPPERCASE and should have ("PATIENT_ID" & "SUBTYPE") colnames

# Once all the above steps are done, select everything from this line till the end and hit run.
clinical_df <- read.delim2(clinical_location, header = TRUE, skip = 4, stringsAsFactors = FALSE)
# log a summary
write(paste(print_statement, "in clinical dataframe are:", length(unique(clinical_df$PATIENT_ID))), file = details_file, append = TRUE)

# clinical data with only cancer
cancer_df <- subset(clinical_df, clinical_df$SUBTYPE == subtype)
cancer_patient_ID <- cancer_df$PATIENT_ID #list of patient IDs with only cancer
# log a summary
write(paste(print_statement, "in clinical dataframe having ", subtype, "are:", length(unique(cancer_df$PATIENT_ID))), file = details_file, append = TRUE)

cancer_mut_df <- subset(mut_clean_df, mut_clean_df$PATIENT_ID %in% cancer_patient_ID) #ONLY cancer patients

#df of cancer & gene_name +ve patients
gene <- subset(cancer_mut_df, grepl(gene_name, cancer_mut_df$Hugo_Symbol))
# log a summary
write(paste(print_statement, "with ", gene_name, " mutation are:", list(unique(gene$PATIENT_ID)), ". Total No. of patients: ", length(unique(gene$PATIENT_ID))), file = details_file, append = TRUE)

rm_silent_gene <- subset(gene, gene$Variant_Classification != "Silent")
final_gene_patient_id <- rm_silent_gene$PATIENT_ID #list of gene & cancer patient ids
silent_gene_patients <- subset(gene, gene$Variant_Classification == "Silent")
# log a summary
write(paste(print_statement, "in silent", gene_name, "mutation are:", list(unique(silent_gene_patients$PATIENT_ID)), ". Total No. of 'non_silent' patients: ", length(unique(rm_silent_gene$PATIENT_ID))), file = details_file, append = TRUE)
rm(silent_gene_patients)

#all the gene mut data of gene & cancer +ve patients
all_gene_df <- subset(mut_clean_df, mut_clean_df$PATIENT_ID %in% final_gene_patient_id)

#all the gene mut data of gene & cancer +ve patients with ref & alt == 1
snp_df <- subset(all_gene_df, all_gene_df$Variant_Type == "SNP")

#-Preparing the df for palimpsest in vcf format----
vcf_df <- snp_df[, c(ncol(snp_df), 11, 5, 6, 12, 14, 34, 35, 36, 1)]

colnames(vcf_df) <- c("Sample", "Type", "CHROM", "POS", "REF", "ALT", "Tumor_Varcount", "Tumor_Depth", "Normal_Depth", "Driver")
vcf_df$Type <- gsub("SNP", "SNV", vcf_df$Type)
vcf_df$CHROM <- paste0("chr", vcf_df$CHROM)
View(vcf_df)


#-Palimpsest analysis----
# Scanning for de novo mutations

#-STEP: 1 - Load Palimpsest & reference genome packages----

#-STEP: 2 - Saving the vcf file----
write.csv(vcf_df, file = paste0(result_dir, "/", "vcf_", gene_name, "_", subtype, "_.csv"))

#-STEP: 3 - Load and annotate mutation data----
vcf <- annotate_VCF(vcf = vcf_df, ref_genome = BSgenome.Hsapiens.UCSC.hg19)

#-STEP: 4 - De Novo Extraction of Single Base Substitution (SBS) Signatures----
SBS_result_dir <- file.path(result_dir,"SBS_denovo"); if(!file.exists(SBS_result_dir)) dir.create(SBS_result_dir)
SBS_input <- palimpsest_input(vcf = vcf, Type = "SBS")

#NMF extraction
NMF_data <- file.path(SBS_result_dir, "NMF_extraction"); if(!file.exists(NMF_data)) dir.create(NMF_data)
#to automatically generate the optimal number of denovo signatures
# SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input, range_of_sigs = 1:10, nrun = 10, resdir = NMF_data)
# OR
SBS_denovo_sigs <- NMF_Extraction(input_matrices = SBS_input, num_of_sigs = num_of_denovo_sings, nrun = 10, resdir = NMF_data)

# Getting proportions of SBS signatures
proportion_SBS_denovo <- SBS_denovo_sigs
proportion_SBS_denovo <- proportion_SBS_denovo * 100
proportion_SBS_denovo <- round(proportion_SBS_denovo, 2)
write.csv(proportion_SBS_denovo, paste0(result_dir, "/denovo_SBS_proportion.csv"), row.names = TRUE)

#-STEP: 5 - Compare Extracted Signatures with known COSMIC Signatures----
compare_tab <- compare_results(reference_sigs = SBS_cosmic, extraction_1 = SBS_denovo_sigs);
compare_tab <- compare_tab[order(compare_tab$Palimp_Equivalent), ]
compare_tab$Ref_Palimp_cosine_score <- round(compare_tab$Ref_Palimp_cosine_score, 2)
compare_tab <- compare_tab[!is.na(compare_tab$Palimp_Equivalent), ]
readr::write_delim(compare_tab, file = file.path(SBS_result_dir,"Comparison_table.txt"))

# Cosine_Similarity_Heatmap.pdf
pdf(file.path(SBS_result_dir, "Cosine_Similarity_Heatmap.pdf"), width = 11, height = 10)
SBS_cosine_similarities <- deconvolution_compare(SBS_denovo_sigs, SBS_cosmic)
dev.off()

save(SBS_cosine_similarities, file = file.path(SBS_result_dir, "Cosine_Similarity_matrix.RData"))

#-STEP: 6 - Calculate and Visualize Signature Exposures----
SBS_col <- signature_colour_generator(rownames(SBS_denovo_sigs))
SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_denovo_sigs, resdir = SBS_result_dir, signature_colours = SBS_col, input_vcf = vcf)

# signature_content_plot.pdf
pdf(file.path(SBS_result_dir, "signature_content_plot.pdf"), width=15, height=10)
deconvolution_exposure(signature_colours = SBS_col, signature_contribution = SBS_signatures_exp)
dev.off()

#-Scanning for all known SBS signatures----
# # Adding custom colors
sig_cols <- c("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD", "#8C8C8C", "#E8C547", "#56B4E9", "#D55E00", "#F0E442", "#009E73",
                       "#AA4499", "#999999", "#117A65", "#D2B48C", "#CC79A7", "#77AB43","#E1A692", "#82BC9F", "#A593E0", "#B3B3B3", "#BC80BD", "#CCEBC5",
                       "#DECBE4", "#FED9A6", "#FBB4AE", "#B3CDE3", "#CAB2D6", "#FFED6F","#B15928", "#FF6F61", "#6A3D9A", "#FFFFB3", "#E31A1C", "#FDBF6F",
                       "#1F78B4", "#33A02C", "#FB9A99", "#FF7F00", "#6A3D9A", "#B15928","#F0F0F0", "#A6D854", "#A3A3A3", "#8DD3C7", "#BC80BD", "#FFDB58",
                       "#FFFF99", "#BEBADA", "#CCEBC5", "#DC267F", "#648FFF", "#FE6100","#FFB000", "#F4A582", "#92C5DE", "#7BAFDE", "#CA0020", "#0571B0",
                       "#8B3A62", "#B2ABD2", "#E78AC3", "#A6CEE3", "#33A02C", "#B3DE69","#FFB6C1", "#9E0142", "#D73027", "#FDAE61", "#FEE08B", "#A6D96A","#1A9850"
)

all_known_sbs <- c(rownames(SBS_cosmic))
SBS_cancer_sigs_all <- SBS_cosmic[all_known_sbs, ]

# calculate and plot the exposure of the signatures across the series
all_SBS_result_dir <- file.path(SBS_result_dir,"all_SBS_result_dir"); if(!file.exists(all_SBS_result_dir)) dir.create(all_SBS_result_dir)
all_SBS_signatures_exp <- deconvolution_fit(input_matrices = SBS_input, input_signatures = SBS_cancer_sigs_all, signature_colours = sig_cols, resdir = all_SBS_result_dir)

graphics.off()
pdf(file.path(all_SBS_result_dir, "all_known_signature_plot.pdf"), width=15, height=10)
deconvolution_exposure(signature_contribution = all_SBS_signatures_exp, signature_colours = sig_cols)
dev.off()

#-ISSUE - Pie charts in the samples are not picking up colors.----
# FIX - Need to generate the pie chart from raw data

all_known_sbs_prop_df <- all_SBS_signatures_exp$sig_props

graphics.off()
pdf(file.path(all_SBS_result_dir, "all_known_signature_piecharts.pdf"), width=15, height=10)
j <- 1
while (j <= length(rownames(all_known_sbs_prop_df))) {
  
  all_known_sbs_prop_df_values <- as.numeric(all_known_sbs_prop_df[j, ])
  all_known_sbs_prop_df_lables <- colnames(all_known_sbs_prop_df)
  non_zero_all_known_sbs_prop_df_values <- all_known_sbs_prop_df_values[all_known_sbs_prop_df_values != 0]
  non_zero_all_known_sbs_prop_df_lables <- all_known_sbs_prop_df_lables[all_known_sbs_prop_df_values != 0]
  
  # Percentage
  all_known_sbs_value_percentage <- round(100 * non_zero_all_known_sbs_prop_df_values / sum(non_zero_all_known_sbs_prop_df_values), 1)
  # Generate labels with names and percentages
  all_known_sbs_labels_with_percentages <- paste(non_zero_all_known_sbs_prop_df_lables, all_known_sbs_value_percentage, "%", sep = " ")
  
  
  pie(non_zero_all_known_sbs_prop_df_values, labels = all_known_sbs_labels_with_percentages, main = rownames(all_known_sbs_prop_df)[j], col = sig_cols)
  j <- j+1
  
}
dev.off()

#-Checking for known DNA Repair Gene mutations----
dna_rep_genes <- c("APEX1", "UNG", "LIG1", "POLB", "RPA", "DDB1", "MLH1", "MSH2", "MSH6", "PMS2", "EXO1", "PARP1", "XRCC1", "XRCC4","53BP1", "ERCC1",
                   "XPA", "BRCA1", "BRCA2", "RAD51", "RAD52", "ATM", "ATR", "CHEK2", "PALB2", "BARD1", "BRIP1", "BLM", "ERCC5", "POLQ", "TREX1", "POLK", "PRKDC")

check_df <- subset(all_gene_df, toupper(Hugo_Symbol) %in% dna_rep_genes) # taking from the df that have all the gene mutation data for patient IDs who ONLY & gene +ve.
write.csv(check_df, file = paste0(SBS_result_dir, "/check_genes.csv"), quote = FALSE, row.names = FALSE)

#-RESULT COMPILATION----
# gene ~~
names(rm_silent_gene) <- toupper(names(rm_silent_gene))
gene_final <- rm_silent_gene[, c(ncol(rm_silent_gene), 1, 10, 40)]
gene_final$Target_Gene <- paste(gene_final$HUGO_SYMBOL, gene_final$VARIANT_CLASSIFICATION, gene_final$HGVSP_SHORT, sep = "; ")
gene_final <- gene_final[, c(1, ncol(gene_final))]
# MERGE
gene_grouped <- group_by(gene_final, PATIENT_ID)
gene_merged <- summarise(gene_grouped, Target_Gene = paste(Target_Gene, collapse = "  ||  "))
rm(gene_grouped)
# View(gene_merged) #df 1

# dna_rep_genes ~~
names(check_df) <- toupper(names(check_df))
dna_rep_final <- check_df[, c(ncol(check_df), 1)]
colnames(dna_rep_final)[2] <- "DNA_Repair_Genes"
# MERGE
dna_rep_grouped <- group_by(dna_rep_final, PATIENT_ID)  # Group by PALIMP_EQUIVALENT
dna_rep_merged <- summarise(dna_rep_grouped, DNA_Rep = paste(DNA_Repair_Genes, collapse = ",  "))
rm(dna_rep_grouped)
# View(dna_rep_merged) #df 2

# denovo count ~~
denovo_count <- SBS_signatures_exp[["sig_nums"]]
denovo_count$PATIENT_ID <- rownames(denovo_count)
denovo_count <- pivot_longer(denovo_count, cols = c(1:ncol(denovo_count)-1), names_to = "PALIMP_EQUIVALENT", values_to = "Count")
denovo_count <- subset(denovo_count, denovo_count$Count != 0)
denovo_count <- denovo_count[, c(1,2)]
# View(denovo_count)  #df 3
# MERGE
# no need

# compare ~~
names(compare_tab) <- toupper(names(compare_tab))
compare_final <- compare_tab[, c(1, 2, 3)]
compare_final$PALIMP_EQUIVALENT <- gsub("Palimp_", "SBS_denovo_", compare_final$PALIMP_EQUIVALENT)
compare_final$SBS_Corr <- paste(compare_final$REF_SIGNATURE, compare_final$REF_PALIMP_COSINE_SCORE, sep = "; ")
compare_final <- compare_final[, c(2, ncol(compare_final))]
# MERGE
compare_grouped <- group_by(compare_final, PALIMP_EQUIVALENT)  # Group by PALIMP_EQUIVALENT
compare_merged <- summarise(compare_grouped, SBS = paste(SBS_Corr, collapse = ",   "))  # Summarize SBS
rm(compare_grouped)
# View(compare_merged) #df 4

# count ~~
count_df <- vcf_df[, c(1, ncol(vcf_df))]
# MERGE
count_merged <- count_df %>% group_by(Sample) %>% summarise(Mutation_Count = n())
colnames(count_merged)[1] <- "PATIENT_ID"
# View(count_merged) #df 5

# proportion ~~
prop_df <- SBS_signatures_exp[["sig_props"]]
prop_df <- prop_df * 100
prop_df <- round(prop_df, digits = 2)
prop_df$PATIENT_ID <- rownames(prop_df)
prop_df <- pivot_longer(prop_df, cols = c(1:ncol(prop_df)-1), names_to = "PALIMP_EQUIVALENT", values_to = "Proportion")
prop_df <- subset(prop_df, prop_df$Proportion != 0)
prop_df$Proportion <- paste0(prop_df$Proportion, "%")
# View(prop_df) #df6
# MERGE
# merge not needed

#-final merge----
df_pnt_ids_list <- list(count_merged, gene_merged, dna_rep_merged, denovo_count)
df_pnt_ids_merged_df <- reduce( df_pnt_ids_list,full_join, by = "PATIENT_ID")
# View(df_pnt_ids_merged_df)

final_merged <- left_join(df_pnt_ids_merged_df, prop_df, by = c("PALIMP_EQUIVALENT", "PATIENT_ID"))

final_merged <- left_join(final_merged, compare_merged, by = "PALIMP_EQUIVALENT")
final_merged <- final_merged[order(final_merged$Mutation_Count), ]
View(final_merged)

#logging
write_tsv(final_merged, file = paste0(result_dir, "/", "final_merged_", gene_name, "_", subtype, "_.tsv"))
#done
