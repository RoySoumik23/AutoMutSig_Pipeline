# INSTRUCTIONS
# 1. This pipeline is designed for mutation signature analysis using Palimpsest and works best with data from cBioPortal.
# 2. Ensure that mutation and clinical data are based on the **hg19** reference genome, as hg38 may cause errors in annotation steps.
# 3. Prepare mutation data to include a `PATIENT_ID` column that matches with clinical data.
# 4. Clinical data should include `PATIENT_ID` and `SUBTYPE` columns, and all column names should be in uppercase.
# 5. Update all placeholder paths and variables (gene_name, subtype, file paths, etc.) before running the script.
# 6. Run each section sequentially to ensure proper setup and avoid runtime errors.

# Requried libraries
# List of required packages
required_packages <- c(
  "dplyr", 
  "readr", 
  "purrr", 
  "tidyr", 
  "BSgenome.Hsapiens.UCSC.hg19", # "BSgenome.Hsapiens.UCSC.hg38"
  "Palimpsest"
)

# Install any missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)


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

# Source the main analysis
source("https://raw.githubusercontent.com/RoySoumik23/AutoMutSig_Pipeline/refs/heads/main/main_analysis.R")

# For supplementary visualisation
source ("https://raw.githubusercontent.com/RoySoumik23/AutoMutSig_Pipeline/refs/heads/main/supplementary_visualisation.R
