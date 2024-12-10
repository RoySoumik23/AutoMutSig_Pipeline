# List of required packages
required_packages <- c(
  "dplyr", 
  "readr", 
  "purrr", 
  "tidyr", 
  "BSgenome.Hsapiens.UCSC.hg19", # "BSgenome.Hsapiens.UCSC.hg38"
  "Palimpsest",
  "tidyr", 
  "pheatmap", 
  "tidyverse", 
  "corrgram"
)

# Install any missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)