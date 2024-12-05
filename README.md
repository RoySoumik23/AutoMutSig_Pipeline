# Palimpsest Automation Pipeline  

## Overview  
This **Mutation Signature Analysis Automation Pipeline** is designed to streamline the preparation, execution, and analysis of mutation signatures using the **Palimpsest R package**. It is tailored for cancer genomics research, automating tasks from data preprocessing to de novo mutation signature extraction and comparison with **COSMIC (Catalogue Of Somatic Mutations In Cancer)** signatures.  

The pipeline processes mutation and clinical data, aligning and cleaning the inputs before performing a detailed analysis. It enables visualization and interpretation of mutation signatures with a focus on specific genes and cancer subtypes. By automating repetitive tasks, this pipeline accelerates research into the mutational processes underlying cancer.  

---

## Features  

### Data Preparation  
- Processes mutation data to ensure compatibility with Palimpsest.  
- Filters clinical data for specific cancer subtypes.  

### Patient Filtering  
- Filters patients based on cancer subtype or specific gene mutations.  

### De Novo Signature Analysis  
- Extracts mutation signatures using **Non-Negative Matrix Factorization (NMF)**.  
- Determines proportions of de novo signatures.  

### COSMIC Signature Comparison  
- Compares extracted signatures to known COSMIC SBS signatures.  
- Visualizes similarity using **cosine similarity heatmaps**.  

### Visualization  
- Generates pie charts to display signature contributions.  
- Produces signature content plots for enhanced interpretation.  

### DNA Repair Gene Analysis  
- Summarizes mutations in DNA repair genes for selected patients.  

### Patient-Level Summary  
- Generates detailed reports for individual patients, including:  
  - Mutation counts.  
  - Key gene mutations.  
  - Signature contributions.  

---

## Prerequisites  

### Software Requirements  
- **R version 4.0 or higher**  

### Required R Packages  
- **Palimpsest**  
- **BSgenome.Hsapiens.UCSC.hg19** (or hg38 if applicable)  
- `dplyr`, `readr`, `purrr`, `ggplot2`, `tidyr`  

### Input Files  
- **Mutation Data**: CSV/TSV file with columns like `Hugo_Symbol`, `Variant_Classification`, `Variant_Type`, and `PATIENT_ID` (if not present, add it during data cleaning step).  
- **Clinical Data**: CSV/TSV file with columns like `PATIENT_ID` and `SUBTYPE`.  

---

## Installation  

### Clone the Repository  
```bash
git clone https://github.com/RoySoumik23/AutoMutSig_Pipeline.git
cd AutoMutSig_Pipeline
```  

### Install Required R Packages  
```R
install.packages(c("Palimpsest", "BSgenome.Hsapiens.UCSC.hg19", "dplyr", "readr", "purrr", "ggplot2", "tidyr"))  
```  

---

## Usage  

### Step 1: Configure Variables  
Update the following variables in the script:  
- **`print_statement`**: Message describing the unique patient count.  
- **`gene_name`**: Name or pattern of the target gene (e.g., `"ALKBH"`).  
- **`subtype`**: Cancer subtype for analysis (e.g., `"TNBC"`).  
- **`mut_location`**: File path to the mutation data.  
- **`clinical_location`**: File path to the clinical data.  
- **`num_of_denovo_sings`**: Number of de novo signatures to extract.  
- **`parent_subfile_name`**: Name of the subfile for this analysis.  
- **`Parent_file_name`**: Location of the parent file.  

### Step 2: Execute the Pipeline  
Run the script in R:  
```R
source("palimpsest_pipeline.R")  
```  

---

## Outputs  

1. **Summary File**: `details.txt` summarizing patient counts at each stage.  
2. **VCF Files**: Mutation data formatted for Palimpsest analysis.  
3. **Comparison Table**: Cosine similarity scores between de novo and COSMIC signatures.  
4. **Visualizations**:  
   - Heatmaps of cosine similarities.  
   - Pie charts showing SBS signature contributions (`all_known_signature_piecharts.pdf`).  
5. **Final Merged Table**: Contains:  
   - Patient IDs.  
   - SBS signatures (de novo and known COSMIC).  
   - Correlation and gene mutation types.  
6. **De Novo SBS Proportions**: `denovo_SBS_proportion.csv` detailing the percentage contribution of each of the 96 SBS conversions.  

---

## Customization  

- Ensure that column names in mutation and clinical data files match the script’s requirements.  
- Modify filtering logic for datasets with unique requirements (e.g., additional subtypes or clinical conditions).  

---

## Known Issues  

1. **Compatibility with hg38**:  
   - If using `BSgenome.Hsapiens.UCSC.hg38`, ensure the mutation data aligns with genome boundaries.  

2. **Pie Chart Colors**:  
   - Colors in pie charts may not render as intended. Adjust the color mapping manually if needed.  

---
## Visualisation
1. This is a MS PowerBI dashboard for palimpsest analysis of ALKBH gene mutation in Triple Negative Breast Cancer (TNBC) cancer. (Need to be made separately)
<img width="722" alt="image" src="https://github.com/user-attachments/assets/0d1c57ea-4bf8-4a66-ae16-12d68685bc7f">


## Additional Instructions  

1. **Data Format Requirements**:  
   - Mutation data must include `Hugo_Symbol`, `Variant_Classification`, `Variant_Type`, and `PATIENT_ID`.  
   - Clinical data must include `PATIENT_ID` and `SUBTYPE`.  

2. **File Paths**:  
   - Use absolute file paths for `mut_location` and `clinical_location` to avoid path-related issues.  

3. **Handling Missing Data**:  
   - Clean datasets containing `NA` values using `na.omit()` or similar methods.  

4. **Runtime Considerations**:  
   - Use a machine with at least **8 GB of RAM** for optimal performance.  

5. **Testing with Sample Data**:  
   - Validate the workflow using a small, curated dataset before running it on large datasets.  

6. **Version Control**:  
   - Track changes to the pipeline using Git for reproducibility.  

---

## Future Improvements  

- Automating input file validation.  
- Adding advanced visualization options.  
- Enhancing compatibility with hg38.  

---

## Contributions  

Contributions are welcome! Please fork the repository, make changes, and submit a pull request. Bug reports and feature suggestions are encouraged.  

---

## Acknowledgements  

Special thanks to the developers of:  
- The **Palimpsest R package** for mutation signature analysis.  
- The **BSgenome project** for genomic reference data integration.  
