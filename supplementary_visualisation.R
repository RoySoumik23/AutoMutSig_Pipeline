# Instructions for Users:
# 1. `result_parent_dir` will be fetched once 'palimpsest_pipeline.R' is executed
          # =OR= Specify the file location of the parent directory of all the outputs (e.g., "/home/shivalik2/Soumik/my_codes/Results/Discovery/")
# 2. Set `cancer_type` to the specific folder name for the cancer type of interest.
# 3. Ensure your data files are properly formatted (i.e., tab-separated or CSV) as expected by the script.
# 4. Run the script to generate the heatmap and bar plot PDF outputs.
# 5. Customize the break values, colors, and plot dimensions as needed.
# 6. The generated PDFs will be saved in the specified output directory (`SBS_denovo` folder).

# Source the libraries setup
source("/home/shivalik/Downloads/R codes/libraries_and_setup.R")

result_parent_dir <- "/home/shivalik2/Soumik/my_codes/Results/Discovery/"
file_type <- "SBS_denovo/Comparison_table.txt" # The end location of the comparison file created by the pipeline

# CHANGE THESE VARIABLES
cancer_type <- "^ALKBH_BRCA_Basal_Palimpsest/" # Specify the cancer type folder (e.g., "Non^ALKBH_BRCA_Basal_Palimpsest/")
final_merged_file <- "final_merged_^ALKBH_BRCA_Basal_.tsv" # Specify the final merged folder (e.g., "final_merged_^ALKBH_BRCA_Basal_.tsv")

# Read the data
cancer_df <- read.delim2(paste0(result_parent_dir, cancer_type, file_type), sep = " ")

# Clean and transform the data
df <- cancer_df[, -4]
# Pivot the data wider
wide_data <- pivot_wider(df, names_from = Ref_Signature, values_from = Ref_Palimp_cosine_score)
# Remove the null column(s)
wide_data <- wide_data[, colnames(wide_data) != "NA"]
# Convert all columns (except the first one, which contains row names) to numeric
wide_data[, -1] <- lapply(wide_data[, -1], function(x) as.numeric(as.character(x)))
wide_data[is.na(wide_data)] <- 0
numeric_matrix <- as.matrix(wide_data[, -1])
rownames(numeric_matrix) <- gsub("^Palimp", "SBS_denovo", wide_data$Palimp_Equivalent)

# VISUALISATION
# Output directory
output_dir <- paste0(result_parent_dir, cancer_type, "SBS_denovo");
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# HEATMAP
# Custom breaks for heatmap
breaks <- c(seq(0.0, 0.30, length.out = 10), seq(0.31, 0.60, length.out = 10), seq(0.61, 1, length.out = 30))

# Generate and save heatmap
pdf(file.path(output_dir, "Correlation_Heatmap.pdf"), width = 11, height = 10)
pheatmap(numeric_matrix, color = colorRampPalette(c("white", "orange", "red", "darkred"))(length(breaks) - 1), breaks = breaks,
         main = "Cosine Similarity Heatmap", cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10)
dev.off()

# BAR PLOT
# Process the merged data
df_merged <- read.delim2(paste0(result_parent_dir, cancer_type, final_merged_file), sep = "\t")
df_merged <- df_merged[, c(5, 2, 6)]
df_merged$Proportion <- as.numeric(gsub("%", "", df_merged$Proportion)) / 100
df_merged[, -1] <- data.frame(lapply(df_merged[, -1], as.numeric))
df_merged$count <- df_merged$Mutation_Count * df_merged$Proportion

# Summarize and create bar plot
totals <- df_merged %>%
  group_by(PALIMP_EQUIVALENT) %>%
  summarise(Total = sum(count))

bar_plot <- ggplot(totals, aes(x = PALIMP_EQUIVALENT, y = Total, fill = PALIMP_EQUIVALENT)) +
  geom_bar(stat = "identity") + 
  labs(title = "", x = "Signature", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.title = element_blank())

# Save the bar plot
ggsave(paste0(result_parent_dir, cancer_type, "SBS_denovo/bar_plot_sbs_denovo.pdf"), plot = bar_plot, width = 8, height = 6)
