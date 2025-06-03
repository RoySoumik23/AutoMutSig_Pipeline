# Install and load necessary packages
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr")
library(ggplot2)
library(ggpubr)

df <- data.frame(denovo_signatures)

# Step 2: Set rownames and remove ID column
rownames(df) <- df$Signatures
df$Signatures <- NULL

# Step 3: Split into groups
group1 <- df[, 1:48]
group2 <- df[, 49:ncol(df)]

# Step 4: Flatten values for unpaired test
group1_values <- unlist(group1)
group2_values <- unlist(group2)

# Step 5: Perform t-test
test_result <- t.test(group1_values, group2_values, alternative = "two.sided")

# Step 6: Create data frame for plotting
plot_df <- data.frame(
  Value = c(group1_values, group2_values),
  Group = factor(c(rep("C > A/G/T", length(group1_values)),
                   rep("T > A/C/G", length(group2_values))))
)

# Step 7: Save as JPEG with comparison bar
jpeg("codon_substitution_test_2.jpg", width = 3000, height = 2000, res = 300)

p <- ggplot(plot_df, aes(x = Group, y = Value)) +
  geom_boxplot(fill = c("#66c2a5", "#fc8d62")) +
  theme_minimal(base_size = 14) +
  ggtitle("Denovo Signature: Non-TNBC") +
  stat_compare_means(method = "t.test", label.y = max(plot_df$Value) * 1.1)

print(p)
dev.off()
