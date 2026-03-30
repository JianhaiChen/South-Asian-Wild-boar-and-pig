#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

cat("\n============================================================\n")
cat(" 🧬 True XCI Skew Background Extractor (Filtering 0 and Inf)\n")
cat("============================================================\n")

# 1. Load data
df <- fread("uncorrected_chrX_1_125M_ratio_merged1.csv", data.table = FALSE)
ratio_cols <- grep("^ratio_", names(df), value = TRUE)

# 2. Define the "real median" function
get_real_median_af <- function(row_vec) {
  vals <- as.numeric(row_vec)
  # Key filter: exclude 0 (no difference or missing), Inf (monoallelic expression), and NA
  clean_vals <- vals[vals != 0 & !is.infinite(vals) & !is.na(vals)]
  
  if(length(clean_vals) < 5) return(0.5) # If there are too few valid genes, fallback to 0.5
  
  # Calculate the median AF of these valid genes
  af_values <- (2^clean_vals) / (1 + 2^clean_vals)
  return(median(af_values))
}

# 3. Calculate per sample
cat(">>> Filtering 67% of zero values, calculating based on valid signals only...\n")
results <- data.frame(
  sample = df$sample,
  Background_XCI_Skew = apply(df[, ratio_cols], 1, get_real_median_af)
)

# 4. Print the real distribution
cat("\n--- [Filtered Real XCI Distribution Statistics] ---\n")
print(summary(results$Background_XCI_Skew))

# Find the truly extreme (most skewed) samples
results$deviation <- abs(results$Background_XCI_Skew - 0.5)
top_skewed <- results %>% arrange(desc(deviation)) %>% head(5)

cat("\n--- [True values of the top 5 most skewed samples] ---\n")
print(top_skewed[, c("sample", "Background_XCI_Skew")])

# 5. Save results
final_input <- df %>% 
  select(-any_of(ratio_cols)) %>%
  inner_join(results %>% select(sample, Background_XCI_Skew), by = "sample") %>%
  mutate(Skew_State = ifelse(Background_XCI_Skew <= 0.25 | Background_XCI_Skew >= 0.75, 1, 0))

write.csv(final_input, "Strict_LMM_Input.csv", row.names = FALSE)
cat("\n✅ True background covariates generated: Strict_LMM_Input.csv\n")