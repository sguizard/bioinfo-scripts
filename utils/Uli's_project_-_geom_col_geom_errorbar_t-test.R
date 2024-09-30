# Data (ZNF141.csv)
# ,MOLM13,MOLM13,MOLM13,NOMO-1,NOMO-1,NOMO-1,THP-1,THP-1,THP-1,MV4-11,MV4-11,MV4-11
# ZNF141,0.95847,1.04124,0.62779,4.92125,5.84829,6.73001,24.32614,24.24827,26.86165,13.31491,17.47254,18.66192

# Libraries
library(tidyverse)
library(ggpubr)

# Set colors
colors <- c(
  "MOLM13" = "#ffb74d", 
  "NOMO-1" = "#64b5f6", 
  "THP-1"  = "#bdbdbd", 
  "MV4-11" = "#aed581")

# Load csv + pivot wider
data_long <-
  read_csv("ZNF141.csv") %>%
  rename(Gene = ...1) %>%
  pivot_longer(
    -Gene, 
    names_to  = "CellLine", 
    values_to = "TPM") %>%
  mutate(CellLine = factor(
    str_replace(CellLine, "\\.{3}\\d+$", ""), 
    levels = c("MOLM13", "NOMO-1", "THP-1", "MV4-11")))

# Summarize data
data_summary <- 
  data_long %>%
  group_by(CellLine) %>%
  summarise(
    mean_TPM = mean(TPM), 
    sd_TPM   = sd(TPM), 
    ymin     = mean_TPM - sd_TPM,
    ymax     = mean_TPM + sd_TPM)

  data_long2 <- left_join(data_long, data_summary)

# Plot col/bar
data_long2 %>% 
  ggplot(aes(x = CellLine, y = TPM)) +
  geom_col(
    mapping = aes(
      CellLine, 
      mean_TPM, 
      fill = CellLine), 
    data = data_summary) + 
  geom_errorbar(
      mapping = aes(
        ymin = ymin, 
        ymax = ymax), 
      width = 0.2) +
  geom_jitter(
    width = 0.2,
    size  = 2) + 
  stat_compare_means(
    comparisons = list(
      c("MOLM13", "MV4-11"), 
      c("NOMO-1", "MV4-11"), 
      c("THP-1", "MV4-11")), 
    method = "t.test", 
    method.args = list(var.equal = TRUE), # Welch's t-test
    label       = "p.signif", 
    size        = 5, 
    label.y     = c(30, 32, 34), # Adjust p-value star positioning
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
      symbols   = c("***", "**", "*", "ns"))) + 
  scale_fill_manual(values = colors) +  # Use the defined color palette
  labs(
    title = "ZNF141", 
    y     = "TPM", 
    x     = "") +
  theme_classic() +
  theme(
    legend.position = "none", 
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 12), 
    axis.title.y = element_text(size = 14))
