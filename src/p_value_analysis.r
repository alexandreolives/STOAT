library(ggplot2)

output_dir <- "output/"  # Replace with your desired directory
file_path <- "output/quantitative_output.tsv"  # Replace with the path to your file
data <- read.table(file_path, header = TRUE, sep = "\t")
data$P_value <- as.numeric(data$P_value) 

# --- P-value distribution plot ---
p_value_plot <- ggplot(data, aes(x = P_value)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "blue") +
  theme_minimal() +
  xlab("P-Value") +
  ylab("Frequency") +
  ggtitle("Distribution of P-values")

ggsave(filename = file.path(output_dir, "p_value_distribution.png"), plot = p_value_plot, width = 8, height = 6)
data$CHR <- as.numeric(as.factor(data$Snarl))
data$BP <- seq(1, nrow(data))
