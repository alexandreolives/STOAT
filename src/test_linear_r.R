# Load necessary libraries
args <- commandArgs(trailingOnly = TRUE)

# Function to fit an OLS model
sm_ols <- function(x, y) {
  print(x)
  print(y)

  # Perform OLS regression (y ~ x)
  model <- lm(y ~ x)
  
  # Print the summary of the regression model
  print(summary(model))
  
  # Extract p-values from the model
  p_values <- summary(model)$coefficients[, 4]
  formatted_p_values <- ifelse(is.na(p_values), "N/A", p_values)
  
  return(formatted_p_values)
}

# Function to parse phenotype file and reshape y as a column vector (matrix)
parse_pheno_file <- function(pheno_file_path) {
  pheno_data <- as.matrix(read.table(pheno_file_path, sep = " ", header = TRUE, colClasses = c("NULL", "numeric")))
  pheno_data <- matrix(pheno_data, ncol = 1)
  return(pheno_data)
}

# Function to parse data file (independent variables)
parse_data_file <- function(data_file_path) {
  data <- as.matrix(read.table(data_file_path, sep = " ", header = TRUE, colClasses = c("NULL", "numeric", "numeric")))
  return(data)
}

# Main function to run OLS regression
main <- function(pheno_file_path, data_file_path) {
  y <- parse_pheno_file(pheno_file_path)
  x <- parse_data_file(data_file_path)
  
  # Perform OLS regression and get p-values
  p_values <- sm_ols(x, y)
  print("P-values:")
  print(p_values)
}

# Get the command-line arguments for file paths
if (!interactive()) {
  main(args[1], args[2])
}

# Rscript src/test_linear_r.R pheno.txt table.twt
