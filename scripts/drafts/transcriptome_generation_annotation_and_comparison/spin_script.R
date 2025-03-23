#!/usr/bin/env Rscript
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a script name is provided
if (length(args) == 0) {
  stop("Please provide the script name as a command-line argument.")
}

script_name <- args[1]

# Run knitr::spin on the provided script name
cat("Running knitr::spin on:", script_name, "\n")

tryCatch({
  knitr::spin(script_name)
  cat("Successfully processed", script_name, "\n")
}, error = function(e) {
  cat("Error processing script:", e$message, "\n")
  quit(status = 1)
})
