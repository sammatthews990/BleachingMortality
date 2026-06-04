# Extract and run code up to scale_params definition
code_lines <- readLines("scratch/temp_code.R")

# Run up to line 1170 to execute data prep and scale_params definition
run_lines <- code_lines[1:1170]
run_lines <- c(run_lines, "
print('=== scale_params ===')
print(scale_params)
")

writeLines(run_lines, "scratch/temp_run.R")
source("scratch/temp_run.R")
