args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]
output_dir <- args[2]
log_file <- args[3]

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

library(Biostrings)
library(dplyr)
library(ggplot2)

tryCatch({
  # Read the fasta file
  reads <- readAAStringSet(fasta_file)
  
  # Extract sequence names and types
  extract_info <- function(header) {
    # Extract the sequence name
    seq_name <- strsplit(header, " ")[[1]][1]
    # Extract the type
    type_detail <- sub(".*type:([^ ]+).*", "\\1", header)
    type_info <- sub(".*type:([^ -]+).*", "\\1", header)
    return(data.frame(SequenceName = seq_name, Type = type_info, Type_detailed = type_detail, stringsAsFactors = FALSE))
  }
  
  # Apply the extraction function to all headers
  info_list <- lapply(names(reads), extract_info)
  info_df <- bind_rows(info_list)
  
  # Generate Plots 
  types <- ggplot(info_df, aes(x = Type, fill=Type)) +
              geom_histogram(stat = "count", color = "black", binwidth = 50) +
              labs(title = "Distribution of sequence types", x = "Sequence Type", y = "Number of Reads") +
              theme_minimal()

  types_detailed <- ggplot(info_df, aes(x = Type, fill=Type_detailed)) +
                      geom_histogram(stat = "count", color = "black") +
                      labs(title = "Distribution of sequence types (detailed)", x = "Sequence Type (detailed)", y = "Number of Reads") +
                      theme_minimal() 
  
  # Save plots as PNG files
  png(paste0(output_dir, "/protein_completeness_plot.png"), width = 800, height = 600, res = 120)
  print(types)
  dev.off()
  
  png(paste0(output_dir, "/protein_completeness_plot_detailed.png"), width = 800, height = 600, res = 120)
  print(types_detailed)
  dev.off()
  
  # Write the dataframe to a CSV file
  write.csv(info_df, paste0(output_dir, "/protein_completeness.csv"), row.names = FALSE)
 
  cat("Plot created successfully.\n")
}, error = function(e) {
  cat("Error occurred: ", e$message, "\n")
}, finally = {
  sink(type = "message")
  close(log_con)
})

