args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]
png_file <- args[2]
log_file <- args[3]


log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

library(Biostrings)

tryCatch({
  # Read the fasta file
  reads <- readAAStringSet(fasta_file)
  read_lengths <- width(reads)

  # Calculate N50
  sorted_lengths <- sort(read_lengths, decreasing = TRUE)
  cumsum_lengths <- cumsum(sorted_lengths)
  total_length <- sum(sorted_lengths)
  n50 <- sorted_lengths[which(cumsum_lengths >= total_length * 0.5)[1]]

  # Plot histogram
  png(png_file)
  hist(read_lengths, breaks=50, main="Read Length Distribution", xlab="Read Length", ylab="Number of Reads", col="lightblue", border="black")
  abline(v=n50, col="black", lwd=2)
  dev.off()
  
  cat("Plot created successfully.\n")
}, error = function(e) {
  cat("Error occurred: ", e$message, "\n")
}, finally = {
  sink(type = "message")
  close(log_con)
})



