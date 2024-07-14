# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
outfmt=args[1]
sample_name=args[2]
fmt=args[3]
outdir=args[4]

# Load required packages
library(tidyverse)

# Required functions
log_message <- function(message) {
  cat(sprintf("[%s %s\n]", Sys.time(), message))
}

save_plot <- function(plot, filename) {
  if (fmt == "pdf") {
    ggsave(filename = filename, plot = plot, width = 15, height = 7)
  } else {
    ggsave(filename = filename, plot = plot, width = 15, height = 7)
  }
  log_message(paste("Saved plot:", filename))
}

# Log the start of the script
log_message("#------------| Species analysis started |-----------#")

# Log the received arguments
log_message(paste("NR Blastx OUTFMT File:", outfmt))
log_message(paste("Sample name:", sample_name))
log_message(paste("Format:", fmt))
log_message(paste("Output directory:", outdir)) 


# Read the data into a data frame
log_message("Reading outfmt file")
data <- read.table(outfmt, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Extract the last column which contains the species information
last_column <- data[, ncol(data)]

# Extract species names using a regular expression
species <- sub(".*\\[(.*)\\]", "\\1", last_column)

# Create a data frame with species counts
species_count <- data.frame(table(species))
log_message("Top 10 species extracted from outfmt")
top10_blastx_species <- species_count %>% arrange(desc(Freq)) %>% head(., 10)


# Plot the histogram
species_plot <- top10_blastx_species %>% 
  ggplot(aes(reorder(species, -Freq), Freq)) + 
  geom_bar(stat = "identity", colour = "black", fill = "firebrick" ,width=.8) +
  theme_bw() +
  labs(title = paste0(sample_name, " Top 10 Species NR BLASTX"), x = "Species", y = "Number of Hits") +
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1))

# Save 
save_plot(species_plot, paste(outdir, sample_name, "_top10_species_dist_NR_blastx.", fmt, sep = ""))
#save_plot(cog_plot, paste(outdir, sample_name, "_COG_dist.", fmt, sep = ""))


log_message("#------------| Species analysis finished |------------#")