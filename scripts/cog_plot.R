# Get the command line arguments
args = commandArgs(trailingOnly=TRUE)
annotation=args[1]
cog=args[2]
sample_name=args[3]
fmt=args[4]
outdir=args[5]

# Load necessary libraries
library(tidyverse)
library(ggthemes)
library(svglite)

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
log_message("#------------| COG analysis started |-----------#")

# Log the received arguments
log_message(paste("Annotation File:", annotation))
log_message(paste("COG File:", cog))
log_message(paste("Sample name:", sample_name))
log_message(paste("Format:", fmt))
log_message(paste("Output directory:", outdir))

# 
log_message("Reading annotation and cog files")
annotation = read_delim(annotation, delim = "\t", col_names = TRUE) %>% rename(COG = EggNM.COG_category)
cog = read_delim(cog, delim = "\t", col_names = c("category", "X2", "description")) %>%  select(category, description)

# Filter and select relevant columns
annot_cog <- annotation %>%
  filter(COG != ".") %>%
  select(COG)

# Count COG categories
cog_table <- data.frame(category = unlist(strsplit(annot_cog$COG, ""))) %>% 
  group_by(category) %>% 
  summarise(count = n())

# Merge and create legend
merged_df <- cog %>%
  left_join(cog_table, by = "category") %>%
  mutate(legend = paste(category, description, sep = ": "))

# Plot
cog_plot <- merged_df %>%
  ggplot(aes(category, count, fill = legend)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "seashell") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  labs(title = paste0(sample_name, " COG Distribution"), x = "COG category", y = "Number of Sequences") +
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

# Save 
save_plot(cog_plot, paste(outdir, sample_name, "_COG_dist.", fmt, sep = ""))

log_message("#------------| COG analysis finished |------------#")