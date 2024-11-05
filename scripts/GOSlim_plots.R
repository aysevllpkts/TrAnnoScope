# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
file=args[1]
sample_name=args[2]
fmt=args[3]
outdir=args[4]

# Load necessary libraries
library(tidyverse)
library(ggthemes)
library(svglite)

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
log_message("#---------| Starting GOSlim analysis |----------#")

# Log the received arguments
log_message(paste("File:", file))
log_message(paste("Sample name:", sample_name))
log_message(paste("Format:", fmt))
log_message(paste("Output directory:", outdir))

# Upload file
log_message("Reading GOSlims data file...")
dataGOSlims <- read.delim(file, header = FALSE, sep = "\t")

# Tidy dataset
log_message("Tidying dataset")
dataGOSlims <- dataGOSlims %>%
  select(V1, V3, V4) %>%
  filter(!(V3 == V1)) %>% # Remove sub-terms with the same name as main group
  mutate(V1 = factor(V1, levels = c("cellular_component", "molecular_function", "biological_process"),
                    labels = c("Cellular Component", "Molecular Function", "Biological Process"))) %>%
  mutate(V3 = str_to_title(V3)) %>%
  group_by(V1) %>%
  slice_max(order_by = V4, n = 20) # Select top 20 terms for each main group

names(dataGOSlims) <- c("category", "goslim", "count")
head(dataGOSlims)

# Create plots
log_message("Creating and saving plots")
custom_colors <- c("green", "blue", "red")

pALL <- dataGOSlims %>%
  ggplot(aes(x = reorder(goslim, count), y = count, fill = category)) +
  facet_wrap(~category, scales = "free", ncol = 3) +
  scale_fill_manual(values = custom_colors[1:length(unique(dataGOSlims$category))]) +
  coord_flip() +
  geom_bar(stat = "identity", width=.5) +
  labs(x="GOSlim Classification", y="Number of Sequences") +
  theme(axis.text=element_text(size=10), text = element_text(size = 15)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=10, angle=0),
        axis.title=element_text(size=15, face="bold")) +
  ggtitle(paste(sample_name, "GOSlims Distribution", sep=" ")) +
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"),
        legend.position = "none")

# Save
save_plot(pALL, paste(outdir, sample_name, "_ALL_GOSlims.", fmt, sep = ""))

# Plot for each category
plot_category <- function(category_name, color, filename_suffix) {
  plot <- dataGOSlims %>% filter(category == category_name) %>%
    ggplot(aes(x = reorder(goslim, count), y = count)) +
    geom_bar(stat = "identity", fill = color, width = .5) +
    coord_flip() + labs(x = "Classification", y = "Number of Sequences") +
    theme_bw() +
    ggtitle(paste(sample_name, category_name, "GOslims", sep=" ")) +
    theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))
  
  save_plot(plot, paste(outdir, sample_name, filename_suffix, fmt, sep = ""))
}

plot_category("Cellular Component", "green", "_Cellular_Component_GOSlims.")
plot_category("Molecular Function", "blue", "_Molecular_Function_GOSlims.")
plot_category("Biological Process", "red", "_Biological_Process_GOSlims.")

# Log the end of the script
log_message("#-----------| GOSlims analysis finished |-----------#")
