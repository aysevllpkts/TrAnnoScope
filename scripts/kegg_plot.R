# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
file=args[1]
sample_name=args[2]
kegg=args[3]
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
log_message("#------------| Species analysis started |-----------#")

# Log the received arguments
log_message(paste("Annotation File:", file))
log_message(paste("KEGG File:", kegg))
log_message(paste("Sample name:", sample_name))
log_message(paste("Format:", fmt))
log_message(paste("Output directory:", outdir)) 


# Read files
log_message("Reading annotation and KEGG file")
annotation = read_delim(file, delim = "\t", col_names = TRUE)
kegg_df <- read_csv(kegg)

#unique(kegg_df$a_class)


kegg_groups_extra <- c("09100 Metabolism", "09150 Organismal Systems", "09130 Environmental Information Processing", "09140 Cellular Processes", "09120 Genetic Information Processing", "09160 Human Diseases")
kegg_data_extra <- kegg_df %>% filter(a_class %in% kegg_groups_extra)

# Filter and select relevant columns
annot_ko <- annotation %>%
  filter(EggNM.KEGG_ko != ".") %>%
  select(EggNM.KEGG_ko) 
#head(annot_ko)

# Count KO categories
ko_table <- data.frame(category = unlist(strsplit(annot_ko$EggNM.KEGG_ko, ","))) %>% 
  group_by(category) %>% 
  summarise(count = n())
ko_table <- ko_table %>% separate(category, c("x", "reaction_ko"), sep = ":") %>% select(reaction_ko, count)

#head(ko_table)

# Merge and create legend
merged_df <- ko_table %>%
  left_join(kegg_data_extra, by = "reaction_ko") 

#unique(merged_df$a_class)

# Main groups Plot
main_plot <- merged_df %>%
  filter(!is.na(a_class) & !is.na(b_class)) %>%
  group_by(a_class) %>%
  summarize(count = n()) %>%
  ggplot(aes(a_class, count, fill = a_class)) + 
  geom_bar(stat = "identity", position = position_dodge(), colour = "seashell", label = NULL) + # Set label = NULL
  theme_bw() +
  labs(title = paste0(sample_name, " KEGG Distribution"), fill = "KEGG orthology (KO)", x ="KO category", y = "Number of Sequences") +  # Empty x-axis label
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

# Save plot
save_plot(main_plot, paste(outdir, sample_name, "_KEGG_slim_dist.", fmt, sep = ""))


# Subgroups Plot
subgroup_plot <- merged_df %>%
  filter(!is.na(a_class) & !is.na(b_class)) %>%
  group_by(a_class, b_class) %>%
  summarize(count = n()) %>%
  mutate(b_class = fct_reorder(b_class, count)) %>%
  ggplot(aes(b_class, count, fill = a_class)) + 
  geom_bar(stat = "identity", position = position_dodge(), colour = "seashell") +
  theme_bw() +
  labs(title = paste0(sample_name, " KEGG Distribution"), fill = "KEGG orthology (KO)", x = "KEGG category", y = "Number of Sequences") +
  coord_flip() + 
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

# Save plot
save_plot(subgroup_plot, paste(outdir, sample_name, "_KEGG_dist.", fmt, sep = ""))

log_message("#------------| KEGG analysis finished |------------#")