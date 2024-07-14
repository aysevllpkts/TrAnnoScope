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
log_message("#---------| Starting GOSlim analysıs |----------#")

# Log the received arguments
log_message(paste("File:", file))
log_message(paste("Sample name:", sample_name))
log_message(paste("Format:", fmt))
log_message(paste("Output directory:", outdir))

# Upload file
log_message("Reading GOSlims data file...")
dataGOSlims <- read.delim(file, header = FALSE, sep = "\t")

# tidy dataset
log_message("Tidying dataset")
dataGOSlims <- dataGOSlims %>% 
  select(V1, V3, V4) %>%
  mutate(V1 = factor(V1, levels = c("cellular_component", "molecular_function", "biological_process"),
                             labels = c("Cellular Component", "Molecular Function", "Biological Process"))) %>%
  mutate(V3 = str_to_title(V3))
names(dataGOSlims) <- c("category", "goslim", "count")
head(dataGOSlims)


# Create plots
log_message("Creating and saving plots")

custom_colors <- c("green", "blue", "red")

pALL <- dataGOSlims %>% filter(count != 0) %>% 
  ggplot(aes(x = reorder(goslim,count), y = count, fill = category)) +
  facet_wrap(~category, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = custom_colors[1:length(unique(dataGOSlims$category))]) + 
  coord_flip() +
  geom_bar(stat = "identity", width=.5) +  # Create bar plots +
  labs(x="GOSlim Classification",y="Number of Sequences")+
  theme(axis.text=element_text(size=10))+theme(text = element_text(size = 15))+
  theme_bw() +
  theme(axis.text.x=element_text(size=10,angle=0))+theme(axis.title=element_text(size=15,face="bold"))+
  ggtitle(paste(sample_name,"GOSlims Distribution",sep=" "))+
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"), legend.position = "none")

# Save
save_plot(pALL, paste(outdir, sample_name, "_ALL_GOSlims.", fmt, sep = ""))


pCC <- dataGOSlims %>% filter(count != 0 & category == 'Cellular Component') %>%
    ggplot(aes(x=reorder(goslim,count), y=count))+
    geom_bar(stat="identity", fill="green", width=.5)+
    coord_flip()+labs(x="Classification",y="Number of Sequences")+
    theme_bw() +
    ggtitle(paste(sample_name,"Cellular Component GOslims",sep=" "))+
    theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

# Save
save_plot(pCC, paste(outdir, sample_name, "_Cellular_Component_GOSlims.", fmt, sep = ""))



pMF <- dataGOSlims %>% filter(count != 0 & category == 'Molecular Function')  %>%
    ggplot(aes(x=reorder(goslim,count), y=count))+
    geom_bar(stat="identity", fill="blue", width=.5)+
    theme_bw() +
    coord_flip()+labs(x="Classification",y="Number of Sequences")+
    ggtitle(paste(sample_name,"Molecular Function GOSlims",sep=" "))+
    theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

# Save
save_plot(pMF, paste(outdir, sample_name, "_Molecular_Function_GOSlims.", fmt, sep = ""))


pBP <- dataGOSlims %>% filter(count != 0 & category == 'Biological Process')  %>%
    ggplot(aes(x=reorder(goslim,count), y=count))+
    geom_bar(stat="identity", fill="red", width=.5)+
    theme_bw() +
    coord_flip()+labs(x="Classification",y="Number of Sequences")+
    ggtitle(paste(sample_name,"Biological Process GOSlims",sep=" "))+
    theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

# Save
save_plot(pBP, paste(outdir, sample_name, "_Biological_Process_GOSlims.", fmt, sep = ""))

# Log the end of the script
log_message("#-----------| GOSlims analysis finished |-----------#")