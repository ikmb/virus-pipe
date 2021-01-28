#!/usr/bin/env Rscript

# Import library
library(tidyverse)
# Import argument. Set variables
args = commandArgs(trailingOnly=TRUE) 
fin <- args[1] 
basename <- args[2] 
fout <- paste0(basename, ".jpg")

#fin <- "/home/lsilva/IKMB/projects/random/covid.coverage.plot/21Ord340-L1_S1_L001.coverage.samtools.txt" basename <- "Patient0"
# Load data
cov <-
  fin %>%
  read.delim(stringsAsFactors = F, header = F) %>%
  rename("Position" = "V2") %>%
  rename("Coverage" = "V3")
# Plot
jpeg(filename = fout, units = "cm", height = 7.5, width = 18, res = 300) 
cov %>%
  ggplot(aes(x = Position, y = Coverage)) +
  geom_area() +
  scale_y_continuous(limits = c(0,200)) +
  scale_x_continuous(limits = c(0,30000)) +
  theme(text = element_text(size = 12,family="sans"),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "lightgrey"),
        axis.line = element_line(colour = "black")) +
  labs(title = basename)
  
dev.off()
