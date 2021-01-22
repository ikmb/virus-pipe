#!/usr/bin/env Rscript
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# Load data
cov <- as.tbl(read.table(args[1])) %>% 
	rename("Position" = V2) %>% rename("Coverage" = V3)
limit <- args[2]

# Plot
pdf(args[3] + ".pdf")
cov %>% select(Position, Coverage) %>% 
	ggplot(aes(Position, Coverage)) + 
	geom_area() + ggtitle(args[1]) + ylim(0,200)

jpg(args[3] + ".jpg")
cov %>% select(Position, Coverage) %>%
        ggplot(aes(Position, Coverage)) +
        geom_area() + ggtitle(args[1]) + ylim(0,200)


dev.off()
