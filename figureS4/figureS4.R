
rm(list = ls())
library(tidyverse)
library(stats)
library(multcomp)

# load New Gene Name file
Gene_name <- read_tsv("../PF_GeneName.tsv")

tajima <- read_tsv("data/CountryPop_tajD.txt") %>%
    inner_join(Gene_name) %>%
    separate(`Population Gene`, c("Population", "Gene"), " ") %>%
    relocate(Chr, .after = Chrom) %>% 
    dplyr::select(-c(11,12))


tajima <- tajima %>% 
    group_by(Population) %>% 
    group_by(Chr) %>% 
    summarise(chr_len = max(Start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(tajima, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, Start) %>%
    mutate( BPcum = Start + tot) %>% 
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(Tajima >= 1, "yes", "no")) %>%
    mutate( is_annotate=ifelse(Tajima >= 1, "yes", "no"))

axis_set <- tajima %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

ylimits <- c(floor(min(tajima$Tajima)), abs(floor(max(tajima$Tajima))) + 2)

colors <- c("gray30", "firebrick") 
# Plot
ggplot(tajima, aes(x = BPcum, y = Tajima, color = as_factor(Chr))) +
    
    # Show all points
    geom_point(size = 0.8) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center, expand = c(0.02,0.05)) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +
    scale_color_manual(values = rep(colors, unique(length(axis_set$Chr)))) +   
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Tajima's D") +
    
    # Custom the theme:
    theme_bw() +
    theme( 
        legend.position = "none",
        axis.line = element_line(color = "black"),
        axis.text = element_text(face = "bold", color = "black", size = 6),
        axis.title.x = element_text(vjust = 0, size = 9),
        axis.title.y = element_text(vjust = 2, size = 9),
        axis.title = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold.italic", hjust = 0.5, size = 10)) +
    facet_wrap(~factor(Population)) 

ggsave("figureS4/Figure.S4.jpeg", width = 30, height = 20, units = "cm")

#=======================
# Test of significance
#=======================
# convert Pop to a factor
tajima$Population <- as.factor(tajima$Population)

# Multiple comparison test
# ANOVA
model <- lm(Tajima~Population, data=tajima)
anova(model)

summary(glht(model, mcp(Population = "Tukey")))
