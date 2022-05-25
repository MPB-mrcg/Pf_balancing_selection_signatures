
rm(list = ls())
library(tidyverse)
library(stats)
library(multcomp)

# load New Gene Name file
Gene_name <- read_tsv("../PF_GeneName.tsv")

# 1. Tajima's D
tajima <- "../Analysis/Selection_Analysis/Selection_indices/FuLi/Population_clusters/combined.txt"

# Format the output file
tajima <- read_tsv(tajima) %>%
    inner_join(Gene_name) %>% 
    select(-c(2,9,11,13)) %>% 
    relocate(Pop, .before = Gene) %>% 
    relocate(Chr, .after = Gene) %>% 
    filter(Pop != "POPx")

tajima <- tajima %>% 
    group_by(Pop) %>% 
    group_by(Chr) %>% 
    summarise(chr_len = max(Start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(tajima, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, Start) %>%
    mutate( BPcum = Start + tot)

axis_set <- tajima %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

ylimits <- c(floor(min(tajima$TajimasD)), abs(floor(max(tajima$TajimasD))) + 2)

colors <- c("gray30", "firebrick") 
# Plot

ggplot(tajima, aes(x = BPcum, y = TajimasD, color = as_factor(Chr))) +
    
    # Show all points
    geom_point(size = 1) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center, expand = c(0.02, 0.05)) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +
    scale_color_manual(values = rep(colors, unique(length(axis_set$Chr)))) +   
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Tajima's D") +
    
    # Custom the theme:
    theme_bw() +
    theme( 
        legend.position = "none",
        axis.line = element_line(color = "black"),
        axis.text = element_text(face = "bold", color = "black", size = 8),
        axis.title.x = element_text(vjust = 0, size = 10),
        axis.title.y = element_text(vjust = 2, size = 10),
        axis.title = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold.italic", hjust = 0.5, size = 10)) +
    facet_wrap(~factor(Pop, levels = c("WAF", "CAF", "SCAF", "EAF", "SEAF", "HAF"))) 

ggsave("figureS6/Figure.S6.jpeg", width = 30, height = 20, units = "cm")

#=======================
# Test of significance
#=======================
# convert Pop to a factor
tajima$Pop <- as.factor(tajima$Pop)

# Multiple comparison test
# ANOVA
model <- lm(TajimasD~Pop, data=tajima)
anova(model)

# Tukey multiple comparisons
test_significance <- summary(glht(model, mcp(Pop = "Tukey")))

#----------------------------------------

model1 <- aov(TajimasD~Pop, data=tajima)
summary(model1)
TukeyHSD(model1, conf.level=.95)
