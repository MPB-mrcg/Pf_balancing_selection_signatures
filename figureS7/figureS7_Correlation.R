
library(tidyverse)
library(rlang)
library(ggrepel)
library(stringr)
theme_set(theme_bw())

#=================
# Load Tajima's D
#=================
tajima <- read_tsv("../balancing_selection_signatures/data/tajima_for_plot_correlation.txt")

#======================
# Identified top genes
#======================
top <- tajima %>% 
    arrange(desc(Tajima)) %>% 
    group_by(Population) %>% slice(1:5)

#================================================
## Correlation between Tajima's D and Fu&Li D*
#================================================
plot <- tajima %>% 
    ggplot(aes(x=Tajima, y=Fu.Li_D)) + 
    geom_point(size=1, alpha = 3/5) +
    geom_label_repel(data=top, 
                     aes(x=Tajima, y = Fu.Li_D, label = GeneName, fill = factor(GeneName)), force = 2,
                     size = 2, force_pull = 2, max.overlaps = Inf) + 
    geom_point(size=1, data = tajima[tajima$Tajima  >1 & tajima$Fu.Li_D > 1,], color = "red") +
    labs( x = "Tajima's D", y = "Fu & Li D*") + theme_bw() + 
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    facet_wrap(~Population, nrow = 3)


print(plot)

#==============
# save to file
#==============
file_ext <- c("jpeg", "pdf", "png")
for (i in seq_along(file_ext)) {
    ggsave(sprintf("figure1_Correlation/figure1_Correlation.%s", file_ext[i]),
           plot = plot, device = file_ext[i], width = 179, height = 180, units = "mm")
}   

#================================================
## Extract genes (Table) with both tajima's D 
## and Fu&Li F*/D* > 1 for all single populations
#================================================

#=================
# Load Tajima's D
#=================
tajima <- read_tsv("../balancing_selection_signatures/tajima_all_populations.txt")

#======================
# Identified top genes
#======================
top <- tajima %>% 
    arrange(desc(Tajima)) %>% 
    group_by(Population) %>% slice(1:5)

#====================================================
## PLOT CORRELATION PLOT FOR ALL SINGLE POPULATIONS 
#====================================================

#============================
# Identifying outliers genes
#============================
outliers <- tajima %>% 
    arrange(desc(Tajima)) %>% 
    group_by(Population) %>% slice(1:5)

#===============================================
## Correlation between Tajima's D and Fu&Li F*
#===============================================
plot1 <- tajima %>% 
    ggplot(aes(x = Tajima, y = Fu.Li_F)) + 
    geom_point(size = 1, alpha = 3/5) +
    geom_label_repel(data = outliers, 
                     aes(x = Tajima, y = Fu.Li_F, label = GeneName, fill = factor(GeneName)), 
                     size = 1.3, force = 2, max.overlaps = Inf, force_pull = 2) + 
    geom_point(size=1, data = outliers, color = "red") +
    labs( x="Tajima's D", y="Fu & Li F*") + theme_bw() + 
    scale_x_continuous(limits = c(-4, 5)) +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    facet_wrap(~Population, nrow = 4)

print(plot1)

#==============
# save to file
#==============
file_ext <- c("jpeg", "pdf", "png")
for (i in seq_along(file_ext)) {
    ggsave(sprintf("tajD_FuLiF_Correlation.%s", file_ext[i]),
           plot = plot1, device = file_ext[i], width = 180, height = 180, units = "mm")
}

#===============================================
## Correlation between Tajima's D and Fu&Li D*
#===============================================

plot2 <- tajima %>% 
    ggplot(aes(x = Tajima, y = Fu.Li_D)) + 
    geom_point(size = 1, alpha = 3/5) +
    geom_label_repel(data = outliers, 
                     aes(x = Tajima, y = Fu.Li_F, label = GeneName, fill = factor(GeneName)), 
                     size = 1.3, force = 2, max.overlaps = Inf, force_pull = 2) + 
    geom_point(size=1, data = outliers, color = "red") +
    labs( x="Tajima's D", y="Fu & Li D*") + theme_bw() + 
    scale_x_continuous(limits = c(-4, 5)) +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=8.5, face="bold"),
          axis.title=element_text(size=10,face="bold")) +
    facet_wrap(~Population, nrow = 4)


print(plot2)

#==============
# save to file
#==============
file_ext <- c("jpeg", "pdf", "png")
for (i in seq_along(file_ext)) {
    ggsave(sprintf("tajD_FuLiD_Correlation.%s", file_ext[i]),
           plot = plot2, device = file_ext[i], width = 180, height = 180, units = "mm")
}
