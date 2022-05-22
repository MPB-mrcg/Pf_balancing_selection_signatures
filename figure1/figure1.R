
# figure1.R
#
# (c) Mouhamadou Fadel DIOP (MRC at LSHTM)
# Date: 2022-05-21
#
# Purpose:
# Generate distribution, correlation and Manhattan plots.
#
# ------------------------------------------------------

rm(list = ls())

# NOTE - Uncomment these lines to install packages as needed
# install.packages("patchwork")
# install.packages("ggridges")
# install.packages("ggrepel")
# install.packages("tidyverse")
# install.packages("forcats")

# Load library
library(ggridges)
library(tidyverse)
library(forcats)
library(patchwork)
library(ggrepel)

# Load annotation file
Gene_name <- read_tsv("../../balancing_selection_signatures/data/PF_GeneName.tsv")

# Load raw data
tajima <- readRDS("../../balancing_selection_signatures/data/CountryPop_tajD.rds") %>%
    inner_join(Gene_name, by = c("Gene" = "Gene_ID")) %>% 
    select(-c("Chromosome", "Start.y", "End.y", "GeneName" ))

# Define color legend
couleurs <- c("#FFFF00","#EEC591","#912CEE","#00FFFF","#87CEFF","#006400","#7CFC00","#40E0D0",
              "#828282","#000000","#EE7AE9","#0000CD","#FF3399","#EE0000", "#8B3636","#FF7F00")

#========================
# RIDGELINE DENSITY
#========================
myorder <- tajima %>% 
    pull(Population) %>% 
    unique()

plot1 <- tajima %>%
    mutate(Order = factor(Population, levels = rev(myorder))) %>%
    ggplot( aes(y = Order, x = Tajima,  fill = Order)) +
    geom_density_ridges_gradient(scale = 4) +
    geom_vline(xintercept = 0, linetype="dashed", color = "#828282", size=1) +
    scale_fill_manual(values = couleurs) +
    labs(x = "Tajima's D", y = "") + ggtitle('A') +
    theme_ridges() +
    theme(
        legend.position="none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
        plot.title.position = "plot",
        panel.spacing = unit(0.1, "lines"),
        axis.text = element_text(size = 8, color = 'black', face = "bold.italic"),
        axis.title = element_text(size = 7, color = 'black', face = "bold.italic"))

#==========================
# CORRELATION PLOT 
#=========================
# Extract genes with highest values of Tajima's D
top <- tajima %>%
    filter(Population != "Ethiopia" & (Tajima > 2 & Fu.Li_F > 2)) %>%
    group_by(Gene) %>%
    arrange(Gene) %>%
    mutate(NewGeneName = recode(NewGeneName,
                                "apoptosis-inducing_factor" = "AIF",
                                "RNA-binding_protein" = "RNABP",
                                "ADP/ATP_CP" = "ADP",
                                "debranching_enzyme-associated_ribonuclease" = "DEAR"))

# Genes to label
highlight <- top %>% 
    filter(Tajima == max(Tajima))

## Zoomed on tail genes
# Identifying top last genes on Population
tail <- tajima %>%
    filter(Population != "Ethiopia" & (Tajima < -2 & Fu.Li_F < -5.5)) %>%
    mutate(NewGeneName = recode(NewGeneName,
                                "TatD-like_deoxyribonuclease" = "TatD",
                                "LETM1-like_protein" = "LETM1",
                                "DnaJ_protein" = "DnaJ",
                                "V-type_proton_ATPase_catalytic_subunit_A" = "ACA",
                                "importin_subunit_beta" = "beta",
                                "alpha_tubulin_1" = "AT1",
                                "proteasome_subunit_beta_type-5" = "beta5",
                                "60S_rp_P0" = "SPO",
                                "peroxisome_assembly_protein_22" = "PA22",
                                "cAMP-dependent_protein_kinase_regulatory_subunit" = "cAMP",
                                "regulator_of_nonsense_transcripts_3B" = "RT3B",
                                "histone-lysine_N-methyltransferase" = "HLN",
                                "60S_rp_L1" = "SL1",
                                "U3_small_nucleolar_ribonucleoprotein_protein_IMP3" = "IMP3",
                                "heat_shock_protein_101" = "H101",
                                "ATP-dependent_zinc_metalloprotease_FTSH" = "FTSH",
                                "repetitive_organellar_protein" = "ROP",
                                "calcium-binding_protein" = "CBP",
                                "AP2_domain_transcription_factor_AP2-O" = "AP2O",
                                "peptide_chain_release_factor_subunit_1" = "PCF1",
                                "NOT_family_protein" = "NFP",
                                "zinc_finger_protein" = "ZFP"))

#================================================
## Correlation between Tajima's D and Fu & Li F*
#================================================
plot2 <- tajima %>% 
    filter(Population != "Ethiopia") %>% 
    ggplot(aes(x = Tajima, y = Fu.Li_F)) + 
    geom_point(size = 0.8, alpha = 3/5) +
    geom_label_repel(data=highlight, aes(x=Tajima, y = Fu.Li_F, label = NewGeneName, 
                                         fill = factor(NewGeneName), fontface = "bold"),
                     force = 5, size = 1.65, force_pull = 5, max.overlaps = Inf) +
    geom_point(data = tail, color = "red") +
    geom_point(data = top , color = "red") +
    labs( x = "Tajima's D", y = "Fu & Li F*") + 
    theme_classic() + 
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
          plot.title.position = "plot",
          axis.line = element_line(colour = "black"),
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 8, face = "bold"),
          axis.title = element_text(size = 7, face = "bold.italic")) +
    ggtitle('B')

# Create zoomed plot in the correlation plot
p1 <- tail %>% 
    ggplot(aes(x = Tajima, y = Fu.Li_F)) + 
    geom_point(size = 1, color = "red") +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data = tail, aes(label=NewGeneName, fontface = "bold"), size=1.2, color = "black") +
    labs( x = "Tajima's D", y = "Fu & Li F*") + 
    theme_bw() + 
    theme(legend.position = "none",
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 4, face = "bold"),
          axis.title = element_text(size = 4, face = "bold.italic")) 

#===================
# MANHATTAN PLOT 
#===================

tajima1 <- tajima %>% 
    
    # Compute chromosome size
    group_by(Chrom) %>% 
    summarise(chr_len = max(Start.x)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(tajima, ., by=c("Chrom"="Chrom")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chrom, Start.x) %>%
    mutate( BPcum = Start.x + tot) %>% 
    
    # Add highlight and annotation information
    mutate( is_highlight = ifelse(Tajima >= 1, "yes", "no")) %>%
    mutate( is_annotate = ifelse(Tajima >= 1, "yes", "no")) %>% 
    ungroup() %>% 
    
    # Remove Ethiopia
    filter(Population != "Ethiopia")


to_be_highlited <- tajima1 %>% 
    filter(Tajima >= 2) %>%
    group_by(Gene) %>% 
    filter(Tajima == max(Tajima)) %>% 
    arrange(Gene) %>% 
    mutate(NewGeneName = recode(NewGeneName, 
                                "apoptosis-inducing_factor" = "AIF",
                                "RNA-binding_protein" = "RNABP",
                                "ADP/ATP_CP" = "ADP",
                                "debranching_enzyme-associated_ribonuclease" = "DEAR"))

to_be_highlited <- to_be_highlited %>%
    add_column(Gene_Code = seq(1:nrow(to_be_highlited)))

axis_set <- tajima1 %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum))

ylimits <- c(floor(min(tajima1$Tajima)), abs(floor(min(tajima1$Tajima))) + 2)

# Ready to make the plot using ggplot2:
plot3 <- ggplot(tajima1, aes(x = BPcum, y = Tajima, color = as_factor(Chr), size = Tajima)) +
    # Show all points
    geom_point(alpha = 0.75, size = 0.7) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = 2, color = "red", linetype = "dashed") +
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +
    scale_color_manual(values = rep(c("gray30", "firebrick"), unique(length(axis_set$Chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "Tajima's D") +
    ggtitle('C') +
    
    # Add highlighted points
    geom_point(data = to_be_highlited, size=1.6) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data = to_be_highlited, aes(label=Gene_Code, fontface = "bold"), size=1.8, color = "black") +
    
    # Custom the theme:
    theme_minimal() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
        plot.title.position = "plot",
        axis.line = element_line(colour = "black"),
        axis.title = element_text(face = "bold.italic", size = 7),
        axis.text = element_text(colour = "black", face = "bold", size = 8))

# Combine plots
(plot1 | (plot2 + inset_element(p1, 0.45, 0, 1, 0.55)) ) / 
    plot3 + plot_layout(guides = 'collect')


# save to file
ggsave("figure1/figure1.jpeg", width = 250, height = 180, units = "mm")
    