#
# figure2.R
#
# (c) Mouhamadou Fadel DIOP (MRC at LSHTM)
# Date: 2022-03-14
#
# Purpose:
# Generate Beta Scores and HKA Manhattan plots
# across African populations.
#=====================================================
rm(list = ls())

# Note: Uncomment these lines to install packages as needed
# install.packages("tidyverse")
# install.packages("stringr")
# install.packages("ggrepel")

## Load packages
library(tidyverse)
library(stringr)
library(ggrepel)
library(patchwork)

# load annotation file
Gene_name <- read_tsv("data/PF_GeneName.tsv")

#===========
# Beta
#===========
beta_data <- readRDS("data/CountryPop_beta.rds")

beta <- beta_data %>%
    # mutate(Gene = gsub("\\.1", "", Gene)) %>% 
    inner_join(Gene_name, by = c("Gene" = "Gene_ID")) %>% 
    dplyr::select(-c(7, 8, 11)) %>% 
    filter(Beta > 0) %>% 
    relocate(Gene, .before = Chr.x) %>% 
    relocate(Beta, .after = End) %>% 
    relocate(Position, .after = End)

## Top Beta Scores
top_scores <- beta %>% 
    filter(Beta >= 5) %>%
    group_by(Gene) %>% 
    summarise(Max = max(Beta))

beta_formatted <- beta %>% 
    # Compute chromosome size
    group_by(Chrom) %>% 
    summarise(chr_len = max(Start)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(beta, ., by=c("Chrom"="Chrom")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chrom, Start) %>%
    mutate( BPcum = Start + tot) %>% 
    mutate(is_annotate = "no")


for (i in 1:nrow(top_scores)) {
    beta_formatted[which(beta_formatted$Gene == top_scores$Gene[i] & 
                             beta_formatted$Beta == top_scores$Max[i]), "is_annotate"] <- "yes"
}

highlight <- beta_formatted %>%
    filter(is_annotate == "yes") %>% 
    group_by(Gene) %>% 
    filter(Beta == max(Beta)) %>% 
    distinct(Population, GeneName, .keep_all = TRUE) %>% 
    mutate(NewGeneName = recode(GeneName, "ABC_transporter_B_family_member_4" = "ABC4",
                                "structural_maintenance_of_chromosomes_protein_3" = "SMCP3",
                                "mitochondrial_rp_S5_precursor" = "MTS5",
                                "RNA-binding_protein" = "RNABP",
                                "anaphase-promoting_complex_subunit_1" = "APCS1",
                                "alpha/beta_hydrolase" = "ABH",
                                "asparagine-rich_antigen_Pfa55-14" = "Pfa55",
                                "peptidase_family_C50" = "C50",
                                "Snf2-related_CBP_activator" = "Snf2",
                                "crystalloid-specific_PH" = "CS",
                                "HECT-like_E3_ubiquitin_ligase" = "HECT",
                                "ATP-dependent_6-phosphofructokinase" = "ATP6-",
                                "perforin-like_protein_3" = "PLP3",
                                "protein_kinase" = "PK",
                                "ADP/ATP_CP" = "ADP_CP",
                                "U2_snRNA/tRNA_PS" = "snRNA",
                                "probable_protein" = "PP",
                                "pre-mRNA-processing_factor_6" = "mRNA6",
                                "protein_phosphatase_2C" = "PP2C",
                                "antigen_332,_DBL-like_protein" = "DBL",
                                "telomere_repeat-binding_zinc_finger_protein" = "TRZF",
                                "translation_initiation_factor_IF-2" = "IF2",
                                "tyrosine_kinase-like_protein"= "TKLP",
                                "DNA_polymerase_epsilon_catalytic_subunit_A" = "DPECA",
                                "reticulocyte_binding_protein_2_homologue_b" = "RBP2",
                                "hypothetical_protein" = "HP",
                                "potassium_channel_K1" = "K1",
                                "Cg1_protein" = "Cg1",
                                "KELT_protein" = "KELT",
                                "ABC_transporter_B_family_member_3" = "ABC3",
                                "ABC_transporter_B_family_member_6" = "ABC6",
                                "acetyl-CoA_carboxylase" = "ACoA",
                                "phosducin-like_protein_1" = "PLP1",
                                "dynein-related_AAA-type_ATPase" = "DRA"))

# set x-axis
axis_set <- beta_formatted %>%
    group_by(Chrom) %>%
    summarize(center = mean(BPcum))

# set y-axis 
ylimits <- c(floor(min(beta_formatted$Beta)), floor(max(beta_formatted$Beta))+1)

# define color legend
col_pal <- c("gray30", "firebrick")

# Create Beta Scores Manhattan plot
plot1 <- beta_formatted %>% 
    
    ggplot(aes(x=BPcum, y=Beta, color = as_factor(Chrom))) +
    # Show all points
    geom_point(alpha = 0.8, size = 0.2) +
    geom_hline(yintercept = 5, color = "gray40", linetype = "dashed") +
    
    # custom X axis:
    scale_x_continuous(label = axis_set$Chrom, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = ylimits) +
    scale_color_manual(values = rep(col_pal, unique(length(axis_set$Chrom)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", 
         y = "Beta Scores") + 
    
    # Add highlighted points
    geom_point(data = highlight, size = 1.5) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data = highlight, aes(label=GeneName, fontface = "bold"), 
                      size=1.5, max.overlaps = Inf, colour = "black", hjust = 0) +
    
    # Custom the theme:
    theme_classic() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
        plot.title.position = "plot",
        axis.line = element_line(colour = "black"),
        panel.grid  = element_blank(),
        axis.title = element_text(face = "bold.italic"),
        axis.text = element_text(face = "bold", size = 6)) +
    ggtitle('A')

#==================================================
# LOAD AND FORMAT HKA RESULTS FOR ALL POPULATIONS
#==================================================

## Load data
ncd.hka_data <- readRDS("data/Combined_NCD_HKA.rds")

ncd.hka <- ncd.hka_data %>% 
    mutate(CHR = recode(CHROM, 
                        "Pf3D7_01_v3" = 1, "Pf3D7_02_v3" = 2, "Pf3D7_03_v3" = 3, "Pf3D7_04_v3" = 4, 
                        "Pf3D7_05_v3" = 5, "Pf3D7_06_v3" = 6, "Pf3D7_07_v3" = 7, "Pf3D7_08_v3" = 8, 
                        "Pf3D7_09_v3" = 9, "Pf3D7_10_v3" = 10, "Pf3D7_11_v3" = 11, "Pf3D7_12_v3" = 12, 
                        "Pf3D7_13_v3" = 13, "Pf3D7_14_v3" = 14))

## Add Gene Name Column
ncd.hka <- ncd.hka %>% 
    add_column(GeneName = NA)

for (i in 1:dim(Gene_name)[1]) {
    target <- Gene_name[i,]
    
    index <- which(ncd.hka$CHROM == target$Chromosome & 
                       (ncd.hka$sitePos >= target$Start & 
                            ncd.hka$sitePos <= target$End))
    ncd.hka$GeneName[index] <- target$NewGeneName
}


## Extract HKA below 300
ncd.hka_low <- ncd.hka %>% 
    filter(HKA < 200 & HKA >0) %>% 
    dplyr::select(-c("NCDopt", "NCDsub", "optF")) %>% 
    arrange(CHROM)

## Get Outliers
top_scores2 <- ncd.hka_low %>%
    filter(HKA >100) %>% 
    # group_by(CHROM, GeneName) %>%
    group_by(CHROM) %>%
    summarise(Max = max(HKA))

ncd.hka_formatted <- ncd.hka_low %>% 
    # Compute chromosome size
    group_by(CHROM) %>% 
    summarise(chr_len = max(sitePos)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(ncd.hka_low, ., by=c("CHROM"="CHROM")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHROM, sitePos) %>%
    mutate( BPcum = sitePos + tot) %>% 
    mutate(is_annotate = "no")

for (i in 1:nrow(top_scores2)) {
    ncd.hka_formatted[which(ncd.hka_formatted$CHROM == top_scores2$CHROM[i] &
                                ncd.hka_formatted$HKA == top_scores2$Max[i]), "is_annotate"] <- "yes"
}

highlight2 <- ncd.hka_formatted %>%
    filter(is_annotate == "yes") %>% 
    group_by(GeneName) %>%
    filter(HKA == max(HKA)) %>% 
    arrange(CHROM) %>% 
    distinct(GeneName, .keep_all = TRUE) %>% 
    mutate(NewGeneName = recode(GeneName, "serine/threonine_protein_kinase" = "STPK",
                                "CX3CL1-binding_protein_1" = "CX3CL1",
                                "GPCR-like_receptor_SR25" = "SR25",
                                "FIKK_family" = "FIKK",
                                "KELT_protein" = "KELT",
                                "1-cys-glutaredoxin-like_protein-1" = "CGLP1",
                                "oocyst_capsule_protein_Cap380" = "Cap380",
                                "nucleoporin_NUP221" = "NUP221",
                                "zinc_finger_protein" = "zinc_finger",
                                "ubiquitin-conjugating_enzyme_E2"= "UCE2",
                                "ATP-dependent_6-phosphofructokinase" = "ADP",
                                "ATP-dependent_helicase" = "ADH",
                                "major_facilitator_superfamily" = "MFS",
                                "U2_snRNA/tRNA_PS" = "U2PS",
                                "RESA-like_protein_with_PHIST_and_DnaJ_domains" = "RESA",
                                "osmiophilic_body_protein_G377" = "G377",
                                "cytochrome_c_oxidase_subunit_2" = "CCO2",
                                "LCCL/lectin_adhesive-like_protein_1" = "LCCL",
                                "glutamate_synthase_[NADH]" = "NADH",
                                "plasmepsin_I" = "PLI",
                                "lysophospholipase"= "LPase",
                                "probable_protein"= "PP",
                                "hypothetical_protein"= "HP",
                                "pre-mRNA-processing-splicing_factor_8" = "mRNA8",
                                "RNA_polymerase_I" = "RNAI",
                                "DNA_polymerase_epsilon_catalytic_subunit_A" = "DNA",
                                "cysteine_repeat_modular_protein_2" = "CRMP2",
                                "dynein_heavy_chain" = "DHC",
                                "antigen_332,_DBL-like_protein" = "DBL",
                                "DnaJ_protein" = "DnaJ",
                                "lysine-rich_membrane-associated_PHISTb_protein" = "PHISTb",
                                "TRAP-like_protein" = "TRAP",
                                "regulator_of_chromosome_condensation" = "RCC",
                                "inner_membrane_complex_protein_1k" = "IMCP",
                                "glideosome-associated_connector" = "GAC",
                                "secreted_ookinete_protein" = "SOP",
                                "rhoptry_neck_protein_2" = "RNP",
                                "von_Willebrand_factor_A_domain-related_protein" = "von",
                                "EMP1-trafficking_protein" = "EMP1"))

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead.

axis_set2 <- ncd.hka_formatted %>%
    group_by(CHR) %>%
    summarize(center = mean(BPcum))

ylimits2 <- c(floor(min(ncd.hka_formatted$HKA)), 
              floor(max(ncd.hka_formatted$HKA))+5)

plot2 <- ncd.hka_formatted %>%  
    
    ggplot(aes(x = BPcum, y = HKA, color = as_factor(CHR))) + 
    
    # Show all points
    geom_point(alpha = 0.75, size = 0.3) +
    geom_hline(yintercept = 100, color = "gray40", linetype = "dashed") +
    
    # custom X axis:
    scale_x_continuous(label = axis_set2$CHR, breaks = axis_set2$center, expand = c(0.02,0.05)) +
    scale_y_continuous(expand = c(0,0), limits = ylimits2) +
    scale_color_manual(values = rep(col_pal, unique(length(axis_set2$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "HKA") + 
    
    # Add highlighted points
    geom_point(data = highlight2, size = 1.5) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data = highlight2, aes(label = NewGeneName, fontface = "bold"), 
                     size=1.5, max.overlaps = 50, colour = "black") +
    
    # Custom the theme:
    theme_classic() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
        plot.title.position = "plot",
        panel.grid  = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(face = "bold.italic"),
        axis.text = element_text(face = "bold", size = 6))

# #============================================================
# ## Extract HKA above 300
# ncd.hka_high <- ncd.hka %>% 
#     filter(HKA > 200) %>% 
#     select(-c("NCDopt", "NCDsub", "optF")) %>% 
#     arrange(CHROM)
# 
# 
# ## Get Outliers
# top_scores3 <- ncd.hka_high %>%
#     group_by(CHROM) %>% 
#     arrange(desc(HKA)) %>%
#     slice(1) %>% 
#     summarise(Max = max(HKA))
# 
# ncd.hka_high_format <- ncd.hka_high %>% 
#     add_row(Population= NA, CHROM = "Pf3D7_06_v3", sitePos = 1291787,   
#             HKA =210, NCD=0.1, CHR = 6, GeneName = NA) %>% 
#     arrange(CHROM)
# 
# ncd.hka_high_format <- ncd.hka_high_format %>% 
#     
#     # Compute chromosome size
#     group_by(CHROM) %>% 
#     summarise(chr_len = max(sitePos)) %>% 
#     
#     # Calculate cumulative position of each chromosome
#     mutate(tot = cumsum(chr_len) - chr_len) %>%
#     select(-chr_len) %>%
#     ungroup() %>% 
#     # Add this info to the initial dataset
#     left_join(ncd.hka_high_format, ., by=c("CHROM"="CHROM")) %>%
#     
#     # Add a cumulative position of each SNP
#     arrange(CHROM, sitePos) %>%
#     mutate( BPcum = sitePos + tot) %>% 
#     mutate(is_annotate = "no")
# 
# for (i in 1:nrow(top_scores3)) {
#     ncd.hka_high_format[which( ncd.hka_high_format$CHROM == top_scores3$CHROM[i] & 
#                                   ncd.hka_high_format$HKA == top_scores3$Max[i]), "is_annotate"] <- "yes"
# }
# 
# highlight3 <- ncd.hka_high_format %>%
#     filter(is_annotate == "yes") %>% 
#     group_by(GeneName) %>% 
#     filter(HKA == max(HKA)) %>% 
#     arrange(CHROM) %>% 
#     distinct(Population, GeneName, .keep_all = TRUE) %>% 
#     mutate(NewGeneName = recode(GeneName,
#                                 "pre-mRNA-processing-splicing_factor_8" = "mRNA8",
#                                 "2-oxoisovalerate_dehydrogenase_subunit_beta,_mitochondrial" = "beta",
#                                 "cysteine_repeat_modular_protein_2" = "CRP2",
#                                 "dynein_heavy_chain" = "DHC",
#                                 "EMP1-trafficking_protein" = "EMP1",
#                                 "nucleoporin_NUP221" = "NUP221",
#                                 "glutamate-rich_protein_GLURP" = "GLURP",
#                                 "osmiophilic_body_protein_G377" = "G377",
#                                 "FIKK_family" = "FIKK",
#                                 "LCCL/lectin_adhesive-like_protein_1" = "LCCL",
#                                 "glutamate_synthase_[NADH]" = "NADH",
#                                 "1-cys-glutaredoxin-like_protein-1" = "CGLP1",
#                                 "zinc_finger_protein" = "zinc_finger",
#                                 "ubiquitin-conjugating_enzyme_E2"= "UCE2",
#                                 "lysophospholipase"= "LP",
#                                 "probable_protein"= "PP",
#                                 "hypothetical_protein"= "HP",
#                                 "pre-mRNA-processing-splicing_factor_8" = "mRNA8",
#                                 "RNA_polymerase_I" = "RNAI",
#                                 "DNA_polymerase_epsilon_catalytic_subunit_A" = "DNA",
#                                 "TRAP-like_protein" = "TRAP",
#                                 "antigen_332,_DBL-like_protein" = "DBL",
#                                 "DnaJ_protein" = "DnaJ",
#                                 "lysine-rich_membrane-associated_PHISTb_protein" = "PHISTb",
#                                 "nuclear_polyadenylated_RNA-binding_protein_NAB2" = "NAB2",
#                                 "regulator_of_chromosome_condensation" = "RCC",
#                                 "KELT_protein" = "KELT",
#                                 "glideosome-associated_connector" = "GAC",
#                                 "secreted_ookinete_protein" = "SOP",
#                                 "rhoptry_neck_protein_2" = "RNP",
#                                 "1-deoxy-D-xylulose_5-phosphate_synthase" = "DDPS",
#                                 "phosducin-like_protein_1" = "PLP1"))
# 
# # Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, 
# # but just show the chromosome name instead.
# 
# axis_set3 <- ncd.hka_high_format %>%
#     group_by(CHR) %>%
#     summarize(center = mean(BPcum))
# 
# ylimits3 <- c(floor(min(ncd.hka_high_format$HKA)), 
#               floor(max(ncd.hka_high_format$HKA))+5)
# 
# plot3 <- ncd.hka_high_format %>%  
#     
#     ggplot(aes(x = BPcum, y = HKA, color = as_factor(CHR))) + 
#     
#     # Show all points
#     geom_point(alpha = 0.75, size = 0.7, na.rm = TRUE) +
#     
#     # custom X axis:
#     scale_x_continuous(label = axis_set3$CHR, breaks = axis_set3$center, expand = c(0.02,0.05)) +
#     scale_y_continuous(expand = c(0,0), limits = ylimits3) +
#     scale_color_manual(values = rep(col_pal, unique(length(axis_set3$CHR)))) +
#     scale_size_continuous(range = c(0.5,3)) +
#     labs(x = "Chromosomes", y = "HKA") + 
#     
#     # Add highlighted points
#     geom_point(data = highlight3, size = 0.5) +
#     
#     # Add label using ggrepel to avoid overlapping
#     geom_label_repel( data = highlight3, aes(label=NewGeneName, fontface = "bold"), 
#                       size=1.4, max.overlaps = 50, colour = "black") +  
#     
#     # Custom the theme:
#     theme_classic() +
#     theme( 
#         legend.position = "none",
#         plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
#         plot.title.position = "plot",
#         panel.grid = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         axis.text.y = element_text(face = "bold", size = 6)) +
#     ggtitle('B')
# 
# plot1 / plot3 / plot2 +
#     plot_layout(widths = 2, heights = unit(c(6, 2, 6), c('cm', 'cm', 'cm')))


#============================================================
## Extract HKA above 300
ncd.hka_high_test <- ncd.hka %>% 
    filter(Population %in% c("Ethiopia", "Gabon", "Madagascar")) %>% 
    dplyr::select(-c("NCDopt", "NCDsub", "optF")) %>% 
    arrange(CHROM)


## Get Outliers
top_scores <- ncd.hka_high_test %>%
    group_by(CHROM) %>% 
    arrange(desc(HKA)) %>%
    slice(1) %>% 
    summarise(Max = max(HKA))

# ncd.hka_high_format <- ncd.hka_high %>% 
#     add_row(Population= NA, CHROM = "Pf3D7_06_v3", sitePos = 1291787,   
#             HKA =210, NCD=0.1, CHR = 6, GeneName = NA) %>% 
#     arrange(CHROM)

ncd.hka_high_test <- ncd.hka_high_test %>% 
    
    # Compute chromosome size
    group_by(CHROM) %>% 
    summarise(chr_len = max(sitePos)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    ungroup() %>% 
    # Add this info to the initial dataset
    left_join(ncd.hka_high_test, ., by=c("CHROM"="CHROM")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHROM, sitePos) %>%
    mutate( BPcum = sitePos + tot) %>% 
    mutate(is_annotate = "no")

for (i in 1:nrow(top_scores)) {
    ncd.hka_high_test[which( ncd.hka_high_test$CHROM == top_scores$CHROM[i] & 
                                 ncd.hka_high_test$HKA == top_scores$Max[i]), "is_annotate"] <- "yes"
}

highlight4 <- ncd.hka_high_test %>%
    filter(is_annotate == "yes") %>% 
    group_by(CHROM, GeneName) %>% 
    filter(HKA == max(HKA)) %>% 
    arrange(CHROM) %>% 
    distinct(Population, GeneName, .keep_all = TRUE) %>% 
    mutate(NewGeneName = recode(GeneName,
                                "pre-mRNA-processing-splicing_factor_8" = "mRNA8",
                                "2-oxoisovalerate_dehydrogenase_subunit_beta,_mitochondrial" = "beta",
                                "cysteine_repeat_modular_protein_2" = "CRP2",
                                "dynein_heavy_chain" = "DHC",
                                "EMP1-trafficking_protein" = "EMP1",
                                "nucleoporin_NUP221" = "NUP221",
                                "glutamate-rich_protein_GLURP" = "GLURP",
                                "osmiophilic_body_protein_G377" = "G377",
                                "FIKK_family" = "FIKK",
                                "LCCL/lectin_adhesive-like_protein_1" = "LCCL",
                                "glutamate_synthase_[NADH]" = "NADH",
                                "1-cys-glutaredoxin-like_protein-1" = "CGLP1",
                                "zinc_finger_protein" = "zinc_finger",
                                "ubiquitin-conjugating_enzyme_E2"= "UCE2",
                                "lysophospholipase"= "LP",
                                "probable_protein"= "PP",
                                "hypothetical_protein"= "HP",
                                "pre-mRNA-processing-splicing_factor_8" = "mRNA8",
                                "RNA_polymerase_I" = "RNAI",
                                "DNA_polymerase_epsilon_catalytic_subunit_A" = "DNA",
                                "TRAP-like_protein" = "TRAP",
                                "antigen_332,_DBL-like_protein" = "DBL",
                                "DnaJ_protein" = "DnaJ",
                                "lysine-rich_membrane-associated_PHISTb_protein" = "PHISTb",
                                "nuclear_polyadenylated_RNA-binding_protein_NAB2" = "NAB2",
                                "regulator_of_chromosome_condensation" = "RCC",
                                "KELT_protein" = "KELT",
                                "glideosome-associated_connector" = "GAC",
                                "secreted_ookinete_protein" = "SOP",
                                "rhoptry_neck_protein_2" = "RNP",
                                "1-deoxy-D-xylulose_5-phosphate_synthase" = "DDPS",
                                "phosducin-like_protein_1" = "PLP1"))

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead.

axis_set3 <- ncd.hka_high_test %>%
    group_by(CHR) %>%
    summarize(center = mean(BPcum))

ylimits3 <- c(100, floor(max(ncd.hka_high_test$HKA))+5)

plot3 <- ncd.hka_high_test %>%  
    
    ggplot(aes(x = BPcum, y = HKA, color = as_factor(CHR))) + 
    
    # Show all points
    geom_point(alpha = 0.75, size = 0.3, na.rm = TRUE) +
    
    # custom X axis:
    scale_x_continuous(label = axis_set3$CHR, breaks = axis_set3$center, expand = c(0.02,0.05)) +
    scale_y_continuous(expand = c(0,0), limits = ylimits3) +
    scale_color_manual(values = rep(col_pal, unique(length(axis_set3$CHR)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosomes", y = "HKA") + 
    
    # Add highlighted points
    geom_point(data = highlight4, size = 1.5) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data = highlight4, aes(label=NewGeneName, fontface = "bold"), 
                      size=1.5, max.overlaps = 50, colour = "black") +  
    
    # Custom the theme:
    theme_classic() +
    theme( 
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0),
        plot.title.position = "plot",
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 6)) +
    ggtitle('B')

plot1 / plot3 / plot2 +
    plot_layout(widths = 2.5, heights = unit(c(5.5, 2.3, 6), c('cm', 'cm', 'cm')))

## Save the plot
ggsave("figure2/figure2.jpeg",
       width = 250, height = 180, units = "mm")

