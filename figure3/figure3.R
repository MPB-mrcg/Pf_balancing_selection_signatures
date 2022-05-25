#
# figure3.R
#
# (c) Mouhamadou Fadel DIOP (MRC at LSHTM)
# Date: 2022-03-14
#
# Purpose:
# Generate histogramme of Gene Ontology Terms
# between balnced and not balanced genes.
#=====================================================
rm(list = ls())
library(tidyverse)
library(patchwork)

# Load tajima's D results
tajima <- read_tsv("balancing_selection_signatures/data/tajima_all_populations.txt")

# Extract top first 100 genes
head_genes <- tajima %>% 
    filter(Population != "Ethiopia") %>% 
    group_by(Population) %>% 
    filter(Tajima >=0 & Fu.Li_F>= 0) %>% 
    arrange(Gene, .by_group = TRUE) %>% 
    ungroup()

# Transform data from long to wide
head_wide <- as_tibble(reshape2::dcast(head_genes, Gene+Chrom+Start+End+GeneName ~ Population, value.var="Tajima")) %>% 
    rename(DR_Congo = `DR Congo`) %>% 
    arrange(Gene)

## Fill NA values with their exact values

for (i in 6:ncol(head_wide)) {
    cat(paste0("Processing ", colnames(head_wide)[i], "\n"))
    
    for (j in 1:nrow(head_wide)) {
        if(is.na(head_wide[j,i])) {
            gene <- head_wide[j,1] %>% pull(); 
            pop <- colnames(head_wide)[i]
            values <- tajima %>% 
                filter(Population == pop) %>%
                filter(Gene == gene) %>% 
                select(Tajima) %>% pull()
            
            if (!identical(values,numeric(0))) {
                head_wide[j,i] <- round(values, 3)
            }
            
        }
    }
}

head_wide %>% 
    arrange(Gene) %>% 
    mutate_if(is.numeric, round, digits=3) %>% 
    write.table(paste0("balancing_selection_signatures/figure3/balanced_genes.xlsx"), 
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#============================
# Extract top last 100 genes
#============================
tail_genes <- tajima %>% 
    filter(Population != "Ethiopia") %>% 
    group_by(Population) %>% 
    arrange(Tajima, .by_group = TRUE) %>% 
    slice(1:150) %>% 
    arrange(Gene, .by_group = TRUE) %>% 
    ungroup()


# Transform data from long to wide
tail_wide <- as_tibble(reshape2::dcast(tail_genes, Gene+Chrom+Start+End+GeneName ~ Population, value.var="Tajima")) %>% 
    rename(DR_Congo = `DR Congo`) %>% 
    arrange(Gene)

## Remove rows with more than 50% NA
tail_wide <- tail_wide[which(rowMeans(!is.na(tail_wide)) > 0.3), ]

for (i in 6:ncol(tail_wide)) {
    cat(paste0("Processing ", colnames(tail_wide)[i], "\n"))
    
    for (j in 1:nrow(tail_wide)) {
        if(is.na(tail_wide[j,i])) {
            gene <- tail_wide[j,1] %>% pull(); 
            pop <- colnames(tail_wide)[i]
            values <- tajima %>% 
                filter(Population == pop) %>%
                filter(Gene == gene) %>% 
                select(Tajima) %>% pull()
            
            if (!identical(values,numeric(0))) {
                tail_wide[j,i] <- round(values, 3)
            }
            
        }
    }
}


tail_wide %>% 
    arrange(Gene) %>% 
    mutate_if(is.numeric, round, digits=3) %>% 
    write.table(paste0("balancing_selection_signatures/figure3/coldSpots.xlsx"), 
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ==================================================================
## FUNCTION TO ADD DOTS AFTER A LONG NAMES
tidy_name <- function(name, n_char) {
    ifelse(nchar(name) > (n_char - 2), 
           {substr(name, 1, n_char) %>% paste0(., "..")},
           name)
}

#====================
# Biological Process
#====================
cold.data <- readxl::read_xlsx("balancing_selection_signatures/figure3/David_GO/coldSpot_genes_BF.xlsx", sheet = 2) %>%
    add_column(Condition = "Not balanced")

hot.data <- readxl::read_xlsx("balancing_selection_signatures/figure3/David_GO/balanced_genes_BF.xlsx", sheet = 2) %>%
    add_column(Condition = "Balanced")

biological <- tibble(rbind(cold.data, hot.data))
biological$Term <- str_replace(biological$Term, "^\\w{1}", toupper)
biological$Description <- biological$Term %>% tidy_name(20)

p1 <- ggplot(biological, aes(fill = Condition, y = Count, x = Description)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    scale_fill_manual(values=c("black", "white")) +
    ggtitle("Biological Process") +
    labs(x = "", y = "# of genes") +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 7.5, colour = "black"),
          axis.text = element_text(face = "bold", size = 7, color = 'black'),
          axis.title = element_text(size = 8.5, face = "bold.italic"))

#====================
# Cellular Component
#====================
cold.data <- readxl::read_xlsx("balancing_selection_signatures/figure3/David_GO/coldSpot_genes_CC.xlsx", sheet = 2) %>%
    add_column(Condition = "Not balanced")

hot.data <- readxl::read_xlsx("balancing_selection_signatures/figure3/David_GO/balanced_genes_CC.xlsx", sheet = 2) %>%
    add_column(Condition = "Balanced")


cellular <- tibble(rbind(cold.data, hot.data))
cellular$Term <- str_replace(cellular$Term, "^\\w{1}", toupper)
cellular$Description <- cellular$Term %>% tidy_name(20)

p2 <- ggplot(cellular, aes(fill = Condition, y = Count, x = Description)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    scale_fill_manual(values=c("black", "white")) +
    ggtitle("Cellular Component") +
    labs(x = "", y = "# of genes") +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 7.5, colour = "black"),
          axis.text = element_text(face = "bold", size = 7, color = 'black'),
          axis.title = element_text(size = 8.5, face = "bold.italic"))

#====================
# Molecular Function
#====================
cold.data <- readxl::read_xlsx("balancing_selection_signatures/figure3/David_GO/coldSpot_genes_MF.xlsx", sheet = 2) %>%
    add_column(Condition = "Not balanced")

hot.data <- readxl::read_xlsx("balancing_selection_signatures/figure3/David_GO/balanced_genes_MF.xlsx", sheet = 2) %>%
    add_column(Condition = "Balanced")


molecular <- tibble(rbind(cold.data, hot.data))
molecular$Term <- str_replace(molecular$Term, "^\\w{1}", toupper)
molecular$Description <- molecular$Term %>% tidy_name(20)

p3 <- ggplot(molecular, aes(fill = Condition, y = Count, x = Description)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    scale_fill_manual(values=c("black", "white")) +
    ggtitle("Molecular Function") +
    labs(x = "", y = "# of genes") +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 7.5, colour = "black"),
          axis.text = element_text(face = "bold", size = 7, color = 'black'),
          axis.title = element_text(size = 8.5, face = "bold.italic"))

#================
# combined.plot
#================
p1 / p2 / p3 + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 8, face = 'bold'),
          legend.position = 'right',
          legend.title = element_blank())

ggsave("balancing_selection_signatures/figure3/figure3_v1.jpeg", 
       width = 250, height = 210, units = "mm")

#=====================
#   SECOND VERSION
#====================

p1_1 <- ggplot(biological, aes(fill = Condition, y = Count, x = Term)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    scale_fill_manual(values=c("black", "white")) +
    aes(stringr::str_wrap(Term, 20), Count) +
    ggtitle("Biological Process") +
    labs(x = "", y = "# of genes") +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 7.5, colour = "black"),
          axis.text = element_text(face = "bold", size = 7, color = 'black'),
          axis.title = element_text(size = 8.5, face = "bold.italic"))

p2_1 <- ggplot(cellular, aes(fill = Condition, y = Count, x = Term)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    scale_fill_manual(values=c("black", "white")) +
    aes(stringr::str_wrap(Term, 20), Count) +
    ggtitle("Cellular Component") +
    labs(x = "", y = "# of genes") +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 7.5, colour = "black"),
          axis.text = element_text(face = "bold", size = 7, color = 'black'),
          axis.title = element_text(size = 8.5, face = "bold.italic"))

p3_1 <- ggplot(molecular, aes(fill = Condition, y = Count, x = Term)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    scale_fill_manual(values=c("black", "white")) +
    aes(stringr::str_wrap(Term, 35), Count) +
    ggtitle("Molecular Function") +
    labs(x = "", y = "# of genes") +
    theme_classic() +
    theme(legend.position = "none",
          legend.text = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 7.5, colour = "black"),
          axis.text = element_text(face = "bold", size = 7, color = 'black'),
          axis.title = element_text(size = 8.5, face = "bold.italic"))

# combined.plot
p1_1 / p2_1 / p3_1 + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 8, face = 'bold'),
          legend.position = 'right',
          legend.title = element_blank())

ggsave("balancing_selection_signatures/figure3/figure3_v2.jpeg", 
       width = 280, height = 210, units = "mm")
