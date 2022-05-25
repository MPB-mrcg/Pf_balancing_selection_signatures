
setwd('/home/fadel/BalancingSelection/Analysis')

library(SeqArray)
seqVCF2GDS("CodingRegions.vcf.gz", "AfricanPopComplexity.gds")


my_vcf <-seqOpen("AfricanPopComplexity.gds") 
seqSummary(my_vcf)

# save sample identifiers
sample.id <- seqGetData(my_vcf, "sample.id")

library(moimix)

# get genomic coordinates of all variants
coords <- getCoordinates(my_vcf)

# filter variant.ids not on apicoplast
seqSetFilter(my_vcf, variant.id = coords$variant.id[coords$chromosome != "Pf3D7_API_v3"])

fws_all <- getFws(my_vcf)

Fws <- cbind(sample.id, fws_all)
write.table(Fws, 'IsolatesComplexity.txt', col.names=T, row.names=F, quote=F)

Fws=as.matrix(fread('/home/fadel/BalancingSelection/Analysis/Fws/Isolates_Complexity-Fws.txt', header = T))

cols=c('brown', 'darkgreen', 'cadetblue1', 'darkorchid', 'hotpink', 'green3', 'cornflowerblue', 'yellow', 'lightskyblue',
       'gold', 'gray49', 'orange', 'blue', 'darkblue', 'black', 'red')

Colors <- c('#FFFF00', '#FFCC99','#4B0082','#00FFFF','#87CEFA','#008000','#32CD32','#66CDAA','#696969','#000000',
            '#FF66B2','#000080','#238b45', '#FF0000','#663300','#FF800C')
#******************* PLOT *******************
library(data.table)
library(tidyverse)

My_Theme = theme(
    legend.position = 'none',
    axis.line = element_line(color = 'black', size=0.5),
    axis.title.x = element_text(size = 6,color = 'black',face = "bold"),
    axis.title.y = element_text(size = 6,color = 'black',face = "bold"),
    legend.text = element_text(size = 5,face="bold",margin = margin(r = 5, unit = "pt")),
    legend.key.size = unit(0.8, 'lines'),
    legend.spacing.x = unit(0, 'cm'))

fws <- readxl::read_xlsx('data/Isolates_Complexity_Fws.xlsx', )

couleurs <- c("#FFFF00","#EEC591","#912CEE","#00FFFF","#87CEFF","#006400","#7CFC00","#40E0D0","#828282",
              "#000000","#EE7AE9","#0000CD","#FF3399","#EE0000", "#8B3636","#FF7F00")

### FWS 
ggplot(data = fws, aes(x = Pop, y = fws_all, col = Pop)) + 
    geom_jitter(position = position_jitter(0.25), cex=0.6) + 
    stat_summary(fun = mean, show.legend = FALSE, geom = "crossbar", width = 0.6, size = 0.3) +
    scale_color_manual(values = couleurs) + 
    scale_y_continuous(breaks=seq(0.0,1.0,0.1)) +
    theme_minimal() + My_Theme +
    labs(y="Fws", x = "Isolates Origin") + 
    theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black',face = "bold", hjust = 1),
          axis.text.y = element_text(size = 5, color = 'black',face = "bold"))

ggsave("figureS1/figureS1_v1.jpeg", width = 180, height = 150, units = "mm")

### FWS 
ggplot(data=fws, aes(x=Pop, y=fws_all, col=Pop)) + 
    geom_jitter(position=position_jitter(0.25), cex=0.6) + 
    scale_color_manual(values = couleurs) +
    theme_bw() + My_Theme +
    scale_y_continuous(breaks=seq(0.0,1.0,0.1)) +
    labs(y="Fws", x = "Isolates Origin") + 
    theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black',face = "bold", hjust = 1),
          axis.text.y = element_text(size = 5, color = 'black',face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

ggsave("figureS1/figureS1_v2.jpeg", width = 180, height = 150, units = "mm")


# geom_jitter
ggplot(data = fws, aes(x = Pop, y = fws_all, fill = Pop)) + 
    geom_jitter(show.legend = FALSE, width = 0.25, shape = 21, color = 'black') +
    stat_summary(fun = mean, show.legend = FALSE, geom = "crossbar", width = 0.6, size = 0.3) +
    scale_fill_manual(values = couleurs) +
    scale_y_continuous(breaks=seq(0.0,1.0,0.1)) +
    theme_minimal() + My_Theme + 
    labs(y="Fws", x = "Isolates Origin") + 
    theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black',face = "bold", hjust = 1),
          axis.text.y = element_text(size = 5, color = 'black',face = "bold"))

ggsave("figureS1/figureS1_v1.jpeg", width = 180, height = 150, units = "mm")


ggplot(data=fws, aes(x=Pop, y=fws_all, fill=Pop)) + 
    geom_jitter(show.legend = FALSE, width = 0.25, shape = 21, color = 'black') + 
    stat_summary(fun = mean, show.legend = FALSE, geom = "crossbar", width = 0.6, size = 0.3) +
    scale_fill_manual(values = couleurs) +
    theme_bw() + My_Theme +
    scale_y_continuous(breaks=seq(0.0,1.0,0.1)) +
    labs(y="Fws", x = "Isolates Origin") + 
    theme(axis.text.x = element_text(angle = 45, size = 6, color = 'black',face = "bold", hjust = 1),
          axis.text.y = element_text(size = 5, color = 'black',face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

ggsave("figureS1/figureS1_v2.jpeg", width = 180, height = 150, units = "mm")
