
library(ggplot2)
library(ggpubr)
library(rstatix)

fws <- readxl::read_xlsx('data/Isolates_Complexity_Fws.xlsx')

#=================================
# SIGNIFICANT TEST ANALYSIS
#=============================
res.aov <- fws %>%
    anova_test(fws_all ~ Pop)

pwc <- fws %>%
    pairwise_t_test(fws_all ~ Pop, p.adjust.method = "bonferroni")


fws$Pop <- factor(fws$Pop)
pairwise.wilcox.test(fws$fws_all, fws$Pop, p.adjust.method = "bonferroni")

model <- lm(fws_all~Pop, data=fws)
anova(model)

summary(glht(model, mcp(Pop = "Tukey")))

# Show adjusted p-values
pwc <- pwc %>% add_xy_position(x = "Pop")

ggboxplot(fws, x = "Pop", y = "fws_all") +
    stat_pvalue_manual(pwc, label = "p.adj", tip.length = 0, step.increase = 0.1) +
    labs(
        subtitle = get_test_label(res.aov, detailed = TRUE),
        caption = get_pwc_label(pwc))

couleurs <- c("#FFFF00","#EEC591","#912CEE","#00FFFF","#87CEFF","#006400","#7CFC00","#40E0D0","#828282",
              "#000000","#EE7AE9","#0000CD","#FF3399","#EE0000", "#8B3636","#FF7F00")

My_Theme = theme(
    legend.position = 'none',
    axis.line = element_line(color = 'black', size=0.5),
    axis.title.x = element_text(size = 9,color = 'black',face = "bold"),
    axis.title.y = element_text(size = 9,color = 'black',face = "bold"),
    legend.text = element_text(size = 5,face="bold",margin = margin(r = 5, unit = "pt")),
    legend.key.size = unit(0.8, 'lines'),
    legend.spacing.x = unit(0, 'cm'))

# Create scatter plot
# geom_jitter
ggplot(data = fws, aes(x = Pop, y = fws_all, fill = Pop)) + 
    geom_jitter(show.legend = FALSE, width = 0.35, shape = 21, color = 'black') +
    stat_summary(fun = median, show.legend = FALSE, geom = "crossbar", width = 0.8, size = 0.35) +
    scale_fill_manual(values = couleurs) +
    scale_y_continuous(breaks=seq(0.0,1.0,0.1)) +
    theme_minimal() + My_Theme + 
    labs(y="Fws", x = "Isolates Origin") + 
    theme(axis.text.x = element_text(angle = 45, size = 8, color = 'black',face = "bold", hjust = 1),
          axis.text.y = element_text(size = 7, color = 'black',face = "bold"))

# Save plot
ggsave("balancing_selection_signatures/figureS1/figureS1_v1.png", width = 180, height = 150, units = "mm")


ggplot(data=fws, aes(x=Pop, y=fws_all, fill=Pop)) + 
    geom_jitter(show.legend = FALSE, width = 0.25, shape = 21, color = 'black') + 
    stat_summary(fun = mean, show.legend = FALSE, geom = "crossbar", width = 0.8, size = 0.35) +
    scale_fill_manual(values = couleurs) +
    theme_bw() + My_Theme +
    scale_y_continuous(breaks=seq(0.0,1.0,0.1)) +
    labs(y="Fws", x = "Isolates Origin") + # , subtitle = get_test_label(res.aov, detailed = TRUE), caption = get_pwc_label(pwc)
    theme(axis.text.x = element_text(angle = 45, size = 8, color = 'black',face = "bold", hjust = 1),
          axis.text.y = element_text(size = 7, color = 'black',face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

ggsave("balancing_selection_signatures/figureS1/figureS1_v2.png", width = 180, height = 150, units = "mm")


