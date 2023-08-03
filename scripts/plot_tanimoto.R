#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# specify our desired options in a list
option_list <- list( 
  make_option("--min-jacc-inf", type='character', dest='min_jacc_inf',
              help="Infile with minimum jaccard distances, as produced by chemical_space_UMAP.py"),
  make_option("--plot-prefix", type='character', dest='plot_prefix',
              help="Prefix for plots"))

opt <- parse_args(OptionParser(option_list=option_list))

# Read in the minimum jaccard distances
jacc_dist <- read.delim(opt$min_jacc_inf)

# Plot the distances as an ECDF, with y-axis flipped
p <- jacc_dist %>%
  ggplot(aes(1-min_jacc, colour=comparison)) +
  
  scale_colour_manual(values=c("#f0e442", "#359B73", 'black', 'grey'),
                      name='Screen compound set') +
  geom_line(aes(y = 100*(1 - ..y..)), stat='ecdf') +
  xlab('Tanimoto similarity threshold') +
  ylab('PubChem compounds (%)') +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size =20) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        aspect.ratio = 1,
        axis.line = element_line(colour = "black", linewidth = 0.5),
        legend.position=c(1,1), legend.justification = c("right", "top"))

ggsave(sprintf('%s_tanimoto_background_vs_screen.png', opt$plot_prefix),
       width=6, height=6)
ggsave(sprintf('%s_tanimoto_background_vs_screen.pdf', opt$plot_prefix),
       width=6, height=6)

perc_above_0.5 <- jacc_dist %>% group_by(comparison) %>%
  summarise(above_0.5=100*mean(min_jacc<0.5)) %>%
  mutate(x=0.5)

p2 <- p +
  coord_cartesian(ylim=c(0,5), xlim=c(0.25,0.75)) +
  geom_point(data=perc_above_0.5, aes(x, above_0.5), size=3) +
  geom_text(data=perc_above_0.5,
            aes(0.4, above_0.5,
                label=format(round(above_0.5, 2), nsmall = 2)),
            size=10, hjust=1) +
  geom_vline(xintercept=0.5, linetype=2, colour='grey') +
  theme(legend.position='none') +
  xlab('') +
  ylab('')


ggsave(sprintf('%s_tanimoto_background_vs_screen_zoom.png', opt$plot_prefix),
       width=4, height=4)
ggsave(sprintf('%s_tanimoto_background_vs_screen_zoom.pdf', opt$plot_prefix),
       width=4, height=4)