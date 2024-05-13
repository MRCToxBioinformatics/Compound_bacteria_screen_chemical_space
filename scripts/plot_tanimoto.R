#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# specify our desired options in a list
option_list <- list( 
  make_option("--min-jacc-inf", type='character', dest='min_jacc_inf',
              help="Infile with minimum jaccard distances, as produced by chemical_space_UMAP.py"),
  make_option("--min-jacc-pol-drug-inf", type='character', dest='min_jacc_pol_drug_inf',
              help="Infile with minimum jaccard distances between pollutants and drugs, as produced by chemical_space_UMAP.py"),
  make_option("--plot-prefix", type='character', dest='plot_prefix',
              help="Prefix for plots"))

opt <- parse_args(OptionParser(option_list=option_list))

# Read in the minimum jaccard distances
jacc_dist <- read.delim(opt$min_jacc_inf)
jacc_dist_pol_drug <- read.delim(opt$min_jacc_pol_drug_inf)

p <- jacc_dist_pol_drug %>%
  ggplot(aes(Tanimoto)) +
  geom_histogram() +
  geom_vline(xintercept=0.75, linetype=2, colour='grey') +
  xlab('Pollutants vs Pharmaceutical drugs\nmaximum Tanimoto similarity') +
  ylab('Count') +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size =20) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        aspect.ratio = 1,
        axis.line = element_line(colour = "black", linewidth = 0.5),
        legend.position=c(1,1), legend.justification = c("right", "top"))

ggsave(sprintf('%s_tanimoto_pollutant_vs_drug.png', opt$plot_prefix),
       width=6, height=6, plot=p)
ggsave(sprintf('%s_tanimoto_pollutant_vs_drug.pdf', opt$plot_prefix),
       width=6, height=6, plot=p)


# Plot the distances as an ECDF, with y-axis flipped
p <- jacc_dist %>%
  ggplot(aes(1-min_jacc, colour=comparison)) +
  
  scale_colour_manual(values=c("#f0e442", "#359B73", 'black', 'grey'),
                      name='Screen compound set') +
  ggrastr::rasterise(geom_line(aes(y = 100*(1 - ..y..)), stat='ecdf'),dpi=400) +
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


# from https://stackoverflow.com/questions/27233738/finding-location-of-maximum-d-statistic-from-ks-test
compare <- function(x, y) {
  n <- length(x); m <- length(y)
  w <- c(x, y)
  o <- order(w)
  z <- cumsum(ifelse(o <= n, m, -n))
  i <- which.max(abs(z))
  w[o[i]]
}

all_jacc <- jacc_dist %>% filter(comparison=='All') %>% pull(min_jacc) * 100
drugs_jacc <- jacc_dist %>% filter(comparison=='Pharmaceutical drugs') %>% pull(min_jacc) * 100

u <- compare(all_jacc, drugs_jacc)
e.x <- ecdf(all_jacc)
e.y <- ecdf(drugs_jacc)

D_yend = 100*(e.y(u))
D_ystart = 100*(e.x(u))
D_x = 1-(u/100)

ks_res <- ks.test(all_jacc, drugs_jacc)

p1 <- p + geom_segment(x=D_x, xend=D_x, y=D_ystart, yend=D_yend,
                       colour="#d55e00", linewidth=0.25) +
  annotate(geom='text', label=sprintf('KS-test:\nD = %s\np%s',
                                      round(ks_res$statistic,3),
                                      ifelse(ks_res$p.value>2.2e-16, paste0(' = ', ks_res$p.value), ' < 2.2e-16')),
           x=0.5,
           y=mean(c(D_ystart, D_yend)),
           vjust=0.5, hjust=0,
           colour="#d55e00",
           size=5)

ggsave(sprintf('%s_tanimoto_background_vs_screen.png', opt$plot_prefix),
       width=6, height=6, plot=p1)
ggsave(sprintf('%s_tanimoto_background_vs_screen.pdf', opt$plot_prefix),
       width=6, height=6, plot=p1)

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
       width=4, height=4, plot=p2)
ggsave(sprintf('%s_tanimoto_background_vs_screen_zoom.pdf', opt$plot_prefix),
       width=4, height=4, plot=p2)