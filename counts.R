library(readr)
library(dplyr)
library(ggplot2)
library(grid)

# Read data containing the number of matches in each replicate
df <- read_tsv("results/summary.tsv") %>%
  mutate(bio_rep=as.factor((replicate-1)%%6 + 1),
         tech_rep=as.factor(floor((replicate-1)/6) + 1)) %>%
  mutate(matches_per_mil=matches/total*1e6)

# Set global plotting parameters
theme_set(theme_bw())
theme_update(strip.placement="outside", strip.background=element_blank(),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank())

total_cutoff <- 1e5
df_agg <-
  df %>%
    group_by(subject) %>%
    filter(total>total_cutoff) %>%
    summarize(median_matches_per_mil=median(matches_per_mil))

total_scale <- function(name) {
  scale_fill_viridis_c(limits=c(1, max(df$total)), trans="log10",
                       breaks=c(1, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 5*(10^6)),
                       labels=c("1", "10", expression(10^2), expression(10^3),
                                expression(10^4), expression(10^5),
                                expression(10^6), expression(5%*%10^6)),
                       name=name)
}

# Plot the overall data
p <- ggplot(df, aes(x=bio_rep)) +
  geom_hline(data=df_agg, aes(yintercept=median_matches_per_mil), color="gray") +
  geom_point(shape=21, position=position_dodge(0.5),
             aes(y=matches_per_mil, fill=total, group=tech_rep)) +
  scale_y_continuous(limits=c(0, 100), n.breaks=5) +
  total_scale("Total number of reads") +
  facet_wrap(vars(subject), nrow=2, labeller=as_labeller(function(x) paste("Subject", x))) +
  xlab("Biological replicate") + 
  ylab("VH1-18 QXXV reads per million") +
  theme(legend.position="bottom", legend.key.width=unit(1.5, "in"),
        plot.title=element_text(size=11, hjust=0.5))
# Add dividers between biological replicates (https://stackoverflow.com/a/61550522)
for (i in 1:6)
  p <- p + annotation_custom(linesGrob(x=unit(c((1/6)*i, (1/6)*i), "npc"),
                                       y=unit(c(0,1), "npc"),
                                       gp=gpar(lty="dashed", lwd=0.5)))
ggsave(paste0("results/counts.pdf"), p, width=10, height=5)

# Compute statistics
print(df %>% filter(total>=total_cutoff) %>% count())
print(
  df %>%
    filter(total>total_cutoff) %>%
    summarize(median_matches_per_mil=median(matches_per_mil),
              iqr_matches_per_mil=quantile(matches_per_mil, 0.75) - quantile(matches_per_mil, 0.25))
)
write_tsv(df_agg, "results/count_medians.tsv")

# Boxplot the overall data together
set.seed(42)
ggplot(df, aes(x=0)) +
  geom_boxplot(data=df %>% filter(total>total_cutoff),
               aes(y=matches_per_mil), outlier.shape=NA) +
  geom_jitter(aes(y=matches_per_mil, fill=total), shape=21) +
  total_scale("Total reads") +
  ylab("VH1-18 QXXV reads per million") +
  guides(fill=guide_colourbar(label.position="left")) +
  theme(legend.key.width=unit(0.25, "in"),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste0("results/counts_boxplot.pdf"), width=3.5, height=3)
