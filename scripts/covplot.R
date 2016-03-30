library(ggplot2)
library(ggthemes)
args = commandArgs(trailingOnly = TRUE)

depth = read.delim(args[1], header=FALSE)

p = ggplot(data=depth, aes(x=V2, y=V3)) +
  geom_bar(stat='identity') +
  xlab('genome position') +
  ylab('coverage [reads]') +
  ggtitle(paste0('Sequence accession number: ', args[2])) +
  scale_y_log10() +
  theme_solarized_2()
ggsave(args[3], width=297, height=210, units='mm')
