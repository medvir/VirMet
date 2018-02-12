library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
depth = read.delim(args[1], header=FALSE)

p = ggplot(data=depth, aes(x=V2, y=V3)) +
  # geom_bar(stat='identity') +
  geom_line() +
  xlab('genome position') +
  ylab('coverage [reads]') +
  ggtitle(paste0('Sequence accession number: ', args[2])) +
  scale_y_log10() +
  expand_limits(x=c(0, as.numeric(args[3])), y=c(0, 10000))

library(ggthemes)
if(is.element("ggthemes", installed.packages()[,1])){
    library(ggthemes)
    p = p + theme_solarized()
}

ggsave(args[4], width=297, height=210, units='mm')
