library(ggplot2)
if(is.element("ggthemr", installed.packages()[,1])){
    library(ggthemr)
    ggthemr('solarized')
}

args = commandArgs(trailingOnly = TRUE)
depth = read.delim(args[1], header=FALSE)

p = ggplot(data=depth, aes(x=V2, y=V3)) +
  geom_bar(stat='identity') +
  xlab('genome position') +
  ylab('coverage [reads]') +
  ggtitle(paste0('Sequence accession number: ', args[2])) +
  scale_y_log10()

ggsave(args[3], width=297, height=210, units='mm')
