library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
depth_file = args[1]
acc = args[2]
seq_len = as.numeric(args[3])
image_name = args[4]
n_reads = as.numeric(args[5])

cov_exp = round(seq_len * (1 - exp(-n_reads * 151 / seq_len)), 0)

depth = read.delim(depth_file, header=FALSE)

cov_obs = sum(depth$V3 > 0)

p = ggplot(data=depth, aes(x=V2, y=V3)) +
  # geom_bar(stat='identity') +
  geom_line() +
  xlab('genome position') +
  ylab('coverage [reads]') +
  ggtitle(paste0('Sequence accession number: ', acc, '\nCovered bases: expected ', cov_exp, ' - observed ', cov_obs)) +
  scale_y_log10() +
  expand_limits(x=c(0, seq_len), y=c(0, 10000))

if(is.element("ggthemes", installed.packages()[,1])){
    library(ggthemes)
    p = p + theme_economist()
}

ggsave(image_name, width=297, height=210, units='mm')
