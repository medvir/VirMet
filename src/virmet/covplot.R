#!/usr/bin/env Rscript

library(ggplot2)
oldw <- getOption("warn")
options(warn = -1)

args = commandArgs(trailingOnly = TRUE)
depth_file = args[1]
acc = args[2]
seq_len = as.numeric(args[3])
image_name = args[4]
n_reads = as.numeric(args[5])

cov_exp = round(seq_len * (1 - exp(-n_reads * 151 / seq_len)), 0)
perc_exp = round(100 * cov_exp / seq_len, 0)

depth = tryCatch(
  {
    read.delim(depth_file, header=FALSE)
  },
  error = function(cond) {
    message("File seems empty")
    message("Here's the original error message:")
    message(cond)
    message('')
    # Choose a return value in case of error
    quit(status=0)
    return(NA)
  }
)
cov_obs = tryCatch(
  {
    sum(depth$V3 > 0)
  },
  error = function(cond){
    message("Error here")
    message(cond)
    quit(stats=0)
    return(0)
  }
)
perc_obs = round(100 * cov_obs / seq_len, 0)

line_1 <- sprintf('Reference: %s length: %d', acc, seq_len)
line_2 <- sprintf('Covered bases: expected %d (%d%% of the reference) - observed %d (%d%% of the reference)',
                  cov_exp, perc_exp, cov_obs, perc_obs)

p = ggplot(data=depth, aes(x=V2, y=V3)) +
  # geom_bar(stat='identity') +
  geom_line() +
  xlab('genome position') +
  ylab('coverage [reads]') +
  ggtitle(label = line_1, subtitle = line_2) +
  scale_y_log10() +
  expand_limits(x=c(0, seq_len), y=c(0, 10000))

if(is.element("ggthemes", installed.packages()[,1])){
    library(ggthemes)
    p = p + theme_economist()
}

ggsave(image_name, width=297, height=210, units='mm')

options(warn = oldw)
print(perc_obs)
