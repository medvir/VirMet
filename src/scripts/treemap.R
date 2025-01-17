#!/usr/bin/env Rscript

library("treemap")
library("grid")
library("colorspace")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

cmargs = commandArgs(TRUE)
if(identical(cmargs, character(0))){
  cmargs = c("hier_stats.csv")
}

hier <- read.csv(cmargs[1])
tm <- treemap(hier, index=c('category', 'id'), vSize='reads',
        algorithm='squarified', vColor='id', type='index', palette=cbPalette,
        align.labels = list(c("right", "center"), c("left", "bottom")),
        overlap.labels=0.05,
        title="Origin of the reads", fontsize.title=18,
        fontface.labels=2, fontsize.labels=14,
        aspRatio=1)

# This does not work in scripting
#lapply(tail(grid.ls(print=FALSE)$name, 2), grid.remove)
