library(numbat)
library(ggplot2)

# JS2
obj = Numbat$new("data/JS2/output", i=1, gtf=gtf_hg38, verbose=TRUE)
heatmap = obj$plot_phylo_heatmap(line_width=1, tvn_line=FALSE)
ggsave("data/JS2/output/heatmap.png", width=9, height=5)
ggsave("data/JS2/output/heatmap.pdf", width=9, height=5)

# AH2
obj = Numbat$new("data/AH2/output", i=1, gtf=gtf_hg38, verbose=TRUE)
heatmap = obj$plot_phylo_heatmap(line_width=1, tvn_line=FALSE)
ggsave("data/AH2/output/heatmap.png", width=9, height=5)
ggsave("data/AH2/output/heatmap.pdf", width=9, height=5)

# AH1
obj = Numbat$new("data/AH1/output", i=1, gtf=gtf_hg38, verbose=TRUE)
heatmap = obj$plot_phylo_heatmap(line_width=2, tvn_line=FALSE)
ggsave("data/AH1/output/heatmap.png", width=9, height=5)
ggsave("data/AH1/output/heatmap.pdf", width=9, height=5)
