library(ggplot2)

setwd("/Users/ianbrennan/Documents/GitHub/Pygopodidae")

astats <- read.csv("Sampling/Pygopodidae_SampleInfo.csv")

#ingroup <- dplyr::filter(astats, Source == "AusARG_Agamidae")

testo <- tidyr::pivot_longer(data=astats, cols=c("AHE","gene","uce"))

ggplot(testo, aes(fill=name, y=value, x=lineage)) + 
  geom_bar(position="stack", stat="identity") + scale_x_discrete(limits=rev) + 
  coord_flip() + theme_classic(base_size = 6)
