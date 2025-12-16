library(hypertrapsct)
library(hypermk)
library(ggpubr)
library(ggplot2)

kp.dir = "~/Dropbox/klebevo/kp-evolution-inference-curate-wip/kleborate-analysis/"
w.dir = "/Users/iain/Dropbox/Documents/2026_Projects/Hyper/hypermk-main/"

setwd(kp.dir)

countries <- gsub("\\.nwk", "", list.files("clean", pattern=".*\\.nwk"))

c.data = data.frame()
for(country in countries) {
  tree.path <- paste0("clean/",country,".nwk")
  if (!file.exists(tree.path)) {
    stop(paste("No newick-tree for", country, "in clean directory"))
  }
  
  if (!file.exists("clean/kleborate-dichotomized.csv")) {
    stop("Run preprocess_kleborate.R first!")
  }
  
  resistance.df <- read.csv("clean/kleborate-dichotomized.csv")
  featurenames <- setdiff(colnames(resistance.df), "id")
  
  ctree <- curate.tree(tree.path, "clean/kleborate-dichotomized.csv")
  
  data.df = ctree$data[ctree$data$label %in% ctree$tree$tip.label,]
  colsums = colSums(data.df[,2:ncol(data.df)])
  nonzeros = which(colsums != 0 & colsums != nrow(data.df))
  if(length(nonzeros) > 1) {
    prune.data.df = data.df[,c(1,1+nonzeros)]
    c.data = rbind(c.data, data.frame(country=country, features=ncol(prune.data.df)-1))
  }
}

c.interest = c.data$country[c.data$features <= 6]

for(country in c.interest) {
  message(country)
  tree.path <- paste0("clean/",country,".nwk")
  if (!file.exists(tree.path)) {
    stop(paste("No newick-tree for", country, "in clean directory"))
  }
  
  if (!file.exists("clean/kleborate-dichotomized.csv")) {
    stop("Run preprocess_kleborate.R first!")
  }
  
  resistance.df <- read.csv("clean/kleborate-dichotomized.csv")
  featurenames <- setdiff(colnames(resistance.df), "id")
  
  ctree <- curate.tree(tree.path, "clean/kleborate-dichotomized.csv")
  
  plotHypercube.curated.tree(ctree)
  data.df = ctree$data[ctree$data$label %in% ctree$tree$tip.label,]
  colsums = colSums(data.df[,2:ncol(data.df)])
  nonzeros = which(colsums != 0 & colsums != nrow(data.df))
  if(length(nonzeros) > 1) {
    prune.data.df = data.df[,c(1,1+nonzeros)]
    c.data = rbind(c.data, data.frame(country=country, features=ncol(prune.data.df)-1))
  }
  
  df_ordered <- prune.data.df[match(ctree$tree$tip.label, prune.data.df$label), ]
  
  ctree$tree$tip.label == df_ordered$label
  mk.data = as.matrix(df_ordered[,2:ncol(df_ordered)])
  rownames(mk.data) <- NULL
  colnames(mk.data) <- NULL
  
  fit.ht <- HyperTraPS(ctree$dests, # set to one thousand walkers
                       initialstates=ctree$srcs,
                       length = 5,
                       walkers = 200,
                       penalty = 1,
                       seed = 1,
                       featurenames=featurenames)
  
  fit.mk = mk_prune_model(mk_infer_phylogenetic(mk.data, ctree$tree))
  fit.mk.irrev = mk_prune_model(mk_infer_phylogenetic(mk.data, ctree$tree, reversible = FALSE))
  
  presences = colSums(df_ordered[2:ncol(df_ordered)])
  disp.plot = data.frame(names=names(presences), counts=as.vector(presences))
  
  sf = 2
  png(paste0(w.dir, "output-", country, ".png", collapse=""), width=800*sf, height=600*sf, res=72*sf)
  print(ggarrange(
    ggarrange(plotHypercube.curated.tree(ctree, font.size = 2),
              plotHypercube.bubbles(fit.ht)+theme(legend.position="none"),
              nrow=1),
    ggarrange(
      ggplot(disp.plot, aes(x=names, y=counts)) + geom_col() + theme(axis.text.x = element_text(hjust=1,angle=45)),
      #plotHypercube.sampledgraph2(fit.ht, node.labels = FALSE, no.times = TRUE)
      mk.inference.plot(fit.mk),
      mk.inference.plot(fit.mk.irrev), nrow=1),
    nrow=2
  ))
  dev.off()
}
  