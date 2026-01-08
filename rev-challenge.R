library(hyperinf)
library(ggpubr)
library(hypertrapsct)
library(ape)

#### Fig 2 -- comparing reversible and irreversible fits
L = 4
# data involving a clear pathway with a reversible features
m = matrix(c(1,0,0,0,
             1,1,0,0,
             1,1,1,0,
             1,1,1,1,
             1,0,0,0,
             1,1,0,0,
             1,1,1,0,
             1,1,1,1,
             0,1,0,0,
             0,1,1,0), ncol=4, byrow = TRUE)
# HyperHMM and HyperMk fits
fit = hyperinf(m)   
plot_hyperinf(fit)
fit.rev = hyperinf(m, reversible=TRUE)
plot_hyperinf(fit.rev)

# plot the comparison
ggarrange(plot_hyperinf_data(m), 
          plot_hyperinf(fit.rev), 
          plot_hyperinf(fit), 
          nrow=1, widths=c(1,2,2),
          labels=c("A", "B", "C"))

###### Fig 6 looking at inferred interactions

# simple case, single pathway
m = matrix(rep(c(1,0,0,0,
              1,1,0,0,
              1,1,1,0,
              1,1,1,1), 10), ncol=4, byrow = TRUE)
fit.ht = hyperinf(m, method="hypertraps", length=4, penalty=0.5)
plot_hyperinf(fit.ht)
plotHypercube.influencegraph(fit.ht)

# produce dataframe of posterior means and plot
post.df = data.frame(i = 0:(L*L-1), m = colMeans(fit.ht$posterior.samples))
l1 = floor(post.df$i / L)
l2 = post.df$i %% L
post.df$label = paste(l1+1, "->", l2+1, sep="")
#ggplot(post.df, aes(x=label, y=m)) + geom_point()

# harder case, single pathway plus reversible feature 1
m.rev = matrix(rep(c(1,0,0,0,
                 1,1,0,0,
                 1,1,1,0,
                 1,1,1,1,
                 0,1,0,0,
                 0,1,1,0), 10), ncol=4, byrow = TRUE)
fit.ht.rev = hyperinf(m.rev, method="hypertraps", length=4, penalty=0.5)
plot_hyperinf(fit.ht.rev)
plotHypercube.influencegraph(fit.ht.rev)

# produce dataframe of posterior means and plot
post.df.rev = data.frame(i = 0:(L*L-1), m = colMeans(fit.ht.rev$posterior.samples))
l1.rev = floor(post.df.rev$i / L)
l2.rev = post.df.rev$i %% L
post.df.rev$label = post.df$label
#ggplot(post.df.rev, aes(x=label, y=m)) + geom_point()

# set up plot comparing the outputs
post.df$expt = "i"
post.df.rev$expt = "ii"
both.df = rbind(post.df, post.df.rev)
ggarrange(
  ggarrange(plot_hyperinf_data(m[1:4,]), 
            plot_hyperinf_data(m.rev[1:6,]), nrow=2,
            labels=c("   i", "ii")),
  ggarrange(plot_hyperinf(fit.ht),
            plot_hyperinf(fit.ht.rev),
            nrow=2),
ggplot(both.df, aes(x=label, y=m, color=expt)) + geom_point() +
  labs(x = "Parameter", y = "Posterior mean", color="Dataset") + 
  theme_minimal(),
nrow=1, widths=c(0.4,0.4,1), labels=c("A", "B", "C")
)


######## Fig 5 -- accounting for vs ignoring phylogeny

# explicitly encode a tree
tree <- read.tree(text =
                    "(Deep_lineage:0.2,(t1:0.1,t2:0.1,t3:0.1,t4:0.1,t5:0.1,t6:0.1,t7:0.1,t8:0.1):0.1);"
)

plot(tree, show.tip.label = TRUE)

# explicit encode this case studies' data
m = matrix(rep(c(1,0,0,0), 8), ncol=4, byrow = TRUE)
m = rbind(matrix(c(0,0,0,1), ncol=4, byrow = TRUE), m)
colnames(m) = 1:4

# do inference with and without phylogeny
tree.yes = hyperinf(m, tree)
tree.no = hyperinf(m)

# comparison plot
ggarrange(plot_hyperinf_data(m, tree),
          plot_hyperinf(tree.no), 
          plot_hyperinf(tree.yes),
          nrow=1, widths=c(1,2,2),
          labels = c("A", "B", "C"))


######## Fig 1 -- illustrating concepts

# simple observation matrix
m1 = matrix(c(1,0,0,0,
             1,1,0,0,
             1,1,1,0,
             1,1,1,1), ncol=4, byrow = TRUE)
colnames(m1) = 1:4
# two pathway observations
m2 = rbind(m1, 1-m1)
# random tree
tree = rtree(4)
# it's easier to restyle HyperTraPS native plots
# so use HyperTraPS for the first example plot
tmp = hyperinf(m1, tree, method="hypertraps")
# produce plot with big node labels for illustration
demo.graph = plotHypercube.sampledgraph2(tmp, thresh=0.2,
                                         no.times = TRUE,
                                         node.label.size = 5,
                                         edge.label.size = 3,
                                         featurenames = 1:4)
# compile collection of illustrations
ggarrange(plot_hyperinf_data(m1), 
          plot_hyperinf_data(m2), 
          plot_hyperinf_data(m1, tree), 
          demo.graph,
          plot_hyperinf(hyperinf(m2)),
          plot_hyperinf(hyperinf(m1, tree)),
          nrow=2, ncol=3, heights=c(0.5,1),
          labels=c("A i", "B i", "C i", "ii", "ii", "ii"))


##########


