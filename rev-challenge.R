library(hyperinf)
library(ggpubr)

L = 4
m = matrix(0, ncol=L, nrow=L)
m[lower.tri(m)] = 1
m = rbind(m, c(0, rep(1, L-2), 0))
m = rbind(m, c(0, 0, rep(1, L-3), 0))
#m = rbind(m, c(0, 0, 0, rep(1, L-3)))

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
fit = hyperinf(m)   
plot_hyperinf(fit)
fit.rev = hyperinf(m, reversible=TRUE)
plot_hyperinf(fit.rev)

ggarrange(plot_hyperinf_data(m), 
          plot_hyperinf(fit.rev), 
          plot_hyperinf(fit), 
          nrow=1, widths=c(1,2,2),
          labels=c("A", "B", "C"))

########

library(ape)

tree <- read.tree(text =
                    "(Deep_lineage:0.2,(t1:0.1,t2:0.1,t3:0.1,t4:0.1,t5:0.1,t6:0.1,t7:0.1,t8:0.1):0.1);"
)

plot(tree, show.tip.label = TRUE)

m = matrix(rep(c(1,0,0,0), 8), ncol=4, byrow = TRUE)
m = rbind(matrix(c(0,0,0,1), ncol=4, byrow = TRUE), m)
tree.yes = hyperinf(m, tree)
tree.no = hyperinf(m)

ggarrange(plot_hyperinf_data(m, tree),
          plot_hyperinf(tree.no), 
          plot_hyperinf(tree.yes),
          nrow=1, widths=c(1,2,2),
          labels = c("A", "B", "C"))
