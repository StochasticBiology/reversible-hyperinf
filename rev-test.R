library(phytools)
library(hypertrapsct)
library(hypermk)

set.seed(2)

pull_orders = function(fluxes, L, order=1, w=10000, rev=FALSE) {
  fwds = which(fluxes$To > fluxes$From)
  revs = which(fluxes$From > fluxes$To)
  fluxes$change = 0
  fluxes$change[fwds] = L-log(fluxes$To[fwds]-fluxes$From[fwds], base=2)
  fluxes$change[revs] = L-log(fluxes$From[revs]-fluxes$To[revs], base=2)
  fluxes$states = unlist(lapply(fluxes$From, DecToBin, len=L))
  
  m = matrix(0, ncol = L, nrow=L)
  for(i in 1:nrow(fluxes)) {
    s = as.numeric(strsplit(fluxes$state[i], split="")[[1]])
    if(i %in% fwds) {
      m[fluxes$change[i], which(s==0)] = m[fluxes$change[i], which(s==0)]+fluxes$Flux[i]/w
    }
  }
  return(m)
}
# m_ij = how often i changes before j

pull_orders_hct = function(routes) {
  L = ncol(routes)
  m = matrix(0, ncol = L, nrow=L)
  for(i in 1:nrow(routes)) {
    for(j in 1:L) {
      m[routes[i,j]+1, setdiff((1:L), (routes[i,1:j]+1))] = 
        m[routes[i,j]+1, setdiff((1:L), (routes[i,1:j]+1))] + 1
    }
  }
  m = m / nrow(routes)
  return(m)
}

mk.from.ct = function(my.ct, reversible=FALSE) {
  L = ncol(my.ct$srcs)
  states = 1+apply(my.ct$dests, 1, function(row)
    sum(row * 2^((length(row)-1):0))
  )
  tip.states = states[1:length(my.ct$tree$tip.label)] 
  return(mk.inference( my.ct$tree, L, 
                       use.priors=FALSE, tip.states, 
                       reversible = reversible,
                       Ntrials = 1))
}

model.fits = list()

for(expt in 18:35) {
  
  L = 5
  # parameterisation for tree construction
  tree.size = 32
  birth.rate = 1
  death.rate = 0.1
  expt.type = expt %% 7

  # accumulation rate for features (and loss rate, for reversible setup)
  if(expt.type == 0) {
    set.paths = TRUE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0, L)
    start.loci = c(1,1)
  } else if(expt.type == 1) {
    set.paths = TRUE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0, L)
    start.loci = c(1,L)
  } else if(expt.type == 2) {
    set.paths = TRUE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0.1, L)
    start.loci = c(1,1)
  } else if(expt.type == 3) {
    set.paths = TRUE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0.1, L)
    start.loci = c(1,L)
  } else if(expt.type == 4) {
    set.paths = FALSE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0, L)
    start.loci = 1:L
  } else if(expt.type == 5) {
    set.paths = FALSE
    accumulation.rate = (1:L)/L
    loss.rate = rep(0, L)
    start.loci = 1:L
  } else if(expt.type == 6) {
    set.paths = FALSE
    accumulation.rate = rep(1, L)
    loss.rate = (1:L)/L
    start.loci = 1:L
  }
  #c(L, 1, 2, 3, L)/L
  #accumulation.rate = c(0.05, 0.2, 2, 1)
  #loss.rate = c(0.1, 0.1, 0.1, 0.1)*0
  
  # create random phylogeny with tree.size nodes from birth-death process parameterised as above
  my.tree = ape::rphylo(tree.size, birth=birth.rate, death=death.rate)
  my.tree$node.label = as.character(1:my.tree$Nnode)
  tree.labels = c(my.tree$tip.label, my.tree$node.label)
  
  # generate state for all nodes traversing tree breadth-first
  # and setting the state of the child nodes according to
  # accumulation (and, if reversible, loss.rate) and branch length
  my.root = phangorn::getRoot(my.tree)
  to.do = c(my.root)
  # initialise state list
  x = list()
  x[[my.root]] = rep(0,L)
  # while we still have vertices to simulate
  while(length(to.do) > 0) {
    # initialise a new to-do list for next iteration
    new.to.do = c()
    # loop through each node in current to-do list
    for(i in to.do) {
      this.outgoing.edges = which(my.tree$edge[,1] == i)
      # loop over this node's children
      for(j in this.outgoing.edges) {
        this.child = my.tree$edge[j,2]
        this.branch.length = my.tree$edge.length[j]
        # construct state for this child based on its parent
        x[[this.child]] = x[[i]]
        
        if(set.paths == TRUE) {
        if(sum(x[[this.child]]) == 0) {
          ref = sample(start.loci, 1)
        } else if(x[[this.child]][L] == 0) {
          ref = which(x[[this.child]] == 0)[1]
        } else if(x[[this.child]][1] == 0) {
          ref = rev(which(x[[this.child]] == 0))[1]
        } else {
          ref = 1
        }
        } else {
          poss = which(x[[this.child]] == 0)
          if(length(poss) == 1) {
            ref = poss
          } else {
            ref = sample(which(x[[this.child]] == 0), 1, prob=accumulation.rate[which(x[[this.child]] == 0)])
          }
        }
        # 00001. 10000
        if(runif(1) < accumulation.rate[ref]*this.branch.length) { x[[this.child]][ref] = 1 } 
        for(ref in 1:L) {
          if(runif(1) < loss.rate[ref]*this.branch.length) { x[[this.child]][ref] = 0 } 
        }
        # in the reversible case, allow the leftmost feature ("first feature" in the ms.:
        # second paragraph of "Synthetic case studies") to revert with some probability
        
        # XXXX
        
        # add this child to to state list, and to next iteration's to-do
        new.to.do = c(new.to.do, this.child)
      }
    }
    # update to-do list
    to.do = new.to.do
  }
  
  my.tree$tip.label = x[1:length(my.tree$tip.label)]
  tip.states = unlist(lapply(my.tree$tip.label, BinToDec))+1
  
  ddf = cbind(data.frame(Species=as.character(1:tree.size)),
              as.data.frame(matrix(unlist(x[1:length(my.tree$tip.label)]), ncol=L, byrow = TRUE)))
  
  my.tree2 = my.tree 
  my.tree2$tip.label = as.character(1:tree.size)
  my.ct = curate.tree(my.tree2, ddf)
  
  my.hct.out = HyperTraPS(my.ct$dests, my.ct$srcs,
                          length = 4, kernel = 4)
  
  phct= plotHypercube.sampledgraph2(my.hct.out, node.labels = FALSE, no.time=TRUE, edge.label.size = 3) +
    theme(legend.position = "none")
  
  my.hct.cs.out = HyperTraPS(my.ct$dests, 
                          length = 4, kernel = 4)
  
  phct.cs= plotHypercube.sampledgraph2(my.hct.cs.out, node.labels = FALSE, no.time=TRUE, edge.label.size = 3) +
    theme(legend.position = "none")
  
  #plotHypercube.bubbles(my.hct.out)
 # barcodes = unlist(lapply(tip.states-1, DecToBin, L))
#  my.tree2$tip.label = barcodes
  
#  data.plot = ggtree(my.tree2, layout="circular") + geom_tiplab2(size=3)
#  data.plot
  
  mk.out.irrev = mk.inference(my.tree, L, 
                              use.priors=FALSE, tip.states, 
                              reversible = FALSE,
                              Ntrials = 1)
  mk.out.rev = mk.inference(my.tree, L, 
                            use.priors=FALSE, tip.states, 
                            reversible = TRUE,
                            Ntrials = 1)
  
  m.irrev = pull_orders(mk.out.irrev$mk_fluxes, L)
  m.rev = pull_orders(mk.out.rev$mk_fluxes, L)
  m.rev = apply(m.rev, c(1,2), min, 1)
  m.hct = pull_orders_hct(my.hct.out$routes)
  m.hct.cs = pull_orders_hct(my.hct.cs.out$routes)
  diag(m.rev) = NA
  diag(m.irrev) = NA
  diag(m.hct) = NA
  diag(m.hct.cs) = NA
  
  model.fits[[length(model.fits)+1]] = list(my.ct,
                                            my.hct.out, my.hct.cs.out,
                                            mk.out.irrev, mk.out.rev)
  df.irrev <- as.data.frame(as.table(m.irrev))
  df.rev = as.data.frame(as.table(m.rev))
  df.hct = as.data.frame(as.table(m.hct))
  df.hct.cs = as.data.frame(as.table(m.hct.cs))
  
  p.irrev = ggplot(df.irrev, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile() 
  p.rev = ggplot(df.rev, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile()
  p.hct = ggplot(df.hct, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile()
  p.hct.cs = ggplot(df.hct.cs, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile()
  
  sf = 2
  fname = paste0("tester-", L, "-", tree.size, "-", expt, ".png", collapse="")
  png(fname, width=1200*sf, height=600*sf, res=72*sf)
  print(ggarrange(plotHypercube.curated.tree(my.ct), ggarrange(mk.inference.plot(mk.out.irrev),
                                 mk.inference.plot(mk.out.rev),
                                 phct, phct.cs,
                                 p.irrev, p.rev, p.hct, p.hct.cs, ncol=4, nrow=2), widths=c(1,2.5)))
  dev.off()
}

save(model.fits, file="model-fits.Rdata")

m.list = list()
m.names = c()
for(i in 1:length(model.fits)) {
  m.list[[4*(i-1)+1]] = pull_orders(model.fits[[i]][[4]]$mk_fluxes, L)
  m.list[[4*(i-1)+2]] = pull_orders(model.fits[[i]][[5]]$mk_fluxes, L)
  m.list[[4*(i-1)+2]] = apply(m.list[[4*(i-1)+2]], c(1,2), min, 1)
  m.list[[4*(i-1)+3]] = pull_orders_hct(model.fits[[i]][[2]]$routes)
  m.list[[4*(i-1)+4]] = pull_orders_hct(model.fits[[i]][[3]]$routes)
  m.names = c(m.names, paste(i-1, c("ht", "htcs", "mkirr", "mkrev")))
  diag(m.list[[4*(i-1)+1]]) = NA
  diag(m.list[[4*(i-1)+2]]) = NA
  diag(m.list[[4*(i-1)+3]]) = NA
  diag(m.list[[4*(i-1)+4]]) = NA
}

refs = 1:(length(model.fits)*4)
refs = refs[refs%%4 != 2]

X <- do.call(
  rbind,
  lapply(m.list[refs], function(m) as.vector(m))
)

rownames(X) <- m.names[refs]
keep_cols <- !is.na(X[1, ])
X_clean <- X[, keep_cols]
pca <- prcomp(X_clean, center = TRUE, scale. = FALSE)

df_pca <- data.frame(
  name = m.names[refs],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

df_pca$group = rep(1:length(model.fits), each = 3)-1

hulls <- df_pca %>%
  group_by(group) %>%
 # filter(n() >= 3) %>%          # need ≥ 3 points for a polygon
  slice(chull(PC1, PC2))

library(ggrepel)
pca.plot = ggplot(df_pca, aes(PC1, PC2, colour = factor(group))) +
  geom_point(size = 1) +
  geom_text_repel(aes(x=PC1, y=PC2, label=name)) +
  geom_polygon(
    data = hulls,
    aes(fill = factor(group)),
    alpha = 0.2,
    colour = NA
  ) +
  theme_minimal() +
  labs(
    colour = "Group",
    fill   = "Group"
  )


library(dplyr)

spread_df <- df_pca %>%
  group_by(group) %>%
  summarise(
    n = n(),
    spread = mean(
      sqrt((PC1 - mean(PC1))^2 + (PC2 - mean(PC2))^2)
    )
  )

expt.names = c("single path", "two path", "one path\nloss", "two path\nloss", "random", "random\nordered", "random\nordered loss")
spread.plot = ggplot(spread_df, aes(x = expt.names[(group%%7)+1], y = spread, label=group)) +
  geom_boxplot() + geom_point() + geom_text_repel() +
  labs(
    x = "",
    y = "Mean distance to centroid",
    title = "Within-group dispersion in PCA space"
  ) +
  theme_minimal()

sf = 2
png("outputs.png", width=800*sf, height=800*sf, res=72*sf)
ggarrange(pca.plot + theme(legend.position = "none"), spread.plot, nrow=2)
dev.off()


