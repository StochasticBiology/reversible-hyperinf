library(phytools)
library(hyperhmm)
library(hypertrapsct)
library(hypermk)
library(hyperinf)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(viridis)

# current Fig 3 is the collection of megaplots around l296 (megaplots from l269, in loop)
# current Fig 4 is (B) the PCA plot around l479 and (A) the errors plot around l541

set.seed(2)

# binary to decimal function
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

# decimal to binary function
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

# extract information on relative feature ordering from a fitted model
# w is the scaling factor for fluxes (as some methods report walker counts, not probabilities)
pull_orders = function(fluxes, L, order=1, w=10000, rev=FALSE) {
  # get fwd and rev transitions
  fwds = which(fluxes$To > fluxes$From)
  revs = which(fluxes$From > fluxes$To)
  # add states and locus changed to dataframe
  fluxes$change = 0
  fluxes$change[fwds] = L-log(fluxes$To[fwds]-fluxes$From[fwds], base=2)
  fluxes$change[revs] = L-log(fluxes$From[revs]-fluxes$To[revs], base=2)
  fluxes$states = unlist(lapply(fluxes$From, DecToBin, len=L))
  
  # initialise matrix
  m = matrix(0, ncol = L, nrow=L)
  # go through dataframe (inefficiently)
  for(i in 1:nrow(fluxes)) {
    # add to matrix elements where this change occurs and other features are zero
    s = as.numeric(strsplit(fluxes$state[i], split="")[[1]])
    if(i %in% fwds) {
      m[fluxes$change[i], which(s==0)] = m[fluxes$change[i], which(s==0)]+fluxes$Flux[i]/w
    }
  }
  return(m)
}
# m_ij = how often i changes before j

# extract ordering information from HyperTraPS output
pull_orders_hct = function(routes) {
  L = ncol(routes)
  # initialise matrix
  m = matrix(0, ncol = L, nrow=L)
  # go through list of simulated routes
  for(i in 1:nrow(routes)) {
    for(j in 1:L) {
      # increase matrix elements corresponding to absence at point of acquisition
      m[routes[i,j]+1, setdiff((1:L), (routes[i,1:j]+1))] = 
        m[routes[i,j]+1, setdiff((1:L), (routes[i,1:j]+1))] + 1
    }
  }
  m = m / nrow(routes)
  return(m)
}

# wrapper function doing HyperMk inference from a curated.tree output
# perhaps obsolete?
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

# initialise lists for output
model.fits = list()
m.list = list()
m.names = c()
plot.list = list()

# key degrees of freedom in this case study: "cs" vs "phy" data structure; tree.size = 40 vs 400
cs.str = "phy"
L = 12
# parameterisation for tree construction
tree.size = 100
birth.rate = 10
sf = 2 # for plotting

# we have 8 different types of generative dynamics
# looping through a higher number repeats each with a different seed, getting more samples
for(expt in 0:35) {

  expt.type = expt %% 8
  
  # accumulation rate for features (and loss rate, for reversible setup)
  # different choices for different generative mechanisms
  # 0 -- hard single path; 1 -- hard two paths; 2 -- hard single path with loss
  # 3 -- hard two paths with loss; 4 -- uniform; 5 -- soft single path
  # 6 -- soft single pathwith loss; 7 -- soft single path FROM loss
  if(expt.type == 0) {
    set.paths = TRUE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0, L)
    start.loci = c(1,1)
  } else if(expt.type == 1) {
    set.paths = TRUE
    accumulation.rate = rep(0.5, L)
    loss.rate = rep(0, L)
    start.loci = c(1,L)
  } else if(expt.type == 2) {
    set.paths = TRUE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0.02, L)
    start.loci = c(1,1)
  } else if(expt.type == 3) {
    set.paths = TRUE
    accumulation.rate = rep(0.5, L)
    loss.rate = rep(0.02, L)
    start.loci = c(1,L)
  } else if(expt.type == 4) {
    set.paths = FALSE
    accumulation.rate = rep(1, L)
    loss.rate = rep(0, L)
    start.loci = 1:L
  } else if(expt.type == 5) {
    set.paths = FALSE
    accumulation.rate = 1-(1:L)/L
    loss.rate = rep(0, L)
    start.loci = 1:L
  } else if(expt.type == 6) {
    set.paths = FALSE
    accumulation.rate = 1-(1:L)/L
    loss.rate = rep(0.02, L)
    start.loci = 1:L
  } else if(expt.type == 7) {
    set.paths = FALSE
    accumulation.rate = rep(1, L)
    loss.rate = (1:L)/L * 0.05
    start.loci = 1:L
  }
  #c(L, 1, 2, 3, L)/L
  #accumulation.rate = c(0.05, 0.2, 2, 1)
  #loss.rate = c(0.1, 0.1, 0.1, 0.1)*0
  
  # create random phylogeny with tree.size nodes from birth-death process parameterised as above
  # my.tree = ape::rphylo(tree.size, birth=birth.rate, death=death.rate)
  my.tree = ape::rtree(tree.size, br=function(...) { return(birth.rate*runif(...)) } )
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
          # for hard paths
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
          # several error handling cases to make sure we're sensible
          poss = which(x[[this.child]] == 0)
          if(length(poss) == 1) {
            ref = poss
          } else {
            if(sum(x[[this.child]]) == L) {
              ref = 1
            } else {
              ref = sample(which(x[[this.child]] == 0), 1, prob=accumulation.rate[which(x[[this.child]] == 0)])
            }
          }
        }
        # 00001. 10000
        if(runif(1) < accumulation.rate[ref]*this.branch.length) { x[[this.child]][ref] = 1 } 
        for(ref in 1:L) {
          if(runif(1) < loss.rate[ref]*this.branch.length) { x[[this.child]][ref] = 0 } 
        }
        # in the reversible case, allow the leftmost feature ("first feature" in the ms.:
        # second paragraph of "Synthetic case studies") to revert with some probability
        
        # add this child to to state list, and to next iteration's to-do
        new.to.do = c(new.to.do, this.child)
      }
    }
    # update to-do list
    to.do = new.to.do
  }
  
  # put the results of the simulation in a consistent form
  my.tree$tip.label = x[1:length(my.tree$tip.label)]
  tip.states = unlist(lapply(my.tree$tip.label, BinToDec))+1
  
  ddf = cbind(data.frame(Species=as.character(1:tree.size)),
              as.data.frame(matrix(unlist(x[1:length(my.tree$tip.label)]), ncol=L, byrow = TRUE)))
  
  # produce a curated tree object
  my.tree2 = my.tree 
  my.tree2$tip.label = as.character(1:tree.size)
  my.ct = curate.tree(my.tree2, ddf)
  colnames(my.ct$data) = c("Species", 1:L)
  
  # do inference, either cross-sectional or phylogenetic
  if(cs.str != "cs") {
    my.hhmm.out = hyperinf(ddf, my.tree2, method="hyperhmm", nboot=0)
  } else {
    my.hhmm.out = hyperinf(ddf, method="hyperhmm", nboot=0)
  }
  
  # plot summary output
  phhmm = plot_hyperinf(my.hhmm.out)
  
  # pull ordering matrix from output of inference
  m.hhmm = pull_orders(my.hhmm.out$transitions, L, w=1)
  diag(m.hhmm) = NA
  
  # populate our output structures from this experiment
  model.fits[[length(model.fits)+1]] = list(my.hhmm.out)
  m.list[[length(m.list)+1]] = m.hhmm
  m.names = c(m.names, paste(i-1, c("hhmm")))
  
  # pull the ordering information into a plottable format
  tab.hhmm = as.table(m.hhmm)
  colnames(tab.hhmm) = 1:L
  rownames(tab.hhmm) = 1:L
  df.hhmm = as.data.frame(tab.hhmm)
  
  # produce ordering plot
  p.hhmm = ggplot(df.hhmm, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile() + scale_fill_viridis(option="magma") + 
    labs(x = "Acquired earlier", y = "Acquired\nlater", fill = "Prob")
  
  # add a megaplot to our growing list
  plot.list[[expt+1]] = 
    ggarrange(
      plotHypercube.curated.tree(my.ct, scale.fn = NULL, hjust = 1, font.size = 2.5) +
        coord_cartesian(clip = "off") +
        theme(
          plot.margin = margin(t = 0, r = 0, b = 20, l = 0)
        ), 
      ggarrange(
        phhmm, p.hhmm, ncol=1, nrow=2,
        heights=c(2,1),
        labels = c("ii", "iii")), 
      widths=c(0.6,1),
      labels=c("     i", "")
    )      
  
  if(FALSE) {
    sf = 2
    fname = paste0("tester-large-hhmm-", cs.str, "-", L, "-", tree.size, "-", expt, ".png", collapse="")
    png(fname, width=500*sf, height=350*sf, res=72*sf)
    #  print(ggarrange(plotHypercube.curated.tree(my.ct), phct, nrow=1, widths=c(1,2.5)))
    print(plot.list[[expt+1]])
    dev.off()
  }
}

# output megaplot set to file
fname = paste0("tester-all-hhmm-", cs.str, "-", L, "-", tree.size, ".png", collapse="")
png(fname, width=1000*sf, height=700*sf, res=72*sf)
print(ggarrange(plot.list[[3]],
                plot.list[[4]],
                plot.list[[5]],
                plot.list[[7]], 
                nrow=2, ncol=2, labels=c("A", "B", "C", "D")))
dev.off()

# save outputs
save(model.fits, file=paste0("model-large-fits-hhmm-", cs.str, "-", L, "-", tree.size, ".Rdata", collapse=""))

if(FALSE) {
  m.list = list()
  m.names = c()
  for(i in 1:length(model.fits)) {
    m.list[[i]] = pull_orders(model.fits[[i]][[1]]$transitions, L, w=1)
    m.names = c(m.names, paste(i-1, c("hhmm")))
    diag(m.list[[i]]) = NA
  }
}


# pull together a collection of vectors summarising ordering matrices through experiments
refs = 1:(length(model.fits))
X <- do.call(
  rbind,
  lapply(m.list[refs], function(m) as.vector(m))
)

# get rid of the diagonal elements
rownames(X) <- m.names[refs]
keep_cols <- !is.na(X[1, ])
X_clean <- X[, keep_cols]

# do PCA on our collection of orderings
pca <- prcomp(X_clean, center = TRUE, scale. = FALSE)

# store output in a dataframe with a useful group structure
df_pca <- data.frame(
  name = m.names[refs],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
df_pca$group = 1:length(model.fits)-1

# produce convex hulls for each group -- somewhat obsolete
hulls <- df_pca %>%
  group_by(group) %>%
  # filter(n() >= 3) %>%          # need ≥ 3 points for a polygon
  slice(chull(PC1, PC2))

# produce initial PCA plot
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

# summarise statistics of spread -- obsolete
spread_df <- df_pca %>%
  group_by(group) %>%
  summarise(
    n = n(),
    spread = mean(
      sqrt((PC1 - mean(PC1))^2 + (PC2 - mean(PC2))^2)
    )
  )

expt.names = c("single path", "two path", "one path\nloss", "two path\nloss", "random", "random\nordered", "random\nordered loss")
spread.plot = ggplot(spread_df, aes(x = expt.names[(group%%8)+1], y = spread, label=group)) +
  geom_boxplot() + geom_point() + geom_text_repel() +
  labs(
    x = "",
    y = "Mean distance to centroid",
    title = "Within-group dispersion in PCA space"
  ) +
  theme_minimal()

sf = 2
png(paste0("large-outputs-hhmm-", cs.str, "-", tree.size, "-", L, ".png", collapse=""), width=800*sf, height=800*sf, res=72*sf)
ggarrange(pca.plot + theme(legend.position = "none"), spread.plot, nrow=2)
dev.off()

######## second (used) PCA plot -- comparing to known truth cases

# produce ordering matrices for known cases
true.m = list()
true.m[[1]] = matrix(0, nrow=L, ncol=L)
true.m[[1]][upper.tri(true.m[[1]])] = 1
true.m[[2]] = matrix(0, nrow=L, ncol=L)
true.m[[2]][lower.tri(true.m[[2]])] = 1
true.m[[3]] = matrix(0.5, nrow=L, ncol=L)

# make a random matrix obeying the right properties for orderings
make_matrix <- function(n) {
  M <- matrix(0, n, n)
  
  # diagonal must be 1/2
  diag(M) <- 0.5
  
  # fill upper triangle randomly
  inds <- which(upper.tri(M), arr.ind = TRUE)
  vals <- runif(nrow(inds), 0, 1)
  
  M[inds] <- vals
  M[cbind(inds[,2], inds[,1])] <- 1 - vals
  
  M
}

# example
M <- make_matrix(L)

# add these truth and sample cases to our collection of ordering matrix vectors
X2 = rbind(X, as.vector(true.m[[1]]))
X2 = rbind(X2, as.vector(true.m[[2]]))
X2 = rbind(X2, as.vector(true.m[[3]]))
# random ordering matrices
for(i in 1:5) {
  X2 = rbind(X2, as.vector(make_matrix(L)))
}
# random rearrangements of the single-path solution
for(i in 1:5) {
  new.order = sample(1:ncol(M), ncol(M))
  M = true.m[[1]][new.order,new.order]
  X2 = rbind(X2, as.vector(M))
}

# remove the diagonals
rownames(X2) <- c(m.names[refs], "one-a", "one-b", "two", paste("rnd", 1:5), paste("ord", 1:5))
keep_cols <- !is.na(X2[1, ])
X2_clean <- X2[, keep_cols]

# do the PCA with truths included
pca2 <- prcomp(X2_clean, center = TRUE, scale. = FALSE)

# set up dataframe for output
df_pca2 <- data.frame(
  name = rownames(X2),
  PC1 = pca2$x[, 1],
  PC2 = pca2$x[, 2]
)

# styling choices
quick.labs = c("1", "2", "1-", "2-", "u", "1s", "1s-", "l")
df_pca2$group = c(1:length(model.fits), 100, 101, 102, 103:107, 108:112)
df_pca2$rgroup = c(quick.labs[1+(1:length(model.fits)-1) %% 8], "[1a]", "[1b]", "[2]", paste("R", 1:5, sep=""), paste("O", 1:5, sep=""))

# hulls -- obsolete
hulls2 <- df_pca2 %>%
  group_by(group) %>%
  # filter(n() >= 3) %>%          # need ≥ 3 points for a polygon
  slice(chull(PC1, PC2))

# produce plot
pca2.plot = ggplot(df_pca2, aes(PC1, PC2, colour = factor(rgroup))) +
  geom_point(size = 1) +
  geom_point(data=df_pca2[df_pca2$group>=100,], size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=rgroup), max.overlaps = 20) +
  geom_polygon(
    data = hulls2,
    aes(fill = factor(group)),
    alpha = 0.2,
    colour = NA
  ) +
  theme_minimal() +
  labs(
    colour = "Group",
    fill   = "Group"
  ) + 
  theme(legend.position="none")

# output plot
png(paste0("large-outputs-hhmm-truth-", cs.str, "-", tree.size, "-", L, ".png", collapse=""), width=800*sf, height=250*sf, res=72*sf)
pca2.plot
dev.off()

########### quantification of errors

# initialise dataframe for ground-truth fits
fit.df = model.fits[[1]][[1]]$transitions 

# dataframe for ground-truth single pathway
df1 <- fit.df
df1$Flux <- 0
state = 0
for(i in 1:L) {
  ref = which(df1$From == state & df1$To == state+2**(L-i))
  df1$Flux[ref] = 1
  state = state+2**(L-i)
}

# dataframe for ground-truth two pathway
df2 <- fit.df
df2$Flux <- 0
state = 0
for(i in 1:L) {
  ref = which(df2$From == state & df2$To == state+2**(L-i))
  df2$Flux[ref] = 0.5
  state = state+2**(L-i)
}
state = 0
for(i in 1:L) {
  ref = which(df2$From == state & df2$To == state+2**(i-1))
  df2$Flux[ref] = 0.5
  state = state+2**(i-1)
}

# function to report various quality statistics from a model fit compared to a ground truth result
fit.stats = function(fitted, truth, thresh = 0.05) {
  
  rmse = sqrt( sum((fitted$Flux - truth$Flux)**2) / nrow(truth) )
  nzs = which(fitted$Flux > thresh | truth$Flux > thresh)
  nzs = which(truth$Flux > thresh)
  
  nzfitted = fitted[nzs,]
  nztruth = truth[nzs,]
  nzrmse = sqrt( sum((nzfitted$Flux - nztruth$Flux)**2) / nrow(nztruth) )
  
  tp = length(which(truth$Flux > thresh & fitted$Flux > thresh)) / length(which(truth$Flux > thresh))
  fp = length(which(truth$Flux <= thresh & fitted$Flux > thresh)) / length(which(truth$Flux <= thresh))
  tn = length(which(truth$Flux < thresh & fitted$Flux < thresh)) / length(which(truth$Flux <= thresh))
  fn = length(which(truth$Flux > thresh & fitted$Flux <= thresh)) / length(which(truth$Flux > thresh))
  return(data.frame(rmse=rmse, nzrmse=nzrmse, tp=tp, fp=fp, tn=tn, fn=fn))
}

# apply statistics function to all fitted models
res.df = data.frame()
for(i in 1:length(model.fits)) {
  fit.df = model.fits[[i]][[1]]$transitions
  res.df = rbind(res.df, data.frame(expt=i-1, gen=1, fit.stats(fit.df, df1)))
  res.df = rbind(res.df, data.frame(expt=i-1, gen=2, fit.stats(fit.df, df2)))
}

# plot the proportion of shared edges as a summary statistic
ggplot(res.df[res.df$gen==1 | res.df$expt%%8 %in% c(1,3),], 
       aes(x=quick.labs[expt%%8 + 1], y=tp, color=factor(gen))) + 
  stat_summary(
    fun = mean,
    fun.min = \(x) mean(x) - sd(x),
    fun.max = \(x) mean(x) + sd(x),
    geom = "pointrange",
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(position = position_dodge(width = 0.6)) + 
  labs(x = "Generative model",
       y = "Proportion of edges shared\nwith ground-truth generator",
       color = "True\nmodel") +
  theme_minimal()

# other plots
if(FALSE) {
ggplot(res.df, aes(x=expt%%8, y=nzrmse, color=factor(gen))) + geom_jitter(width=0.2, height=0)

ggarrange(
  ggplot(res.df, aes(x=expt%%8, y=tp, color=factor(gen))) + geom_jitter(width=0.2),
  ggplot(res.df, aes(x=expt%%8, y=fp, color=factor(gen))) + geom_jitter(width=0.2),
  ggplot(res.df, aes(x=expt%%8, y=tn, color=factor(gen))) + geom_jitter(width=0.2),
  ggplot(res.df, aes(x=expt%%8, y=fn, color=factor(gen))) + geom_jitter(width=0.2),
  labels = c("TP","FP","TN","FN")
)
}


