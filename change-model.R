library(hyperinf)

# simply returns a binary (character string) of length len from a decimal 
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}

# simply converts a binary to a decimal
BinToDec <- function(state) {
  this.ref = 0
  for(j in 1:length(state)) {
    this.ref = this.ref + state[j]*(2**(length(state)-j))
  }
  return(this.ref)
}

if(FALSE) {
  L = 5
expt = 2
if(expt == 1) {
m = matrix(c(0,0,0,0,1,
             0,0,0,1,1,
             0,0,1,1,1,
             0,1,1,1,1), ncol=L, byrow = TRUE)
} else {
  m = matrix(c(0,0,0,0,1,
               0,0,0,1,1,
               0,0,1,1,1,
               0,1,1,1,1,
               1,0,0,0,0,
               1,1,0,0,0,
               1,1,1,0,0,
               1,1,1,1,0), ncol=L, byrow = TRUE)
}
fit = hyperinf(m)
plot_hyperinf(fit)

df = fit$transitions
}

fitted.m = model.fits[[1]][[1]]

df = fitted.m$transitions
L = fitted.m$L

d.mat = replicate(
  L,
  matrix(numeric(0), nrow = 0, ncol = L),
  simplify = FALSE
)
weights = resp = replicate(
  L,
  matrix(numeric(0), nrow = 0, ncol = 1),
  simplify = FALSE
)

for(i in 1:nrow(df)) {
  src = DecToBin(df$From[i], L)
  dest = DecToBin(df$To[i], L)
  change = which(dest != src)
  present = which(src != 0)
  coeffs = rep(0, L)
  coeffs[present] = 1
  d.mat[[change]] = rbind(d.mat[[change]], coeffs)
  resp[[change]] = rbind(resp[[change]], df$Probability[i])
  weights[[change]] = rbind(weights[[change]], df$Flux[i])
}

l2 = matrix(0, nrow=L, ncol=L)
for(i in 1:L) {
  this.d.mat = d.mat[[i]][,-i]
  this.lm = lm(log(resp[[i]]) ~ this.d.mat, weights=weights[[i]])
  l2[i,-i] = this.lm$coefficients[2:L]
  l2[i,i] = this.lm$coefficients[1]
}
l2[is.na(l2)] = 0
l2v = as.vector(t(l2))

fitht = hyperinf(m, method="hypertraps")
est = list(posterior.samples=t(as.matrix(l2v)),
           model=2,
           L=L)

library(hypertrapsct)
est.fit = PosteriorAnalysis(est)
plotHypercube.sampledgraph2(est.fit)
