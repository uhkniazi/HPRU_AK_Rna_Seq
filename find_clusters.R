# Name: find_clusters.R
# Auth: u.niazi@imperial.ac.uk
# Date: 9/3/2016
# Desc: find possible groupings in the data


## functions used in the script

## function to simulate posterior predictive for missing data
simpost = function(ivDat){
  #set.seed(123)
  # calculate using non-informative prior parameters  
  sigma.0 = 0
  k.0 = 0
  v.0 = k.0 - 1
  mu.0 = 0
  
  ## look at page 68 of Bayesian Data Analysis (Gelman) for formula
  sim.post = function(dat.grp){
    # calculate conjugate posterior
    n = length(dat.grp)
    k.n = k.0 + n
    v.n = v.0 + n
    y.bar = mean(dat.grp)
    s = sd(dat.grp)
    mu.n = (k.0/k.n * mu.0) + (n/k.n * y.bar)
    sigma.n = (( v.0*sigma.0 ) + ( (n-1)*(s^2) ) + ( (k.0*n/k.n)*((y.bar-mu.0)^2) )) / v.n
    #post.scale = ((prior.dof * prior.scale) + (var(dat.grp) * (length(dat.grp) - 1))) / post.dof
    ## simulate variance
    sigma = (sigma.n * v.n)/rchisq(1, v.n)
    mu = rnorm(1, mu.n, sqrt(sigma)/sqrt(k.n))
    return(list(mu=mu, var=sigma))
  }
  
  # get a sample of the posterior variance for each group
  s = sim.post(ivDat)
  return(rnorm(1, s$mu, sqrt(s$var)))
}


plot.scale = function(x, ...){
  # plot the z-scaled data vs original values
  x = na.omit(x)
  z = scale(x)
  plot(z, x, ...)
}


##### read the data
dfData = read.csv('Data_external/MSD_data.csv', header=T, row.names=1)

dfData = na.omit(dfData)

par(mfrow=c(2,2))

sapply(seq_along(1:nrow(dfData)), function(x){
  d = as.numeric(dfData[x,])
  names(d) = colnames(dfData)
  plot.scale(d, main=rownames(dfData)[x], xlab='Z scaled', ylab='Raw')
})

## any clusters in the data?
par(mfrow=c(1,1), mar=c(2,2,1,1))
m1 = t(as.matrix(dfData))
# scale the data across columns i.e each mediator
mCounts = scale(m1)

pr.out = prcomp(mCounts, scale=F)
biplot(pr.out)

d = dist(mCounts)
hc = hclust(d)
plot(hc)



