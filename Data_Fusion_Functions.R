# Functions for data fusion
# Elea McDonnell Feit, eleafeit@gmail.com
# 11 March 2016

data.mvn.split <- function(K1=2, K2=2, Kb=3, N1=100, N2=100,
                           mu=rep(0, K1+K2+Kb), Sigma=diag(1, K1+K2+Kb)) 
{
  y <- mvrnorm(n=N1+N2, mu=mu, Sigma=Sigma)
  list(data=list(K1=K1, K2=K2, Kb=Kb, N1=N1, N2=N2, 
                 y1=as.matrix(y[1:N1, 1:K1], col=K1), 
                 y2=as.matrix(y[N1+1:N2, K1+1:K2], col=K2),
                 yb=as.matrix(y[,K1+K2+1:Kb], col=Kb)), 
       true=list(mu=mu, Sigma=Sigma, 
                 y1mis=y[1:N1, K1+1:K2], 
                 y2mis=y[N1+1:N2, 1:K1]))
}   

data.mvp.split <- function(K1=2, K2=2, Kb=3, N1=100, N2=100,
                           mu=rep(0, K1+K2+Kb), Sigma=diag(1, K1+K2+Kb)) 
{
  z <- mvrnorm(n=N1+N2, mu=mu, Sigma=Sigma)
  y <- z
  y[y>0] <- 1
  y[y<0] <- 0
  y1mis <- y[1:N1, K1+1:K2]
  y2mis <- y[N1+1:N2, 1:K1]
  y[1:N1, K1+1:K2] <- NA
  y[N1+1:N2, 1:K1] <- NA
  true=list(mu=mu, Sigma=Sigma, z=z, y=y, y1mis=y1mis, y2mis=y2mis)
  y[is.na(y)] <- 0
  data=list(K1=K1, K2=K2, Kb=Kb, N1=N1, N2=N2, y=y) 
  list(data=data, true=true)
} 

plot.post.density <- function(m.stan, pars, true, prefix=NULL){
  for (i in 1:length(pars)) {
    draws <- As.mcmc.list(m.stan, pars=pars[i])
    if (!is.null(prefix)) {
      filename <- paste(prefix, "Post",  pars[i], ".png", sep="")
      png(filename=filename, width=600, height=400)
    }
    beanplot(data.frame(draws[[1]]), 
             horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
             main=paste("Posterior Density of", pars[[i]]))
    if (!is.null(prefix)) dev.off()
  }
}

plot.true.v.est <- function(m.stan, pars, true, prefix=NULL){
  for (i in 1:length(pars)) {
    draws <- As.mcmc.list(m.stan, pars=pars[i]) 
    est <- summary(draws)
    if (!is.null(prefix)) {
      filename <- paste(prefix, "TrueVEst", pars[i], ".png", sep="")
      png(filename=filename, width=600, height=400)
    }
    plot(true[[i]], est$quantiles[,3], col="blue", 
         xlab=paste("True", pars[i]), 
         ylab=paste("Estiamted", pars[i], "(posterior median)"))
    abline(a=0, b=1)
    arrows(true[[i]], est$quantiles[,3], true[[i]], est$quantiles[,1], 
           col="gray90", length=0)
    arrows(true[[i]], est$quantiles[,3], true[[i]], est$quantiles[,5], 
           col="gray90", length=0)
    points(true[[i]], est$quantiles[,3], col="blue")
    if (!is.null(prefix)) dev.off()
  }
}