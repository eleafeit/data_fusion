# Data fusion example in Stan
# Elea McDonnell Feit, eleafeit@gmail.com
# 11 March 2016

library(MASS)
library(coda)
library(beanplot)
library(rstan)
source("Data_Fusion_Functions.R")

# Example 1a: MVN ================================================================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(c(1, 0.3, -0.2, 0.7, 0.3, 1, -0.6, 0.4, 
                  -0.2, -0.6, 1, 0.1, 0.7, 0.4, 0.1, 1), nrow=4)
d1 <- data.mvn.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
str(d1$data)
# Call to Stan to generate posterior draws
m1 <- stan(file="Data_Fusion_MVN.stan", data=d1$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
# Summaries of posterior draws for population-level parameters
summary(m1, par=c("mu"))  
summary(m1, par=c("tau"))
summary(m1, par=c("Omega"))
plot.post.density(m1, pars=c("mu", "tau"), prefix="Ex1",
                  true=list(d1$true$mu, sqrt(diag(d1$true$Sigma)), 
                            cov2cor(d1$true$Sigma)))
draws <- As.mcmc.list(m1, pars=c("Omega"))
png(filename="./images/Ex1PostOmega.png", width=600, height=600)
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega (correlations)", log=""), 
         cex.axis=0.5)
dev.off()
# Summaries of posterior draws for missing data
summary(extract(m1, par=c("y1mis"))$y1mis[,3,]) 
png("./images/Ex1y13mis.png")
plot(density(extract(m1, par=c("y1mis"))$y1mis[,3,]), 
     main="Posterior of Unobserved y_1", xlab="y_1")
dev.off()
summary(m1, par=c("y"))  # posteriors of observed data place a point mass at the observed value
plot.true.v.est(m1, pars=c("y1mis", "y2mis"), prefix="Ex1",
                true=list(d1$true$y1mis, d1$true$y2mis))

# Example 1b: MVN with zero correlations =======================================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(0, nrow=4, ncol=4)
diag(Sigma) <- 1
# Call to Stan to generate posterior draws
d2 <- data.mvn.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m2 <- stan(file="Data_Fusion_MVN.stan", data=d2$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
# Summarize posteriors of population-level parameters
summary(m2, par=c("mu"))  
summary(m2, par=c("tau"))
summary(m2, par=c("Omega"))
plot.post.density(m2, pars=c("mu", "tau"), prefix="Ex2",
                  true=list(d1$true$mu, sqrt(diag(d1$true$Sigma)), 
                            cov2cor(d1$true$Sigma)))
draws <- As.mcmc.list(m2, pars=c("Omega"))
png(filename="./images/Ex2PostOmega.png", width=600, height=400)
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""), cex.axis=0.5)
dev.off()
# Summaries of posterior draws for missing data
plot.true.v.est(m2, pars=c("y1mis", "y2mis"), prefix="Ex2",
                true=list(d2$true$y1mis, d2$true$y2mis))

# Example 1c: MVN with strong positive correlations ===========================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(0.9, nrow=4, ncol=4)
diag(Sigma) <- 1
# Call to Stan to generate posterior draws
d3 <- data.mvn.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m3 <- stan(file="Data_Fusion_MVN.stan", data=d3$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
# Summaries of population-level parameters
summary(m3, par=c("mu"))  
summary(m3, par=c("tau"))
summary(m3, par=c("Omega"))
plot.post.density(m3, pars=c("mu", "tau"), prefix="Ex3",
                  true=list(d1$true$mu, sqrt(diag(d1$true$Sigma))))
draws <- As.mcmc.list(m3, pars=c("Omega"))
png(filename="./images/Ex3PostOmega.png", width=600, height=400)
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))
dev.off()
# Summaries of posterior draws for missing data
plot.true.v.est(m3, pars=c("y1mis", "y2mis"), prefix="Ex3",
                true=list(d3$true$y1mis, d3$true$y2mis))

# Example 2: MVP =================================================================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(c(1, 0.3, -0.2, 0.7, 0.3, 1, -0.6, 0.4, 
                  -0.2, -0.6, 1, 0.1, 0.7, 0.4, 0.1, 1), nrow=4)
d1 <- data.mvp.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
# Call to Stan to generate posterior draws
m1 <- stan(file="Data_Fusion_MVP.stan", data=d1$data, 
           iter=10000, warmup=2000, chains=1, seed=35)
# Summaries of posteriors of population-level parameters
summary(m1, par=c("mu", "Omega"))
plot.post.density(m1, pars=c("mu"), prefix="Ex1MVP", true=list(d1$true$mu))
png(filename="./images/Ex1MVPPostOmega.png", width=600, height=400)
draws <- As.mcmc.list(m1, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))
dev.off()
# Summarize posteriors for one of missing values 
y1mis.draws <- extract(m1, par=c("y1mis"))[[1]][,1,1]  # draws for third respondent
mean(y1mis.draws > 0)
# Confusion matrix for missing data
y1mis.est <- summary(m1, par=c("y1mis"))$summary[, "50%"]>0
xtabs(~y1mis.est + (d1$true$y1mis>0))
y2mis.est <- summary(m1, par=c("y1mis"))$summary[, "50%"]>0
xtabs(~y2mis.est + (d1$true$y2mis>0))
z.est <- data.frame(z.true=as.vector(t(d1$true$z)), y=as.vector(t(d1$true$y)), 
                    z.postmed=summary(m1, pars=c("z"))$summary[,"50%"])
png(filename="./images/Ex1MVPTrueVEstz.png", width=600, height=400)
plot(z.est[,c(1,3)], xlab="True Latent Variable", ylab="Posterior Mean of Latent Variable")
points(z.est[is.na(z.est$y), c(1,3)], col="red")
abline(h=0, v=0)
dev.off()

Sigma <- matrix(0, nrow=4, ncol=4)
diag(Sigma) <- 1
d2 <- data.mvp.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m2 <- stan(file="Data_Fusion_MVP.stan", data=d2$data, 
           iter=10000, warmup=2000, chains=1, seed=35)
print(m2, par=c("mu", "Omega"))
plot.post.density(m2, pars=c("mu"), prefix="Ex2MVP", true=list(d2$true$mu))
png(filename="./images/Ex2MVPPostOmega.png", width=600, height=400)
draws <- As.mcmc.list(m2, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))
dev.off()
y1mis.est <- summary(m2, par=c("y1mis"))$summary
xtabs(~y1mis.est[,"50%"]+d2$true$y1mis)
y2mis.est <- summary(m2, par=c("y1mis"))$summary
xtabs(~y2mis.est[,"50%"]+d2$true$y2mis)
png(filename="./images/Ex2MVPTrueVEstz.png", width=600, height=400)
z.est <- data.frame(z.true=as.vector(t(d2$true$z)), y=as.vector(t(d2$true$y)), 
                    z.postmed=summary(m2, pars=c("z"))$summary[,"50%"])
plot(z.est[,c(1,3)], xlab="True Latent Variable", ylab="Posterior Mean of Latent Variable")
points(z.est[is.na(z.est$y), c(1,3)], col="red")
abline(h=0, v=0)
dev.off()

Sigma <- matrix(0.9, nrow=4, ncol=4)
diag(Sigma) <- 1
d3 <- data.mvp.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m3 <- stan(file="Data_Fusion_MVP.stan", data=d3$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
print(m3, par=c("mu", "Omega"))
plot.post.density(m3, pars=c("mu"), prefix="Ex3MVP", true=list(d2$true$mu))
png(filename="./images/Ex3MVPPostOmega.png", width=600, height=400)
draws <- As.mcmc.list(m3, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))
dev.off()
y1mis.est <- summary(m3, par=c("y1mis"))$summary[, "50%"]
xtabs(~y1mis.est + (d3$true$y1mis>0))
y2mis.est <- summary(m3, par=c("y1mis"))$summary[,"50%"]
xtabs(~y2mis.est + (d3$true$y2mis>0))
png(filename="./images/Ex3MVPTrueVEstz.png", width=600, height=400)
z.est <- data.frame(z.true=as.vector(t(d3$true$z)), y=as.vector(t(d3$true$y)), 
                    z.postmed=summary(m3, pars=c("z"))$summary[,"50%"])
plot(z.est[,c(1,3)], xlab="True Latent Variable", ylab="Posterior Mean of Latent Variable")
points(z.est[is.na(z.est$y), c(1,3)], col="red")
abline(h=0, v=0)
dev.off()

Sigma <- matrix(0.5, nrow=7, ncol=7)
diag(Sigma) <- 1
d4 <- data.mvp.split(K1=2, K2=2, Kb=3, N1=100, N2=100, mu=rep(0,7), Sigma=Sigma)
m4 <- stan(file="Data_Fusion_MVP.stan", data=d4$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
print(m4, pars=c("mu", "Omega"))
