library(tidyverse)

pops <- read_csv("data/initial-density-dependence.csv")
glimpse(pops)

fit <- read_csv("data/NBBg-samples/NBBg-samples-combined.csv")
fit_lwr <- apply(fit, 2, quantile, prob = 0.025)
fit_mean <- colMeans(fit)
fit_upr <- apply(fit, 2, quantile, prob = 0.975)


N0 <- seq(0, 200)
model_pops <- data.frame(Nt = N0)

for (i in seq_along(N0)) {
  Ntp1 <- N0[i] * fit$R0 * exp(-1 * fit$alpha * N0[i])
  model_pops$Ntp1[i] <- list(Ntp1)
  model_pops$lwr[i] <- quantile(Ntp1, probs = 0.025)
  model_pops$mean[i] <- mean(Ntp1)
  model_pops$upr[i] <- quantile(Ntp1, probs = 0.975)
}

glimpse(model_pops)

K <- log(fit$R0) / fit$alpha

K_lwr <- quantile(K, 0.025)
K_mean <- log(fit_mean["R0"]) / fit_mean["alpha"]
K_upr <- quantile(K, 0.975)

pdf("figures/density-dependence.pdf", height = 3.14961, width = 3.14961)
# postscript("figures/density-dependence.eps", height = 3.14961, width = 3.14961)
par(mar = c(3.5, 3.5, 1, 1), family = "Helvetica", mgp = c(2.25, 1, 0))

plot(x = pops$N0, y = pops$census, xlab = expression(N[t]), ylab = expression(N[t+1]), las = 1, type = "n")

# shaded 95% credible interval around the mean
polygon(x = c(model_pops$Nt, rev(model_pops$Nt)), y = c(model_pops$lwr, rev(model_pops$upr)), col = adjustcolor("lightgrey", alpha.f = 0.5), border = NA)

# data points from the 125 populations
points(x = pops$N0, y = pops$census, pch = 19, cex = 0.5)

# deterministic Ricker function using mean values of R0 and alpha
lines(x = model_pops$Nt, y = model_pops$mean, col = "red", lwd = 1)


# lines(x = model_pops$Nt, y = model_pops$lwr, col = "blue", lwd = 2)
# lines(x = model_pops$Nt, y = model_pops$upr, col = "blue", lwd = 2)
# abline(v = K_mean, col = "blue", lwd = 2)
# abline(v = K_lwr)
# abline(v = K_upr)
dev.off()