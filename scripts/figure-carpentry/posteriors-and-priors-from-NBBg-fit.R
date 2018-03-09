library(tidyverse)

# These samples are already burned in
# Originally 3 separate chains of 50,000 samples each and 10,000 sample burn in
# Now all three chains are included together, so ((50000 - 10000) * 3 chains) = 120,000 samples
post_samps <- read_csv("data/NBBg-samples/NBBg-samples-combined.csv")

# All priors are gamma distributions with shape and scale parameters as given
priors <- list(R0 = c(shape = 2.6, scale = 1),
               alpha = c(shape = 0.0037, scale = 1),
               kE = c(shape = 17.6, scale = 1),
               kD = c(shape = 1.07, scale = 1))

glimpse(posterior_samps)

pdf("figures/posteriors-and-priors-from-NBBg-fit.pdf", height = 4.5, width = 6)
# postscript("figures/posteriors-and-priors-from-NBBg-fit.eps", height = 5, width = 5)
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = c(2, 2, 0, 0), family = "Helvetica")
# Plot R0, the density independent per capita growth rate in an average environment
plot(density(post_samps$R0), main = expression(R[0]), xlab = NA, ylab = NA, las = 1)
curve(dgamma(x, shape = priors$R0["shape"], scale = priors$R0["scale"]), add = TRUE, col = "red")

# Plot alpha, the cannibalism rate
plot(density(post_samps$alpha), main = expression(alpha), xlab = NA, ylab = NA, las = 1)
curve(dgamma(x, shape = priors$alpha["shape"], scale = priors$alpha["scale"]), add = TRUE, col = "red")

# Plot kE
plot(density(post_samps$kE), main = expression(k[E]), xlab = NA, ylab = NA, las = 1)
curve(dgamma(x, shape = priors$kE["shape"], scale = priors$kE["scale"]), add = TRUE, col = "red")

# Plot kD
plot(density(post_samps$kD), main = expression(k[D]), xlab = NA, ylab = NA, las = 1)
curve(dgamma(x, shape = priors$kD["shape"], scale = priors$kD["scale"]), add = TRUE, col = "red")

legend("topright", legend = c("Posterior", "Prior"), col = c("black", "red"), lty = 1, bty = "n")

mtext(side = 1, text = "Parameter value", outer = TRUE)
mtext(side = 2, text = "Probability density", outer = TRUE)

dev.off()