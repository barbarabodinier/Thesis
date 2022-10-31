rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter3/")

n <- 100000
mean1 <- 1
mean2 <- 3

rGaussianMixture <- function(n, mu_set, sd_set, pz) {
  # Step 1: drawing the latent variable Z
  z <- rmultinom(n = n, size = 1, prob = pz)
  z <- sharp:::DummyToCategories(t(z))

  # Step 2: drawing the variable X conditionally on the latent variable Z
  x1 <- rnorm(n, mean = mu_set[1], sd = sd_set[1])
  x2 <- rnorm(n, mean = mu_set[2], sd = sd_set[2])
  x <- ifelse(z == 1, yes = x1, no = x2)

  return(x)
}

pdf(paste0("Working_figures/Gaussian_mixture_examples.pdf"),
  width = 12, height = 4
)
par(mar = c(5, 5, 1, 1), mfrow = c(1, 3))
set.seed(1)

# Example 1
x <- rGaussianMixture(
  n = n,
  mu_set = c(mean1, mean2),
  sd_set = c(0.1, 0.1),
  pz = c(0.5, 0.5)
)
plot(density(x),
  las = 1, lwd = 2, col = "navy",
  xlab = "x", cex.lab = 1.5, main = ""
)
abline(v = mean1, col = "forestgreen", lty = 2, lwd = 2)
abline(v = mean2, col = "tomato", lty = 2, lwd = 2)

# Example 2
x <- rGaussianMixture(
  n = n,
  mu_set = c(mean1, mean2),
  sd_set = c(0.5, 0.5),
  pz = c(0.5, 0.5)
)
plot(density(x),
  las = 1, lwd = 2, col = "navy",
  xlab = "x", ylab = "", cex.lab = 1.5, main = ""
)
abline(v = mean1, col = "forestgreen", lty = 2, lwd = 2)
abline(v = mean2, col = "tomato", lty = 2, lwd = 2)

# Example 3
x <- rGaussianMixture(
  n = n,
  mu_set = c(mean1, mean2),
  sd_set = c(0.5, 0.5),
  pz = c(0.2, 0.8)
)
plot(density(x),
  las = 1, lwd = 2, col = "navy",
  xlab = "x", ylab = "", cex.lab = 1.5, main = ""
)
abline(v = mean1, col = "forestgreen", lty = 2, lwd = 2)
abline(v = mean2, col = "tomato", lty = 2, lwd = 2)
dev.off()
