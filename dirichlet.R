rdirichlet <- function(p, phi) {
  alpha <- p * phi
  obs <- rgamma(length(p), alpha, 1)
  obs / sum(obs)
}

rdirichlet_linear <- function(p, theta, N) {
  alpha <- p * theta * N
  obs <- rgamma(length(p), alpha, 1)
  obs / sum(obs)
}

out1 <- sapply(1:10, function(i) rdirichlet(c(0.1, 0.1, 0.4, 0.4), 1e6))
rowMeans(out1)

# lower phi = more noise on probabilities p:
out2 <- sapply(1:10, function(i) rdirichlet(c(0.1, 0.1, 0.4, 0.4), 1))
rowMeans(out2)

rdirmultinom <- function(N, p, phi, type = "saturating") {
  if (type[[1]] != "saturating") {
    dp <- rdirichlet_linear(p, phi, N)
  } else {
    dp <- rdirichlet(p, phi) # note bug in WHAM simulate(); should be done once per multinomial
  }
  obs <- rmultinom(n = 1, prob = dp, size = N)
  obs
}

p <- c(0.1, 0.1, 0.4, 0.4)
p <- p / sum(p)

sim_rdirmult <- function(N, phi) {
  out <- matrix(ncol = length(p), nrow = 100, data = 0)
  for (i in 1:nrow(out)) out[i, ] <- rdirmultinom(N, p, phi)
  plot(1, ylim = c(1, ncol(out)), xlim = c(1, nrow(out)), type = "n")
  out <- out / N
  for (i in seq(1, nrow(out))) {
    symbols(rep(i, ncol(out)), seq_along(out[1, ]),
      circles = out[i, ] * 4, inches = FALSE, add = TRUE
    )
  }
}
sim_rdirmult(100, 0.1)
sim_rdirmult(100, 100)

get_neff <- function(N, phi) {
  (N + N * phi) / (N + phi)
}
ns <- seq(1, 100, length.out = 100)
plot(ns, get_neff(ns, 1000), type = "l", col = "red", ylab = "Neff", xlab = "N")
lines(ns, get_neff(ns, 100), type = "l", col = "blue")
lines(ns, get_neff(ns, 10), type = "l", col = "green")

# so these should be about the same:
sim_rdirmult(300, 40)
sim_rdirmult(get_neff(300, 40), 1e9)

sim_rdirmult(300, 5)
sim_rdirmult(get_neff(300, 5), 1e9)

sim_rdirmult(3000, 50)
get_neff(3000, 50)
sim_rdirmult(get_neff(3000, 50), 1e9)

get_neff_linear <- function(N, theta) {
  N * (theta / (1 + theta))
}
get_neff_linear(300, 1)

# linear
ns <- seq(1, 100, length.out = 100)
plot(ns, get_neff_linear(ns, 100), type = "l", col = "red", ylab = "Neff", xlab = "N")
lines(ns, get_neff_linear(ns, 1), type = "l", col = "blue")
lines(ns, get_neff_linear(ns, 0.1), type = "l", col = "green")

library(TMB)
compile("dirichlet.cpp")
dyn.load(dynlib("dirichlet"))

set.seed(10293)
p <- runif(50)
p <- p / sum(p)
phi <- 200
N <- 200
get_neff(N, phi)

naa <- matrix(ncol = length(p), nrow = 40, data = 0)
for (i in 1:nrow(naa)) naa[i, ] <- rdirmultinom(N, p, phi = phi)
plot(1, ylim = c(1, ncol(naa)), xlim = c(1, nrow(naa)), type = "n")
paa <- naa / N
for (i in seq(1, nrow(naa))) {
  symbols(rep(i, ncol(naa)), seq_along(naa[1, ]),
    circles = paa[i, ] * 4, inches = FALSE, add = TRUE
  )
}

data <- list(paa_obs = paa, Neff = rep(N, nrow(paa)))
parameters <- list(log_phi = rnorm(1, log(phi), 0.5), p = rep(1, ncol(paa)))

obj <- MakeADFun(data, parameters, DLL="dirichlet")
opt <- nlminb(obj$par, obj$fn, obj$gr)
# sdr <- sdreport(obj)

p_hat <- as.numeric(opt$par[-1])
paa_pred <- exp(p_hat) / sum(exp(p_hat))

phi
exp(as.numeric(opt$par[1]))

plot(paa_pred, p);abline(0, 1)

paa_pred_matrix <- matrix(paa_pred, ncol = ncol(paa), nrow = nrow(paa), byrow = TRUE)

plot(1, ylim = c(1, ncol(paa)), xlim = c(1, nrow(paa)),
  type = "n", xlab = "Year", ylab = "Age bin")
for (i in seq(1, nrow(paa))) {
  symbols(rep(i, ncol(paa)), seq_along(paa[1, ]),
    circles = paa[i, ] * 4, inches = FALSE, add = TRUE
  )
  symbols(rep(i, ncol(paa)), seq_along(paa[1, ]),
    circles = paa_pred_matrix[i, ] * 4, inches = FALSE, add = TRUE, bg = "red"
  )
}
