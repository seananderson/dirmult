set.seed(10293)
p <- runif(50)
p <- p / sum(p)
# paa as calculated from this p in dirichlet.R was manually pasted into
# admb/dirmult.dat and is called paa_obs in the admb/dirmult.tpl code

# https://github.com/pbs-assess/gfiscamutils/blob/de8336de4e04c0d3033db904c5963b19fbceee84/R/load-models.R#L685
rep <- gfiscamutils::read.report.file("admb/dirmult.rep")
p_hat <- rep$p
paa_pred_admb <- exp(p_hat) / sum(exp(p_hat))

phi <- exp(rep$log_phi)
plot(paa_pred_admb, p);abline(0, 1)

# Add TMB estimates - source("dirichlet.R") first to make opt object
p_hat <- as.numeric(opt$par[-1])
paa_pred <- exp(p_hat) / sum(exp(p_hat))
points(paa_pred, p, col = "red")
