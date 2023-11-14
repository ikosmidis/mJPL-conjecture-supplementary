library("parallel")
library("brglm2")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ## Number of cores to use for determining whether MLE exists across settings
    ncores <- 10
    ## Number of observations per dataset
    nobs <- 3000
    ## beta star setting
    beta_star_setting <- "s1"
    ## seed
    seed <- 0
}
maxit <- 250
tolerance <- 1e-03
ns <- 200000

image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

out_path <- file.path(image_path, paste0("estimates-n-", nobs,
                                         "-beta-", beta_star_setting,
                                         "-test",
                                         ".rda"))

## kappa-gamma combinations to consider in the experiment
kappa_gamma <- data.frame(kappa = c(0.01, 0.01, 0.01,
                                    0.05, 0.05, 0.05,
                                    0.15, 0.15, 0.15,
                                    0.22, 0.22,
                                    0.25, 0.25, 0.25,
                                    0.30, 0.30,
                                    0.35,
                                    0.35, 0.35, 0.35,
                                    0.40, 0.40,
                                    0.45, 0.45, 0.45,
                                    0.50, 0.50,
                                    0.55, 0.55, 0.55
                                    ),
                          gamma = c(1, 8, 15,
                                    4.5, 11.5, 18.5,
                                    1, 8, 15,
                                    8, 15,
                                    4.5, 11.5, 18.5,
                                    8, 15,
                                    1,
                                    4.5, 11.5, 18.5,
                                    8, 15,
                                    4.5, 11.5, 18.5,
                                    8, 15,
                                    4.5, 11.5, 18.5
                                    ))

rhosq_psi <- expand.grid(rhosq = seq.int(0, 9, 1) / 10,
                         psi = c(0, 0.3, 0.6, 0.9))

e_settings <- NULL
for (s in 1:nrow(rhosq_psi)) {
    e_settings <- rbind(e_settings, cbind(psi = rhosq_psi$psi[s], rhosq = rhosq_psi$rhosq[s], kappa_gamma))
}
n_settings <- nrow(e_settings)


## Add distinct seeds to experiment settings
set.seed(seed)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
e_settings$mle_exists <- mclapply(seq.int(nrow(e_settings)), function(j) {
    rhosq <- e_settings$rhosq[j]
    gamma <- e_settings$gamma[j]
    kappa <- e_settings$kappa[j]
    ## MLE exists if below the phase transition curve or if the kappa
    ## we get for each gamma in the e_settings is larger than the kappa in
    ## the e_settings
    compute_pt(beta0 = sqrt(rhosq) * gamma, gamma, ncores = 1, XZU =  xzu)$kappa > kappa
}, mc.cores = ncores) |> unlist()


dup_check <- TRUE
while (dup_check) {
    e_settings$seed <- round(runif(n_settings) * 1000000)
    dup_check <- length(unique(e_settings$seed)) != n_settings
}

## Get estimates
results <- mclapply(1:n_settings, function(s) {
    setting <- e_settings[s, ]
    compute_estimates(setting,
                      beta_star_setting = beta_star_setting,
                      n = nobs,
                      maxit = maxit,
                      tolerance = tolerance,
                      verbose = 0,
                      current_setting = s,
                      max_setting = n_settings)
}, mc.cores = ncores)

estimates <- do.call("rbind", results)
performance <- do.call("rbind", lapply(results, function(x) attr(x, "performance")))

## Save results
save(estimates, performance, file = out_path)
