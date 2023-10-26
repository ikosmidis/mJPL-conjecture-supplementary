library("parallel")
library("brglm2")

if (interactive()) {
    project_path <- "~/Desktop/mJPL-conjecture"
    ## Number of cores to use for the computation of phase transition curves
    ncores <- 10
    ## Number of observations per dataset
    nobs <- 2000
    ## psi value
    psi <- 0.5
    ## beta star setting
    beta_star_setting <- "u2"
    ## Number of independent repetitions per rhosq-kappa-gamma setting
    repetitions <- 50
    ## seed for generating unique seeds
    seed <- 0

}
maxit <- 250
tolerance <- 1e-03

image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

out_path <- file.path(image_path, paste0("estimates-n-", nobs,
                                         "-beta-", beta_star_setting,
                                         "-psi-", psi,
                                         "-training",
                                         ".rda"))

## Get settings
load(file.path(image_path, paste0("design-training.rda")))

design <- cbind(psi, design)
e_settings <- design[rep(1:nrow(design), each = repetitions), ]
n_settings <- nrow(e_settings)

## Add distinct seeds to experiment settings
set.seed(seed)
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
