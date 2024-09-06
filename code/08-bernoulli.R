library("parallel")
library("brglm2")
library("detectseparation")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ## Number of cores to use for determining whether MLE exists across settings
    ncores <- 10
    ## Number of observations per dataset
    nobs <- 2000
    ## beta star setting
    beta_star_setting <- "s1"
    ## seed
    seed <- 0
    rsqs <- 0.3
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
}

maxit <- 250
tolerance <- 1e-03

image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

## ## kappa-gamma combinations to consider in the experiment
## kappa_gamma <- data.frame(kappa = c(0.01, 0.01, 0.01,
##                                     0.05, 0.05, 0.05,
##                                     0.15, 0.15, 0.15,
##                                     0.22, 0.22,
##                                     0.25, 0.25, 0.25,
##                                     0.30, 0.30,
##                                     0.35,
##                                     0.35, 0.35, 0.35,
##                                     0.40, 0.40,
##                                     0.45, 0.45, 0.45,
##                                     0.50, 0.50,
##                                     0.55, 0.55, 0.55
##                                     ),
##                           gamma = c(1, 8, 15,
##                                     4.5, 11.5, 18.5,
##                                     1, 8, 15,
##                                     8, 15,
##                                     4.5, 11.5, 18.5,
##                                     8, 15,
##                                     1,
##                                     4.5, 11.5, 18.5,
##                                     8, 15,
##                                     4.5, 11.5, 18.5,
##                                     8, 15,
##                                     4.5, 11.5, 18.5
##                                     ))

kappa_gamma <- expand.grid(gamma = c(2, 4, 8, 12, 16, 20),
                           kappa = c(0.05, 0.1, 0.3, 0.5, 0.7))
rhosq_prob <- expand.grid(rhosq = rsqs,
                          prob = c(0.1))

e_settings <- NULL
for (s in 1:nrow(rhosq_prob)) {
    e_settings <- rbind(e_settings, cbind(prob = rhosq_prob$prob[s], rhosq = rhosq_prob$rhosq[s], kappa_gamma))
}
n_settings <- nrow(e_settings)

set.seed(seed)
dup_check <- TRUE
while (dup_check) {
    e_settings$seed <- round(runif(n_settings) * 1000000)
    dup_check <- length(unique(e_settings$seed)) != n_settings
}

## Get estimates
results <- mclapply(1:n_settings, function(s) {
    setting <- e_settings[s, ]
    compute_estimates_bern(setting,
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

## ## Save results
save(estimates, performance, file = file.path(image_path, "bernoulli.rda"))

load(file.path(image_path, conjecture_model))
## Get conjecture constants
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]


library("ggplot2")
library("patchwork")
library("dplyr")
## "#E69F00", , "#0072B2", , "#CC79A7", "#56B4E9"

okabe <- c("#D55E00", "#009E73", "#F0E442")
mjpls <- rmjpls <- as.list(rep(0, length(rsqs)))
for (j in seq_along(rsqs)) {
    crhosq <- rsqs[j]
    ss <- estimates |> filter(method == "mJPL", rhosq == crhosq, parameter > 0)
    ss1 <- ss |> correct_mJPL_estimates(pk = po_k, pg = po_g, pr = po_r)
    ss1$method <- "rmJPL"
    ss <- rbind(ss, ss1)
    ss <- ss |> mutate(kappa_fac = paste("kappa == ", kappa),
                       gamma_fac = paste("gamma == ", gamma))
    ss$kappa_fac <- factor(ss$kappa_fac, levels = unique(ss$kappa_fac), ordered = TRUE)
    ss$gamma_fac <- factor(ss$gamma_fac, levels = rev(unique(ss$gamma_fac)), ordered = TRUE)

    mjpls[[j]] <-  ggplot(ss |> filter(method == "mJPL", parameter > 0)) +
        geom_point(aes(x = truth, y = estimate, col = mle_exists), size = 0.3, alpha = 0.5) +
        geom_abline(aes(intercept = 0, slope = 1)) +
        geom_smooth(aes(x = truth, y = estimate), col = okabe[3], method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
        facet_grid(gamma_fac ~ kappa_fac, labeller = label_parsed) +
        coord_cartesian(x = c(-10, 10), y = c(-10, 10)) +
        theme_minimal() +
        labs(x = NULL) +
        theme(legend.position = "top") +
        scale_color_manual(values = okabe, name = "ML estimate exists")
        ## labs(title = substitute(paste("mJPL, ", rho^2 == r), list(r = crhosq)))

    rmjpls[[j]] <- ggplot(ss |> filter(method == "rmJPL", parameter > 0)) +
        geom_point(aes(x = truth, y = estimate, col = mle_exists), size = 0.3, alpha = 0.5) +
        geom_abline(aes(intercept = 0, slope = 1)) +
        geom_smooth(aes(x = truth, y = estimate), col = okabe[3], method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
        facet_grid(gamma_fac ~ kappa_fac, labeller = label_parsed) +
        coord_cartesian(x = c(-10, 10), y = c(-10, 10)) +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_color_manual(values = okabe)
        ## labs(title = substitute(paste("rescaled mJPL, ", rho^2 == r), list(r = crhosq)))
}


pdf(file.path(project_path, "figures/bernoulli.pdf"), width = 6.5, height = 9)
## pdf(file.path("~/Repositories/mJPL-conjecture/figures", "bernoulli.pdf"), width = 6.5, height = 9)
print(mjpls[[1]] / rmjpls[[1]])
dev.off()
