library("parallel")
library("brglm2")
library("dplyr")


if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    image_path <- file.path(project_path, "images")
    code_path <- file.path(project_path, "code")
    figure_path <- file.path(project_path, "figures")
    nobs <- 1000
    repetitions <- 20
    ncores <- 10
    ## seed for generating unique seeds
    seed <- 201
}
maxit <- 300
tolerance <- 1e-03

source(file.path(code_path, "helper-functions.R"))

kappas <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
gammas <- c(1, 5, 10, 15.0)

settings0 <- expand.grid(kappa = kappas, gamma = gammas)

set.seed(123)
ns <- 200000
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
settings0$kappa_diff <- mclapply(seq.int(nrow(settings0)), function(j) {
    gamma <- settings0$gamma[j]
    kappa <- settings0$kappa[j]
    ## MLE exists if below the phase transition curve or if the kappa
    ## we get for each gamma in the design is larger than the kappa in
    ## the design
    compute_pt(beta0 = 0, gamma, ncores = 1, XZU =  xzu)$kappa - kappa
}, mc.cores = ncores) |> unlist()

settings0$mle_exists <- settings0$kappa_diff > 0


settings <- settings0[rep(1:nrow(settings0), each = repetitions), ]
n_settings <- nrow(settings)
set.seed(seed)
dup_check <- TRUE
while (dup_check) {
    settings$seed <- round(runif(n_settings) * 100000000)
    dup_check <- length(unique(settings$seed)) != n_settings
}

settings <- settings |> transform(p = nobs * kappa,
                                  psi = 0,
                                  rhosq = 0)

## base_beta <- c(-10, 0, 0, 0)
results <- mclapply(1:n_settings, function(s) {
    setting <- settings[s, ]
    beta_s <- rep(c(1, 2), each = settings[s, "p"] / 10)
    beta_s <- c(beta_s, rep(0, settings[s, "p"] - length(beta_s)))
    out <- compute_estimates(setting,
                             n = nobs,
                             beta_star = beta_s,
                             maxit = maxit,
                             tolerance = tolerance,
                             verbose = 0,
                             current_setting = s,
                             max_setting = n_settings)
    out$sample <- s
    out
}, mc.cores = ncores)


estimates <- do.call("rbind", results)

out_file <- file.path(image_path, "ridge-mDYPL-comparisons.rda")
save(base_beta, maxit, tolerance, n_settings, settings, estimates, file = out_file)















######## Summaries


library("ggplot2")
library("dplyr")
library("tidyr")

if (interactive()) {
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    image_path <- file.path(project_path, "images")
    code_path <- file.path(project_path, "code")
    figure_path <- file.path(project_path, "figures")
}
load(file.path(image_path, "ridge-mDYPL-comparisons.rda"))
load(file.path(image_path, conjecture_model))
## Get conjecture constants
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]


source(file.path(code_path, "helper-functions.R"))


mses_mJPL <- estimates |>
    correct_mJPL_estimates(pk = po_k, pg = po_g, pr = po_r) |>
    ## Rescale estimates to match covariate setting (i.e. x_i ~ N(0, 1/ p I))
    mutate(estimate = sqrt(p) * estimate, truth = sqrt(p) * truth) |>
    group_by(method, gamma, kappa, p, sample) |>
    summarize(mse = mean((estimate - truth)^2),
              bias = mean(estimate - truth)) |>
    filter(method == "mJPL")

mses_mJPL <- mses_mJPL |>
    group_by(method, gamma, kappa, p) |>
    summarize(sd = sd(mse) / sqrt(n()),
              mse = mean(mse),
              bias = mean(bias))


J_mse <- rbind(data.frame(mses_mJPL[c("gamma", "kappa")],
                          alpha = 0,
                          mses_mJPL["mse"],
                          method = "mJPL",
                          mses_mJPL["sd"]),
               data.frame(mses_mJPL[c("gamma", "kappa")],
                          alpha = 1,
                          mses_mJPL["mse"],
                          method = "mJPL",
                          mses_mJPL["sd"]))

DY_mus <- read.csv(file.path(image_path, "DY_mus.csv"))
DY_sigmas <- read.csv(file.path(image_path, "DY_sigmas.csv"))
ridge_mus <- read.csv(file.path(image_path, "ridge_mus.csv"))
ridge_sigmas <- read.csv(file.path(image_path, "ridge_sigmas.csv"))
MLE_pars <- read.csv(file.path(image_path, "ML_pars.csv"))

MLE_pars <- MLE_pars |>
    transform(mse = kappa * sigma^2 / mu^2,
              method = "ML",
              sd = NA,
              alpha = 0)
MLE_pars <- rbind(MLE_pars,
                  MLE_pars |> transform(alpha = 1))

DY_pars <- data.frame(DY_mus, sigma = DY_sigmas$sigma) |>
    transform(mse = kappa * sigma^2 / mu^2,
              method = "mDYPL",
              sd = NA)
ridge_pars <- data.frame(ridge_mus, sigma = ridge_sigmas$sigma) |>
    transform(mse = sigma^2 / mu^2,
              method = "Ridge",
              alpha = 1 - lambda / max(lambda),
              sd = NA)

MLE_pars <- MLE_pars[c("gamma", "kappa", "alpha", "mse", "method", "sd")]
DY_pars <- DY_pars[c("gamma", "kappa", "alpha", "mse", "method", "sd")]
ridge_pars <- ridge_pars[c("gamma", "kappa", "alpha", "mse", "method", "sd")]


ggplot(rbind(MLE_pars, DY_pars, ridge_pars, J_mse)) +
    geom_line(aes(alpha, mse, col = method)) +
    geom_ribbon(aes(x = alpha, ymin = mse - 3 * sd, ymax = mse + 3 * sd, fill = method),
                alpha = 0.2) +
    ## geom_line(aes(alpha, mse + 3 * sd, col = method), lty = 3) +
    ## geom_line(aes(alpha, mse - 3 * sd, col = method), lty = 3) +
    facet_grid(gamma ~ kappa, scales= "free") +
    scale_x_reverse() +
    theme_minimal()

mses_other <- rbind(DY_pars, ridge_pars)
mses_other <- mses_other |>
    group_by(kappa, gamma, method) |>
    summarise(mse = min(mse, na.rm = TRUE))

mses_all <- rbind(data.frame(mses_mJPL[c("kappa", "gamma")],
                             method = "mJPL",
                             mse = mses_mJPL["mse"]),
                  mses_other) |>
    arrange(kappa, gamma, method)


rel_change <- function(mse, min, digits = 2) {
    round(mse / min, digits)
}


mses_all |>
    pivot_wider(names_from = method, values_from = mse) |>
    mutate(min_mse = pmin(Ridge, mDYPL, mJPL)) |>
    mutate(Ridge = rel_change(Ridge, min_mse),
           mDYPL = rel_change(mDYPL, min_mse),
           mJPL = rel_change(mJPL, min_mse))

















if (FALSE) {
    j <- 1
    oo <- compute_estimates(setting = settings[j, ],
                            n = nobs,
                            beta_star = rep(c(-10, 10, 0, 0, 0, 0, 0, 0), each = settings[j, "p"] / 8))

    kappa <- kappas[2]
    gamma <- gammas[3]
    beta_star <- rep(c(-10, 10, 0, 0, 0, 0, 0, 0), each = nobs * kappa / 8)
    beta0 <- 0
    dd <- simulate_ZSC2022(n = nobs,
                           kappa = kappa,
                           gamma = gamma,
                           beta0 = beta0,
                           beta_star = beta_star,
                           R = diag(length(beta_star)))
    true_betas <- c(attr(dd, "beta0"), attr(dd, "beta"))
    dd <- data.frame(Y = dd$Y, X = dd$X)
    dd$Y <- 0.5 * (dd$Y + 1)
    p <- length(beta_star)
    form <- as.formula(paste("Y ~ -1 + ", paste("X", 1:(p + 1), sep = ".", collapse = " + ")))









seeds <- unique(round(runif(nrow(settings), 0, 1000 * nrow(settings))))
settings$seed <- seeds
## Just set to FALSE to avoid computing the MLE
settings$mle_exists <- FALSE

time_br <- system.time(
    fit_br <- glm(form, family = binomial(), data = dd,
                  type = "MPL_Jeffreys",
                  method = "brglm_fit",
                  max_step_factor = 1,
                  check_aliasing = FALSE,
                  start = rep(0, p + 1),
                  epsilon = 1e-03,
                  maxit = 250,
                  trace = 0)
)


correct_mJPL_estimates(pk = po_k, pg = po_g, pr = po_r)
}
