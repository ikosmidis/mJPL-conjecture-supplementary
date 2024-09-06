library("parallel")
library("brglm2")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    image_path <- file.path(project_path, "images")
    code_path <- file.path(project_path, "code")
    nobs <- 2000
    repetitions <- 50
    ncores <- 10
    ## seed for generating unique seeds
    seed <- 207
}
maxit <- 100
tolerance <- 1e-03

source(file.path(code_path, "helper-functions.R"))

kappas <- (1:6) / 10
gammas <- c(1, seq(2.5, 15, 2.5))

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

settings <- settings |> transform(p = nobs * kappa)
settings$rhosq <- 0

results <- mclapply(1:n_settings, function(s) {
    setting <- settings[s, ]
    p <- settings[s, "p"]
    g <- settings[s, "gamma"]
    k <- settings[s, "kappa"]
    mle_exists <- settings[s, "mle_exists"]
    set.seed(setting$seed)
    beta_s <- seq(-10, 10, length = p)
    beta <- c(0, g * beta_s / sqrt(sum(beta_s^2)))
    X <- cbind(1, matrix(rnorm(nobs * p, 0, 1), nrow = nobs))
    probs <- plogis(X %*% beta)
    y <- rbinom(nobs, 1, probs)
    start <- rep(0, p + 1)
    fit_br <- glm(y ~ -1 + X, method = brglm_fit, family = binomial(),
                  max_step_factor = 1,
                  start = start,
                  check_aliasing = FALSE,
                  epsilon = 1e-03,
                  maxit = maxit,
                  trace = 0)
    libeta_br <- max(abs(vcov(fit_br) %*% fit_br$grad[1:(p + 1)]), na.rm = TRUE)
    diagnostics <- paste0("|delta_BR|_inf = ", round(libeta_br, 5), paste0("(", fit_br$iter,")"))
    cat(s, "/", n_settings, " : ",
        "κ = ", sprintf("%.2f", k),
        ", γ = ", sprintf("%2.2f", g),
        ", p = ",  sprintf("%2.2f", p),
        ", n = ", nobs,
        " : ", diagnostics,  "\n", sep = "")
    res <- data.frame(method = "mJPL",
                      estimate = coef(fit_br),
                      parameter = 0:p,
                      truth = beta,
                      kappa = k,
                      gamma = g,
                      n = nobs,
                      p = p,
                      rhosq = settings[s, "rhosq"],
                      sample = s,
                      mle_exists = settings[s, "mle_exists"])
}, mc.cores = ncores)

res <- do.call("rbind", results)

save(res, file = file.path(image_path, "mJPL_ridge_comparison.rda"))

