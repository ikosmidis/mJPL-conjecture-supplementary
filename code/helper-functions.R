## Estimate probability of separation
estimate_inf_prob <- function(setting,
                              beta_star_setting = "u2",
                              n = 2000,
                              n_simu = 50,
                              ncores = 10) {
    kappa <- setting$kappa
    gamma <- setting$gamma
    seed <- setting$seed
    rhosq <- setting$rhosq
    beta0 <- gamma * sqrt(rhosq)
    psi <- setting$psi
    p <- ceiling(n * kappa)
    R <- chol(AR1cor(p, psi))
    nz <- ceiling(p * 0.2)
    beta_star <- switch(beta_star_setting,
                        "s1" = seq(-10, 10, length.out = p),
                        "s2" = c(rep(-10, nz), rep(10, nz), rep(0, p - 2 * nz)),
                        "u1" = c(rep(-3, nz), rep(-1, nz), rep(0, p - 3 * nz), rep(1, nz)),
                        "u2" = seq(1, 10, length.out = p),
                        stop("invalid beta star setting"))
    set.seed(seed)

    datasets <- lapply(1:n_simu, function(j) {
        dd <- simulate_ZSC2022(n = n,
                               kappa = kappa,
                               gamma = gamma,
                               beta0 = beta0,
                               R = R,
                               beta_star = beta_star)
        true_betas <- c(attr(dd, "beta0"), attr(dd, "beta"))
        dd <- data.frame(Y = dd$Y, X = dd$X)
        dd$Y <- 0.5 * (dd$Y + 1)
        dd
    })

    form <- as.formula(paste("Y ~ -1 + ", paste("X", 1:(p + 1), sep = ".", collapse = " + ")))
    res_sep <- mclapply(1:n_simu, function(j) {
        glm(form, family = binomial(), data = datasets[[j]],
            method = detectseparation::detect_infinite_estimates)$outcome
    }, mc.cores = ncores)
    mean(unlist(res_sep))
}

## Phase transition curves
compute_pt <- function(beta0 = 0,
                       gamma_grid = seq(0, 20, length = 100),
                       nsimu = 1000000,
                       ncores = 10,
                       XZU = NULL) {
    ## Theorem 2.1 of DOI: 10.1214/18-AOS1789
    if (is.null(XZU)) {
        X <- rnorm(nsimu)
        U <- runif(nsimu)
        Z <- rnorm(nsimu)
    } else {
        X <- XZU$X
        U <- XZU$U
        Z <- XZU$Z
    }
    kappa <- function(gamma) {
        gamma0 <- sqrt(gamma^2 - beta0^2)
        Y <- -1 + 2 * (plogis(beta0 + gamma0 * X) > U)
        V <- X * Y
        obj <- function(ts) {
            mean(pmax(ts[1] * Y + ts[2] * V - Z, 0)^2)
        }
        optim(c(0, 0), obj, method = "BFGS")$value
    }
    kappa <- unlist(mclapply(1:length(gamma_grid), function(k) {
        kappa(gamma = gamma_grid[k])
    }, mc.cores = ncores))
    data.frame(kappa = kappa, gamma = gamma_grid)
}


plot_pt <- function(PT, max_kappa = NULL, cols = c("#F1F1F1", "#FFFFFF")) {
    PT <- PT[order(PT$kappa), ]
    max_gamma <- max(PT$gamma)
    max_kappa <- ifelse(is.null(max_kappa), max(PT$kappa), max_kappa)
    polygon1 <- data.frame(kappa = c(PT$kappa, max_kappa, max_kappa),
                           gamma = c(PT$gamma, 0, max_gamma))
    polygon2 <- data.frame(kappa = c(0, 0, PT$kappa),
                           gamma = c(0, max_gamma, PT$gamma))
    ggplot(PT) +
        geom_polygon(data = polygon1, aes(kappa, gamma), fill = cols[1]) +
        geom_polygon(data = polygon2, aes(kappa, gamma), fill = cols[2]) +
        geom_line(aes(kappa, gamma), col = "grey") +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = expression(kappa), y = expression(gamma))
}

AR1cor <- function(p, rho) {
    sapply(1:p, \(k) rho^abs(k - 1:p))
}

simulate_ZSC2022_Bern <- function(n,
                                  kappa = 0.2,
                                  gamma = 1,
                                  beta0 = 0,
                                  prob = 0.5,
                                  beta_star = rnorm(p)) {
    ## Section 3.2.1. of DOI: 10.3150/21-BEJ1401
    p <- ceiling(kappa * n)
    gamma0 <- sqrt(gamma^2 - beta0^2)
    X <- cbind(1, matrix(rbinom(n * p, 1, prob), nrow = n))
    beta <- gamma0 * beta_star / sqrt(sum(beta_star^2) * prob * (1 - prob))
    eta <- drop(X %*% c(beta0, beta))
    Y <- 2 * (plogis(eta) > runif(n)) - 1
    out <- list(Y = Y, X = X)
    attr(out, "gamma0") <- gamma0
    attr(out, "beta") <- beta
    attr(out, "beta0") <- beta0
    attr(out, "kappa") <- kappa
    out
}


simulate_ZSC2022 <- function(n,
                             kappa = 0.2,
                             gamma = 1,
                             R = NULL, ## The upper triangular factor of the Cholesky deco, i.e. R'R = x.
                             beta0 = 0,
                             beta_star = rnorm(p)) {
    ## Section 3.2.1. of DOI: 10.3150/21-BEJ1401
    p <- ceiling(kappa * n)
    gamma0 <- sqrt(gamma^2 - beta0^2)
    if (is.null(R)) {
        R <- chol(AR1cor(p, 0.5))
    }
    Sigma <- crossprod(R)
    X <- cbind(1, MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma))
    v_star <- R %*% beta_star
    beta <- gamma0 * beta_star / sqrt(sum(v_star^2))
    eta <- drop(X %*% c(beta0, beta))
    Y <- 2 * (plogis(eta) > runif(n)) - 1
    out <- list(Y = Y, X = X)
    attr(out, "gamma0") <- gamma0
    attr(out, "beta") <- beta
    attr(out, "beta0") <- beta0
    attr(out, "kappa") <- kappa
    out
}


compute_estimates <- function(setting,
                              beta_star_setting = "u1",
                              n = 1000,
                              maxit = 250,
                              tolerance = 1e-03,
                              verbose = 0,
                              current_setting = NULL,
                              max_setting = NULL) {
    kappa <- setting$kappa
    gamma <- setting$gamma
    seed <- setting$seed
    rhosq <- setting$rhosq
    beta0 <- gamma * sqrt(rhosq)
    psi <- setting$psi
    p <- ceiling(n * kappa)
    R <- chol(AR1cor(p, psi))
    nz <- ceiling(p * 0.2) ## nz_perc -x and nz_perc +x
    beta_star <- switch(beta_star_setting,
                        "s1" = seq(-10, 10, length.out = p),
                        "s2" = c(rep(-10, nz), rep(10, nz), rep(0, p - 2 * nz)),
                        "u1" = c(rep(-3, nz), rep(-1, nz), rep(0, p - 3 * nz), rep(1, nz)),
                        "u2" = seq(1, 10, length.out = p),
                        stop("invalid beta star setting"))
    set.seed(seed)

    dd <- simulate_ZSC2022(n = n,
                           kappa = kappa,
                           gamma = gamma,
                           beta0 = beta0,
                           R = R,
                           beta_star = beta_star)
    true_betas <- c(attr(dd, "beta0"), attr(dd, "beta"))
    dd <- data.frame(Y = dd$Y, X = dd$X)
    dd$Y <- 0.5 * (dd$Y + 1)


    form <- as.formula(paste("Y ~ -1 + ", paste("X", 1:(p + 1), sep = ".", collapse = " + ")))
    time_br <- system.time(
        fit_br <- glm(form, family = binomial(), data = dd,
                      type = "MPL_Jeffreys",
                      method = "brglm_fit",
                      max_step_factor = 1,
                      check_aliasing = FALSE,
                      start = rep(0, p + 1),
                      epsilon = tolerance,
                      maxit = maxit,
                      trace = verbose)
    )
    coefs_br <- coef(fit_br)
    ses_br <- summary(fit_br)$coef[, "Std. Error"]
    elapsed_br <- time_br[["elapsed"]]
    iter_br <- fit_br$iter

    libeta_br <- max(abs(vcov(fit_br) %*% fit_br$grad[1:(p + 1)]), na.rm = TRUE)
    if (setting$mle_exists) {
        time_ml <- system.time(
            fit_ml <- glm(form, family = binomial(), data = dd,
                          type = "ML",
                          method = "brglm_fit",
                          max_step_factor = 1,
                          check_aliasing = FALSE,
                          start = coefs_br,
                          epsilon = tolerance,
                          maxit = maxit,
                          trace = verbose)
        )

        elapsed_ml <- time_ml[["elapsed"]]
        iter_ml <- fit_ml$iter
        libeta_ml <- max(abs(vcov(fit_ml) %*% fit_ml$grad[1:(p + 1)]), na.rm = TRUE)
        ## If ML did not converge then set estimates and ses to NA
        if (libeta_ml > 10) {
            coefs_ml <- ses_ml <- rep(NA, length(true_betas))
        } else {
            coefs_ml <- coef(fit_ml)
            ses_ml <- summary(fit_ml)$coef[, "Std. Error"]
        }
    } else {
        coefs_ml <- ses_ml <- rep(NA, length(true_betas))
        elapsed_ml <- iter_ml <- libeta_ml <- NA
    }
    if (!is.null(current_setting) & !is.null(max_setting)) {
        cat(current_setting, "/", max_setting, "|| ")
    }
    diagnostics <- paste0("|delta_ML|_inf = ", round(libeta_ml, 5), paste0("(", iter_ml,")"),
                         ", |delta_BR|_inf = ", round(libeta_br, 5), paste0("(", iter_br,")"))
    cat("β* = ", beta_star_setting,
        ", ψ = ", sprintf("%.2f", psi),
        ", ρ^2 = ", sprintf("%.2f", rhosq),
        ", κ = ", sprintf("%.2f", kappa),
        ", γ = ", sprintf("%2.2f", gamma),
        ", p = ",  sprintf("%2.2f", p),
        ", n = ", nobs,
        " : ", diagnostics,  "\n", sep = "")

    results <- data.frame(method = rep(c("ML", "mJPL"), each = length(true_betas)),
                          estimate = c(coefs_ml, coefs_br),
                          se = c(ses_ml, ses_br),
                          truth = c(true_betas, true_betas),
                          parameter = rep(seq_along(true_betas) - 1, 2),
                          psi = psi,
                          rhosq = rhosq,
                          kappa = kappa,
                          gamma = gamma,
                          mle_exists = setting$mle_exists,
                          n = n,
                          p = p)
    performance <- data.frame(method = c("ML", "mJPL"),
                              elapsed = c(elapsed_ml, elapsed_br),
                              iter = c(iter_ml, iter_br),
                              libeta = c(libeta_ml, libeta_br),
                              psi = psi,
                              rhosq = rhosq,
                              kappa = kappa,
                              gamma = gamma,
                              mle_exits = setting$mle_exists,
                              n = n,
                              p = p)
    rownames(results) <- rownames(performance) <- NULL
    attr(results, "performance") <- performance
    attr(results, "seed") <- seed
    results
}

compute_estimates_bern <- function(setting,
                                   beta_star_setting = "u1",
                                   n = 1000,
                                   maxit = 250,
                                   tolerance = 1e-03,
                                   verbose = 0,
                                   current_setting = NULL,
                                   max_setting = NULL) {
    kappa <- setting$kappa
    gamma <- setting$gamma
    seed <- setting$seed
    rhosq <- setting$rhosq
    beta0 <- gamma * sqrt(rhosq)
    prob <- setting$prob
    p <- ceiling(n * kappa)
    nz <- ceiling(p * 0.2) ## nz_perc -x and nz_perc +x
    beta_star <- switch(beta_star_setting,
                        "s1" = seq(-10, 10, length.out = p),
                        "s2" = c(rep(-10, nz), rep(10, nz), rep(0, p - 2 * nz)),
                        "u1" = c(rep(-3, nz), rep(-1, nz), rep(0, p - 3 * nz), rep(1, nz)),
                        "u2" = seq(1, 10, length.out = p),
                        stop("invalid beta star setting"))
    set.seed(seed)

    dd <- simulate_ZSC2022_Bern(n = n,
                                kappa = kappa,
                                gamma = gamma,
                                beta0 = beta0,
                                prob = prob,
                                beta_star = beta_star)
    true_betas <- c(attr(dd, "beta0"), attr(dd, "beta"))
    dd <- data.frame(Y = dd$Y, X = dd$X)
    dd$Y <- 0.5 * (dd$Y + 1)

    form <- as.formula(paste("Y ~ -1 + ", paste("X", 1:(p + 1), sep = ".", collapse = " + ")))
    time_br <- system.time(
        fit_br <- glm(form, family = binomial(), data = dd,
                      type = "MPL_Jeffreys",
                      method = "brglm_fit",
                      max_step_factor = 1,
                      check_aliasing = FALSE,
                      start = rep(0, p + 1),
                      epsilon = tolerance,
                      maxit = maxit,
                      trace = verbose)
    )
    coefs_br <- coef(fit_br)
    ses_br <- summary(fit_br)$coef[, "Std. Error"]
    elapsed_br <- time_br[["elapsed"]]
    iter_br <- fit_br$iter

    libeta_br <- max(abs(vcov(fit_br) %*% fit_br$grad[1:(p + 1)]), na.rm = TRUE)
    if (!is.null(current_setting) & !is.null(max_setting)) {
        cat(current_setting, "/", max_setting, "|| ")
    }
    diagnostics <- paste0("|delta_BR|_inf = ", round(libeta_br, 5), paste0("(", iter_br,")"))
    cat("β* = ", beta_star_setting,
        ", π = ", sprintf("%.2f", prob),
        ", ρ^2 = ", sprintf("%.2f", rhosq),
        ", κ = ", sprintf("%.2f", kappa),
        ", γ = ", sprintf("%2.2f", gamma),
        ", p = ",  sprintf("%2.2f", p),
        ", n = ", nobs,
        " : ", diagnostics,  "\n", sep = "")

    results <- data.frame(method = rep("mJPL", each = length(true_betas)),
                          estimate = c(coefs_br),
                          se = c(ses_br),
                          truth = c(true_betas),
                          parameter = seq_along(true_betas) - 1,
                          prob = prob,
                          rhosq = rhosq,
                          kappa = kappa,
                          gamma = gamma,
                          n = n,
                          p = p)
    performance <- data.frame(method = "mJPL",
                              elapsed = elapsed_br,
                              iter = iter_br,
                              libeta = libeta_br,
                              prob = prob,
                              rhosq = rhosq,
                              kappa = kappa,
                              gamma = gamma,
                              n = n,
                              p = p)
    rownames(results) <- rownames(performance) <- NULL
    attr(results, "performance") <- performance
    attr(results, "seed") <- seed
    results
}


correct_mJPL_estimates <- function(results, pk = 1.2, pg = 1.1, pr = 1/3) {
    results |> mutate(estimate = ifelse(method == "mJPL" & !mle_exists,
                                        estimate / (kappa^pk * gamma^pg * (1 - rhosq)^pr),
                                        estimate))
}

## For a single kappa/gamma combination
plot_results <- function(summaries, cols = c("#CF4446", "#00AD9A", "#FB9A06"),
                         p_alpha = 0.2, p_size = 1,
                         type = c("estimate_vs_truth", "estimate_and_truth")) {
    type <- match.arg(type)
    summaries$method <- factor(summaries$method, levels = c("ML", "mJPL"), ordered = TRUE)
    ## If the ML estimates do not exist plot nothing
    if (all(!summaries$mle_exists)) {
        summaries <- summaries |> subset(method != "ML")
    }
    lims <- with(summaries, range(c(truth, estimate)))
    if (type == "estimate_vs_truth") {
        p1 <- ggplot(summaries) +
            geom_point(aes(x = truth, y = estimate, col = method), alpha = p_alpha, size = p_size) +
            geom_abline(aes(intercept = 0, slope = 1), col = "black", lwd = 0.5) +
            geom_smooth(aes(x = truth, y = estimate),
                        method = "lm", formula = "y ~ x",
                        se = FALSE, col = cols[2], lwd = 0.5) +
            coord_cartesian(x = lims, y = lims) +
            scale_colour_manual(values = c("ML" = cols[1], "mJPL" = cols[3]))
    }
    if (type == "estimate_and_truth") {
        sd <- summaries |> group_by(kappa, gamma, method, truth) |>
            summarize(e_mean = mean(estimate)) |>
            inner_join(summaries, c("kappa", "gamma", "method", "truth"))
        p1 <- ggplot(summaries) +
            geom_point(aes(x = parameter, y = estimate, col = method), alpha = p_alpha, size = p_size) +
            geom_line(aes(x = parameter, y = truth), col = "black", lwd = 0.5) +
            geom_line(data = sd, aes(x = parameter, y = e_mean), col = cols[2], lwd = 0.5) +
           coord_cartesian(y = lims) +
            scale_colour_manual(values = c("ML" = cols[1], "mJPL" = cols[3]))
    }
    p1 + theme_minimal() +
        facet_grid(~ method) +
        theme(legend.position = "none") +
        theme(
            plot.background = element_rect(fill= 'transparent', color = "grey"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_blank()
        ) +
        labs(x = NULL, y = NULL)
}




