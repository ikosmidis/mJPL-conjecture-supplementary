library("parallel")
library("ggplot2")
library("ggpp")
library("patchwork")
library("dplyr")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ## Number of observations per dataset
    nobs <- 2000
    ## beta star setting
    beta_star_setting <- "u1"
    ncores <- 10
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.5.rda"
}
ns <- 200000
exp_h <- 5 * 200 * 1.5
vp_h <- 0.145
exp_w <- exp_h * sqrt(2)
exp_ratio <- exp_h / exp_w
vp_w <- vp_h * exp_ratio
## transparency of the points
p_size <- 0.5
alpha_fac <- 1000

plot_type <- switch(beta_star_setting,
                    "s1" = "estimate_vs_truth",
                    "u2" = "estimate_vs_truth",
                    "s2" = "estimate_and_truth",
                    "u1" = "estimate_and_truth")

image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

results_path <- file.path(image_path, paste0("estimates-n-", nobs,
                                             "-beta-", beta_star_setting,
                                             "-test",
                                             ".rda"))

load(results_path)
load(file.path(image_path, conjecture_model))
## Get conjecture constants
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]

prkg <- unique(estimates[c("psi", "rhosq", "kappa", "gamma", "mle_exists")])
rhosq_grid <- unique(prkg$rhosq)
psi_grid <- unique(prkg$psi)

## Get pt curves
load(file.path(image_path, "pt-curves-test.rda"))

## Plots
for (ind_psi in seq_along(psi_grid)) {
    for (ind_rhosq in seq_along(rhosq_grid)) {
        cpsi <- psi_grid[ind_psi]
        crhosq <- rhosq_grid[ind_rhosq]
        tt <- substitute(paste(beta, " setting: "*b0, ", ", n == x0, ", ", psi == z0, ", ", rho^2 == y0),
                         list(x0 = nobs, y0 = crhosq, z0 = cpsi, b0 = beta_star_setting))

        kappa_gamma <- prkg |> filter(psi == cpsi, rhosq == crhosq)

        insets_estimates <- as.list(rep(NA, nrow(kappa_gamma)))
        for (wh in 1:nrow(kappa_gamma)) {
            ckappa <- kappa_gamma[wh, "kappa"]
            cgamma <- kappa_gamma[wh, "gamma"]
            p_alpha <- alpha_fac * (1 - ckappa) / nobs
            ep <- estimates |>
                dplyr::filter(parameter > 0, psi == cpsi, rhosq == crhosq, kappa == ckappa, gamma == cgamma) |>
                plot_results(p_alpha = p_alpha, p_size = p_size, type = plot_type)
            insets_estimates[[wh]] <- tibble(x = ckappa,
                                             y = cgamma,
                                             plot = list(ep))
        }
        out <- plot_pt(pt_curve[[ind_rhosq]], max_kappa = 0.6) +
            geom_point(data = kappa_gamma, aes(x = kappa, y = gamma), pch = 18, size = 1, col = "grey")
        for (wh in 1:nrow(kappa_gamma)) {
            ckappa <- kappa_gamma[wh, "kappa"]
            cgamma <- kappa_gamma[wh, "gamma"]
            out <- out +  geom_plot(data = insets_estimates[[wh]],
                                    aes(x = x, y = y, label = plot),
                                    vp.width = ifelse(kappa_gamma[wh, "mle_exists"], 2 * vp_w, vp_w),
                                    vp.height = vp_h,
                                    vjust = "bottom", hjust = "left")
        }
        out_u <- out + labs(title = tt)

        ## Plot mJPL and ML estimates (the latter only when they exist),
        ## adjusting the former when the ML estimates do not exist
        insets_estimates <- as.list(rep(NA, nrow(kappa_gamma)))
        for (wh in 1:nrow(kappa_gamma)) {
            ckappa <- kappa_gamma[wh, "kappa"]
            cgamma <- kappa_gamma[wh, "gamma"]
            p_alpha <- alpha_fac * (1 - ckappa) / nobs
            ep <- estimates |>
                filter(parameter > 0, psi == cpsi, rhosq == crhosq, kappa == ckappa, gamma == cgamma) |>
                correct_mJPL_estimates(pk = po_k, pg = po_g, pr = po_r) |>
                plot_results(p_alpha = p_alpha, p_size = p_size, type = plot_type)
            insets_estimates[[wh]] <- tibble(x = ckappa,
                                             y = cgamma,
                                             plot = list(ep))
        }
        out <- plot_pt(pt_curve[[ind_rhosq]], max_kappa = 0.6) +
            geom_point(data = kappa_gamma, aes(x = kappa, y = gamma), pch = 18, size = 1,
                       col = "grey")
        for (wh in 1:nrow(kappa_gamma)) {
            ckappa <- kappa_gamma[wh, "kappa"]
            cgamma <- kappa_gamma[wh, "gamma"]
            out <- out +  geom_plot(data = insets_estimates[[wh]],
                                    aes(x = x, y = y, label = plot),
                                    vp.width = ifelse(kappa_gamma[wh, "mle_exists"], 2 * vp_w, vp_w),
                                    vp.height = vp_h,
                                    vjust = "bottom", hjust = "left")
        }
        out_c <- out + labs(title = tt)

        f_u <- out_u +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            labs(x = NULL)
        f_c <- out_c + labs(title = NULL)
        base_name <- paste0("validation-n-", nobs, "-beta-", beta_star_setting, "-psi-", cpsi, "-rhosq-", crhosq)
        print(base_name)
        ## pdf(file.path(figure_path, paste0(base_name, ".pdf")),
        ##     width = exp_w / 300, height = exp_h * 1.8 / 300)
        print(f_u / f_c)
        ## dev.off()
   }
}








if (FALSE) {

    library("ggplot2")
    
    coefs <- estimates |>
        filter(parameter > 0) |>
        group_by(kappa, gamma, rhosq, psi, method, mle_exists) |>
        summarize(fit = {
            xmat <- cbind(1, truth)
            colnames(xmat) <- c("a0", "a1")
            if (all(is.na(estimate))) {
                list(rep(NA, 3))
            } else {
                obj <- lm.fit(x = xmat, y = estimate)
                list(c(coef(obj)[1],
                       coef(obj)[2],
                       sigma = sqrt(sum(obj$residuals^2) / (obj$df.residual))))
            }
        }) |> ## Unlist estimates
        mutate(
            a0 = sapply(fit, "[", "a0"),
            a1 = sapply(fit, "[", "a1"),
            sigma = sapply(fit, "[", "sigma")) |>
        mutate(psi_fac = paste("psi ==", psi),
               method = factor(method, levels = c("ML", "mJPL"), ordered = TRUE))


    Rsq <- function(slopes, kappa, gamma, rhosq, pk, pg, pr) {
        m_la1 <- mean(ll <- log(slopes))
        tss <- sum((ll - m_la1)^2)
        resid <- ll - pk * log(kappa) - pg * log(gamma) - pr * log(1 - rhosq)
        rss <- sum(resid^2)
        1 - rss / tss
    }

    rsqs <- coefs |>
        filter(method == "mJPL", !mle_exists) |>
        group_by(psi_fac, rhosq) |>
        summarize(rsq = Rsq(a1, kappa, gamma, rhosq, po_k, po_g, po_r))

    p_rsq_ne <- ggplot(rsqs, aes(rhosq, rsq)) +
        geom_line() +
        geom_point() +
        geom_hline(aes(yintercept = 1), col = "grey", lty = 2) +
        facet_grid(~ psi_fac, labeller = label_parsed) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = expression(rho^2), y = expression(R["test"]^2))

   p_hist_e <- ggplot(coefs |> filter(mle_exists, a1 < 4)) +
        geom_histogram(aes(y = a1)) +
        geom_rug(aes(y = a1)) +
        geom_hline(aes(yintercept = 1), col = "grey", lty = 2) +
        facet_grid(method ~ psi) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(y = expression(delta[1]^"*"))



    bias <- function(estimates, pk, pg, pr) {
        dev_un <- estimates$estimate - estimates$truth
        dev_ad <- estimates$estimate / (estimates$kappa^pk * estimates$gamma^pg * (1 - estimates$rhosq)^pr) - estimates$truth
        b_un <- mean(dev_un)
        b_ad <- mean(dev_ad)
        mse_un <- mean(dev_un^2)
        mse_ad <- mean(dev_ad^2)
        data.frame(bias = c(b_un, b_ad), mse = c(mse_un, mse_ad), type = c("vanilla", "rescaled"))
    }

}

if (FALSE) {

    coefs_mJPL <- estimates |>
        filter(parameter > 0, method == "mJPL") |>
        group_by(kappa, gamma, rhosq, psi, method, mle_exists) |>
        summarize(fit = {
            xmat <- cbind(1, truth)
            colnames(xmat) <- c("a0", "a1")
            obj <- lm.fit(x = xmat, y = estimate)
            list(c(coef(obj)[1],
                   coef(obj)[2],
                   sigma = sqrt(sum(obj$residuals^2) / (obj$df.residual))))
        }) |> ## Unlist estimates
        mutate(
            a0 = sapply(fit, "[", "a0"),
            a1 = sapply(fit, "[", "a1"),
            sigma = sapply(fit, "[", "sigma"))




    library("ggplot2")

    coefs_mJPL$rhosq_fac <- ifelse(coefs_mJPL$rhosq <= cut_point,
                                   paste("rhosq <=", cut_point),
                                   paste("rhosq >", cut_point))



    ggplot(coefs_mJPL |> filter(method == "mJPL")) +
        geom_histogram(aes(a0)) +
        geom_vline(aes(xintercept = 0), col = "red") +
        lims(x = c(-0.1, 0.1)) +
        facet_grid(mle_exists ~ rhosq_fac)

    ggplot(coefs_mJPL |> filter(mle_exists, method == "mJPL")) +
        geom_point(aes(x = seq_along(a1), y = a1, col = kappa)) +
        geom_hline(aes(yintercept = 1), col = "red") +
        facet_grid( ~ rhosq_fac) +
        lims(y = c(0, 2))


    ggplot(coefs_mJPL |> filter(!mle_exists, method == "mJPL")) +
        geom_point(aes(x = log(kappa^po_k * gamma^po_g * (1 - rhosq)^po_r),
                       y = log(a1), col = kappa)) +
        geom_abline(aes(intercept = 0, slope = 1), col = "red") +
        facet_grid(psi ~ rhosq_fac)


    rhosq0 <- 0
    psi0 <- 0.9
    ggplot(estimates |>
           correct_mJPL_estimates(pk = po_k, pg = po_g, pr = po_r) |>
           filter(parameter > 0, method == "mJPL", psi == psi0, rhosq == rhosq0),
           aes(x = truth, y = estimate)) +
        geom_point(aes(col = mle_exists)) +
        geom_abline(aes(intercept = 0, slope = 1), col = "red") +
        geom_smooth(method = "lm", formula = y ~ x) +
        facet_wrap(kappa ~ gamma, scales = "free")

    ggplot(estimates |>
           filter(parameter > 0, method == "mJPL", psi == psi0, rhosq == rhosq0),
           aes(x = truth, y = estimate)) +
        geom_point(aes(col = mle_exists)) +
        geom_abline(aes(intercept = 0, slope = 1), col = "red") +
        geom_smooth(method = "lm", formula = y ~ x) +
        facet_wrap(kappa ~ gamma, scales = "free")


}


## Standard errors
if (FALSE) {

    get_tausq <- function(p, psi, parameter) {
        case_when(
            parameter == 0 ~ NA,
            parameter == 1 ~ (1 - psi^2),
            parameter == p ~ (1 - psi^2),
            .default = (1 - psi^2) / (1 + psi^2)
        )
    }

    ests <- estimates |>
        filter(method == "mJPL", parameter > 0) |>
        group_by(kappa, gamma, rhosq, psi, method, mle_exists, p) |>
        mutate(tausq = get_tausq(p, psi, parameter),
               scaled_truth = ifelse(mle_exists, truth, truth * (kappa^po_k * gamma^po_g * (1 - rhosq)^po_r))) |>
        summarize(estv = var(sqrt(n) * (estimate - scaled_truth) * sqrt(tausq), na.rm = TRUE),
                  estimated_sd = sd(estimate),
                  mean_se = mean(se))

    ## Estimate q() based on the slopes of mJPL estimates vs truth when
    ## the MLE exists and when it does not
    mod_var <- lm(log(estv) ~ log(gamma) + log(kappa) + log(1 - rhosq),
                   data = ests |> filter(!mle_exists, rhosq < 0.6, gamma < 15, kappa < 0.5))


    ggplot(ests |> filter(!mle_exists, rhosq < 0.6, gamma < 15, kappa < 0.4),
           aes(coef(mod_var)["(Intercept)"] +
               log(kappa^coef(mod_var)["log(kappa)"] *
                   gamma^coef(mod_var)["log(gamma)"] *
                   (1 - rhosq)^coef(mod_var)["log(1 - rhosq)"]),
           log(estv))) +
        geom_point(aes(col = gamma)) +
        geom_abline(aes(intercept = 0,  slope = 1))


    ggplot(ests |> filter(!mle_exists, rhosq < 0.6, gamma < 15, kappa < 0.4),
           aes(coef(mod_var)["(Intercept)"] +
               log(kappa^coef(mod_var)["log(kappa)"] *
                   gamma^coef(mod_var)["log(gamma)"] *
                   (1 - rhosq)^coef(mod_var)["log(1 - rhosq)"]),
           log(estv))) +
        geom_point(aes(col = gamma)) +
        geom_abline(aes(intercept = 0,  slope = 1))


    ggplot(ests |> filter(!mle_exists),
           aes(coef(mod_var)["(Intercept)"] +
               log(kappa^coef(mod_var)["log(kappa)"] *
                   gamma^coef(mod_var)["log(gamma)"] *
                   (1 - rhosq)^coef(mod_var)["log(1 - rhosq)"]),
               log(estv))) +
        geom_point(aes(col = gamma)) +
        geom_abline(aes(intercept = 0,  slope = 1)) +
        facet_grid(gamma ~rhosq)


    ggplot(ests |> filter(mle_exists),
           aes(log(estimated_sd),
               log(mean_se))) +
        geom_point(aes(col = gamma)) +
        geom_abline(aes(intercept = 0,  slope = 1))



}
