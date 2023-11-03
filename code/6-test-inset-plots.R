library("parallel")
library("ggplot2")
library("ggpp")
library("patchwork")
library("dplyr")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ## Number of observations per dataset
    nobs <- c(3000, 1000)
    beta_star_setting <- c("u1", "s1")
    psi <- c(0.3, 0.6)
    rhosq <- c(0.6, 0.1)
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
}

## Some plotting constants
exp_h <- 5 * 200 * 1.5
vp_h <- 0.145
exp_w <- exp_h * sqrt(2)
exp_ratio <- exp_h / exp_w
vp_w <- vp_h * exp_ratio
p_size <- 0.5
alpha_fac <- 1000

image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

load(file.path(image_path, conjecture_model))
## Get conjecture constants
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]

## Get pt curves
load(file.path(image_path, "pt-curves-test.rda"))



settings <- data.frame(
    nobs = nobs,
    beta_star_setting = beta_star_setting,
    psi = psi,
    rhosq = rhosq)
plots <- as.list(rep(0, nrow(settings)))

for (j in seq.int(nrow(settings))) {
    bss <- settings[j, "beta_star_setting"]
    cpsi <- settings[j, "psi"]
    crhosq <- settings[j, "rhosq"]
    cnobs <- settings[j, "nobs"]

    plot_type <- switch(bss,
                        "s1" = "estimate_vs_truth",
                        "u2" = "estimate_vs_truth",
                        "s2" = "estimate_and_truth",
                        "u1" = "estimate_and_truth")

    results_path <- file.path(image_path, paste0("estimates-n-", cnobs,
                                                 "-beta-", bss,
                                                 "-test",
                                                 ".rda"))

    load(results_path)

    prkg <- unique(estimates[c("psi", "rhosq", "kappa", "gamma", "mle_exists")])
    rhosq_grid <- unique(prkg$rhosq)


    ind_rhosq <- which(crhosq == rhosq_grid)
    tt <- substitute(paste(beta, " configuration: "*b0, ", ",
                           n == x0, ", ",
                           psi == z0, ", ",
                           rho^2 == y0),
                     list(x0 = cnobs, y0 = crhosq, z0 = cpsi, b0 = bss))

    kappa_gamma <- prkg |> filter(psi == cpsi, rhosq == crhosq)

    insets_estimates <- as.list(rep(NA, nrow(kappa_gamma)))
    for (wh in 1:nrow(kappa_gamma)) {
        ckappa <- kappa_gamma[wh, "kappa"]
        cgamma <- kappa_gamma[wh, "gamma"]
        p_alpha <- alpha_fac * (1 - ckappa) / cnobs
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
        p_alpha <- alpha_fac * (1 - ckappa) / cnobs
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
    plots[[j]] <- f_u / f_c
    cat(j, "\n")
}


pdf(file.path(project_path, "figures/illustration.pdf"),
    width = 9, height = 9 / sqrt(2))
print(plots[[1]] | plots[[2]])
dev.off()

