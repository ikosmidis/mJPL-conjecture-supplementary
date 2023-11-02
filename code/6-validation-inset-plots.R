library("parallel")
library("ggplot2")
library("ggpp")
library("patchwork")
library("dplyr")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ## Number of observations per dataset
    nobs <- 3000
    beta_star_setting <- "u1"
    psi <- 0.3
    rhosq <- 0.5
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.5.rda"
}

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

## Get pt curves
load(file.path(image_path, "pt-curves-test.rda"))

ind_rhosq <- which(rhosq == rhosq_grid)
cpsi <- psi
crhosq <- rhosq
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

pdf(file.path("~/Repositories/mJPL-conjecture/figures", "illustration.pdf"),
    width = exp_w / 300, height = exp_h * 1.8 / 300)
print(f_u / f_c)
dev.off()

