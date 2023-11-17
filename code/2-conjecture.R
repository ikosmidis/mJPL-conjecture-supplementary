library("dplyr")
library("ggplot2")
library("patchwork")
library("colorspace")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    nobs <- 2000
    psi <- 0
    beta_star_setting <- "u2"
}
image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

results_file <- file.path(image_path, paste0("estimates-n-", nobs,
                                         "-beta-", beta_star_setting,
                                         "-psi-", psi,
                                         "-training",
                                         ".rda"))
load(results_file)
load(file.path(image_path, "design-training.rda"))

## Extract mJPL estimates, regress estimates vs truth across
## repetitions for each rhosq-kappa-gamma combination, and collect
## intercepts, slopes and residual standard errors
coefs_mJPL <- estimates |>
    filter(parameter > 0, method == "mJPL", gamma > 0.5) |>
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
        sigma = sapply(fit, "[", "sigma"),
        gamma0 = gamma * sqrt(1 - rhosq))

coefs_mJPL <- coefs_mJPL |> inner_join(design, by = c("kappa", "gamma", "rhosq", "mle_exists"))

## Set a cutpoint that determines small/moderate rhosq values
cut_point <- 0.7


## Estimate q() based on the slopes of mJPL estimates vs truth when
## the MLE exists and when it does not

## Linear regression log(a1) ~ log(gamma) + log(kappa) + log(1 - rhosq)
lm_noMLE_rhosq <- lm(log(a1) ~ log(kappa) + log(gamma) + log(1 - rhosq),
                     data = coefs_mJPL |> filter(!mle_exists, rhosq < cut_point))

lm_noMLE_gamma0 <- lm(log(a1) ~ log(kappa) + log(gamma) +  log(gamma0),
                      data = coefs_mJPL |> filter(!mle_exists, rhosq < cut_point))

## Gamma-response GLM with log link a1 ~ log(gamma) + log(kappa) + log(1 - rhosq)
glm_noMLE_rhosq <- glm(
    a1 ~ log(kappa) + log(gamma) + log(1 - rhosq),
    data = filter(coefs_mJPL, !mle_exists, rhosq < cut_point),
    family = Gamma(link = "log"))

glm_noMLE_gamma0 <- glm(
    a1 ~ log(kappa) + log(gamma) + log(gamma0),
    data = filter(coefs_mJPL, !mle_exists, rhosq < cut_point),
    family = Gamma(link = "log"))


out_file <- file.path(image_path, paste0("conjecture-n-", nobs,
                                         "-beta-", beta_star_setting,
                                         "-psi-", psi, ".rda"))

save(coefs_mJPL, cut_point,
     lm_noMLE_rhosq, lm_noMLE_gamma0,
     glm_noMLE_rhosq, glm_noMLE_gamma0,
     file = out_file)




## Plots and table for main text
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


coefs_mJPL <- coefs_mJPL |>
    mutate(rhosq_fac = ifelse(rhosq <= cut_point,
                              paste("rho^2 <=", cut_point),
                              paste("rho^2 >", cut_point)))

p_hi_int <- ggplot(coefs_mJPL |> filter(method == "mJPL")) +
    geom_histogram(aes(x = a0), fill = okabe[1]) +
    geom_rug(aes(x = a0)) +
    geom_vline(aes(xintercept = 0), col = okabe[5], lty = 2) +
    theme_minimal() +
    labs(x = expression(delta[0]^"*"))

midp <- 1 + abs(min((coefs_mJPL |> filter(mle_exists, method == "mJPL"))$a1) - 1)
## Slopes when MLE exists
p_kg <- ggplot(coefs_mJPL |> filter(mle_exists, method == "mJPL")) +
    geom_point(aes(kappa, gamma, color = a1, fill = a1), shape = 24) +
    geom_point(data = coefs_mJPL |> filter(!mle_exists, method == "mJPL"),
               aes(kappa, gamma), col = "grey", alpha = 0.2) +
    lims(x = c(0, 0.6), y = c(0, 20)) +
    scale_color_continuous_divergingx(
        palette = "PuOr",
        mid = midp,
        l1 = 40, l2 = 80, l3 = 40,
        c1 = 80, c2 = 50, c3 = 80,
        alpha = 0.8,
        name = expression(delta[1]^"*"),
        n.breaks = 4) +
    scale_fill_continuous_divergingx(
        palette = "PuOr",
        mid = midp,
        l1 = 40, l2 = 80, l3 = 40,
        c1 = 80, c2 = 50, c3 = 80,
        alpha = 0.8,
        name = expression(delta[1]^"*"),
        n.breaks = 4) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = expression(kappa), y = expression(gamma), title = "MLE exists")

## Slopes when MLE does not exists
po_int <- coef(glm_noMLE_rhosq)["(Intercept)"]
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]


xlims <- with(coefs_mJPL |> filter(!mle_exists, method == "mJPL"), {
    r1 <- range(log(kappa^po_k * gamma^po_g * (1 - rhosq)^po_r))
    r2 <- range(log(a1))
    c(min(r1[1], r2[1]), max(r1[2], r2[2]))
})

p_s_noMLE <-
    ggplot(coefs_mJPL |> filter(!mle_exists, method == "mJPL")) +
    geom_point(aes(x = log(kappa^po_k * gamma^po_g * (1 - rhosq)^po_r),
                   y = log(a1), col = rhosq)) +
    geom_abline(aes(intercept = 0, slope = 1), col = "darkgrey") +
    facet_grid(~ rhosq_fac, labeller = label_parsed) +
    scale_color_continuous_divergingx(
        palette = "PiYG",
        mid = 0.5,
        l1 = 30, l2 = 80, l3 = 30,
        c1 = 80, c2 = 50, c3 = 80,
        alpha = 0.8,
        name = expression(rho^2)) +
    lims(x = xlims, y = xlims) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = expression(log * " " * q(kappa, gamma, gamma[0] * " ; " * b^"*")),
         y = expression(log * " " * delta[1]^"*"),
         title = "MLE does not exist")

pdf(file.path(figure_path, "conjecture.pdf"), width = 6, height = 6)
## pdf(file.path("~/Repositories/mJPL-conjecture/figures", "conjecture.pdf"), width = 6, height = 6)
print(((p_s_noMLE + p_kg) + plot_layout(widths = c(2, 1))) / p_hi_int)
dev.off()

## Estimates and confidence intervals
library("car")
set.seed(111)
bglm <- Boot(glm_noMLE_gamma0, method = "case", R = 9999, ncores = 1)
print(memisc:::toLatex.default(Confint(bglm), digits = 3))
## deviance explained
with(summary(glm_noMLE_gamma0), 1 - deviance/null.deviance)


