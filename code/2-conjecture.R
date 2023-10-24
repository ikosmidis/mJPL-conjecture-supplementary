library("dplyr")

if (interactive()) {
    project_path <- "~/Desktop/mJPL-conjecture"
    nobs <- 2000
    psi <- 0.5
    beta_star_setting <- "u2"
}
image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

results_file <- file.path(image_path, paste0("estimates-n-", nobs,
                                         "-beta-", beta_star_setting,
                                         "-psi-", psi,
                                         "-", "estimate",
                                         ".rda"))
load(results_file)

## Extract mJPL estimates, regress estimates vs truth across
## repetitions for each rhosq-kappa-gamma combination, and collect
## intercepts, slopes and residual standard errors
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
            sigma = sapply(fit, "[", "sigma"),
            gamma0 = gamma * sqrt(1 - rhosq))

## Set a cutpoint that determines small/moderate rhosq values
cut_point <- 0.7

## Estimate q() based on the slopes of mJPL estimates vs truth when
## the MLE exists and when it does not
lm_noMLE <- lm(log(a1) ~ log(gamma) + log(kappa) + log(1 - rhosq),
               data = coefs_mJPL |> filter(!mle_exists, rhosq < cut_point))

lm_MLE <- lm(log(a1) ~ log(gamma) + log(kappa) + log(1 - rhosq),
             data = coefs_mJPL |> filter(mle_exists, rhosq < cut_point))

out_file <- file.path(image_path, paste0("conjecture-n-", nobs,
                                         "-beta-", beta_star_setting,
                                         "-psi-", psi, ".rda"))

save(coefs_mJPL, lm_noMLE, lm_MLE, file = out_file)




if (FALSE) {

    library("ggplot2")

    coefs_mJPL <- coefs_mJPL |>
        mutate(rhosq_fac = ifelse(rhosq <= cut_point,
                                  paste("rho^2 <=", cut_point),
                                  paste("rho^2 >", cut_point)),
               `MLE exists` = ifelse(mle_exists, "yes", "no"))

    ## Check that all intercepts are around zero
    ggplot(coefs_mJPL |> filter(method == "mJPL")) +
        geom_boxplot(aes(y = a0)) +
        geom_hline(aes(yintercept = 0), col = "grey", lty = 2) +
        lims(y = c(-0.01, 0.01)) +
        theme_minimal() +
        labs(y = expression(delta[0])) +
        theme(axis.text.x = element_blank())

    ## Slopes when MLE exists
    p_kg <- ggplot(coefs_mJPL |> filter(mle_exists, method == "mJPL")) +
        geom_point(aes(kappa, gamma, color = a1, fill = a1), shape = 24) +
        labs(x = expression(kappa), y = expression(gamma)) +
        lims(x = c(0, 0.6), y = c(0, 20)) +
        scale_colour_continuous(type= "viridis") +
        scale_fill_continuous(type= "viridis") +
        theme_minimal()

    p_gr <- ggplot(coefs_mJPL |> filter(mle_exists, method == "mJPL")) +
        geom_point(aes(gamma, rhosq, color = a1, fill = a1), shape = 24) +
        labs(x = expression(gamma), y = expression(rhosq)) +
        lims(x = c(0, 20), y = c(0, 1)) +
        scale_colour_continuous(type= "viridis") +
        scale_fill_continuous(type= "viridis") +
        theme_minimal()

    p_rk <- ggplot(coefs_mJPL |> filter(mle_exists, method == "mJPL")) +
        geom_point(aes(rhosq, kappa, color = a1, fill = a1), shape = 24) +
        labs(x = expression(rho^2), y = expression(kappa)) +
        lims(x = c(0, 1), y = c(0, 0.6)) +
        scale_color_continuous(type = "viridis") +
        scale_fill_continuous(type= "viridis") +
        theme_minimal()


    ## Slopes when MLE does not exists
    po_int <- coef(lm_noMLE)["(Intercept)"]
    po_k <- coef(lm_noMLE)["log(kappa)"]
    po_g <- coef(lm_noMLE)["log(gamma)"]
    po_r <- coef(lm_noMLE)["log(1 - rhosq)"]

    ggplot(coefs_mJPL |> filter(!mle_exists, method == "mJPL")) +
        geom_point(aes(x = log(kappa^po_k * gamma^po_g * (1 - rhosq)^po_r),
                       y = log(a1), col = rhosq)) +
        geom_abline(aes(intercept = 0, slope = 1), col = "red") +
        facet_grid(~ rhosq_fac, labeller = label_parsed)

    ## Residual error
    ggplot(coefs_mJPL |> filter(!mle_exists, method == "mJPL")) +
        geom_histogram(aes(sigma)) +
        facet_grid( ~ rhosq_fac)


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
        summarize(sest = var(sqrt(n) * (estimate - scaled_truth) * sqrt(tausq), na.rm = TRUE))

    ## Estimate q() based on the slopes of mJPL estimates vs truth when
    ## the MLE exists and when it does not
    mod_var <- lm(log(sest) ~ log(gamma) + log(kappa) + log(1 - rhosq),
                   data = ests |> filter(!mle_exists, rhosq < 0.6))


    ggplot(ests |> filter(!mle_exists, rhosq < 0.6),
           aes(coef(mod_var)["(Intercept)"] +
               log(kappa^coef(mod_var)["log(kappa)"] *
                   gamma^coef(mod_var)["log(gamma)"] *
                   (1 - rhosq)^coef(mod_var)["log(1 - rhosq)"]),
           log(sest))) +
        geom_point(aes(col = kappa)) +
        ## geom_smooth() +
        geom_abline(aes(intercept = 0,  slope = 1))

    ggplot(ests |> filter(mle_exists, rhosq < 0.6),
           aes(seq_along(sest), sest)) +
        geom_point(aes(col = kappa)) +
        geom_smooth() +
        geom_hline(aes(yintercept = mean(sest)))


reparam <- function(ckgr) {
    c(ckgr["(Intercept)"],
      ckgr["log(gamma)"] - 2 * ckgr["log(1 - rhosq)"],
      ckgr["log(kappa)"],
      2 * ckgr["log(1 - rhosq)"])
}

}
