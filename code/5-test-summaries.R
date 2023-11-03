library("ggplot2")
library("patchwork")
library("dplyr")

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ns <- c(1000, 2000, 3000)
    beta_star <- c("s1", "s2", "u1", "u2")
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
}

image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

## Get conjecture constants
load(file.path(image_path, conjecture_model))
## Get conjecture constants
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]

## Extract test estimates
all_estimates <- NULL
for (bs in beta_star) {
    for (n in ns) {
        results_path <- file.path(image_path, paste0("estimates-n-", n,
                                                     "-beta-", bs,
                                                     "-test",
                                                     ".rda"))
        if (file.exists(results_path))
            load(results_path)
        else
            next
        cat(results_path,"\n")
        estimates$beta_star <- bs
        all_estimates <- rbind(all_estimates, estimates)
    }
}


## Summary plots
coefs <- all_estimates |>
    filter(parameter > 0) |>
    group_by(kappa, gamma, rhosq, psi, method, mle_exists, beta_star, n) |>
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
           n_fac = factor(n),
           method = factor(method, levels = c("ML", "mJPL"), ordered = TRUE))

Rsq <- function(log_slopes, kappa, gamma, rhosq, pk, pg, pr) {
    m_la1 <- mean(log_slopes, na.rm = TRUE)
    tss <- sum((log_slopes - m_la1)^2, na.rm = TRUE)
    resid <- log_slopes - pk * log(kappa) - pg * log(gamma) - pr * log(1 - rhosq)
    rss <- sum(resid^2, na.rm = TRUE)
    1 - rss / tss
}

rsqs <- coefs |>
    filter(method == "mJPL", !mle_exists) |>
    group_by(n_fac, beta_star, psi_fac, rhosq) |>
    summarize(rsq = Rsq(log(a1), kappa, gamma, rhosq, po_k, po_g, po_r))
    ## summarize(rsq = Rsq(a1, kappa, gamma, rhosq, -1, -1, 0.5))


p_rsq_ne <- ggplot(rsqs, aes(rhosq, rsq)) +
    geom_line(aes(col = n_fac), alpha = 0.5) +
    ## geom_point(aes(col = n_fac), alpha = 0.5) +
    geom_hline(aes(yintercept = 1), col = "darkgrey", lty = 1) +
    facet_grid(psi_fac ~ beta_star, labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "left", axis.text.x = element_text(size = 6)) +
    scale_color_manual(values = okabe, name = "n") +
    coord_cartesian(y = c(0.7, 1)) +
    labs(x = expression(rho^2), y = expression(R["test"]^2), title = "MLE does not exist")

coefs_restr <- coefs |> filter(mle_exists, method == "mJPL")

p_bp_e <- ggplot(coefs_restr) +
    geom_jitter(aes(y = n_fac, x = a1, col = n_fac), position = "identity",
                alpha = 0.1) +
    geom_point(data = coefs_restr |> group_by(n_fac, beta_star, psi_fac) |> summarize(a1 = mean(a1)),
               aes(y = n_fac, x = a1), shape = 1) +
    geom_vline(aes(xintercept = 1), col = "grey", lty = 1) +
    facet_grid(psi_fac ~ beta_star, labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
    scale_color_manual(values = okabe, name = "n") +
    labs(y = "n", x = expression(delta[1]^"*"), title = "MLE exists")

p_bp_int <- ggplot(coefs |> filter(method == "mJPL")) +
    geom_jitter(aes(y = n_fac, x = a0, col = n_fac), position = "identity",
                alpha = 0.1) +
    geom_point(data = coefs |> filter(method == "mJPL") |> group_by(n_fac, beta_star, psi_fac) |> summarize(a0 = mean(a0)),
               aes(y = n_fac, x = a0), shape = 1) +
    geom_vline(aes(xintercept = 0), col = "grey", lty = 1) +
    facet_grid(psi_fac ~ beta_star, labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
    scale_color_manual(values = okabe, name = "n") +
    labs(y = "n", x = expression(delta[0]^"*"))

pdf(file.path(project_path, "figures/test.pdf"), width = 8.5, height = 8.5/sqrt(2))
## pdf(file.path("~/Repositories/mJPL-conjecture/figures", "test.pdf"), width = 9, height = 9/sqrt(2))
layout <- "
AABB
CDD#
"
print((p_rsq_ne + p_bp_e + guide_area() + p_bp_int) +
    plot_layout(design = layout, guides = "collect"))
dev.off()

