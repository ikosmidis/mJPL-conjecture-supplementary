library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")

if (interactive()) {
    conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    image_path <- file.path(project_path, "images")
    code_path <- file.path(project_path, "code")
    figure_path <- file.path(project_path, "figures")
}

source(file.path(code_path, "helper-functions.R"))

## Load conjecture constants
load(file.path(image_path, conjecture_model))
## Get conjecture constants
po_k <- coef(glm_noMLE_rhosq)["log(kappa)"]
po_g <- coef(glm_noMLE_rhosq)["log(gamma)"]
po_r <- coef(glm_noMLE_rhosq)["log(1 - rhosq)"]

load(file.path(image_path, "mJPL_ridge_comparison.rda"))

## Rescale to match experimental setting
res <- res |> mutate(estimate = estimate * sqrt(p),
                     truth = truth * sqrt(p))

## Compute aggregate MSE
mses_mJPL <- res |>
    correct_mJPL_estimates(pk = po_k, pg = po_g, pr = po_r) |>
    group_by(method, gamma, kappa, p, sample) |>
    summarize(mse = mean((estimate - truth)^2),
              bias = mean(estimate - truth),
              mle_exists = unique(mle_exists)) |>
    filter(method == "mJPL")

mses_mJPL <- mses_mJPL |>
    group_by(method, gamma, kappa, p) |>
    summarize(sd = sd(mse) / sqrt(n()),
              mse = mean(mse),
              bias = mean(bias),
              mle_exists = unique(mle_exists),
              R = n())

## Get ridge constats
load(file.path(image_path, "state_evolution_ridge.rda"))
load(file.path(image_path, "state_evolution_ML.rda"))


ridge_mses <- kglmbs |>
    mutate(mse = kappa * sigma^2 / mu^2,
           sd = 0,
           method = "ridge")

ML_mses <- kgmbs |>
    mutate(mse = kappa * sigma^2 / mu^2,
           sd = 0,
           method = "ML")
ML_mses <- rbind(ML_mses |> mutate(lambda = min(kglmbs$lambda)),
                 ML_mses |> mutate(lambda = max(kglmbs$lambda)))

mJPL_mses <- rbind(data.frame(mses_mJPL[c("gamma", "kappa")],
                              lambda = min(kglmbs$lambda),
                              mses_mJPL["mse"],
                              method = "mJPL",
                              mses_mJPL["sd"]),
                   data.frame(mses_mJPL[c("gamma", "kappa")],
                              lambda = max(kglmbs$lambda),
                              mses_mJPL["mse"],
                              method = "mJPL",
                              mses_mJPL["sd"]))

all_mses <- rbind(mJPL_mses,
                  ridge_mses[names(mJPL_mses)],
                  ML_mses[names(mJPL_mses)])


min_mses <- all_mses |>
    group_by(kappa, gamma, method) |>
    summarise(mse = min(mse),
              sd = mean(sd)) |>
    arrange(kappa, gamma, method)

min_mses <- min_mses |>
    pivot_wider(names_from = method, values_from = c(sd, mse)) |>
    mutate(ratio_ridge = mse_mJPL / mse_ridge,
           ratio_ridge_sd = sd_mJPL / mse_ridge,
           ratio_ML = mse_mJPL / mse_ML,
           ratio_ML_sd = sd_mJPL / mse_ML) |>
    data.frame() |>
    left_join(mses_mJPL[c("kappa", "gamma", "mle_exists")], by = c("kappa", "gamma"))

min_mses <- min_mses |> mutate(gamma_lab = factor(paste("gamma ==", gamma),
                                                  levels = paste("gamma ==", unique(gamma)),
                                                  ordered = TRUE))

m_mJPL_ridge <- ggplot(min_mses) +
    geom_col(aes(kappa, ratio_ridge, fill = mle_exists)) +
    geom_point(aes(kappa, ratio_ridge)) +
    geom_segment(aes(x = kappa, y = ratio_ridge - 3 * ratio_ridge_sd,
                     xend = kappa, yend = ratio_ridge + 3 * ratio_ridge_sd)) +
    facet_grid(. ~ gamma_lab, labeller = label_parsed) +
    geom_hline(aes(yintercept = 1)) +
    labs(y = expression(MSE[mJPL] / min[lambda]~MSE[ridge]),
         x = expression(kappa),
         fill = "MLE exists") +
    lims(y = c(0, 3)) +
    scale_x_continuous(breaks = (1:6) / 10) +
    theme_minimal() +
    theme(legend.position = "none", strip.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

m_mJPL_ML <- ggplot(min_mses) +
    geom_col(aes(kappa, ratio_ML, fill = mle_exists)) +
    geom_point(aes(kappa, ratio_ML)) +
    geom_segment(aes(x = kappa, y = ratio_ML - 3 * ratio_ML_sd,
                     xend = kappa, yend = ratio_ML + 3 * ratio_ML_sd)) +
    facet_grid(. ~ gamma_lab, labeller = label_parsed) +
    geom_hline(aes(yintercept = 1)) +
    labs(y = expression(MSE[mJPL] / MSE[ML]),
         x = element_blank(),
         fill = "MLE exists") +
    lims(y = c(0, 3)) +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_blank())

pdf(file.path(project_path, "figures/vs-ridge.pdf"), width = 8, height = 5)
print(m_mJPL_ML / m_mJPL_ridge)
dev.off()

## Range of aggregate bias
range(mses_mJPL$bias)



## ggplot(all_mses |>
##        group_by(kappa, gamma) |>
##        mutate(smse = mse / min(mse, na.rm = TRUE)) |>
##        mutate(kappa_lab = factor(paste("kappa ==", kappa)),
##               gamma_lab = factor(paste("gamma ==", gamma)))) +
##     geom_hline(aes(yintercept = 1), col = "grey") +
##     geom_line(aes(lambda, smse, color = method)) +
##     facet_grid(kappa_lab ~ gamma_lab, labeller = label_parsed) +
##     coord_cartesian(y = c(1, 3), x = c(0, 1)) +
##     labs(x = expression(lambda), y = expression(MSE / minMSE)) +
##     scale_color_manual(values = cols[c(1, 5, 3)]) +
##     theme_minimal() +
##     theme(legend.position = "top")
