library("parallel")
library("minimaxdesign")

if (interactive()) {
    project_path <- "~/Desktop/mJPL-conjecture"
    ## Number of cores to use for determining whether MLE exists across settings
    ncores <- 10
    ## Number of points to generate on the rhosq-kappa-gamma space
    npoints <- 100
    ## seed
    seed <- 0
    ## suffix for image file
    suffix <- "estimate"
}
image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

set.seed(seed)

## rhosq values
max_rhosq <- 1.0
## Maximum kappa value
max_kappa <- 0.6
## Maximum gamma value
max_gamma <- 20
## Number of samples for the estimation of the phase transition curve
ns <- 200000

## Compute a minimax design for rhosq - kappa - gamma using the minimax
## clustering algorithm in Mak and Joseph (2018)
## https://doi.org/10.1080/10618600.2017.1302881
designs <- miniMaxPro(N = npoints, p = 3)
design <- data.frame(designs$miniMaxPro)
names(design) <- c("rhosq", "kappa", "gamma")
design$rhosq <- max_rhosq * design$rhosq
design$kappa <- max_kappa * design$kappa
design$gamma <- max_gamma * design$gamma

n_settings <- nrow(design)

## Check if MLE exists for every point in the design
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
design$kappa_diff <- mclapply(seq.int(nrow(design)), function(j) {
    rhosq <- design$rhosq[j]
    gamma <- design$gamma[j]
    kappa <- design$kappa[j]
    ## MLE exists if below the phase transition curve or if the kappa
    ## we get for each gamma in the design is larger than the kappa in
    ## the design
    compute_pt(beta0 = sqrt(rhosq) * gamma, gamma, ncores = 1, XZU =  xzu)$kappa - kappa
}, mc.cores = ncores) |> unlist()

design$mle_exists <- design$kappa_diff > 0

experiment_settings <- design

save(experiment_settings,
     file = file.path(image_path, paste0("experiment-settings-", suffix, ".rda")))


## Plots for main text
if (FALSE) {

    library("ggplot2")
    library("patchwork")
    library("dplyr")
    load(file.path(image_path, paste0("experiment-settings-", suffix, ".rda")))

    okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    e_set <- experiment_settings |> mutate(`MLE exists` = ifelse(mle_exists, "yes", "no"))

    p_kg <- ggplot(e_set) +
        geom_point(aes(kappa, gamma, col = `MLE exists`, pch = `MLE exists`), alpha = 0.8) +
        lims(x = c(0, 0.6), y = c(0, 20)) +
        theme_minimal() +
        theme(legend.position = "top") +
        labs(x = expression(kappa), y = expression(gamma)) +
        scale_color_discrete(type = okabe)

    p_gr <- ggplot(e_set) +
        geom_point(aes(gamma, rhosq, col = `MLE exists`, pch = `MLE exists`), alpha = 0.8) +
        lims(x = c(0, 20), y = c(0, 1)) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(y = expression(rho^2), x = expression(gamma)) +
        scale_color_discrete(type = okabe)

    p_rk <- ggplot(e_set) +
        geom_point(aes(rhosq, kappa, col = `MLE exists`, pch = `MLE exists`), alpha = 0.8) +
        lims(x = c(0, 1), y = c(0, 0.6)) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(y = expression(kappa), x = expression(rho^2)) +
        scale_color_discrete(type = okabe)

    ## pdf(file.path(figure_path, "design.pdf"), width = 9, height = 4)
    pdf(file.path("~/Repositories/mJPL-conjecture/figures", "design.pdf"), width = 7, height = 3)
    p_kg + p_gr + p_rk
    dev.off()

}
