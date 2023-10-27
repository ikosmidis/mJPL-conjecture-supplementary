library("parallel")
library("minimaxdesign")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ## Number of cores to use for determining whether MLE exists across settings
    ncores <- 10
    ## Number of points to generate on the rhosq-kappa-gamma space
    npoints <- 200
    ## seed
    seed <- 100
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

save(design,
     file = file.path(image_path, paste0("design-training.rda")))


## Plots for main text
if (FALSE) {

    library("ggplot2")
    library("patchwork")
    library("dplyr")
    load(file.path(image_path, paste0("design-training.rda")))

    okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    e_set <- design |> mutate(`MLE exists` = ifelse(mle_exists, "yes", "no"),
                              alpha = gamma * sqrt(rhosq),
                              gamma0 = gamma * sqrt(1 - rhosq))

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

    p_kr <- ggplot(e_set) +
        geom_point(aes(kappa, rhosq, col = `MLE exists`, pch = `MLE exists`), alpha = 0.8) +
        lims(y = c(0, 1), x = c(0, 0.6)) +
        theme_minimal() +
        theme(legend.position = "none") +
        labs(x = expression(kappa), y = expression(rho^2)) +
        scale_color_discrete(type = okabe)

    p_ag0 <- ggplot(e_set) +
        geom_point(aes(x = alpha, y = gamma0, col = rhosq)) +
        labs(y = expression(gamma[0]), x = expression(beta[0])) +
        scale_color_continuous(type = "viridis", name = expression(rho^2)) +
        theme_bw() +
        theme(legend.position = "top")

    ## pdf(file.path(figure_path, "design.pdf"), width = 9, height = 4)
    pdf(file.path("~/Repositories/mJPL-conjecture/figures", "design.pdf"), width = 7, height = 7/sqrt(2))
    (p_kg | p_ag0) / (p_kr | p_gr)
    dev.off()

}
