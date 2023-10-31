library("parallel")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    ncores <- 10
    seed <- 123
}
image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

rhosq_grid <- seq.int(0, 9, 1) / 10

ns <- 200000
set.seed(seed)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))


pt_curve <- as.list(rep(0, length(rhosq_grid)))
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
ga <- c(0.001, 0.01, seq(0, 20, length = 200))
for (j in seq_along(rhosq_grid)) {
    rhosq <- rhosq_grid[j]
    cat("ρ^2", rhosq, "\n")
    pt <- mclapply(1:length(ga), function(k) {
        d <- compute_pt(beta0 = sqrt(rhosq) * ga[k], ga[k], ncores = 1, XZU =  xzu)
        d$rhosq <- rhosq
        d
    }, mc.cores = ncores)
    pt_curve[[j]] <- do.call("rbind", pt)
}

save(pt_curve, file = file.path(image_path, "pt-curves-test.rda"))
