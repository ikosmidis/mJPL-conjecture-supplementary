library("dplyr")
library("ggplot2")

if (interactive()) {
    project_path <- "~/Repositories/mJPL-conjecture-supplementary"
    nobs <- 2000
    beta_star_setting <- "s2"
}
image_path <- file.path(project_path, "images")
code_path <- file.path(project_path, "code")
figure_path <- file.path(project_path, "figures")

source(file.path(code_path, "helper-functions.R"))

results_file <- file.path(image_path, paste0("estimates-n-", nobs,
                                             "-beta-", beta_star_setting,
                                             "-test.rda"))

load(results_file)

mJPL_perf <- performance |> group_by(kappa, gamma, p) |> filter(method == "mJPL") |>
    summarize(mean_elapsed = mean(elapsed / 60),
              min_iter = min(iter),
              mean_iter = mean(iter),
              max_iter = max(iter))

print(range(mJPL_perf$mean_elapsed))
