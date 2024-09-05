# Installing the `minimaxdesign` R package

At the time of producing the current version of the supplementary material, the `minimaxdesign` R package has been archived on CRAN.  Please run the following code in order to install it.

```R
pkg_url <- "https://cran.r-project.org/src/contrib/Archive/minimaxdesign/minimaxdesign_0.1.5.tar.gz"
install.packages(pkg_url, type = "source", repos = NULL)
```
You may have to install additional dependencies.

# Reproduce simulations, tables and figures.

0. Compute 100 $(\rho^2, \kappa, \gamma$) points for the computer
   experiment and check whether MLE asymptotically exists for each
   one, and produce Figure 1 of the manuscript.

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary/";
  ncores <- 10;
  npoints <- 100;
  seed <- 0;
  source(file.path(project_path, "code/0-design.R"))'
```

1. Training phase

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary/";
  ncores <- 10;
  nobs <- 2000;
  psi <- 0;
  beta_star_setting <- "u2";
  repetitions <- 100;
  seed <- 101;
  source(file.path(project_path, "code/1-train-set.R"))'
```

2. Predict conjectured approximation and produce Figure 2 and Table 1
   of the manuscript.

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  nobs <- 2000;
  psi <- 0;
  beta_star_setting <- "u2";
  source(file.path(project_path, "code/2-conjecture.R"))'
```

3. Produce test sets

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "s1";
  seed <- 103;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "s2";
  seed <- 104;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "u1";
  seed <- 105;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "u2";
  seed <- 106;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "s1";
  seed <- 107;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "s2";
  seed <- 108;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "u1";
  seed <- 109;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "u2";
  seed <- 110;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "s1";
  seed <- 111;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "s2";
  seed <- 112;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "u1";
  seed <- 113;
  source(file.path(project_path, "code/3-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "u2";
  seed <- 114;
  source(file.path(project_path, "code/3-test-set.R"))'

```

4. Compute phase transition curves for the setting where the test sets
   have been produced

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  seed <- 123
  source(file.path(project_path, "code/4-test-set-pt.R"))'
```

5. Reproduce Figure 3 of the manuscript
```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ns <- c(1000, 2000, 3000)
  beta_star <- c("s1", "s2", "u1", "u2")
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  source(file.path(project_path, "code/5-test-summaries.R"))'
```

6. Reproduce Figure 4 and Figure 5 in manuscript
```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  nobs <- c(3000, 1000)
  beta_star_setting <- c("u1", "s1")
  psi <- c(0.3, 0.6)
  rhosq <- c(0.6, 0.1)
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  source(file.path(project_path, "code/6-test-inset-plots.R"))'
```

7. Compute reported summaries on elapsed time in Table 2
```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  nobs <- 2000
  beta_star_setting <- "s2"
  source(file.path(project_path, "code/7-computational-performance.R"))'
```

8. Reproduce Figure 6 of the manuscript
```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  nobs <- 2000
  ncores <- 10
  beta_star_setting <- "s1"
  seed <- 111
  rsqs <- 0.3
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  source(file.path(project_path, "code/8-bernoulli.R"))'
```

9. Carry out the simulation study of Section 6.2 for the estimation of
   the aggregate MSE of the rescaled mJPL estimator
```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  image_path <- file.path(project_path, "images")
  code_path <- file.path(project_path, "code")
  nobs <- 2000
  repetitions <- 50
  ncores <- 10
  ## seed for generating unique seeds
  seed <- 207
  source(file.path(project_path, "code/9-ridge.R"))'
```

10. Reproduce Figure 7 of the manuscript
```bash
Rscript --no-init-file -e \
 'conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  image_path <- file.path(project_path, "images")
  code_path <- file.path(project_path, "code")
  figure_path <- file.path(project_path, "figures")
  source(file.path(project_path, "code/10-ridge-summaries.R"))'
```
