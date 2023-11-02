R | 4.3.1
minimaxdesign | 0.1.5
ggplot2 | 3.4.4

Commit 517815a is with fixed psi = 0.5 in training

0. Compute 100 $(\rho^2, \kappa, \gamma$) points for the computer
   experiment and check whether MLE asymptotically exists for each one

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
  psi <- 0.5;
  beta_star_setting <- "u2";
  repetitions <- 100;
  seed <- 101;
  source(file.path(project_path, "code/1-train-set.R"))'
```

2. Predict conjectured approximation

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories//mJPL-conjecture-supplementary";
  nobs <- 2000;
  psi <- 0.5;
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

5. Reproduce Figure X1 in manuscript
```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ns <- c(1000, 2000, 3000)
  beta_star <- c("s1", "s2", "u1", "u2")
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.5.rda"
  source(file.path(project_path, "code/5-test-summaries.R"))'
```




