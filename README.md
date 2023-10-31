R | 4.3.1
minimaxdesign | 0.1.5
ggplot2 | 3.4.4

0. Compute 100 $(\rho^2, \kappa, \gamma$) points for the computer
   experiment and check whether MLE asymptotically exists for each one

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary/";
  ncores <- 10;
  npoints <- 100;
  seed <- 0;
  source(file.path(project_path, "code/0-settings.R"))'
```

1. Simulate 100 data sets per design point and get mJPL and ML
   estimates (whenever the latter exist asymptotically)

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary/";
  ncores <- 10;
  nobs <- 2000;
  psi <- 0.5;
  beta_star_setting <- "u2";
  repetitions <- 100;
  seed <- 101;
  source(file.path(project_path, "code/1-estimates.R"))'
```

2. Estimate conjectured relationship

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Desktop/mJPL-conjecture-supplementary";
  nobs <- 2000;
  psi <- 0.5;
  beta_star_setting <- "u2";
  source(file.path(project_path, "code/2-conjecture.R"))'
```

3. Validate

```bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "s1";
  seed <- 103;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "s2";
  seed <- 104;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "u1";
  seed <- 105;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "u2";
  seed <- 106;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "s1";
  seed <- 107;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "s2";
  seed <- 108;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "u1";
  seed <- 109;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "u2";
  seed <- 110;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "s1";
  seed <- 111;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "s2";
  seed <- 112;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "u1";
  seed <- 113;
  source(file.path(project_path, "code/3-validation-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "u2";
  seed <- 114;
  source(file.path(project_path, "code/3-validation-set.R"))'

```







