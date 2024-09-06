# Supplementary material code for “Jeffreys-prior penalty for high-dimensional logistic regression: A conjecture about aggregate bias”
Ioannis Kosmidis, Patrick Zietkiewicz
September 6, 2024

# Directory structure

The directory `code/` contains the scripts

<table style="width:46%;">
<colgroup>
<col style="width: 45%" />
</colgroup>
<thead>
<tr class="header">
<th>script</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>00-design.R</td>
</tr>
<tr class="even">
<td>01-train-set.R</td>
</tr>
<tr class="odd">
<td>02-conjecture.R</td>
</tr>
<tr class="even">
<td>03-test-set.R</td>
</tr>
<tr class="odd">
<td>04-test-set-pt.R</td>
</tr>
<tr class="even">
<td>05-test-summaries.R</td>
</tr>
<tr class="odd">
<td>06-test-inset-plots.R</td>
</tr>
<tr class="even">
<td>07-computational-performance.R</td>
</tr>
<tr class="odd">
<td>08-bernoulli.R</td>
</tr>
<tr class="even">
<td>09-ridge.R</td>
</tr>
<tr class="odd">
<td>10-ridge-summaries.R</td>
</tr>
<tr class="even">
<td>helper-functions.R</td>
</tr>
</tbody>
</table>

that reproduce all the numerical results and figures in the manuscript

> Kosmidis I, Zietkiewicz P (2023). Jeffreys-prior penalty for
> high-dimensional logistic regression: A conjecture about aggregate
> bias. https://arxiv.org/abs/2311.11290

The directory `images/` is populated by R image files that store the
numerical results the scripts produce, and two images with supporting
results for producing Figure 7 (`state_evolution_ridge.rda`,
`state_evolution_ML.rda`). The directory `figures/` is populated by
graphics that the scripts produce.

# R version and packages

All results are reproducible using R version 4.4.1 (2024-06-14) and the
contributed packages

<table style="width:43%;">
<colgroup>
<col style="width: 26%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th>package</th>
<th>version</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>brglm2</td>
<td>0.9.2</td>
</tr>
<tr class="even">
<td>colorspace</td>
<td>2.1-1</td>
</tr>
<tr class="odd">
<td>detectseparation</td>
<td>0.3</td>
</tr>
<tr class="even">
<td>dplyr</td>
<td>1.1.4</td>
</tr>
<tr class="odd">
<td>ggplot2</td>
<td>3.5.1</td>
</tr>
<tr class="even">
<td>ggpp</td>
<td>0.5.8-1</td>
</tr>
<tr class="odd">
<td>memisc</td>
<td>0.99.31.7</td>
</tr>
<tr class="even">
<td>minimaxdesign</td>
<td>0.1.5</td>
</tr>
<tr class="odd">
<td>patchwork</td>
<td>1.2.0</td>
</tr>
<tr class="even">
<td>tidyr</td>
<td>1.3.1</td>
</tr>
</tbody>
</table>

At the time of producing the current version of the supplementary
material, the `minimaxdesign` R package has been archived on CRAN.
Please run the following code in order to install it.

``` r
pkg_url <- "https://cran.r-project.org/src/contrib/Archive/minimaxdesign/minimaxdesign_0.1.5.tar.gz"
install.packages(pkg_url, type = "source", repos = NULL)
```

You may have to install additional dependencies.

# Reproducing the results

The results can be reproduced in a Unix-based shell environment by
running the following commands in the order they appear. Some of scripts
used below rely on parallel computing, which is implemented through the
`parallel` R package. Those scripts will not work on Windows unless
`ncores <- 1` (which will lead in long compute times and is not
recommended) or they are modified to use a different parallel back-end.
All steps should work directly in macOS and Linux systems.

1.  Compute 100 (*ρ*<sup>2</sup>, *κ*, *γ*) points for the computer
    experiment and check whether MLE asymptotically exists for each one,
    and produce Figure 1 of the manuscript.

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary/";
  ncores <- 10;
  npoints <- 100;
  seed <- 0;
  source(file.path(project_path, "code/00-design.R"))'
```

Outputs: `images/design-training.rda`, `figures/design.pdf` (Figure 1)

1.  Training phase

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary/";
  ncores <- 10;
  nobs <- 2000;
  psi <- 0;
  beta_star_setting <- "u2";
  repetitions <- 100;
  seed <- 101;
  source(file.path(project_path, "code/01-train-set.R"))'
```

Inputs: `images/design-training.rda`

Outputs: `images/estimates-n-2000-beta-u2-psi-0-training.rda`

1.  Predict conjectured approximation and produce Figure 2 and Table 1
    of the manuscript.

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  nobs <- 2000;
  psi <- 0;
  beta_star_setting <- "u2";
  source(file.path(project_path, "code/02-conjecture.R"))'
```

Inputs: `images/design-training.rda`,
`images/estimates-n-2000-beta-u2-psi-0-training.rda`

Outputs: `images/conjecture-n-2000-beta-u2-psi-0.rda`,
`figures/conjecture.pdf` (Figure 2), Table 1

1.  Produce test sets

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "s1";
  seed <- 103;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "s2";
  seed <- 104;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "u1";
  seed <- 105;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 1000;
  beta_star_setting <- "u2";
  seed <- 106;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "s1";
  seed <- 107;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "s2";
  seed <- 108;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "u1";
  seed <- 109;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 2000;
  beta_star_setting <- "u2";
  seed <- 110;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "s1";
  seed <- 111;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "s2";
  seed <- 112;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "u1";
  seed <- 113;
  source(file.path(project_path, "code/03-test-set.R"))'

Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  nobs <- 3000;
  beta_star_setting <- "u2";
  seed <- 114;
  source(file.path(project_path, "code/03-test-set.R"))'
```

Outputs: `images/estimates-n-1000-beta-s1-test.rda`,
`images/estimates-n-1000-beta-s2-test.rda`,
`images/estimates-n-1000-beta-u1-test.rda`,
`images/estimates-n-1000-beta-u2-test.rda`,
`images/estimates-n-2000-beta-s1-test.rda`,
`images/estimates-n-2000-beta-s2-test.rda`,
`images/estimates-n-2000-beta-u1-test.rda`,
`images/estimates-n-2000-beta-u2-test.rda`,
`images/estimates-n-3000-beta-s1-test.rda`,
`images/estimates-n-3000-beta-s2-test.rda`,
`images/estimates-n-3000-beta-u1-test.rda`,
`images/estimates-n-3000-beta-u2-test.rda`

1.  Compute phase transition curves for the setting where the test sets
    have been produced

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ncores <- 10;
  seed <- 123
  source(file.path(project_path, "code/04-test-set-pt.R"))'
```

Outputs: `images/pt-curves-test.rda`

1.  Reproduce Figure 3 of the manuscript

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary";
  ns <- c(1000, 2000, 3000)
  beta_star <- c("s1", "s2", "u1", "u2")
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  source(file.path(project_path, "code/05-test-summaries.R"))'
```

Inputs: `images/conjecture-n-2000-beta-u2-psi-0.rda`,
`images/estimates-n-1000-beta-s1-test.rda`,
`images/estimates-n-1000-beta-s2-test.rda`,
`images/estimates-n-1000-beta-u1-test.rda`,
`images/estimates-n-1000-beta-u2-test.rda`,
`images/estimates-n-2000-beta-s1-test.rda`,
`images/estimates-n-2000-beta-s2-test.rda`,
`images/estimates-n-2000-beta-u1-test.rda`,
`images/estimates-n-2000-beta-u2-test.rda`,
`images/estimates-n-3000-beta-s1-test.rda`,
`images/estimates-n-3000-beta-s2-test.rda`,
`images/estimates-n-3000-beta-u1-test.rda`,
`images/estimates-n-3000-beta-u2-test.rda`

Outputs: `figures/test.pdf` (Figure 3)

1.  Reproduce Figure 4 and Figure 5 in manuscript

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  nobs <- c(3000, 1000)
  beta_star_setting <- c("u1", "s1")
  psi <- c(0.3, 0.6)
  rhosq <- c(0.6, 0.1)
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  source(file.path(project_path, "code/06-test-inset-plots.R"))'
```

Inputs: `images/conjecture-n-2000-beta-u2-psi-0.rda`,
`images/pt-curves-test.rda`, `images/estimates-n-1000-beta-s1-test.rda`,
`images/estimates-n-3000-beta-u1-test.rda`

Outputs: `figures/illustration1.pdf` (Figure 4),
`figures/illustration2.pdf` (Figure 5)

1.  Compute reported summaries on elapsed time in Table 2

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  nobs <- 2000
  beta_star_setting <- "s2"
  source(file.path(project_path, "code/07-computational-performance.R"))'
```

Inputs: `images/estimates-n-2000-beta-s2-test.rda`

Outputs: Table 2

1.  Reproduce Figure 6 of the manuscript

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  nobs <- 2000
  ncores <- 10
  beta_star_setting <- "s1"
  seed <- 111
  rsqs <- 0.3
  conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  source(file.path(project_path, "code/08-bernoulli.R"))'
```

Inputs: `images/conjecture-n-2000-beta-u2-psi-0.rda`

Outputs: `images/bernoulli.rda`, `figures/bernoulli.pdf` (Figure 6)

1.  Carry out the simulation study of Section 6.2 for the estimation of
    the aggregate MSE of the rescaled mJPL estimator

``` bash
Rscript --no-init-file -e \
 'project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  image_path <- file.path(project_path, "images")
  code_path <- file.path(project_path, "code")
  nobs <- 2000
  repetitions <- 50
  ncores <- 10
  ## seed for generating unique seeds
  seed <- 207
  source(file.path(project_path, "code/09-ridge.R"))'
```

Outputs: `images/mJPL_ridge_comparison.rda`

1.  Reproduce Figure 7 of the manuscript

``` bash
Rscript --no-init-file -e \
 'conjecture_model <- "conjecture-n-2000-beta-u2-psi-0.rda"
  project_path <- "~/Repositories/mJPL-conjecture-supplementary"
  image_path <- file.path(project_path, "images")
  code_path <- file.path(project_path, "code")
  figure_path <- file.path(project_path, "figures")
  source(file.path(project_path, "code/10-ridge-summaries.R"))'
```

Inputs: `images/conjecture-n-2000-beta-u2-psi-0.rda`,
`images/mJPL_ridge_comparison.rda`, `images/state_evolution_ridge.rda`,
`images/state_evolution_ML.rda`

Outputs: `figures/vs-ridge.pdf` (Figure 7)
