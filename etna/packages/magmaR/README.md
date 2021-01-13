magmaR is an R client for magma. As such, it provides access to `magma`
functionality directly within R.

## Users
Install with:
```
if (!require("remotes")) install.packages("remotes")
remotes::install_github(
    "mountetna/etna",
    subdir = "packages/magmaR"
```
For further instructions see the rendered version of the vignettes stored at
[TBD location].

## Developers
**RStudio:** It is easiest, and highly recommended due to the many quality of
life tools included, to use RStudio. Here, cloning this repo and then opening
the included .Rproj file is enough to get going. You will then have a `Build`
tab in the upper right of the RStudio window where:

- **Getting package requirements:**
1. Install the package with the below to get everything needed:
```
# R
remotes::install_github(
    "mountetna/etna",
    subdir = "packages/magmaR",
    dependencies = c("Depends", "Imports", "Suggests")
    )
```

- **Working directory** should be the magmaR folder itself, etna/packages/magmaR.
(This is the default for the RStudio project.)
  - Within this directory:
    - Functions (& documentation instructions) are defined in `.R` files of `R/`
    - Unit tests are defined in `.R` files of `tests/testthat/` where `test-*`
    files actually contain testing code and `helper-` files are run as setup.
    - Vignettes are in the `vignettes/` folder.
    - `DESCRIPTION`, base-level file, holds version, dependencies, and more.
    - `NAMESPACE`, & all files in `man/` should never be edited directly. They
    hold package load/unload info & documentation files, respectively. But these
    are managed by Roxygen via Roxygen instructions within `R/*.R` scripts.
    - `clear_stubs_for_true_test.R`, base-level file, is an important helper
    file, with usage described below.
    
- **To fresh install after making changes:**
1. Build source (tar.gz) package: (`Build` tab: More > Build Source Package)
    - Remakes all documentation files, then builds the package in the directory
    one-level up. Runs:
    ```
    # R
    devtools::document(roclets = c('rd', 'collate', 'namespace'))
    devtools::build()
    ```
2. Install and Restart: (`Build` tab: Top-level button)
    - Installs from the last source build, then restarts the R session. Runs:
    ```
    # Terminal
    R CMD INSTALL --no-multiarch --with-keep.source magmaR
    ```

- **To "knit" & view the vignettes:**
  1. Run `usethis::edit_r_environ(scope = "project")` (which is included in the 
`clear_stubs_for_true_test.R` file), and set it as:
```
TOKEN="YOUR-produciton-token"
URL="https://magma.ucsf.edu"
```
  2. ((For a TRUE test, where magma is actually pinged, ensure that you have
proper access to the example project in staging, then actually run the entirety
of `clear_stubs_for_true_test.R`.))
  3. Open the .Rmd file from the `vignettes/` folder & click the `Knit` button in
the window where the file opens. OR run
`rmarkdown::render("vignettes/<which_vignette.Rmd>")`.
  4. Wait for it to process. It will open after.

- **To run unit tests:**
1. Same 1&2 as *To "knit" & view the vignettes* above
2. Same 1&2 as *To "knit" & view the vignettes* above
3. Run them with the automated system: (`Build` tab: More > Test Package)
    - Utilizes the stubs or makes new ones (depending on step 2) and runs
    through the `tests/testthat/` files, starting with `helper-magmaR.R`.
    Actually runs:
    ```
    # R
    devtools::test()
    ```

- **To perform `R CMD Check`, the standard suite of CRAN tests:**
1. Same 1&2 as *To "knit" & view the vignettes* above
2. Same 1&2 as *To "knit" & view the vignettes* above
3. Run with the automated system: (`Build` tab: Top-level button, `Check`)
    - Runs a battery of things, all with short descriptions and results given.
    Includes unit tests, vignette builds, package build and installation, plus
    many more.  Actually runs:
    ```
    # R
    devtools::check()
    ```

# BEFORE any commit:
- Do not commit:
  - Rendered vignettes (`vignettes/<vignette_name>.html`).
  - Updated stubs / fixtures, *unless needed*. Your stubs / fixtures, the .yml
  files in `test/fixtures/` may have changed due to deleting the old ones in
  order to run a full package test. But, please do not commit them unless a
  change has been made which requires a particular stubs to be updated for some
  altered call/behavior. (Clutter reduction)
    
