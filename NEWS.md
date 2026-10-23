# sts 1.6

## Bug fix

* `initialization = "stm"` (including `stmSeed`): the sentiment-discourse
  coefficients `kappa_s` were always initialized to zero, because the STM
  topic coefficients were read where the content-covariate interaction
  coefficients were intended. `kappa_t` and `kappa_s` are now initialized
  from the STM's content-covariate effects (intercept and slope of a
  document-weighted least-squares fit on the centered levels of
  `initializationVar`), which reproduces the STM topic-word distributions
  exactly for a binary `initializationVar`. On the Gadarian example this
  gives a higher final ELBO and faster convergence; fits with
  `initialization = "stm"` will differ from sts 1.5. An informative error is
  now given if `stmSeed` has no content covariate matching
  `initializationVar`.

## Removal of the stm dependency

* The package no longer depends on the stm package, which is scheduled for
  archival on CRAN. Apart from the bug fix above, results are numerically
  identical to sts 1.5 for the same random seed.
  - The spectral (anchor word) initialization used by
    `initialization = "anchor"` (the default) and the STM fit used by
    `initialization = "stm"` are now implemented internally (`R/stmPort.R`,
    `src/stmfuns.cpp`), adapted from stm 1.3.8 (MIT license; the stm authors
    are credited in `Authors@R` and `inst/COPYRIGHTS`).
  - New exported functions `textProcessor()`, `prepDocuments()` and
    `make.heldout()`, ported from stm, so that raw text can be prepared for
    `sts()` and held-out data created without stm.
  - New data set `gadarian` (from stm) used in the examples.
  - A prefit STM model can still be supplied through `stmSeed`.
* `sts()` now gives an informative error for an unknown `initialization`
  value.
* Updated the reference to Chen and Mankad (2025), Management Science 71(7).
* The vignette is now registered through the `R.rsp::asis` engine.
* Added SnowballC (used for stemming in `textProcessor()`) to Imports and
  removed stm, knitr and rmarkdown.

# sts 1.5

* Previous CRAN release.
