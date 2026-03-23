# SIHR 2.1.1

## Bug fixes

- **CVXR 1.8 compatibility**: Updated to work with CVXR 1.8.x which introduced breaking API changes. The package now supports both CVXR 1.0.x and 1.8.x:
  - Replaced deprecated `getValue()` with new `value()` API
  - Replaced `result$status` with `status(prob)` where applicable
  - Treat `optimal_inaccurate` as successful (CVXR 1.8 may return this instead of `optimal`)

- **Robust error handling in `compute_direction()`**: Fixed vignette build failure (`object 'direction' not found`) when CVXR throws errors (e.g., on CRAN check servers). The function now reliably falls back to the diagonal approximation when CVXR fails, using `is.null(direction)` instead of `<<-` assignment which could fail with rlang-based errors.

- **Numerical stability**: Added safeguard against near-zero diagonal elements in the fallback Sigma.hat when using the diagonal approximation.
