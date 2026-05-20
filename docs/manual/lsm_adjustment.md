# Least-Squares Adjustment

!!! note
    This section is a stub. Content will be added in a future release.

## Adjustment Methods

- **LSMDiff** — differential observations (gravity differences between station pairs)
- **LSMNonDiff** — non-differential (absolute setup observations)
- **VGLSM** — vertical gravity gradient estimation
- **BEV Legacy** — multiple linear regression replicating BEV's legacy processing scheme

All methods estimate drift polynomials per survey alongside gravity values at stations.
Results include post-fit residuals and statistical tests (global model test, tau test).
