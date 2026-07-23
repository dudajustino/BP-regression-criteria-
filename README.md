# BP Regression Criteria

This repository contains the R codes used to reproduce the simulation study presented in the article

> Justino, M. E. C., Pereira, T. L., & Souza, T. C. (2026). *Model selection criteria in beta prime regression*. Communications in Statistics – Simulation and Computation. https://doi.org/10.1080/03610918.2026.2695404

## Authors

- Maria Eduarda da Cruz Justino
- Tarciana Liberal Pereira
- Tatiene Correia de Souza

## Programming language

- R

## Repository structure

### `scripts/`

These scripts reproduce the simulation scenarios presented in the paper.

| File | Description |
|------|-------------|
| `fixed_precision.R` | Simulation study under the correctly specified beta prime regression model with fixed precision. |
| `variable_precision.R` | Simulation study under the correctly specified beta prime regression model with variable precision. |
| `error_link_mu_precision.R` | Simulation study under misspecification of the mean link function and the precision submodel. |
| `error_omission_fixed_precision.R` | Simulation study under omission of relevant covariates with fixed precision. |
| `error_omission_variable_precision.R` | Simulation study under omission of relevant covariates with variable precision. |

### `source codes/`

Supporting functions used by the simulation scripts.

| File | Description |
|------|-------------|
| `gamlss_BP.R` | Beta prime family for the `gamlss` framework. |
| `criteria.R` | Functions implementing the model selection criteria studied in the paper. |
| `Residual_H_Log_Like_BP.R` | Functions for residuals, hat matrix, log-likelihood and diagnostic quantities. |

### Root directory

| File | Description |
|------|-------------|
| `README.md` | Repository description and usage information. |
| `sessionInfo.txt` | R session information used to generate the simulation results, including package versions and RNG settings. |

## Requirements

The simulations require R and the packages used in the scripts, including

- `gamlss`
- `extraDistr`

## Citation

If you use this code, please cite

> Justino, M. E. C., Pereira, T. L., & Souza, T. C. (2026). *Model selection criteria in beta prime regression*. Communications in Statistics – Simulation and Computation. https://doi.org/10.1080/03610918.2026.2695404
