# MixtureRegression

A Bayesian inference framework for two-component Gaussian **mixture regression** models, implemented in R. The sampler uses **Metropolis-Hastings within Gibbs** steps combined with **Parallel Tempering** to efficiently explore the posterior distribution of all model parameters.

---

## Overview

Standard linear regression assumes a single Gaussian noise term. Mixture regression relaxes this by allowing each observation to arise from one of *K* latent Gaussian components, each with its own mean and variance. This is useful when the residual distribution is multimodal or when the data come from heterogeneous sub-populations.

This project focuses on the **K = 2** case and introduces a *factored* parameterisation that keeps all parameters identifiable and naturally respects their constraints.

---

## Model

Given covariates **X** (n × k) and a response **Y** (n × 1), the model is:

```
Y_i ~ p · N(μ₁ᵢ, σ₁ᵢ²)  +  (1−p) · N(μ₂ᵢ, σ₂ᵢ²)
```

where the component means and standard deviations are derived from a *factored* parameterisation:

| Parameter | Description | Constraint |
|-----------|-------------|------------|
| **β** (k×1) | Regression coefficients (shared linear predictor μ = Xβ) | — |
| **p** ∈ (0,1) | Mixture weight for component 1 | — |
| **σ** > 0 | Overall standard deviation scale | — |
| **(η₁, η₂, φ)** | Variance-decomposition factors | η₁² + η₂² + φ² = 1 |

From these parameters the per-component quantities are:

```
μ₁ᵢ = μᵢ − σ · φ · √((1−p)/p)
μ₂ᵢ = μᵢ + σ · φ · √(p/(1−p))

σ₁  = σ · η₁ / √p
σ₂  = σ · η₂ / √(1−p)
```

### Priors

| Parameter | Prior |
|-----------|-------|
| p | Beta(τ₁, τ₂) |
| (η₁², η₂², φ²) | Dirichlet(α) |
| σ² | Inverse-Gamma(shape, rate) |
| β | Improper (log-Student-t penalty) |

---

## Algorithm

The sampler is implemented in `FactoredMixReg.R` via the `MixReg_factor()` function.

1. **Metropolis-Hastings updates** for each parameter block in sequence:
   - Mixture weight **p** — Beta proposal
   - Variance-decomposition vector **(η₁², η₂², φ²)** — Dirichlet proposal
   - Standard deviation **σ²** — Inverse-Gamma proposal
   - Regression coefficients **β** — Multivariate Normal (random-walk) proposal

2. **Parallel Tempering** — `N_chains` chains are run simultaneously at a common temperature (or a custom schedule). Every `swap_config$cycle_period` iterations, pairs of chains attempt to exchange their states using a Metropolis acceptance criterion.

3. **Adaptive proposal tuning** — the step-size ε for each parameter is adjusted periodically to keep acceptance rates in the target range [0.30, 0.40].

4. **Gelman-Rubin convergence diagnostics** — the multivariate potential scale reduction factor (PSRF) is computed at regular intervals using the `coda` package and stored in the returned `gelman` list.

---

## Repository Structure

```
MixtureRegression/
├── Model.R            # Log-likelihood and log-posterior functions
├── MH_Samplers.R      # Metropolis-Hastings samplers for each parameter block
│                      #   + chain-swap and state-exchange logic
├── FactoredMixReg.R   # Main MCMC driver: MixReg_factor()
├── utilities.R        # Adaptive tuning, convergence diagnostics,
│                      #   progress bars, and ggplot helpers
├── rnormix.R          # Simulate samples from a two-component Gaussian mixture
├── generate_test.R    # Generate synthetic regression data for testing
├── experiment.Rmd     # R Markdown notebook: end-to-end experiment & diagnostics
└── MixRegTuner.Rmd    # R Markdown notebook: hyperparameter tuning workflow
```

---

## Usage

### 1. Generate synthetic data

```r
source("generate_test.R")

data <- generate_test(n_obs = 3000, k_var = 4, add_intercept = TRUE)
X <- data$X   # Design matrix (n × 5 with intercept)
Y <- data$Y   # Response vector
```

### 2. Set hyperparameters

```r
Tau      <- list(tau1 = 0.5, tau2 = 0.5)   # Beta prior on p
Alpha    <- c(0.5, 0.5, 0.5)               # Dirichlet prior on (η₁², η₂², φ²)
ksi      <- list(shape = 1.5, rate = 2)    # Inverse-Gamma prior on σ²
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% Y  # OLS initialisation for β
```

### 3. Run the sampler

```r
source("FactoredMixReg.R")

T        <- 10000   # Number of MCMC iterations
N_chains <- 10      # Number of parallel chains

set.seed(1999)
L <- MixReg_factor(T, X, Y, Tau, Alpha, ksi, beta_hat,
                   core_temperature = 1,
                   N_chains = N_chains)
```

### 4. Inspect results

The returned list `L` contains:

| Field | Description |
|-------|-------------|
| `L$P` | N_chains × T matrix of sampled mixture weights |
| `L$etaPhi` | N_chains × T × 3 array of (η₁, η₂, φ) samples |
| `L$sigma_sqrd` | N_chains × T matrix of sampled σ² values |
| `L$Beta_chains` | N_chains × T × k array of sampled β values |
| `L$ac_rate_*` | Per-chain cumulative acceptance rates for each block |
| `L$gelman` | Gelman-Rubin diagnostics over iterations |

```r
source("utilities.R")
library(gridExtra)

# Trace plots and histograms for mixture weight p
p1 <- ggplot_overlapping_histograms(L$P, title = "Mixture weight p")
p2 <- ggplot_overlapping_lines(L$P, title = "Trace: p")
grid.arrange(p1, p2, ncol = 2)
```

Full worked examples are available in `experiment.Rmd` and `MixRegTuner.Rmd`.

---

## Dependencies

Install the required R packages before running:

```r
install.packages(c("MASS", "coda", "progress", "invgamma",
                   "mvtnorm", "DirichletReg", "ggplot2",
                   "tidyr", "dplyr", "gridExtra"))
```

---

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file included in the repository.