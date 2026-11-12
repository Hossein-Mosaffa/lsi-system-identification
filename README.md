# System Identification with Least Squares Iterative (LSI) Algorithms

MATLAB implementations and a project report for three papers that apply the **Least Squares Iterative (LSI)** idea to system identification: a multivariable CARMA system, Box-Jenkins models, and Wiener nonlinear systems.

> **Course project** · K. N. Toosi University of Technology, Faculty of Electrical Engineering

> **Instructor:** Dr. Aliyari

> **Authors:** Hossein Mosaffa, Ali Sheikholeslami

> **Date:** February 2020

**Keywords:** LSI, CARMA, Box-Jenkins, Wiener model, parameter estimation

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements and Installation](#requirements-and-installation)
- [Usage](#usage)
- [Method Summary](#method-summary)
- [Results](#results)
- [References](#references)

---

## Overview

Ordinary least squares (LS) fails when the system's internal variables are unavailable, or when the regressor (information) vector contains noise terms that cannot be measured. LSI handles this by replacing the unknown terms with their estimates and refining them over several iterations until the parameter estimate converges within a tolerance.

| Paper | System | Reference |
|---|---|---|
| 1 | Multivariable controlled ARMA (CARMA) | Bao et al. [1] |
| 2 | Box-Jenkins (BJ), recursive vs. iterative | Liu, Wang and Ding [2] |
| 3 | Wiener nonlinear model | Zhou, Li and Pan [3] |

## Repository Structure

```text
.
├── codes/
│   ├── E1/
│   │   └── CARMA_lsi.m                 # Paper 1: multivariable CARMA
│   ├── E2/
│   │   ├── BOXJENKINS_ITERATIVE.m      # Paper 2: iterative BJ (LSI)
│   │   └── BOXJENKINS_RECURSIVE.m      # Paper 2: recursive BJ (BJ-RGELS)
│   └── E3/
│       ├── lsi1.m                      # Paper 3: Wiener model (LSI)
│       ├── lsi2.m
│       └── lsi3.m
└── System_Identification_Report_EN.docx
```

## Requirements and Installation

- MATLAB. The scripts are plain `.m` files, so recent releases should work.
- No other dependencies are documented.

```bash
git clone https://github.com/Hossein-Mosaffa/lsi-system-identification.git
cd lsi-system-identification
```

Then open MATLAB and add the code folders to the path:

```matlab
addpath(genpath('codes'));
```

For the DaISy experiment (heating system and heat exchanger), download the datasets from the [DaISy database](https://homes.esat.kuleuven.be/~smc/daisy/) separately.

## Usage

Run each script from MATLAB. The paths below follow the repository structure above.

```matlab
% Paper 1: multivariable CARMA, L = 1000, 2000, 3000
run('codes/E1/CARMA_lsi.m')

% Paper 2: Box-Jenkins
run('codes/E2/BOXJENKINS_RECURSIVE.m')
run('codes/E2/BOXJENKINS_ITERATIVE.m')

% Paper 3: Wiener model
run('codes/E3/lsi1.m')
```

Settings reported in the project:

| Setting | Value |
|---|---|
| Data length `L` | 1000, 2000, 3000 (Paper 1); 1000, 2000 (Paper 2) |
| Noise variance | 0.25 (Papers 1 and 2), also 1 (Paper 2) |
| Noise std (Paper 3) | 0.1 and 0.2 |
| Iterations (Paper 2 iterative) | fixed at 15 |
| Stop rule | `‖θ_k − θ_(k−1)‖` below a tolerance |

The plots (parameter estimates per iteration, error curves, flowcharts) appear in the original Persian report and are not in the English version. Regenerate them by running the scripts.

## Method Summary

### Paper 1: CARMA

Model: `A(z) y(t) = B(z) u(t) + D(z) v(t)`, with the form `y(t) = θᵀ φ(t) + v(t)`.

1. Collect input/output data for `t = 1..L`.
2. Set `k = 1` and initialize noise estimates `v₀(t)` randomly.
3. Build the information vectors `φ_k(t)`, using `v_{k-1}` in place of the unknown noise.
4. Update `θ_k = (Φ_k Φ_kᵀ)⁻¹ Φ_k Yᵀ`.
5. Update the noise estimate `v_k(t) = y(t) − θ_kᵀ φ_k(t)`.
6. Stop when `‖θ_k − θ_(k−1)‖` is below the tolerance. Otherwise increase `k` and repeat.

### Paper 2: Box-Jenkins

Model: `y(t) = [B(z)/A(z)] u(t) + [D(z)/C(z)] v(t)`, with `na = nb = 2` and `nc = nd = 1`.

- **Recursive (BJ-RGELS):** online, one sample at a time, with covariance and parameter updates.
- **Iterative (LSI):** offline batch over a finite data length. It re-estimates the noise-free output, coloured noise and white noise at each pass.

### Paper 3: Wiener model

A linear dynamic block followed by a static nonlinearity. At each iteration the regressor uses the previous estimate of the unmeasurable `x(t)`, then `θ` is solved by least squares.

If `x` starts at zero or at a constant, the matrix to invert becomes singular. Two fixes are used: random initial conditions, or the regularized form `θ = (Φ_kᵀ Φ_k + μ_k I)⁻¹ Φ_kᵀ Y`.

## Results

### Paper 1: CARMA (noise variance 0.25)

The algorithm reached the optimum in 9 iterations for every `L`. Estimation error (%) falls sharply with iterations:

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| L = 1000 | 29.11 | 9.12 | 4.77 | 3.16 | 2.49 | 2.21 | 2.13 | 2.08 | 2.07 |
| L = 2000 | 30.08 | 7.64 | 2.69 | 1.35 | 1.43 | 1.29 | 1.31 | 1.37 | 1.37 |

### Paper 2: Box-Jenkins (recursive, final estimate)

True θ = `[-0.4, 0.44, -1.8, 0.24, -0.48, 0.85]`.

| Noise variance | L | Error % |
|---|---|---|
| 0.25 | 1000 | 5.78 |
| 0.25 | 2000 | 3.98 |

Higher noise variance degrades the estimate noticeably. The iterative algorithm gave a lower error than the recursive one in every case. See the report for the full per-iteration tables.

### Paper 3: Wiener model

- Fit on separate test data: **92%** for noise std 0.1 and **86.25%** for noise std 0.2.
- The method did not fit well on the DaISy heating-system and heat-exchanger datasets.

> Because the noise is random, individual runs can deviate slightly from the numbers above.

## References

1. Bao, Bo, et al. "Least squares based iterative parameter estimation algorithm for multivariable controlled ARMA system modelling with finite measurement data." *Mathematical and Computer Modelling* 53.9-10 (2011): 1664-1669.
2. Liu, Yanjun, Dongqing Wang, and Feng Ding. "Least squares based iterative algorithms for identifying Box-Jenkins models with finite measurement data." *Digital Signal Processing* 20.5 (2010): 1458-1467.
3. Zhou, Lincheng, Xiangli Li, and Feng Pan. "Least-squares-based iterative identification algorithm for Wiener nonlinear systems." *Journal of Applied Mathematics* 2013 (2013).

[Hossein Mosaffa's personal website](https://hossein-mosaffa.github.io)
