# Bayesian Trans-Dimensional Soil Behaviour Type Inference

This MATLAB package implements a **Bayesian trans-dimensional Markov chain Monte Carlo** (RJ-MCMC) sampler for inferring soil behaviour type (SBT) interfaces and layer-wise Gaussian Process parameters from Cone Penetration Test (CPT) data.  It follows the methodology described in:

> Koch, M. C., Fujisawa, K. & Ray, A. (2025).  
> *Bayesian trans-dimensional soil behaviour type inference using conditional posterior proposals*.  
> Geophysical Prospecting, **73**(5), 1510–1533.  
> https://doi.org/10.1111/1365-2478.70021

---

## 📂 Repository layout

```
bayes_sbt_inference/
│
├── runExample.m             Main script: load data → run sampler → plot results
│
├── data/                    Observation data
│   ├── obs_data_real.mat
│   └── obs_data_synthetic.mat
│
├── src/                     Core sampling routines
│   ├── initializeState.m    Build `data`, `const`, `state`, `stats` structs
│   ├── runMCMC.m            Loop over n_samp and call the three blocks
│   ├── sampleFirstBlock.m   RJ-MCMC birth/death/perturb sampler (k, z, θ)
│   ├── sampleSecondBlock.m  Sample GP variance σ² per layer
│   ├── sampleThirdBlock.m   Sample GP correlation length ℓ per layer
│   ├── trandn.m             Truncated normal sampler (external dependency)
│   ├── updateCorrStr.m      Update correlation‐structure helper
│   ├── computeAvlblInterval.m
│   ├── computeCyInvMult.m
│   ├── computeMarLik.m
│   ├── computeOptFxnCoeff.m
│   ├── computeResiduals.m
│   ├── findInterfaceDepth.m
│   ├── priorProbzk.m
│   └── postProcessing.m     Compute diagnostics and generate plots
│
├── plotting/                Visualization helpers (7–8 figures)
│   ├── plotChainHistk.m
│   ├── plotContourThetaPdfz.m
│   ├── plotLogMarginalLikelihood.m
│   ├── plotMarginalsSigmal.m
│   ├── plotResiduals.m
│   ├── plotSbt.m
│   └── postProcessMaxDensz.m
│
├── README.md                (this file)
└── LICENSE                  MIT License
```

---

## 🔧 Requirements

- **MATLAB** R2019b or later  
- **Statistics and Machine Learning Toolbox** (for `ksdensity`, etc.) 

---

## 🚀 Quick Start

1. **Clone or download** this repository:
   ```bash
   git clone https://github.com/mkoch24/bayes_sbt_inference.git
   cd bayes_sbt_inference
   ```

2. **Open MATLAB**, change folder to the project root, and run:
   ```matlab
   runExample
   ```

   This will:
   1. Add the `src/`, `plotting/`, and `data/` folders to your MATLAB path.
   2. Initialize data, constants, sampler state, and diagnostics.
   3. Run the three-block RJ-MCMC sampler (`n_samp` iterations).
   4. Post-process the chain and generate all figures from the paper.

---

## ⚙️ Configuration

Inside **`runExample.m`** you can adjust:

- **Data file**: switch between real and synthetic  
- **Sampler settings**: number of samples, burn-in, model bounds, priors etc. 

---

## 📈 Outputs

- **Posterior samples** for:
  - Number of layers `k`
  - Interface depths `z`
  - Layer GP means `θ`
  - Layer GP variances `σ²`
  - Layer GP correlation lengths `ℓ`
- **Diagnostic plots** (from `plotting/`), including:
  - Trace histograms of `k`
  - Conditiional posterior pdfs of `z`
  - Contour plots of `θ` 
  - Log-marginal-likelihood evolution
  - Residual distributions
  - SBT classification summary (if needed)

All figures match those presented in the published paper.

---

## 🛠️ Extending or Modifying

- **Plotting**: edit or add functions in `plotting/`. Each plotting helper produces one figure.
- **Diagnostic metrics**: adjust `postProcessing.m` to compute additional statistics.

After any change, rerun:
```matlab
runExample
```
to ensure everything works end-to-end.

---

## 📝 License

This code is released under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## 📚 Citation

If you use this code, please cite:

```bibtex
@article{koch2025bayesian,
  title={Bayesian trans-dimensional soil behaviour type inference using conditional posterior proposals},
  author={Koch, Michael Conrad and Fujisawa, Kazunori and Ray, Anandaroop},
  journal={Geophysical Prospecting},
  volume={73},
  number={5},
  pages={1510--1533},
  year={2025},
  publisher={European Association of Geoscientists \& Engineers}
}
```

---

Enjoy exploring soil behaviour with Bayesian trans-dimensional inference!
