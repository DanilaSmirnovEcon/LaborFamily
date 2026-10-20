# Labor and Family Dynamics in a Joint-Search Framework — model code

MATLAB code that solves, simulates and evaluates the model in Morazzoni and
Smirnov, "Labor and Family Dynamics in a Joint-Search Framework".

## Requirements

MATLAB (R2021a or later) with the Statistics and Machine Learning Toolbox
(`gppdf`, `logncdf`, `logninv`) and the Parallel Computing Toolbox (`parfor`;
without it the simulation runs serially).

## How to run

| Script | What it does | Paper output |
|---|---|---|
| `Main.m` | Solves the baseline model, simulates a panel of households, computes moments, prints them next to the data targets and draws the figures | Tables 3, 4 (except the elasticity of marriage rates to wages), 5 (model with endogenous marriage), 6; Figures 1, 5, 6, 7, 8 |
| `PolicyCounterfactuals.m` | Solves and simulates the baseline and the two policies, with and without endogenous marriage (about 6 times the run time of `Main.m`) | Table 5 (model without endogenous marriage), Tables 7–10 |

Outputs written to the working directory:

- `solutionSP.mat`: equilibrium value functions (baseline)
- `simulation.mat`: simulated panel and moments (baseline)
- `calibration_table.csv`: model vs. data for Tables 3–6 (`calibration_table.m`)
- `Plots/`: figures (`PlotSimul.m`). Figure 1: `Fig1_*.png`; Figure 5:
  `EndogWageDiff_FvM.png`, `EndogUnemploymentDiff_FvM.png`; Figure 6:
  `WageChangeAfterUnemployment_S_vs_M.png`, `WageChangeAfterUnemployment_Kids.png`;
  Figure 7: `EndogMarriedShare.png`, `EndogProductivities.png`; Figure 8:
  `HeatmapAbilitiesMarried.png`
- `policy_tables.csv`, `PolicyCounterfactuals_results.mat`: Tables 7–10

The simulation is reproducible for a given `simParam.seed` regardless of the
number of parallel workers.

## Files

| File | Role |
|---|---|
| `Param_Gen.m` | Parameters common to both genders (Table 2); optional policy hook |
| `Param_F.m`, `Param_M.m` | Gender-specific parameters (Table 2) and single-agent grids |
| `BuildSingleGrid.m`, `BuildDoubleGrid.m` | State grids and finite-difference generators for productivity |
| `SolveSingleProblem.m` | Singles' HJBs without marriage; their stationary distribution is the marriage-market pool |
| `SolveEndogProblem.m` | Singles' and couples' HJBs with endogenous marriage (Section 2.5, Appendix B) |
| `approximateValueFunctions.m` | Interpolants of the value functions used in the simulation |
| `SimulatePanel.m` | Simulated panel with endogenous marriage |
| `SimulatePanelCounter.m` | Simulated panel without endogenous marriage (52% married at birth to a random partner, no further meetings; values solved with μ = 0) |
| `compute_moments.m` | Moments from the simulated panel |
| `calibration_table.m` | Model vs. data table (Tables 3–6) |
| `fiscal_welfare.m` | Fiscal accounts (T, B, credits) and aggregate welfare (Section 5) |
| `policy_tables.m` | Tables 7–10 from the policy runs |
| `PlotSimul.m` | Paper figures; `PlotSimul_Extras.m` adds diagnostic figures |
| `util.m` | CRRA utility |

## Model and code notation

Time is in months. Wages, benefits, child costs and home production are monthly
flows in model consumption units (c.u.); ȳ = 61 c.u. corresponds to $45,000 a
year ($25 an hour) in 2018 dollars.

| Paper | Code | Paper | Code |
|---|---|---|---|
| γ | `param.nu` | μ | `param.m` |
| ρ | `param.rho` | π^s, π^c | `param.pcs`, `param.pcp` |
| ξ | `param.zeta` | κ_min, κ_scale | `param.cmin`, `param.cscale` |
| η1, η2 | `param.lambda_tax`, `param.tau_tax` | h^s_{g,k}, h^c_{g,k} | `pf/pm.hs0, hs1, hc0, hc1` |
| ȳ | `param.averagewage` | θ_w | `param.theta_w` (`pf.wscale`) |
| β^{Y1}_k, β^{Y2}_k | `param.pi1_n0/n1`, `param.pi2_n0/n1` | σ_w | `param.sigma_w` |
| y^{Y1}_k, y^{Y2}_k | `param.E1_n0/n1`, `param.E2_n0/n1` | δ_g | `pf/pm.delta` |
| b_min, b_scale, b_max | `param.bmin`, `param.bscale`, `param.bmax` | λ_g | `pf/pm.lambda_u` |
| φ_e, φ_ne | `param.phi_e`, `param.phi_ne` | ν_g | `pf/pm.lambda_oE` |
| ε, ā | `param.epsilon`, `param.abar` | μ_a0, σ_a0 | `param.mu_a0`, `param.sigma_a0` |

Numerical settings not in Table 2: grids (`param.a`, `param.w`: productivity on
[0, 3], wage offers on [10, 30]), the switching intensity into and out of the
labor force (`pf/pm.lambda_o` = 50) and the trembling-hand parameter
(`param.tremble`), see Appendix B.2.
