# DDCM Duffing Oscillator (ADM + Variational Phase-Space Integrator)

This repository provides a compact MATLAB implementation of a **Data-Driven Computational Mechanics (DDCM)** forward-dynamics solver applied to the **forced Duffing oscillator**. The method combines:

- a **phase-space variational time integrator** (symmetric/"Kane"-type discrete Lagrangian structure), and
- an **Alternating Direction Method (ADM)** that alternates between:
  - a continuous equilibrium/compatibility update (via a KKT solve), and
  - a discrete projection onto an experimental/synthetic constitutive dataset.

The provided example reproduces a canonical Duffing setup (including parameter regimes that can exhibit period-doubling cascades).

---

## Files

- `DDCM_Duffing.m`  
  Executable example script: defines parameters, generates a dataset, runs the solver, and plots results.

- `DDCMSolver.m`  
  Main solver: DDCM forward dynamics with an ADM fixed-point loop at each time step.

- `MyForce.m`  
  External force definitions (cosine, sine with ramp-on, sigmoid ramp, smoothstep).

- `phiESOperator.m`  
  Wrapper that applies the elementwise projection `phiES` to each stage/element.

- `phiES.m`  
  Elementwise data-driven projection: selects the closest discrete point in the dataset.

---

## Requirements

- MATLAB R2018b or newer (earlier versions may work but are not tested).
- No toolboxes required.

---

## Quick start

1. Clone/download the repository.
2. Make sure the folder is on your MATLAB path.
3. Run:

```matlab
DDCM_Duffing
```

You should obtain a constitutive-space plot showing:

- the dataset `(e,f)` samples,
- the **selected** discrete pairs for the two stages `(1-\alpha)` and `\alpha`, and
- the **continuous** mid values returned by the KKT step.

---

## How to modify the example

All user-facing settings are in `DDCM_Duffing.m` under **Problem definition**.

### Duffing parameters

```matlab
DATA.m   = 2.0;   % mass
DATA.c   = 0.5;   % damping
DATA.KL  = 1.0;   % linear stiffness
DATA.KNL = 0.5;   % cubic stiffness
```

### Forcing

```matlab
DATA.F    = 1.0;     % amplitude
DATA.W    = 0.5;     % frequency
DATA.Flag = 0;       % forcing type (see MyForce.m)
```

### Time grid

```matlab
DATA.Dt    = 0.025;
DATA.Tspan = [0 200];
```

### ADM parameters

```matlab
DATA.Iter       = 50;     % max iterations per time step
DATA.Tol        = 1e-12;  % residual tolerance
DATA.TolC       = 1e-6;   % cost change tolerance
DATA.DetectJump = 4;      % heuristic jump detector
```

### Dataset resolution

```matlab
nd   = 101;
qmin = -1.5;
qmax = +1.5;
```

Increase `nd` to densify the dataset.

---

## Output (`Solution` struct)

`DDCMSolver` returns a struct with the most relevant fields:

- `Solution.Time` : time grid
- `Solution.Z`    : phase-space state `[q; p]` over time (`2 x nSteps`)
- `Solution.qmid, Solution.Fmid` : continuous constitutive values from the KKT step (`2 x (nSteps-1)`)
- `Solution.qtil, Solution.Ftil` : discrete pairs selected from the dataset (`2 x (nSteps-1)`)
- `Solution.Iter` : ADM iteration count per time step
- `Solution.Res`  : equilibrium residual per step
- `Solution.Costq, Solution.CostF` : cost contributions
- `Solution.Index` : index history of selected data points (per time step)
- `Solution.Jump` / `Solution.IndxJump` : jump heuristic diagnostics

---

## Notes on the ADM projection (`phiES`)

The projection operator selects the closest data point `(\tilde e, \tilde f)` to a continuous query `(e,f)`.

- For `p = 1` (default in the example), the original code effectively uses **only the elongation distance** `|e-edata|`.
- For `p \ne 1`, `phiES` uses a separable `p/q`-type metric with the Hölder conjugate `q = p/(p-1)`.

If you have an experimental dataset, replace the synthetic generation in `DDCM_Duffing.m` by loading your measured `(e,f)` pairs.

---

## Reproducibility

To reproduce the same results as the example:

- use the provided default values in `DDCM_Duffing.m`,
- keep `DATA.Alpha = 0.5`, and
- keep the dataset symmetric about zero (`qmin = -qmax`).

---

## License

Choose a license appropriate for your journal and institution (e.g., MIT, BSD-3-Clause, Apache-2.0). A permissive license is typically preferred for open-source academic code.

---

## Citation

If you use this code, please cite the associated journal publication:

```bibtex
@article{RocciaDDND2026,
  title   = {Data-driven computational mechanics meets nonlinear dynamics: unlocking bifurcations, limit cycle oscillations, chaos and synchronization},
  author  = {Roccia B.A., and Lind, P., and Gebhardt, C.G.},
  journal = {XXX},
  year    = {2026},
  doi     = {DOI:}
}
```
