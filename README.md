# Data Driven for Nonlinear Dynamics ⚡

This repository contains MATLAB implementations of a **data-driven computational framework** for the simulation of nonlinear dynamic phenomena. The approach replaces traditional constitutive models with **finite sets of data points** (synthetic or experimental), enforcing them through a feedback operator within a **variational time-integration scheme**.

-------------------------------------------------------------------------------

# Repository for the DDND

This software can be used and distributed under the following license:

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

----------------------------------------------------------------------------------------------------
**DDND:** <br />
First version released on February 16, 2026.

**Warning:** <br />
It works on Matlab 7+ <br />
A technical description of the implementation can be found in the following paper:

Roccia, B.A., Lind, P. and Gebhardt, C.G., , "Data-driven computational mechanics meets nonlinear dynamics: unlocking bifurcations, limit cycle oscillations, chaos and synchronization," to be submitted, 2026. 

-------------------------------------------------------------------------------

## ✨ Features
- Variational integrators formulated on QxQ, T*Q, and TQ.
- Discrete Forward Problem for Nonlinear Dynamical Systems for every variational algorithm.
- Data-Driven Discrete forward-dynamics Problem for Nonlinear Dynami-cal Systems.
- Characterization of the associated KKT system.  
- **Alternating Direction Method (ADM)** solver for DDCM applied to nonlinear dynamic systems.  
- DDCM verification.
- Nonlinear systems studied: forced Duffing equation, van der Pol-Duffing equation, Kuramoto oscillator.  
- Tools for post-processing:
  - Poincaré maps
  - Basin of attractions
  - Limit cycle oscillations
  - Chaos
  - Basins of synchronization  
  - Residuum, error cost, and global cost function histories  

------------------------------------------------------------------------------

## 📂 Repository Structure
- `DDEC.m` – Example script that assembles and simulates a circuit using DDCM.  
- `circuit_assembly.m` – Builds circuit matrices (L, C, constraints).  
- `KKT_Time.m` – Constructs KKT matrices for time integration.
- `KKT_initialState.m` - Constructs KKT matrices for time step 0
- `phiESOperator.m` – Feedback operator: projects trial states to closest data points.  
- `DDCM_ADM_Electric.m` – Core solver implementing DDCM with ADM.

---


## Contact

Dr. Bruno Roccia (bruno.roccia@uib.no), University of Bergen, Norway <br />
