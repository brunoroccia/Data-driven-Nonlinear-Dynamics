# Data Driven for Electric Circuits ⚡

This repository contains MATLAB implementations of a **data-driven computational framework** for the simulation of electric circuits. The approach replaces traditional constitutive models with **finite sets of data points** (synthetic or experimental), enforcing them through a feedback operator within a **variational time-integration scheme**.

-------------------------------------------------------------------------------

# Repository for the DDEC

This software can be used and distributed under the following license:

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

----------------------------------------------------------------------------------------------------
**DDEC:** <br />
First version released on September 01, 2025.

**Warning:** <br />
It works on Matlab 7+ <br />
A technical description of the implementation can be found in the following paper:

Gebhardt, C.G., Roccia, B.A., Bossio, J.M., and Bossio, G.R., , "A framework for data-driven simulation of electrical circuits based on discrete-continuous optimization," submitted to IEEE Access, 2025. 

-------------------------------------------------------------------------------

## ✨ Features
- Circuit assembly from **graph representations** with inductors and capacitors.  
- Non-holonomic constraints handling via reduced incidence matrices.  
- **Alternating Direction Method (ADM)** solver for DDCM applied to RLC-type networks.  
- Synthetic capacitor data generation with nonlinear effects.  
- Tools for post-processing:
  - Phase-space plots (charge–voltage pairs vs. data)  
  - Time evolution of charges  
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
