This repository contains the MATLAB and Mathematica codes corresponding to the manuscript  
*"Stability Analysis of a Multistage Compartmental Model of Infection"*, published in the *Journal of Biological Dynamics*.

---

#  Files Included

## 1. `General_Multistage_SIRS.ml`
- Computes the positive coefficients for the candidate Lyapunov functions.
- INPUT: the number of stages (n=1,2,3), Type of equilibrium (DFE or EE)
- OUTPUT: The possitive coefficients for the Lyapunov functions and the derivative of the Lyapunov functions. For further simplification, A_coeff matrix and b_coeff vector can be passed to Mathematica as well.

## 2. `MULTISTAGE_SIRS_SIMPLIFY.nb`
- Mathematica script used to simplify the expressions obtained from MATLAB, particularly for the case n ≥ 3.
- INPUT: A_coeff matrix and b_coeff vector obtained in the MATLAB code above
- OUTPUT: Simplified expresseion for the positive coefficients also conditions to make them positive.

## 3. `Simulation_Multistage_n_2.ml`
- Generates simulation results for the two-stage SIRS model.
- OUTPUT: Simulation results for the two-stage SIRS model. 

## 4. `Simulation_Multistage_n_3.ml`
- Generates simulation results for the three-stage SIRS model.
- OUTPUT: Simulation results for the three-stage SIRS model. 

---

# Notes
- MATLAB Live Scripts (`.mlx`) are used for analytical computations and simulations.
- The Mathematica notebook is used for symbolic simplification.
