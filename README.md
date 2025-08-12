# ModelingW25
ESAM346, other projects
# MATLAB W25 Projects – Numerical Methods for Dynamical Systems

This repository contains two major simulation-based projects developed for Modeling and Computation in Applied Math. Each project focuses on numerical methods applied to different physical models using MATLAB.

## Contents

- [Vortex Interaction Simulation]
- [Fiber Extrusion Simulation]

---

## Vortex Interaction Simulation

**Objective:**  
Simulate the behavior of vortices and a passive particle in 2D flow using numerical integration methods.

**Overview:**  
The vortex simulation analyzes how particles move under the influence of 2D vortices. Several test cases explore the effect of time step size and numerical method accuracy, especially Euler’s method versus 4th-order Runge-Kutta (RK4). The project demonstrates convergence behavior and limitations of explicit methods in chaotic or symmetric systems.

**Key Components:**

- **Part 1a–1c:**  
  Simulates a single vortex flow using both Euler and RK4 methods with varying step sizes. Shows improved accuracy and convergence of RK4:
  - RK4 drastically reduces error from ~1.5 (Euler, N=50) to ~10⁻⁵ (RK4, N=50) and smaller with finer steps.

- **Part 2a–2c:**  
  Simulates multiple interacting vortices, showing emergent group behavior and system sensitivity.
  - Highlights how closely spaced vortices behave like a single vortex.
  - Explores numerical instability due to round-off error and symmetry breaking in tightly packed systems.

**Notable Insight:**  
The RK4 method is far more stable and accurate in preserving expected trajectories in vortex dynamics, but even it is susceptible to error accumulation and chaotic divergence in symmetric systems with small perturbations.

---

## Fiber Extrusion Simulation with and without Noise

**Objective:**  
Model the trajectory of a flexible fiber being extruded in 2D space, both deterministically and stochastically, to investigate convergence, geometric behavior, and noise effects.

**Overview:**  
This project models fiber extrusion using a system of ODEs governed by curvature (via a bending function \( b(r) \)). Two models are implemented:
- **Deterministic model:** Solved with RK4
- **Stochastic model:** Solved using Euler-Maruyama method

**Key Features:**

- `fiber_Anders.m`:  
  Implements RK4 for deterministic simulation of the fiber path given \( b(r) \), initial conditions, and time horizon.

- `fibersde_Anders.m`:  
  Extends the model with stochastic noise using Euler-Maruyama, simulating physical uncertainty in extrusion.

**Highlights:**

- **Problem 1–3:**  
  Confirm expected circular or spiral paths using different \( b(r) \) functions (\( b(r) = r \) and \( b(r) = 1/r \)) and analyze convergence rates.
  
- **Problem 4:**  
  Tests unusual configurations and parameters, analyzing behavior over long extrusions.

- **Problem 5–6:**  
  Examines how added noise (parameter A) affects trajectory accuracy and convergence. Shows how noise can cause deviation even when deterministic paths are well-behaved.

- **Problem 7–8:**  
  Visualizes random sample paths and compares terminal distributions under different bending functions, emphasizing the impact of noise.

**Conclusion:**  
RK4 is shown to be highly accurate for this system, with observed convergence rates near 4. The Euler-Maruyama method captures the effects of stochasticity but is less precise. The combined analysis illustrates the power and limits of deterministic versus stochastic modeling in physical simulations.

---

## Requirements

- MATLAB (tested with R2023+)
- No external libraries required

## Running the Code

Each project is self-contained:
- Run `ESAM346_HW4.m` for Fiber Extrusion (all deterministic and stochastic problems in one script).
- The vortex simulation files can be modularized similarly `ESAM346_HW1.m`.

