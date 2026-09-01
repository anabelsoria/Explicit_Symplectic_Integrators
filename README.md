# Explicit Symplectic Integrators

Symplectic integrators are designed to preserve the geometric properties of Hamiltonian systems over long time horizons, offering improved energy stability compared to traditional methods.

This repository implements a modular framework for symplectic integration of dynamical systems, with a focus on orbital mechanics problems such as the Circular Restricted Three-Body Problem (CR3BP) and the Two-Body Problem (TBP). 

The codebase is organized using MATLAB classes for flexibility and extensibility, and includes examples of orbit propagation and visualization. This framework serves as a foundation for exploring more complex systems, such as the Elliptic Restricted Three-Body Problem (ER3BP), and supports future development of additional features and orbit types.

<p align="center">
    <img src = "./Example_CR3BP/traj_dC.gif" width=80%>
</p>

## Necessary Software

1. MATLAB (tested with version R2023a).

## Directory Structure

- `lib/+astro/`: dynamical-system classes (`CR3BP`, `ER3BP`, `BCR4BP`, `HR4BP`, `TwoBody`) and orbital-element conversion utilities (`+conics`)
- `lib/integrators/`: the integrator classes (`SI`, `RK`, `TimeRegularized`, `Integrator`)
  - `lib/integrators/precision/CR3BP/`, `lib/integrators/precision/ER3BP/`: standalone one-timestep kernels for the compensated-summation/double-double precision ladder (uncompensated baseline, `CompSumScalar`/ICS, `ExtCompSum`/CS, double-double), fixed-step and time-regularized. Not yet wired into the `SI`/`TimeRegularized` classes as a selectable option — that's a follow-up.
- `lib/utils/`: shared numerical helpers (`comp_sum`, double-double arithmetic, plotting)

Example/driver scripts (`Example_TBP/`, `Example_CR3BP/`, etc.) are still at their original top-level locations pending a separate reorganization pass.

## References

[1] Soria-Carro, A., Akella, M., “Long-Duration Explicit Symplectic Approximations and Uncertainty Propagation for Cislunar Regimes,” AAS/AIAA Astrodynamics Specialist Conference, Broomfield, CO, August 2024. AAS Paper Number 24-258. 
