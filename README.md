# Explicit Symplectic Integrators

Symplectic integrators are designed to preserve the geometric properties of Hamiltonian systems over long time horizons, offering improved energy stability compared to traditional methods.

This repository implements a modular MATLAB framework for symplectic integration of dynamical systems, with a focus on orbital mechanics problems. The codebase is organized using MATLAB classes for flexibility and extensibility, and the main branch includes explicit symplectic integrators of different orders (2, 3, 4, 6, 8) for multiple dynamical systems (TBP, CR3BP, ER3BP, and BCR4BP). Each system is accompanied by example scripts that can be run directly to propagate orbits and generate plots of the results.

The codebase is organized using MATLAB classes for flexibility and extensibility, and includes examples of orbit propagation and visualization. This framework serves as a foundation for exploring more complex systems and supports future development of additional features and orbit types.

<p align="center">
    <img src = "./Example_CR3BP/traj_dC.gif" width=80%>
</p>

## Necessary Software

1. MATLAB (tested with version R2023a).

## Directory Structure

The repository contains the following main files:

- `Example_TBP/Example_TBP.m`: The main script with implementation of propagating a Two-Body Orbit
- `Example_CR3BP/Example_CR3BP.m`: The main script with implementation of propagating a Circular Restricted Three-Body orbit
- `Example_ER3BP/Example_ER3BP.m`: The main script with implementation of propagating an Elliptic Restricted Three-Body orbit
- `Example_BCR4BP/Example_BCR4BP.m`: The main script with implementation of propagating a Bicircular Restricted Four-Body orbit


## References

[1] Soria-Carro, A., Akella, M., “Long-Duration Explicit Symplectic Approximations and Uncertainty Propagation for Cislunar Regimes,” AAS/AIAA Astrodynamics Specialist Conference, Broomfield, CO, August 2024. AAS Paper Number 24-258. 
