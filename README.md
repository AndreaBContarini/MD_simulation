# Molecular Dynamics Simulation of a Lennard-Jones System

This project implements a molecular dynamics (MD) simulation in C++ to study the behavior of particles interacting via the Lennard-Jones potential in a canonical ensemble (NTV). The simulation uses advanced numerical algorithms to provide accurate thermodynamic properties and validates results against published literature.

---

## Project Highlights

- **Objective**: Simulate a system of particles under Lennard-Jones interactions and analyze thermodynamic properties like pressure and energy.
- **Key Features**:
  - Canonical ensemble simulation (constant number of particles, volume, and temperature).
  - Implements Tuckerman's reversible integration algorithm for precise motion calculations.
  - Supports observables like potential energy and pressure, with tail corrections for long-range interactions.
- **Algorithms Used**:
  - Velocity Verlet integrator.
  - Nosé-Hoover thermostat for temperature control.
  - Tail corrections for improved accuracy in pressure and energy calculations.

---

## Files in the Project

### 1. Source Code Files
- **`mdsim.cpp`**:
  - Entry point of the simulation program.
  - Initializes the random number generator and the simulation system.
  - Calls functions to prepare initial configurations and execute the simulation.

- **`params.hpp`**:
  - Defines parameters for the simulation, including particle density, temperature, and cutoff radius.
  - Encapsulates configurations for both Monte Carlo and MD simulations.

- **`particle.hpp`**:
  - Defines the `particleLJ` class to handle particle interactions, updates, and calculations using Lennard-Jones forces.

- **`pvector.hpp`**:
  - Provides vector operations such as dot products and vector addition, essential for 3D simulations.

- **`randnumgen.hpp`**:
  - Implements a random number generator used for initializing particle positions and velocities.

- **`sim.hpp`**:
  - Houses simulation classes (`simLJ`, `mcsimLJ`, and `mdsimLJ`) that orchestrate initialization, execution, and data output.

- **`pmatrix.hpp`**:
  - Handles matrix operations, potentially used for advanced calculations or interactions.

---

### 2. Documentation
- **`Article.pdf`**:
  - A detailed report explaining the theoretical foundation of the simulation.
  - Covers:
    - Lennard-Jones potential and its implementation.
    - Tuckerman's algorithm and Velocity Verlet integrator.
    - Observables such as pressure and energy, with formulas for tail corrections.
    - Comparative analysis with literature results.
  - Includes Python scripts for plotting results.

---

## Simulation Details

### 1. Lennard-Jones Potential
The Lennard-Jones potential models interactions between particles:
\[
V(r) = 4\epsilon \left[ \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^6 \right]
\]
where:
- \( \epsilon \): Depth of the potential well.
- \( \sigma \): Finite distance where the potential reaches zero.
- \( r \): Distance between two particles.

### 2. Observables
- **Pressure**:
  \[
  P^* = \rho^* T^* + \frac{1}{3V} \sum_{i<j} \mathbf{f}_{ij} \cdot \mathbf{r}_{ij}
  \]
- **Potential Energy**:
  \[
  u^* = \frac{1}{N_p} \sum_{i<j} 4\epsilon \left[ \left( \frac{\sigma}{r_{ij}} \right)^{12} - \left( \frac{\sigma}{r_{ij}} \right)^6 \right]
  \]

### 3. Algorithms
- **Tuckerman’s Algorithm**: Used for reversible and stable integration of equations of motion.
- **Nosé-Hoover Thermostat**: Ensures realistic temperature fluctuations for the canonical ensemble.

---

## How to Use

### 1. Prerequisites
- **Compiler**: C++11 or later.
- **Libraries**: Standard C++ libraries (`<vector>`, `<cmath>`, etc.).

### 2. Compilation
Run the following command to compile:
```bash
g++ mdsim.cpp -o mdsim -std=c++11
