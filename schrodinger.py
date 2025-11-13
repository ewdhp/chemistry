"""
Schrödinger Equation: Mathematical Theory & 1D Particle-in-a-Box Calculation
----------------------------------------------------------------------------

The Schrödinger equation is the fundamental equation of quantum mechanics, describing how the quantum state of a physical system changes over time.

# Time-Dependent Schrödinger Equation (TDSE)

    iħ ∂Ψ(x, t)/∂t = [ - (ħ²/2m) ∇² + V(x, t) ] Ψ(x, t)

where:
    Ψ(x, t) : wave function (probability amplitude)
    ħ      : reduced Planck constant (h/2π)
    m      : mass of the particle
    ∇²     : Laplacian operator (spatial second derivatives)
    V(x, t): potential energy function

# Time-Independent Schrödinger Equation (TISE)

    [ - (ħ²/2m) ∇² + V(x) ] ψ(x) = E ψ(x)

where:
    ψ(x) : spatial part of the wave function
    E    : energy eigenvalue

# 1D Example: Particle in a Box

For a particle confined in a box of length L (V(x) = 0 inside, ∞ outside):

    - (ħ²/2m) d²ψ(x)/dx² = E ψ(x)
    ψ(0) = ψ(L) = 0

Solution:
    ψ_n(x) = sqrt(2/L) sin(nπx/L),   n = 1, 2, 3, ...
    E_n = (n²π²ħ²)/(2mL²)

# Physical Interpretation

- |Ψ(x, t)|² gives the probability density of finding the particle at position x at time t.
- The equation is linear and allows superposition of solutions.
- Boundary conditions and potentials V(x) determine allowed energy levels and wave functions.

# Key Concepts
- Wave-particle duality
- Quantum states and observables
- Operators: position, momentum, Hamiltonian
- Eigenvalues and eigenfunctions
- Normalization and orthogonality
"""


"""
Schrödinger Equation: Mathematical Theory & 1D Particle-in-a-Box Calculation
----------------------------------------------------------------------------

The Schrödinger equation is the fundamental equation of quantum mechanics, describing how the quantum state of a physical system changes over time.

# Time-Dependent Schrödinger Equation (TDSE)

    iħ ∂Ψ(x, t)/∂t = [ - (ħ²/2m) ∇² + V(x, t) ] Ψ(x, t)

where:
    Ψ(x, t) : wave function (probability amplitude)
    ħ      : reduced Planck constant (h/2π)
    m      : mass of the particle
    ∇²     : Laplacian operator (spatial second derivatives)
    V(x, t): potential energy function

# Time-Independent Schrödinger Equation (TISE)

    [ - (ħ²/2m) ∇² + V(x) ] ψ(x) = E ψ(x)

where:
    ψ(x) : spatial part of the wave function
    E    : energy eigenvalue

# 1D Example: Particle in a Box

For a particle confined in a box of length L (V(x) = 0 inside, ∞ outside):

    - (ħ²/2m) d²ψ(x)/dx² = E ψ(x)
    ψ(0) = ψ(L) = 0

Solution:
    ψ_n(x) = sqrt(2/L) sin(nπx/L),   n = 1, 2, 3, ...
    E_n = (n²π²ħ²)/(2mL²)

# Physical Interpretation

- |Ψ(x, t)|² gives the probability density of finding the particle at position x at time t.
- The equation is linear and allows superposition of solutions.
- Boundary conditions and potentials V(x) determine allowed energy levels and wave functions.

# Key Concepts
- Wave-particle duality
- Quantum states and observables
- Operators: position, momentum, Hamiltonian
- Eigenvalues and eigenfunctions
- Normalization and orthogonality


Quantum Harmonic Oscillator:

Analytic solution for a particle in a quadratic potential.
Shows quantized energy levels and Hermite polynomial wavefunctions.
Finite Square Well:

Particle in a box with finite potential walls.
Demonstrates bound and unbound states, tunneling.
Quantum Tunneling (Barrier Penetration):

Probability of a particle passing through a potential barrier.
Key for understanding alpha decay and scanning tunneling microscopy.
Hydrogen Atom:

3D solution for the electron in a Coulomb potential.
Produces atomic orbitals and quantized energy levels.
Double Well Potential:

Demonstrates quantum superposition and tunneling between wells.
Free Particle:

Plane wave solutions, momentum eigenstates.
Delta Function Potential:

Bound state formation with a singular attractive potential.
Kronig-Penney Model:

Periodic potential, basis for band theory in solids.
"""

import numpy as np
import matplotlib.pyplot as plt

def print_schrodinger_equation_theory():
    """
    Prints the mathematical theory of the Schrödinger equation.
    """
    print(__doc__)

def particle_in_box_demo(L=1.0, m=9.10938356e-31, n_levels=3):
    """
    Numerically demonstrates the 1D particle-in-a-box solution.
    Plots the first n_levels wavefunctions and prints their energies.
    Args:
        L: Length of the box (meters)
        m: Mass of the particle (kg, default is electron)
        n_levels: Number of energy levels to show
    """
    hbar = 1.054571817e-34  # Planck constant / 2pi (J·s)
    x = np.linspace(0, L, 1000)
    print("\nSchrödinger Equation: 1D Particle in a Box")
    print("Potential: V(x) = 0 for 0 < x < L, V(x) = ∞ otherwise")
    print("Boundary conditions: ψ(0) = ψ(L) = 0")
    print("Wavefunctions: ψ_n(x) = sqrt(2/L) * sin(nπx/L)")
    print("Energies: E_n = (n²π²ħ²)/(2mL²)")
    plt.figure(figsize=(8, 6))
    for n in range(1, n_levels+1):
        psi_n = np.sqrt(2/L) * np.sin(n * np.pi * x / L)
        E_n = (n**2 * np.pi**2 * hbar**2) / (2 * m * L**2)
        plt.plot(x, psi_n, label=f"n={n}, E={E_n:.2e} J")
        print(f"n={n}: E_n = {E_n:.4e} J")
    plt.title("Particle in a Box: Wavefunctions ψ_n(x)")
    plt.xlabel("x (meters)")
    plt.ylabel("ψ_n(x)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print_schrodinger_equation_theory()
    particle_in_box_demo(L=1e-9, m=9.10938356e-31, n_levels=3)  # Example: electron in a 1 nm box
