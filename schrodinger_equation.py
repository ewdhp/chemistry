#!/usr/bin/env python3
"""
THE SCHRÖDINGER EQUATION: Quantum Mechanical Wave Equation

================================================================================
THEORY AND FUNDAMENTAL CONCEPTS
================================================================================

1. INTRODUCTION TO THE SCHRÖDINGER EQUATION
--------------------------------------------
The Schrödinger equation is the fundamental equation of quantum mechanics,
describing how the quantum state of a physical system changes with time.
It was formulated by Erwin Schrödinger in 1926 and represents one of the
most important equations in physics.

Historical Context:
• 1900: Planck introduces quantum hypothesis (E = hν)
• 1913: Bohr model of the atom (quantized orbits)
• 1924: de Broglie proposes wave-particle duality (λ = h/p)
• 1925: Heisenberg develops matrix mechanics
• 1926: Schrödinger formulates wave mechanics
• 1926: Born provides probabilistic interpretation

The Schrödinger equation describes particles as wave functions rather than
point particles, fundamentally changing our understanding of nature at the
atomic scale.


2. TIME-DEPENDENT SCHRÖDINGER EQUATION
---------------------------------------
The most general form describing the evolution of a quantum state:

    iℏ ∂Ψ(r,t)/∂t = ĤΨ(r,t)

Where:
• i = imaginary unit (√-1)
• ℏ = h/2π = reduced Planck constant (1.054571817×10⁻³⁴ J·s)
• Ψ(r,t) = wave function (complex-valued probability amplitude)
• Ĥ = Hamiltonian operator (total energy operator)
• r = position vector
• t = time

The Hamiltonian operator:
    Ĥ = -ℏ²/2m ∇² + V(r,t)

Where:
• -ℏ²/2m ∇² = kinetic energy operator
• ∇² = Laplacian operator (∂²/∂x² + ∂²/∂y² + ∂²/∂z²)
• V(r,t) = potential energy function
• m = particle mass


3. TIME-INDEPENDENT SCHRÖDINGER EQUATION
-----------------------------------------
For systems with time-independent potentials V(r), the wave function can be
separated: Ψ(r,t) = ψ(r)φ(t)

This leads to the time-independent Schrödinger equation:

    Ĥψ(r) = Eψ(r)

Or explicitly:
    -ℏ²/2m ∇²ψ(r) + V(r)ψ(r) = Eψ(r)

This is an eigenvalue equation where:
• ψ(r) = spatial wave function (eigenfunction)
• E = energy (eigenvalue)
• Solutions exist only for specific values of E (quantized energy levels)

One-dimensional form:
    -ℏ²/2m d²ψ(x)/dx² + V(x)ψ(x) = Eψ(x)


4. PHYSICAL INTERPRETATION
---------------------------

Born's Probabilistic Interpretation (1926):
The wave function itself has no direct physical meaning, but its square
gives the probability density:

    P(r) = |Ψ(r,t)|² = Ψ*(r,t)Ψ(r,t)

Where Ψ* is the complex conjugate of Ψ.

Normalization Condition:
The total probability of finding the particle anywhere must be 1:

    ∫ |Ψ(r,t)|² d³r = 1  (over all space)

Expectation Values:
The average value of an observable O is:

    ⟨O⟩ = ∫ Ψ* Ô Ψ d³r

Where Ô is the quantum operator corresponding to the observable.


5. CLASSIC QUANTUM SYSTEMS
---------------------------

A. PARTICLE IN A 1D BOX (Infinite Square Well)
-----------------------------------------------
Potential:
    V(x) = 0     for 0 < x < L
    V(x) = ∞     for x ≤ 0 or x ≥ L

Boundary conditions: ψ(0) = ψ(L) = 0

Solutions (Energy Eigenstates):
    ψₙ(x) = √(2/L) sin(nπx/L)    n = 1, 2, 3, ...
    
Energy Levels:
    Eₙ = n²π²ℏ²/(2mL²) = n²h²/(8mL²)

Key Features:
• Quantized energy levels (discrete spectrum)
• Ground state (n=1) has lowest energy: E₁ = π²ℏ²/(2mL²)
• Zero-point energy (E₁ ≠ 0) - quantum mechanical effect
• Energy increases as n²
• Wave functions are standing waves with n-1 nodes


B. HARMONIC OSCILLATOR
-----------------------
Potential (quadratic):
    V(x) = ½kx² = ½mω²x²

Where:
• k = spring constant (force constant)
• ω = √(k/m) = angular frequency
• x = displacement from equilibrium

Energy Levels:
    Eₙ = ℏω(n + ½)    n = 0, 1, 2, 3, ...

Wave Functions:
    ψₙ(x) = (mω/πℏ)^(1/4) × 1/√(2ⁿn!) × Hₙ(√(mω/ℏ)x) × exp(-mωx²/2ℏ)

Where Hₙ is the nth Hermite polynomial.

Key Features:
• Equally spaced energy levels (ΔE = ℏω)
• Zero-point energy: E₀ = ½ℏω (even in ground state)
• Wave functions are Hermite polynomials × Gaussian
• Models molecular vibrations, phonons in solids
• Classical limit: particle oscillates at amplitude A where E = ½kA²


C. HYDROGEN ATOM (3D Central Potential)
----------------------------------------
Potential (Coulomb):
    V(r) = -e²/(4πε₀r) = -ke²/r

Energy Levels (Bohr formula from Schrödinger equation):
    Eₙ = -13.6 eV/n²    n = 1, 2, 3, ...

Wave Functions:
    ψₙₗₘ(r,θ,φ) = Rₙₗ(r) × Yₗₘ(θ,φ)

Where:
• Rₙₗ(r) = radial wave function (depends on n, l)
• Yₗₘ(θ,φ) = spherical harmonic (angular part)
• n = principal quantum number (1, 2, 3, ...)
• l = angular momentum quantum number (0, 1, ..., n-1)
• m = magnetic quantum number (-l, ..., +l)

Key Features:
• Three quantum numbers arise naturally from boundary conditions
• Degeneracy: energy depends only on n, not on l or m
• Wave functions are products of radial and angular parts
• Exact analytical solution (rare in quantum mechanics)


D. QUANTUM TUNNELING
---------------------
Particles can penetrate through potential barriers that would be classically
forbidden (E < V).

Transmission Coefficient (rectangular barrier):
    T ≈ exp(-2κa)

Where:
• κ = √(2m(V-E)/ℏ²) = decay constant in barrier
• a = barrier width
• V = barrier height
• E = particle energy

Applications:
• Alpha decay (radioactivity)
• Scanning tunneling microscope (STM)
• Tunnel diodes
• Quantum dots and resonant tunneling


6. OPERATORS IN QUANTUM MECHANICS
----------------------------------

Position operator:      x̂ = x (multiply by x)
Momentum operator:      p̂ = -iℏ ∂/∂x
Kinetic energy:         T̂ = p̂²/2m = -ℏ²/2m ∂²/∂x²
Potential energy:       V̂ = V(x) (multiply by V(x))
Total energy (Hamiltonian): Ĥ = T̂ + V̂
Angular momentum:       L̂ = r̂ × p̂

Commutation Relations:
    [x̂, p̂] = iℏ (Heisenberg uncertainty relation)


7. HEISENBERG UNCERTAINTY PRINCIPLE
------------------------------------
Fundamental limit on simultaneous measurement precision:

Position-Momentum:
    Δx · Δp ≥ ℏ/2

Energy-Time:
    ΔE · Δt ≥ ℏ/2

Consequences:
• Cannot know both position and momentum exactly
• Zero-point energy (particle cannot be at rest)
• Width of energy levels
• Fundamental to quantum mechanics, not measurement limitation


8. APPLICATIONS
---------------

Chemistry:
• Atomic and molecular structure
• Chemical bonding (molecular orbital theory)
• Spectroscopy (energy level transitions)
• Reaction rates (tunneling through barriers)

Physics:
• Condensed matter physics (band structure)
• Nuclear physics (shell model)
• Particle physics (quantum field theory)
• Quantum optics and photonics

Technology:
• Semiconductors and transistors
• Lasers and LEDs
• Quantum computing (qubits)
• Scanning tunneling microscopy
• Magnetic resonance imaging (MRI)

Materials Science:
• Electronic properties of materials
• Superconductivity
• Quantum dots and nanostructures
• Catalysis

================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.special import hermite, genlaguerre, sph_harm, factorial
from scipy.linalg import eigh_tridiagonal
import matplotlib.patches as mpatches
import math

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

h = 6.62607015e-34  # Planck constant (J·s)
hbar = h / (2 * np.pi)  # Reduced Planck constant (J·s)
m_e = 9.1093837015e-31  # Electron mass (kg)
e = 1.602176634e-19  # Elementary charge (C)
epsilon_0 = 8.8541878128e-12  # Vacuum permittivity (F/m)
c = 299792458  # Speed of light (m/s)
a0 = 5.29177210903e-11  # Bohr radius (m)
eV_to_J = 1.602176634e-19  # Electron volt to Joule conversion

# ============================================================================
# ANALYTICAL SOLUTIONS
# ============================================================================

def particle_in_box_wavefunction(x, n, L):
    """
    Wave function for particle in a 1D box.
    
    ψₙ(x) = √(2/L) sin(nπx/L)
    
    Args:
        x: Position array
        n: Quantum number (1, 2, 3, ...)
        L: Box length
    
    Returns:
        Wave function values
    """
    return np.sqrt(2/L) * np.sin(n * np.pi * x / L)

def particle_in_box_energy(n, L, m=m_e):
    """
    Energy levels for particle in a 1D box.
    
    Eₙ = n²h²/(8mL²)
    
    Args:
        n: Quantum number (1, 2, 3, ...)
        L: Box length (m)
        m: Particle mass (kg)
    
    Returns:
        Energy in Joules
    """
    return n**2 * h**2 / (8 * m * L**2)

def harmonic_oscillator_wavefunction(x, n, omega, m=m_e):
    """
    Wave function for quantum harmonic oscillator.
    
    ψₙ(x) = (mω/πℏ)^(1/4) × 1/√(2ⁿn!) × Hₙ(ξ) × exp(-ξ²/2)
    where ξ = √(mω/ℏ) x
    
    Args:
        x: Position array
        n: Quantum number (0, 1, 2, ...)
        omega: Angular frequency (rad/s)
        m: Particle mass (kg)
    
    Returns:
        Wave function values
    """
    # Dimensionless coordinate
    xi = np.sqrt(m * omega / hbar) * x
    
    # Normalization constant
    norm = (m * omega / (np.pi * hbar))**(1/4) / np.sqrt(2**n * math.factorial(n))
    
    # Hermite polynomial
    H_n = hermite(n)
    
    # Wave function
    psi = norm * H_n(xi) * np.exp(-xi**2 / 2)
    
    return psi

def harmonic_oscillator_energy(n, omega):
    """
    Energy levels for quantum harmonic oscillator.
    
    Eₙ = ℏω(n + 1/2)
    
    Args:
        n: Quantum number (0, 1, 2, ...)
        omega: Angular frequency (rad/s)
    
    Returns:
        Energy in Joules
    """
    return hbar * omega * (n + 0.5)

def hydrogen_radial_wavefunction(r, n, l, Z=1):
    """
    Radial wave function for hydrogen-like atoms.
    
    Simplified for demonstration purposes.
    
    Args:
        r: Radial distance (in units of a0)
        n: Principal quantum number (1, 2, 3, ...)
        l: Angular momentum quantum number (0, 1, ..., n-1)
        Z: Nuclear charge
    
    Returns:
        Radial wave function R_nl(r)
    """
    rho = 2 * Z * r / n
    
    # Normalization (simplified)
    norm = np.sqrt((2*Z/n)**3 * math.factorial(n-l-1) / (2*n*math.factorial(n+l)))
    
    # Laguerre polynomial
    L = genlaguerre(n-l-1, 2*l+1)
    
    # Radial function
    R = norm * np.exp(-rho/2) * rho**l * L(rho)
    
    return R

def hydrogen_energy(n):
    """
    Energy levels for hydrogen atom.
    
    Eₙ = -13.6 eV / n²
    
    Args:
        n: Principal quantum number (1, 2, 3, ...)
    
    Returns:
        Energy in eV
    """
    return -13.6 / n**2

def transmission_coefficient(E, V, a, m=m_e):
    """
    Transmission coefficient for rectangular barrier (tunneling).
    
    T ≈ exp(-2κa) for E < V
    where κ = √(2m(V-E)/ℏ²)
    
    Args:
        E: Particle energy (J)
        V: Barrier height (J)
        a: Barrier width (m)
        m: Particle mass (kg)
    
    Returns:
        Transmission coefficient (0 to 1)
    """
    if E >= V:
        # Classical: particle goes over barrier
        return 1.0
    else:
        # Quantum tunneling
        kappa = np.sqrt(2 * m * (V - E)) / hbar
        T = np.exp(-2 * kappa * a)
        return T

# ============================================================================
# NUMERICAL SOLUTION (FINITE DIFFERENCE METHOD)
# ============================================================================

def solve_schrodinger_1d(V, x_min, x_max, N, m=m_e, num_states=5):
    """
    Numerically solve 1D time-independent Schrödinger equation using
    finite difference method.
    
    Args:
        V: Potential function V(x)
        x_min, x_max: Domain boundaries
        N: Number of grid points
        m: Particle mass
        num_states: Number of energy eigenstates to find
    
    Returns:
        x: Position array
        energies: Energy eigenvalues
        wavefunctions: Energy eigenfunctions
    """
    # Create spatial grid
    x = np.linspace(x_min, x_max, N)
    dx = x[1] - x[0]
    
    # Evaluate potential on grid
    V_grid = V(x)
    
    # Construct Hamiltonian matrix using finite differences
    # -ℏ²/2m d²ψ/dx² + V(x)ψ = Eψ
    
    # Kinetic energy: -ℏ²/2m (ψ[i+1] - 2ψ[i] + ψ[i-1])/dx²
    diagonal = hbar**2 / (m * dx**2) + V_grid[1:-1]
    off_diagonal = -hbar**2 / (2 * m * dx**2) * np.ones(N-3)
    
    # Solve eigenvalue problem (tridiagonal matrix)
    energies, wavefunctions_inner = eigh_tridiagonal(diagonal, off_diagonal, 
                                                      select='i', 
                                                      select_range=(0, num_states-1))
    
    # Add boundary conditions (ψ=0 at boundaries)
    wavefunctions = np.zeros((N, num_states))
    wavefunctions[1:-1, :] = wavefunctions_inner
    
    # Normalize wave functions
    for i in range(num_states):
        norm = np.sqrt(np.trapz(wavefunctions[:, i]**2, x))
        wavefunctions[:, i] /= norm
    
    return x, energies, wavefunctions

# ============================================================================
# TERMINAL OUTPUT - CALCULATIONS
# ============================================================================

print("\n" + "=" * 80)
print("SCHRÖDINGER EQUATION: Solutions and Applications")
print("=" * 80)

# Part 1: The Schrödinger Equation
print("\n" + "─" * 80)
print("PART 1: THE SCHRÖDINGER EQUATION")
print("─" * 80)
print()
print("Time-Dependent Schrödinger Equation:")
print("    iℏ ∂Ψ/∂t = ĤΨ")
print("    where Ĥ = -ℏ²/2m ∇² + V(r,t)")
print()
print("Time-Independent Schrödinger Equation:")
print("    ĤΨ = EΨ")
print("    -ℏ²/2m ∇²ψ + V(r)ψ = Eψ")
print()
print(f"Physical Constants:")
print(f"  • Planck constant (h):     {h:.4e} J·s")
print(f"  • Reduced Planck (ℏ):      {hbar:.4e} J·s")
print(f"  • Electron mass (mₑ):      {m_e:.4e} kg")
print(f"  • Elementary charge (e):   {e:.4e} C")
print(f"  • Bohr radius (a₀):        {a0:.4e} m")

# Part 2: Particle in a 1D Box
print("\n" + "─" * 80)
print("PART 2: PARTICLE IN A 1D BOX (Infinite Square Well)")
print("─" * 80)
print()
print("System: Particle confined to box of length L")
print("Potential: V(x) = 0 for 0 < x < L, V(x) = ∞ elsewhere")
print()
print("Wave functions:  ψₙ(x) = √(2/L) sin(nπx/L)")
print("Energy levels:   Eₙ = n²h²/(8mL²)")
print()

L_box = 1e-9  # Box length: 1 nanometer

print(f"Box length L = {L_box*1e9:.1f} nm")
print()
print(f"{'n':<6} {'Energy (J)':<18} {'Energy (eV)':<18} {'Energy/E₁':<12}")
print("─" * 80)

for n in range(1, 6):
    E_J = particle_in_box_energy(n, L_box)
    E_eV = E_J / eV_to_J
    E_ratio = n**2
    print(f"{n:<6} {E_J:<18.4e} {E_eV:<18.4f} {E_ratio:<12}")

print()
print("Key observations:")
print("  • Quantized energy levels (discrete spectrum)")
print("  • Zero-point energy: E₁ ≠ 0 (quantum effect)")
print(f"  • Ground state energy: E₁ = {particle_in_box_energy(1, L_box)/eV_to_J:.4f} eV")
print("  • Energy scales as n²")
print("  • Number of nodes = n - 1")

# Part 3: Quantum Harmonic Oscillator
print("\n" + "─" * 80)
print("PART 3: QUANTUM HARMONIC OSCILLATOR")
print("─" * 80)
print()
print("System: Particle in parabolic potential well")
print("Potential: V(x) = ½mω²x²")
print()
print("Wave functions:  ψₙ(x) = Nₙ Hₙ(ξ) exp(-ξ²/2)")
print("                 where ξ = √(mω/ℏ)x, Hₙ = Hermite polynomial")
print("Energy levels:   Eₙ = ℏω(n + ½)")
print()

omega = 1e15  # Angular frequency (rad/s)
print(f"Angular frequency ω = {omega:.2e} rad/s")
print(f"Frequency ν = ω/2π = {omega/(2*np.pi):.2e} Hz")
print()
print(f"{'n':<6} {'Energy (J)':<18} {'Energy (eV)':<18} {'E - E₀ (eV)':<15}")
print("─" * 80)

E0 = harmonic_oscillator_energy(0, omega)
for n in range(6):
    E_J = harmonic_oscillator_energy(n, omega)
    E_eV = E_J / eV_to_J
    E_diff = (E_J - E0) / eV_to_J
    print(f"{n:<6} {E_J:<18.4e} {E_eV:<18.4f} {E_diff:<15.4f}")

print()
print("Key observations:")
print("  • Equally spaced energy levels (ΔE = ℏω)")
print(f"  • Zero-point energy: E₀ = ½ℏω = {E0/eV_to_J:.4f} eV")
print("  • Models molecular vibrations")
print("  • Wave functions are Gaussian × Hermite polynomials")

# Part 4: Hydrogen Atom
print("\n" + "─" * 80)
print("PART 4: HYDROGEN ATOM")
print("─" * 80)
print()
print("System: Electron in Coulomb potential of proton")
print("Potential: V(r) = -e²/(4πε₀r)")
print()
print("Energy levels:   Eₙ = -13.6 eV / n²")
print("Wave functions:  ψₙₗₘ(r,θ,φ) = Rₙₗ(r) Yₗₘ(θ,φ)")
print()
print(f"{'n':<6} {'Energy (eV)':<15} {'Orbital(s)':<30} {'Degeneracy':<12}")
print("─" * 80)

orbital_names = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g'}
for n in range(1, 6):
    E = hydrogen_energy(n)
    orbitals = [f"{n}{orbital_names[l]}" for l in range(n)]
    degeneracy = n**2
    print(f"{n:<6} {E:<15.4f} {', '.join(orbitals):<30} {degeneracy:<12}")

print()
print("Key observations:")
print("  • Energy depends only on n (for hydrogen)")
print("  • Ground state: n=1, E₁ = -13.6 eV")
print("  • Ionization energy = 13.6 eV")
print("  • Degeneracy = n² (multiple states with same energy)")
print(f"  • Bohr radius a₀ = {a0*1e10:.4f} Å")

# Part 5: Quantum Tunneling
print("\n" + "─" * 80)
print("PART 5: QUANTUM TUNNELING")
print("─" * 80)
print()
print("Transmission through rectangular barrier:")
print("T ≈ exp(-2κa) where κ = √(2m(V-E)/ℏ²)")
print()

# Example: electron tunneling
V_barrier = 5 * eV_to_J  # 5 eV barrier
a_barrier = 1e-9  # 1 nm width

print(f"Barrier height V = 5.0 eV")
print(f"Barrier width a = {a_barrier*1e9:.1f} nm")
print()
print(f"{'E (eV)':<10} {'E/V':<10} {'Transmission':<20} {'Description':<20}")
print("─" * 80)

energies_tunnel = [0.5, 1.0, 2.0, 3.0, 4.0, 4.5]
for E_eV in energies_tunnel:
    E_J = E_eV * eV_to_J
    E_ratio = E_eV / 5.0
    T = transmission_coefficient(E_J, V_barrier, a_barrier)
    description = "Strong tunneling" if T > 0.1 else ("Weak tunneling" if T > 1e-10 else "Very weak")
    print(f"{E_eV:<10.1f} {E_ratio:<10.2f} {T:<20.4e} {description:<20}")

print()
print("Key observations:")
print("  • Tunneling probability decreases exponentially with barrier width")
print("  • Lower energy → lower transmission (more 'forbidden')")
print("  • Classically impossible (E < V) but quantum mechanically allowed")
print("  • Applications: STM, alpha decay, tunnel diodes")

# Part 6: Uncertainty Principle
print("\n" + "─" * 80)
print("PART 6: HEISENBERG UNCERTAINTY PRINCIPLE")
print("─" * 80)
print()
print("Fundamental limits on measurement precision:")
print()
print("Position-Momentum:  Δx · Δp ≥ ℏ/2")
print("Energy-Time:        ΔE · Δt ≥ ℏ/2")
print()
print(f"Minimum uncertainty product: ℏ/2 = {hbar/2:.4e} J·s")
print()
print("Examples:")
print()

# Example 1: Electron localized to atomic size
Delta_x = 1e-10  # 1 Angstrom
Delta_p_min = hbar / (2 * Delta_x)
Delta_v = Delta_p_min / m_e
print(f"1. Electron localized to Δx = {Delta_x*1e10:.1f} Å:")
print(f"   Minimum momentum uncertainty: Δp ≥ {Delta_p_min:.4e} kg·m/s")
print(f"   Minimum velocity uncertainty: Δv ≥ {Delta_v:.4e} m/s")
print(f"   Minimum kinetic energy: ΔKE ≥ {(Delta_p_min**2)/(2*m_e)/eV_to_J:.4f} eV")

print()

# Example 2: Energy-time for atomic transition
Delta_E = 1 * eV_to_J  # 1 eV transition
Delta_t_min = hbar / (2 * Delta_E)
print(f"2. Atomic transition with ΔE = 1.0 eV:")
print(f"   Minimum time uncertainty: Δt ≥ {Delta_t_min:.4e} s")
print(f"   Natural line width: ΔE ≥ {hbar/(2*Delta_t_min)/eV_to_J:.4e} eV")

print("\n" + "=" * 80)
print("Generating visualizations...")
print("Window 1: Classical Quantum Systems (Particle in Box, Harmonic Oscillator)")
print("Window 2: Hydrogen Atom and Quantum Tunneling")
print("=" * 80)

# ============================================================================
# VISUALIZATIONS - WINDOW 1
# ============================================================================

fig1 = plt.figure(figsize=(18, 12))
gs1 = GridSpec(2, 3, figure=fig1, hspace=0.35, wspace=0.35)

# ============================================================================
# PLOT 1: Particle in a Box - Wave Functions
# ============================================================================

ax1 = fig1.add_subplot(gs1[0, 0])

x_box = np.linspace(0, L_box, 1000)
colors_box = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']

for n in range(1, 6):
    psi = particle_in_box_wavefunction(x_box, n, L_box)
    # Offset by energy level for visualization
    E_n = particle_in_box_energy(n, L_box) / eV_to_J
    ax1.plot(x_box*1e9, psi/np.max(np.abs(psi)) + n, 
            color=colors_box[n-1], linewidth=2, label=f'n={n}')
    ax1.axhline(y=n, color=colors_box[n-1], linestyle='--', alpha=0.3, linewidth=1)

ax1.set_xlabel('Position x (nm)', fontsize=11, fontweight='bold')
ax1.set_ylabel('ψₙ(x) + n (offset)', fontsize=11, fontweight='bold')
ax1.set_title('Particle in a Box\nWave Functions (n=1 to 5)', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, L_box*1e9)

# ============================================================================
# PLOT 2: Particle in a Box - Probability Densities
# ============================================================================

ax2 = fig1.add_subplot(gs1[0, 1])

for n in range(1, 6):
    psi = particle_in_box_wavefunction(x_box, n, L_box)
    prob = psi**2
    # Offset by energy level
    E_n = particle_in_box_energy(n, L_box) / eV_to_J
    ax2.fill_between(x_box*1e9, n, prob/np.max(prob)*0.8 + n, 
                     color=colors_box[n-1], alpha=0.6, label=f'n={n}')
    ax2.axhline(y=n, color=colors_box[n-1], linestyle='--', alpha=0.3, linewidth=1)

ax2.set_xlabel('Position x (nm)', fontsize=11, fontweight='bold')
ax2.set_ylabel('|ψₙ(x)|² + n (offset)', fontsize=11, fontweight='bold')
ax2.set_title('Particle in a Box\nProbability Densities', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9, loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, L_box*1e9)

# ============================================================================
# PLOT 3: Particle in a Box - Energy Levels
# ============================================================================

ax3 = fig1.add_subplot(gs1[0, 2])

n_levels = np.arange(1, 9)
E_levels = [particle_in_box_energy(n, L_box)/eV_to_J for n in n_levels]

for i, (n, E) in enumerate(zip(n_levels, E_levels)):
    color = colors_box[i] if i < len(colors_box) else '#95a5a6'
    ax3.hlines(E, n-0.3, n+0.3, colors=color, linewidth=4)
    ax3.text(n+0.4, E, f'E_{n} = {E:.1f} eV', fontsize=8, va='center')

ax3.set_xlabel('Quantum Number n', fontsize=11, fontweight='bold')
ax3.set_ylabel('Energy (eV)', fontsize=11, fontweight='bold')
ax3.set_title('Particle in a Box\nEnergy Levels (Eₙ ∝ n²)', fontsize=12, fontweight='bold')
ax3.set_xticks(n_levels)
ax3.grid(True, alpha=0.3, axis='y')

# ============================================================================
# PLOT 4: Harmonic Oscillator - Wave Functions
# ============================================================================

ax4 = fig1.add_subplot(gs1[1, 0])

x_ho = np.linspace(-5e-10, 5e-10, 1000)

for n in range(5):
    psi = harmonic_oscillator_wavefunction(x_ho, n, omega)
    E_n = harmonic_oscillator_energy(n, omega) / eV_to_J
    # Normalize for plotting
    psi_norm = psi / np.max(np.abs(psi)) * 0.4
    ax4.plot(x_ho*1e10, psi_norm + E_n, color=colors_box[n], linewidth=2, label=f'n={n}')
    ax4.axhline(y=E_n, color=colors_box[n], linestyle='--', alpha=0.3, linewidth=1)

# Add potential
V_ho = 0.5 * m_e * omega**2 * x_ho**2 / eV_to_J
V_scale = V_ho / np.max(V_ho) * np.max(E_levels[:5])
ax4.plot(x_ho*1e10, V_scale, 'k--', linewidth=2, alpha=0.5, label='V(x)')

ax4.set_xlabel('Position x (Å)', fontsize=11, fontweight='bold')
ax4.set_ylabel('Energy (eV)', fontsize=11, fontweight='bold')
ax4.set_title('Harmonic Oscillator\nWave Functions', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9, loc='upper right')
ax4.grid(True, alpha=0.3)

# ============================================================================
# PLOT 5: Harmonic Oscillator - Probability Densities
# ============================================================================

ax5 = fig1.add_subplot(gs1[1, 1])

for n in range(5):
    psi = harmonic_oscillator_wavefunction(x_ho, n, omega)
    prob = psi**2
    E_n = harmonic_oscillator_energy(n, omega) / eV_to_J
    # Normalize for plotting
    prob_norm = prob / np.max(prob) * 0.4
    ax5.fill_between(x_ho*1e10, E_n, prob_norm + E_n, 
                     color=colors_box[n], alpha=0.6, label=f'n={n}')
    ax5.axhline(y=E_n, color=colors_box[n], linestyle='--', alpha=0.3, linewidth=1)

# Add potential
ax5.plot(x_ho*1e10, V_scale, 'k--', linewidth=2, alpha=0.5, label='V(x)')

ax5.set_xlabel('Position x (Å)', fontsize=11, fontweight='bold')
ax5.set_ylabel('Energy (eV)', fontsize=11, fontweight='bold')
ax5.set_title('Harmonic Oscillator\nProbability Densities', fontsize=12, fontweight='bold')
ax5.legend(fontsize=9, loc='upper right')
ax5.grid(True, alpha=0.3)

# ============================================================================
# PLOT 6: Harmonic Oscillator - Energy Levels
# ============================================================================

ax6 = fig1.add_subplot(gs1[1, 2])

n_ho = np.arange(0, 8)
E_ho = [harmonic_oscillator_energy(n, omega)/eV_to_J for n in n_ho]
E_spacing = hbar * omega / eV_to_J

for i, (n, E) in enumerate(zip(n_ho, E_ho)):
    color = colors_box[i] if i < len(colors_box) else '#95a5a6'
    ax6.hlines(E, n-0.3, n+0.3, colors=color, linewidth=4)
    ax6.text(n+0.4, E, f'n={n}', fontsize=8, va='center')

ax6.set_xlabel('Quantum Number n', fontsize=11, fontweight='bold')
ax6.set_ylabel('Energy (eV)', fontsize=11, fontweight='bold')
ax6.set_title(f'Harmonic Oscillator\nEqually Spaced Levels (ΔE = {E_spacing:.4f} eV)', 
             fontsize=12, fontweight='bold')
ax6.set_xticks(n_ho)
ax6.grid(True, alpha=0.3, axis='y')

# Add arrow showing spacing
if len(E_ho) >= 2:
    ax6.annotate('', xy=(7.5, E_ho[1]), xytext=(7.5, E_ho[0]),
                arrowprops=dict(arrowstyle='<->', color='red', lw=2))
    ax6.text(7.7, (E_ho[0]+E_ho[1])/2, 'ℏω', fontsize=10, color='red', fontweight='bold')

plt.suptitle('Schrödinger Equation: Part 1 - Particle in Box & Harmonic Oscillator', 
            fontsize=16, fontweight='bold', y=0.995)

plt.show(block=False)

# ============================================================================
# VISUALIZATIONS - WINDOW 2
# ============================================================================

fig2 = plt.figure(figsize=(18, 12))
gs2 = GridSpec(2, 3, figure=fig2, hspace=0.35, wspace=0.35)

# ============================================================================
# PLOT 7: Hydrogen Atom - Radial Wave Functions
# ============================================================================

ax7 = fig2.add_subplot(gs2[0, 0])

r = np.linspace(0.01, 20, 500)  # in units of a0

orbitals_H = [(1, 0, '1s'), (2, 0, '2s'), (2, 1, '2p'), (3, 0, '3s'), (3, 1, '3p')]
colors_H = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6']

for (n, l, label), color in zip(orbitals_H, colors_H):
    R = hydrogen_radial_wavefunction(r, n, l)
    ax7.plot(r, R, color=color, linewidth=2, label=label)

ax7.set_xlabel('Radial Distance (a₀)', fontsize=11, fontweight='bold')
ax7.set_ylabel('Radial Wave Function R(r)', fontsize=11, fontweight='bold')
ax7.set_title('Hydrogen Atom\nRadial Wave Functions', fontsize=12, fontweight='bold')
ax7.legend(fontsize=10, loc='upper right')
ax7.grid(True, alpha=0.3)
ax7.set_xlim(0, 20)
ax7.axhline(0, color='black', linewidth=0.5)

# ============================================================================
# PLOT 8: Hydrogen Atom - Radial Probability Density
# ============================================================================

ax8 = fig2.add_subplot(gs2[0, 1])

for (n, l, label), color in zip(orbitals_H, colors_H):
    R = hydrogen_radial_wavefunction(r, n, l)
    P = r**2 * R**2  # Radial probability density
    ax8.plot(r, P, color=color, linewidth=2, label=label)

# Mark most probable radius for 1s
r_max_1s = 1.0  # a0 for 1s orbital
ax8.axvline(r_max_1s, color='#e74c3c', linestyle='--', alpha=0.5, linewidth=1.5)
ax8.text(r_max_1s + 0.3, ax8.get_ylim()[1]*0.9, 'r = a₀\n(1s max)', 
        fontsize=8, color='#e74c3c')

ax8.set_xlabel('Radial Distance (a₀)', fontsize=11, fontweight='bold')
ax8.set_ylabel('Radial Probability Density r²R²(r)', fontsize=11, fontweight='bold')
ax8.set_title('Hydrogen Atom\nRadial Probability Densities', fontsize=12, fontweight='bold')
ax8.legend(fontsize=10, loc='upper right')
ax8.grid(True, alpha=0.3)
ax8.set_xlim(0, 20)

# ============================================================================
# PLOT 9: Hydrogen Atom - Energy Levels
# ============================================================================

ax9 = fig2.add_subplot(gs2[0, 2])

n_H = np.arange(1, 7)
E_H = [hydrogen_energy(n) for n in n_H]

orbital_labels = {
    1: '1s',
    2: '2s, 2p',
    3: '3s, 3p, 3d',
    4: '4s, 4p, 4d, 4f',
    5: '5s, 5p, 5d, 5f',
    6: '6s, 6p, 6d, 6f'
}

for i, (n, E) in enumerate(zip(n_H, E_H)):
    color = colors_H[i] if i < len(colors_H) else '#95a5a6'
    ax9.hlines(E, 0.5, 1.5, colors=color, linewidth=4)
    ax9.text(1.6, E, f'n={n}: {orbital_labels[int(n)]}', fontsize=9, va='center')
    ax9.text(0.3, E, f'{E:.2f} eV', fontsize=8, va='center', ha='right')

ax9.axhline(0, color='red', linestyle='--', linewidth=2, alpha=0.5)
ax9.text(1.6, 0.5, 'Ionization (E=0)', fontsize=9, color='red', fontweight='bold')

ax9.set_xlim(0, 3)
ax9.set_ylabel('Energy (eV)', fontsize=11, fontweight='bold')
ax9.set_title('Hydrogen Atom\nEnergy Levels (Eₙ = -13.6/n² eV)', fontsize=12, fontweight='bold')
ax9.set_xticks([])
ax9.grid(True, alpha=0.3, axis='y')

# ============================================================================
# PLOT 10: Quantum Tunneling - Barrier Penetration
# ============================================================================

ax10 = fig2.add_subplot(gs2[1, 0])

x_tunnel = np.linspace(0, 3e-9, 1000)
V_barrier = 5 * eV_to_J
a_barrier = 1e-9
barrier_start = 1e-9
barrier_end = 2e-9

# Potential
V_plot = np.zeros_like(x_tunnel)
barrier_mask = (x_tunnel >= barrier_start) & (x_tunnel <= barrier_end)
V_plot[barrier_mask] = V_barrier / eV_to_J

ax10.fill_between(x_tunnel*1e9, 0, V_plot, alpha=0.3, color='gray', label='Barrier')

# Plot wave function for different energies
E_values = [2*eV_to_J, 3*eV_to_J, 4*eV_to_J]
colors_tunnel = ['#e74c3c', '#f39c12', '#2ecc71']
labels_tunnel = ['E = 2 eV', 'E = 3 eV', 'E = 4 eV']

for E, color, label in zip(E_values, colors_tunnel, labels_tunnel):
    # Simplified wave function (qualitative)
    k = np.sqrt(2 * m_e * E) / hbar
    psi = np.sin(k * x_tunnel)
    
    # Decay in barrier
    kappa = np.sqrt(2 * m_e * (V_barrier - E)) / hbar
    psi[barrier_mask] = psi[barrier_mask] * np.exp(-kappa * (x_tunnel[barrier_mask] - barrier_start))
    
    # Normalize
    psi = psi / np.max(np.abs(psi)) * 3
    
    ax10.plot(x_tunnel*1e9, psi + E/eV_to_J, color=color, linewidth=2, label=label)

ax10.set_xlabel('Position x (nm)', fontsize=11, fontweight='bold')
ax10.set_ylabel('Energy (eV) / Wave Function', fontsize=11, fontweight='bold')
ax10.set_title('Quantum Tunneling\nWave Function Penetration', fontsize=12, fontweight='bold')
ax10.legend(fontsize=9)
ax10.grid(True, alpha=0.3)

# ============================================================================
# PLOT 11: Transmission Coefficient vs Energy
# ============================================================================

ax11 = fig2.add_subplot(gs2[1, 1])

E_range = np.linspace(0.1, 5.5, 100) * eV_to_J
T_values = [transmission_coefficient(E, V_barrier, a_barrier) for E in E_range]

ax11.semilogy(E_range/eV_to_J, T_values, 'b-', linewidth=2.5)
ax11.axvline(V_barrier/eV_to_J, color='red', linestyle='--', linewidth=2, 
            label=f'Barrier height ({V_barrier/eV_to_J:.0f} eV)')
ax11.axhline(1, color='green', linestyle='--', linewidth=1, alpha=0.5, label='T = 1 (classical)')

ax11.set_xlabel('Particle Energy (eV)', fontsize=11, fontweight='bold')
ax11.set_ylabel('Transmission Coefficient T', fontsize=11, fontweight='bold')
ax11.set_title(f'Quantum Tunneling\nTransmission vs Energy (a = {a_barrier*1e9:.1f} nm)', 
              fontsize=12, fontweight='bold')
ax11.legend(fontsize=9)
ax11.grid(True, alpha=0.3, which='both')
ax11.set_ylim(1e-50, 10)

# ============================================================================
# PLOT 12: Uncertainty Principle Illustration
# ============================================================================

ax12 = fig2.add_subplot(gs2[1, 2])

# Gaussian wave packet
x_unc = np.linspace(-5, 5, 1000)
sigma_values = [0.5, 1.0, 2.0]
colors_unc = ['#e74c3c', '#3498db', '#2ecc71']
labels_unc = ['Narrow (σ=0.5)', 'Medium (σ=1.0)', 'Wide (σ=2.0)']

for sigma, color, label in zip(sigma_values, colors_unc, labels_unc):
    psi_gauss = np.exp(-x_unc**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))
    ax12.plot(x_unc, psi_gauss, color=color, linewidth=2.5, label=label)

ax12.set_xlabel('Position Δx', fontsize=11, fontweight='bold')
ax12.set_ylabel('Wave Function |ψ(x)|', fontsize=11, fontweight='bold')
ax12.set_title('Heisenberg Uncertainty Principle\nΔx · Δp ≥ ℏ/2', fontsize=12, fontweight='bold')
ax12.legend(fontsize=9)
ax12.grid(True, alpha=0.3)

# Add text annotation
ax12.text(0.5, 0.95, 'Narrower position spread\n→ Larger momentum spread', 
         transform=ax12.transAxes, fontsize=9, ha='center', va='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle('Schrödinger Equation: Part 2 - Hydrogen Atom & Quantum Tunneling', 
            fontsize=16, fontweight='bold', y=0.995)

plt.show(block=True)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("KEY INSIGHTS - SCHRÖDINGER EQUATION")
print("=" * 80)
print()
print("FUNDAMENTAL CONCEPTS:")
print("  • Wave function ψ describes quantum state (probability amplitude)")
print("  • |ψ|² gives probability density of finding particle")
print("  • Energy levels are quantized (discrete values)")
print("  • Zero-point energy: ground state E ≠ 0 (quantum effect)")
print("  • Wave-particle duality: particles behave as waves")
print()
print("PARTICLE IN A BOX:")
print("  • Quantized energies: Eₙ = n²h²/(8mL²)")
print("  • Standing wave solutions with n-1 nodes")
print("  • Models electrons in molecules, quantum dots")
print()
print("HARMONIC OSCILLATOR:")
print("  • Equally spaced levels: Eₙ = ℏω(n + ½)")
print("  • Zero-point energy E₀ = ½ℏω")
print("  • Models molecular vibrations, phonons")
print()
print("HYDROGEN ATOM:")
print("  • Energy levels: Eₙ = -13.6 eV/n²")
print("  • Three quantum numbers (n, l, m) from 3D solution")
print("  • Foundation for understanding atomic structure")
print()
print("QUANTUM TUNNELING:")
print("  • Particles penetrate classically forbidden regions")
print("  • Transmission T ≈ exp(-2κa) for E < V")
print("  • Essential for: STM, alpha decay, tunnel diodes")
print()
print("UNCERTAINTY PRINCIPLE:")
print("  • Δx · Δp ≥ ℏ/2 (position-momentum)")
print("  • ΔE · Δt ≥ ℏ/2 (energy-time)")
print("  • Fundamental limit, not measurement error")
print()
print("APPLICATIONS:")
print("  • Chemistry: bonding, spectroscopy, reactions")
print("  • Technology: semiconductors, lasers, quantum computing")
print("  • Materials: electronic properties, superconductivity")
print("  • Medicine: MRI, PET scans, radiation therapy")
print()
print("=" * 80)
print("Visualization complete!")
print("=" * 80)
