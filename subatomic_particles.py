"""
Subatomic Particles: Protons, Neutrons, and Electrons
=====================================================

This script provides a comprehensive exploration of the fundamental particles
that make up atoms: protons, neutrons, and electrons.

INTRODUCTION:
Atoms are not indivisible as once thought. They are composed of smaller
subatomic particles that determine the atom's properties and behavior.
Understanding these particles is fundamental to chemistry, physics, and
our understanding of matter itself.

THEORY:

1. HISTORICAL DISCOVERY:

Timeline of Subatomic Particle Discovery:
- 1897: J.J. Thomson discovers the ELECTRON using cathode ray experiments
  * Determined electron charge-to-mass ratio (e/m)
  * Proposed "plum pudding model" of atom
  
- 1909-1911: Ernest Rutherford discovers the PROTON
  * Gold foil experiment revealed dense, positive nucleus
  * Disproved Thomson's plum pudding model
  * Proposed nuclear model of atom
  
- 1932: James Chadwick discovers the NEUTRON
  * Explained isotopes (same element, different mass)
  * Completed the basic atomic model

2. FUNDAMENTAL PROPERTIES:

A. PROTON (p⁺):
   Symbol: p or p⁺
   Location: Nucleus (center of atom)
   Charge: +1 elementary charge (+1.602 × 10⁻¹⁹ C)
   Mass: 1.6726 × 10⁻²⁷ kg = 1.007276 u (atomic mass units)
   Relative mass: 1 (defined as reference)
   Spin: ½ (fermion)
   Classification: Baryon (made of 3 quarks: uud)
   
   Key Facts:
   - Number of protons = ATOMIC NUMBER (Z)
   - Defines the element (e.g., 6 protons = Carbon, always)
   - Determines positive charge of nucleus
   - Approximately 1836 times heavier than electron
   - Stable particle (does not decay)

B. NEUTRON (n⁰):
   Symbol: n or n⁰
   Location: Nucleus (with protons)
   Charge: 0 (electrically neutral)
   Mass: 1.6749 × 10⁻²⁷ kg = 1.008665 u
   Relative mass: 1.001 (slightly heavier than proton)
   Spin: ½ (fermion)
   Classification: Baryon (made of 3 quarks: udd)
   
   Key Facts:
   - Number of neutrons = Mass number (A) - Atomic number (Z)
   - Isotopes have different numbers of neutrons
   - Adds mass without adding charge
   - Provides nuclear stability (reduces proton-proton repulsion)
   - Free neutron is unstable (half-life ~10.2 minutes, β-decay)
   - Bound neutron in stable nucleus is stable

C. ELECTRON (e⁻):
   Symbol: e or e⁻
   Location: Electron cloud (orbitals around nucleus)
   Charge: -1 elementary charge (-1.602 × 10⁻¹⁹ C)
   Mass: 9.1094 × 10⁻³¹ kg = 0.0005486 u
   Relative mass: 1/1836 (0.000545 compared to proton)
   Spin: ½ (fermion)
   Classification: Lepton (elementary particle, not made of quarks)
   
   Key Facts:
   - Number of electrons = number of protons (in neutral atom)
   - Determines chemical properties and bonding behavior
   - Occupies >99.9% of atom's volume (electron cloud)
   - Contains <0.1% of atom's mass
   - Responsible for electrical conductivity
   - Stable particle (does not decay)
   - Quantum mechanical behavior (wave-particle duality)

3. ATOMIC STRUCTURE:

Nucleus (Nucleons):
- Contains protons and neutrons (collectively called nucleons)
- Extremely dense: ~10¹⁴ g/cm³
- Diameter: ~10⁻¹⁵ m (1 femtometer)
- Held together by STRONG NUCLEAR FORCE
  * Strongest of the four fundamental forces
  * Overcomes electromagnetic repulsion between protons
  * Short-range force (~10⁻¹⁵ m)

Electron Cloud:
- Electrons occupy orbitals (probability distributions)
- Atom diameter: ~10⁻¹⁰ m (100,000 times larger than nucleus)
- If nucleus were a marble, atom would be size of football stadium
- Electrons repel each other (electromagnetic force)
- Quantum mechanics governs electron behavior

4. MATHEMATICAL RELATIONSHIPS:

Atomic Number (Z):
   Z = number of protons = number of electrons (neutral atom)
   Example: Carbon has Z = 6 (6 protons, 6 electrons)

Mass Number (A):
   A = number of protons + number of neutrons
   Example: Carbon-12: A = 12 (6 protons + 6 neutrons)

Number of Neutrons (N):
   N = A - Z
   Example: Carbon-12: N = 12 - 6 = 6 neutrons

Atomic Mass (m):
   m ≈ A × u, where u = 1.66054 × 10⁻²⁷ kg (atomic mass unit)
   More precisely: m = (Z × m_p) + (N × m_n) + (Z × m_e) - Binding Energy
   Note: Binding energy causes mass defect (Einstein's E=mc²)

Charge of Atom:
   Net charge = (# protons) - (# electrons)
   Neutral atom: charge = 0
   Cation (positive ion): fewer electrons than protons
   Anion (negative ion): more electrons than protons

5. ISOTOPES:

Definition: Atoms of the same element with different numbers of neutrons
- Same Z (same element)
- Different N (different mass)
- Different A (different mass number)

Examples:
- Hydrogen isotopes:
  * Protium (¹H): 1p, 0n (99.985% abundant)
  * Deuterium (²H or D): 1p, 1n (0.015% abundant)
  * Tritium (³H or T): 1p, 2n (radioactive, trace amounts)

- Carbon isotopes:
  * Carbon-12 (¹²C): 6p, 6n (98.9% abundant, stable)
  * Carbon-13 (¹³C): 6p, 7n (1.1% abundant, stable)
  * Carbon-14 (¹⁴C): 6p, 8n (radioactive, used for dating)

Notation: ᴬX or X-A, where A is mass number
          ᶻ
Example: ¹²C or Carbon-12
         ₆

6. ATOMIC MASS UNIT (u or amu):

Definition: 1/12 the mass of a Carbon-12 atom
1 u = 1.66054 × 10⁻²⁷ kg

Comparison:
- Proton: 1.007276 u
- Neutron: 1.008665 u
- Electron: 0.0005486 u

Atomic mass ≈ sum of nucleons (protons + neutrons)
Electrons contribute negligible mass (<0.1%)

7. FORCES IN ATOMS:

Electromagnetic Force:
- Proton-proton: repulsion (both positive)
- Electron-electron: repulsion (both negative)
- Proton-electron: attraction (opposite charges)
- Magnitude: F = k(q₁q₂/r²), where k = 8.99 × 10⁹ N·m²/C²

Strong Nuclear Force:
- Binds protons and neutrons in nucleus
- Overcomes electromagnetic repulsion between protons
- Extremely strong but short-range (<10⁻¹⁵ m)
- Attractive between nucleons
- Mediated by gluons and pions

8. QUANTUM PROPERTIES:

Wave-Particle Duality:
- All particles exhibit both wave and particle properties
- de Broglie wavelength: λ = h/p = h/(mv)
  * h = Planck's constant = 6.626 × 10⁻³⁴ J·s
  * For electrons: significant wavelength, quantum effects important
  * For protons/neutrons: much smaller wavelength

Heisenberg Uncertainty Principle:
- Δx·Δp ≥ ℏ/2, where ℏ = h/(2π)
- Cannot simultaneously know exact position and momentum
- Explains why electrons cannot collapse into nucleus
- Explains electron orbital structure

Spin:
- All three particles have spin ½ (fermions)
- Obey Pauli exclusion principle
- Electron spin: basis for magnetism and electron pairing

9. APPLICATIONS AND IMPORTANCE:

Chemistry:
- Protons determine element identity
- Electrons determine chemical bonding and reactivity
- Isotopes used in research and medicine

Nuclear Physics:
- Understanding radioactivity and nuclear reactions
- Nuclear energy (fission and fusion)
- Nuclear medicine (PET scans, cancer treatment)

Technology:
- Electron microscopy
- Particle accelerators
- Semiconductors and electronics
- Mass spectrometry

Medicine:
- Radiation therapy
- Radioactive tracers
- MRI (nuclear magnetic resonance)

Author: ewdhp
Date: November 12, 2025
"""
