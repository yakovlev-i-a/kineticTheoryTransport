# kineticTheoryTransport

An OpenFOAM-v14 library implementing a **mixture-averaged** transport model
based on the standard kinetic theory of gases (Chapman–Enskog /
Lennard–Jones + Stockmayer correction) for reacting multicomponent flows.

A port of the `standardKinetic/of-11` fork to OpenFOAM-v14 as a **standalone
library** (no modification of the native OpenFOAM sources).

## Attribution and license

**The physics (equations) is classical kinetic theory of gases**, not our
invention:

- Chapman–Enskog for μ<sub>k</sub> and D<sub>jk</sub>: Hirschfelder, Curtiss & Bird (1954);
- collision integrals Ω<sup>(1,1)\*</sup>, Ω<sup>(2,2)\*</sup> and their
  approximations: Monchick & Mason (1961), Neufeld et al. (1972);
- Wilke's rule for the mixture viscosity: Wilke (1950);
- combination averaging of the thermal conductivity: Mathur, Tondon & Saxena
  (1967);
- a consolidated presentation of all the equations: **Kee, Coltrin & Glarborg,
  "Chemically Reacting Flow" (2003)**.

**The concrete implementation (class structure, collision-integral
coefficients, Stockmayer correction) is a port** of the GPL fork
`standardKinetic/of-11` from:

> Fei Wang, Liang An, Tat Leung Chan. *An OpenFOAM solver incorporating
> detailed transport model for reacting flow simulations*, Computer Physics
> Communications **309** (2025) 109472.
> https://doi.org/10.1016/j.cpc.2024.109472

The fork is distributed under **GPLv3** (full text in `LICENSE`). This library
is a derivative work and is likewise licensed under GPLv3. If you use the code
in a scientific work, please cite Wang et al. (2025) and Kee et al. (2003).

## What it computes

From the six CHEMKIN molecular parameters (`geometry`, `welldepth`, `diameter`,
`dipole`, `polar`, `relaxation`), the library computes:

- pure-species μ<sub>k</sub>, κ<sub>k</sub>, D<sub>kk</sub>;
- binary D<sub>jk</sub> (Chapman–Enskog, collision integrals
  Ω<sup>(1,1)\*</sup>, Ω<sup>(2,2)\*</sup> with the Stockmayer polar correction);
- mixture viscosity μ — Wilke's rule;
- mixture thermal conductivity κ — combination averaging (Mathur et al.);
- mixture-averaged D<sub>k,mix</sub> = (1 − Y<sub>k</sub>) / Σ<sub>j≠k</sub> (X<sub>j</sub> / D<sub>jk</sub>).

## Layout

```
kineticTheoryTransport/
├── Make/{files,options}
├── specie/transport/standardKineticTransport/   # species transport
├── mixtures/kineticTheoryMixture/               # runtime Djk interface
├── mixtures/kineticTheoryMulticomponentMixture/ # mixture (mu/kappa/Djk)
├── thermophysicalTransport/kineticTheoryFickianFourier/  # transport model
├── thermos.C                    # psi/rho thermo instantiation
├── thermophysicalTransport.C    # kineticTheoryFickianFourier registration
└── validation/                  # test case + validation against Cantera
```

## Building

```bash
source /opt/openfoam14/etc/bashrc
cd kineticTheoryTransport
./Allwmake
```

`Allwmake` builds `libkineticTheoryTransport.so` and installs the user-level
override of the `chemistryModel` template (see below).

## Using in a case

`constant/physicalProperties`:

```cpp
libs ("libkineticTheoryTransport.so");

thermoType
{
    type            hePsiThermo;
    mixture         kineticTheoryMulticomponentMixture;
    transport       standardKinetic;
    thermo          janaf;
    energy          sensibleEnthalpy;
    equationOfState perfectGas;
    specie          specie;
}
```

In `constant/thermo`, each species carries the kinetic parameters instead of
`transport { As ...; Ts ...; }`:

```cpp
transport
{
    geometry        1;      // 0 atom / 1 linear / 2 nonlinear
    welldepth       97.53;  // epsilon/kB [K]
    diameter        3.621;  // sigma [Angstrom]
    dipole          0.0;    // [Debye]
    polar           1.76;   // polarizability [A^3]
    relaxation      4.0;    // Zrot(298) [-]
}
```

`constant/thermophysicalTransport`:

```cpp
laminar
{
    model           kineticTheoryFickianFourier;
}
```

## Important detail: chemistryModel via dynamicCode

In OpenFOAM-v14 the `chemistryModel` class (standard stiff-ODE chemistry) is
templated on the thermoPhysics and compiled **on the fly** (dynamicCode). So
that this dynamic compilation can find `standardKineticTransport.H`,
`Allwmake` installs a user-level template override at

```
~/.OpenFOAM/14/codeTemplates/dynamicCode/chemistryModel
```

which adds `-I<lib>/lnInclude` to `EXE_INC`. On the first run of a case,
OpenFOAM compiles `standard<standardKinetic<...>>` once and caches it.

## Validation

`validation/validate_vs_cantera.py` checks the kinetic-theory equations against
Cantera (gri30, mixture-averaged). Result:

| quantity | T=300 K | T=1800 K |
|---|---|---|
| μ | 1.24% | 0.40% |
| κ | 1.28% | 0.22% |
| D<sub>k,mix</sub> | 0.14–1.35% | 0.01–0.28% |

The `validation/` test case is a 2D counterflow flame (GRI 3.0, 53 species,
325 reactions), run with `foamRun -libs ("libkineticTheoryTransport.so")`.
