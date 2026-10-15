# kineticTheoryTransport

An OpenFOAM-v14 library implementing **detailed multicomponent transport**
based on the standard kinetic theory of gases (Chapman–Enskog /
Lennard–Jones + Stockmayer correction) for reacting flows.

The library provides **three species-diffusion closures**, from the cheap
mixture-averaged approximation up to the strict coupled Stefan–Maxwell
formulation:

| model | diffusion closure | cost |
|---|---|---|
| `kineticTheoryFickianFourier` | mixture-averaged, $\nabla Y_i$ only | low |
| `kineticTheoryHirschfelderCurtiss` | mixture-averaged + $\frac{Y_i}{W_m}\nabla W_m$ + correction velocity | medium |
| `kineticTheoryMaxwellStefanFourier` | strict multicomponent (Stefan–Maxwell → generalised Fick) | high |

A port of the `standardKinetic/of-11` fork to OpenFOAM-v14 as a **standalone
library** (no modification of the native OpenFOAM sources).

## Attribution and license

**The physics (equations) is classical kinetic theory of gases**, not our
invention:

- Chapman–Enskog for $\mu_k$ and $D_{jk}$: Hirschfelder, Curtiss & Bird (1954);
- collision integrals $\Omega^{(1,1)*}$, $\Omega^{(2,2)*}$ and their
  approximations: Monchick & Mason (1961), Neufeld et al. (1972);
- Wilke's rule for the mixture viscosity: Wilke (1950);
- combination averaging of the thermal conductivity: Mathur, Tondon & Saxena
  (1967);
- mixture-averaged and Stefan–Maxwell multicomponent diffusion:
  Taylor & Krishna (1993), *Multicomponent Mass Transfer*;
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

- pure-species $\mu_k$, $\kappa_k$, $D_{kk}$;
- binary $D_{jk}$ (Chapman–Enskog, collision integrals
  $\Omega^{(1,1)*}$, $\Omega^{(2,2)*}$ with the Stockmayer polar correction);
- mixture viscosity $\mu$ — Wilke's rule;
- mixture thermal conductivity $\kappa$ — combination averaging (Mathur et al.).

The three diffusion models then assemble the species mass fluxes $\mathbf{j}_i$
from those binary coefficients.

### Model 1 — `kineticTheoryFickianFourier` (mixture-averaged)

Mixture-averaged mass diffusivity of species $i$:

$$ D_{i,\mathrm{m}} = \frac{1 - Y_i}{W_m \sum_{j \ne i} Y_j/(W_j D_{ij})} = \frac{1 - Y_i}{\sum_{j \ne i} X_j/D_{ij}} $$

Species mass flux:

$$ \mathbf{j}_i = -\rho D_{i,\mathrm{m}} \nabla Y_i $$

### Model 2 — `kineticTheoryHirschfelderCurtiss`

Hirschfelder–Curtiss approximation of the diffusion velocity, including the
molecular-weight gradient term:

$$ \mathbf{V}_i = -D_{i,\mathrm{m}} \nabla Y_i - \frac{D_{i,\mathrm{m}} Y_i}{W_m} \nabla W_m $$

To enforce global mass conservation ($\sum_i Y_i \mathbf{V}_i = 0$) a
correction velocity is added:

$$ \mathbf{V}_c = \sum_k \left( D_{k,\mathrm{m}} \nabla Y_k + \frac{D_{k,\mathrm{m}} Y_k}{W_m} \nabla W_m \right) = -\sum_k \mathbf{V}_k $$

Species mass flux:

$$ \mathbf{j}_i = -\rho D_{i,\mathrm{m}} \nabla Y_i - \rho \frac{D_{i,\mathrm{m}} Y_i}{W_m} \nabla W_m + \rho \mathbf{V}_c Y_i $$

### Model 3 — `kineticTheoryMaxwellStefanFourier` (strict multicomponent)

The coupled Stefan–Maxwell relation between the mole-fraction gradients and the
diffusion velocities of all species:

$$ \nabla X_i = \sum_{j \ne i} \frac{X_i X_j}{D_{ij}} (\mathbf{V}_j - \mathbf{V}_i) $$

This is inverted (Taylor–Krishna transformation, an $(N-1)\times(N-1)$
sub-system excluding the default specie) to the generalised Fick form:

$$ \mathbf{j}_i = -\rho D_{ii}^{\mathrm{Fick}} \nabla Y_i - \sum_{j \ne i} \rho D_{ij}^{\mathrm{Fick}} \nabla Y_j $$

where $D_{ii}^{\mathrm{Fick}}$ is the implicit (diagonal) part and
$D_{ij}^{\mathrm{Fick}}$ the explicit (off-diagonal) part.

### Energy equation

All three models add the differential-diffusion enthalpy flux
$\sum_i \mathbf{j}_i h_i$ to the energy equation, so the Soret-free
multi-component enthalpy transport is closed consistently.

## Layout

```
kineticTheoryTransport/
├── Make/{files,options}
├── specie/transport/standardKineticTransport/   # species transport (mu_k, kappa_k, D_kk)
├── mixtures/kineticTheoryMixture/               # runtime Djk interface
├── mixtures/kineticTheoryMulticomponentMixture/ # mixture (mu/kappa/Djk)
├── thermophysicalTransport/
│   ├── kineticTheoryFickianFourier/             # model 1: mixture-averaged
│   ├── kineticTheoryHirschfelderCurtiss/        # model 2: + grad(W_m) + correction
│   └── kineticTheoryMaxwellStefanFourier/       # model 3: Stefan-Maxwell
├── chemistryModel/standardChemistryModels.C     # static standard chemistry instantiation
├── thermos.C                    # psi/rho thermo instantiation
├── thermophysicalTransport.C    # diffusion-model registration
└── validation/                  # test case + validation against Cantera
```

## Building

```bash
source /opt/openfoam14/etc/bashrc
cd kineticTheoryTransport
./Allwmake
```

`Allwmake` builds `libkineticTheoryTransport.so` and installs a user-level
`chemistryModel` template override (see below).

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

`constant/thermophysicalTransport` — pick one of the three diffusion models:

```cpp
laminar
{
    // model           kineticTheoryFickianFourier;          // mixture-averaged
    // model           kineticTheoryHirschfelderCurtiss;     // + grad(W_m) + correction velocity
    model           kineticTheoryMaxwellStefanFourier;     // strict Stefan-Maxwell
}
```

## Validation

`validation/validate_vs_cantera.py` checks the kinetic-theory equations against
Cantera (gri30, mixture-averaged). Result:

| quantity | T=300 K | T=1800 K |
|---|---|---|
| $\mu$ | 1.24% | 0.40% |
| $\kappa$ | 1.28% | 0.22% |
| $D_{k,\mathrm{mix}}$ | 0.14–1.35% | 0.01–0.28% |

The `validation/` test case is a 2D counterflow flame (GRI 3.0, 53 species,
325 reactions), run with `foamRun -libs ("libkineticTheoryTransport.so")`.

The three diffusion models are exercised on the DRM-19 (21 species) premixed
CH₄/air flame. For CH₄/air the mixture-averaged and Stefan–Maxwell closures
agree to within 1–2% (as expected); the difference becomes significant for
H₂/air or mixtures with heavy hydrocarbons.

## Notes for developers

The `standard` chemistry model is **statically instantiated** in
`chemistryModel/standardChemistryModels.C` for the
`standardKinetic<janaf<perfectGas<specie>>, sensibleEnthalpy>` package and
compiled into `libkineticTheoryTransport.so`, so cases do not create a
`dynamicCode/` directory.

If you change the thermo/transport package (a different EoS, energy form or
species transport) and a `dynamicCode/` directory reappears in a case, it means
the new combination is not covered by `standardChemistryModels.C` — add the
corresponding explicit instantiation there and rebuild.

