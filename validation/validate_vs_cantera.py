#!/usr/bin/env python3
"""
Validate the kinetic-theory transport formulas (Wang 2025 / Kee 2003)
against Cantera's mixture-averaged transport.

Re-implements the paper's equations in Python from scratch, using the SAME
molecular parameters as Cantera (converted to CHEMKIN convention used by the
OpenFOAM standardKineticTransport class), and compares to Cantera.
"""
import numpy as np
import cantera as ct

# Physical constants (SI)
kB = 1.3806488e-23
NA = 6.022140857e23
R = 8.3144621        # J/mol/K
PI = np.pi
DEBYE = 3.33564e-30  # C*m per Debye

# Collision integral coefficients (Kee Table 12.1)
a = np.array([1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5])  # Omega(1,1)*
b = np.array([1.0413, 0.11930, 0.43628, 1.6041, 0.095661, 2.0])  # Omega(2,2)*

GEO = {"atom": 0, "linear": 1, "nonlinear": 2}


def omega(coeff, Tstar, deltastar):
    o = coeff[0]*Tstar**(-coeff[1]) + (Tstar + coeff[2])**(-coeff[3])
    f = 1.0 + (np.exp(coeff[4]/Tstar) - np.exp(-coeff[5]/Tstar))*deltastar**2 \
        / (2.0 + 2.5*deltastar)
    return o*f


class Species:
    """Molecular transport parameters in CHEMKIN convention (K, A, Debye, A^3)."""
    def __init__(self, geometry, welldepth, diameter, dipole, polar, relax, W):
        self.geometry = int(geometry)
        self.welldepth = float(welldepth)   # epsilon/kB [K]
        self.diameter = float(diameter)     # sigma [Angstrom]
        self.dipole = float(dipole)         # [Debye]
        self.polar = float(polar)           # [A^3]
        self.relax = float(relax)
        self.W = float(W)                   # [kg/kmol]

    # SI accessors (mirror the OpenFOAM standardKineticTransport conversions)
    def sigma(self):   return 1e-10*self.diameter
    def eps(self):     return kB*self.welldepth
    def mud(self):     return 1e-21*10**(-7.0/2.0)*self.dipole   # [m^3/2 J^1/2]
    def pol(self):     return 1e-30*self.polar                    # [m^3]


def mu_k(s, T):
    sigma = s.sigma()
    eps = s.eps()
    Tstar = T/s.welldepth
    deltastar = 0.5*s.mud()**2/(eps*sigma**3)
    return (5.0/16.0*np.sqrt((s.W/1000.0)/NA*kB*T/PI)
            / (sigma**2 * omega(b, Tstar, deltastar)))


def Dkk(s, p, T):
    sigma = s.sigma()
    eps = s.eps()
    Tstar = T/s.welldepth
    deltastar = 0.5*s.mud()**2/(eps*sigma**3)
    return (3.0/8.0*np.sqrt(NA*(kB*T)**3/(PI*(s.W/1000.0)))
            / (p*sigma**2*omega(a, Tstar, deltastar)))


def D_jk(s1, s2, p, T):
    Wkj = s1.W*s2.W/(s1.W + s2.W)
    epskj = np.sqrt(s1.eps()*s2.eps())
    sigkj = 0.5*(s1.sigma() + s2.sigma())
    muk = s1.mud()
    muj = s2.mud()
    polk = s1.pol()
    polj = s2.pol()

    deltastar = 0.0
    # Stockmayer polar/nonpolar correction (1:1 port of of-11)
    if polk == 0.0:
        if polj == 0.0:
            deltastar = 0.5*muk*muj/(epskj*sigkj**3)
        else:
            alphajstar = polj/s2.sigma()**3
            miukstar = muk/np.sqrt(s1.eps()*s1.sigma()**3)
            ksi = 1.0 + 0.25*alphajstar*miukstar*np.sqrt(s1.eps()/s2.eps())
            epskj *= ksi**2
            sigkj *= ksi**(-1.0/6.0)
            deltastar = 0.5*muk*muj/(epskj*sigkj**3)
    else:
        if polj == 0.0:
            alphakstar = polk/s1.sigma()**3
            miujstar = muj/np.sqrt(s2.eps()*s2.sigma()**3)
            ksi = 1.0 + 0.25*alphakstar*miujstar*np.sqrt(s2.eps()/s1.eps())
            epskj *= ksi**2
            sigkj *= ksi**(-1.0/6.0)
            deltastar = 0.5*muk*muj/(epskj*sigkj**3)
        else:
            deltastar = 0.0

    Tstar = T*kB/epskj
    o = omega(a, Tstar, deltastar)
    return 3.0/16.0*np.sqrt(2.0*NA*(kB*T)**3/(PI*Wkj/1000.0)) \
        / (p*sigkj**2*o)


def F(s, T):
    w = s.welldepth
    return (1.0 + 0.5*PI**1.5*np.sqrt(w/T) + (PI**2/4.0 + 2.0)*(w/T)
            + PI**1.5*(w/T)**1.5)


def Zrot(s, T):
    return s.relax*F(s, 298.0)/F(s, T)


def kappa_k(s, p, T, cv_molar):
    # cv_molar [J/mol/K]
    Cv_trans = 1.5*R
    Cv_rot = {0: 0.0, 1: R, 2: 1.5*R}[s.geometry]
    Cv_vib = cv_molar - Cv_trans - Cv_rot

    rho = p*(s.W/1000.0)/R/T            # ideal gas [kg/m^3]
    rhoDmu = rho*Dkk(s, p, T)/mu_k(s, T)

    A = 2.5 - rhoDmu
    B = Zrot(s, T) + 2.0/PI*(5.0/3.0*Cv_rot/R + rhoDmu)
    ftrans = 2.5*(1.0 - 2.0/PI*Cv_rot/Cv_trans*A/B)
    frot = rhoDmu*(1.0 + 2.0/PI*A/B)
    fvib = rhoDmu
    return mu_k(s, T)/(s.W/1000.0)*(ftrans*Cv_trans + frot*Cv_rot + fvib*Cv_vib)


def wilke_mu(species, X, p, T):
    N = len(species)
    muk = [mu_k(s, T) for s in species]
    mu = 0.0
    for k in range(N):
        denom = 0.0
        for j in range(N):
            phi = ((1.0 + np.sqrt(muk[k]/muk[j])
                    * (species[j].W/species[k].W)**0.25)**2
                   / np.sqrt(8.0*(1.0 + species[k].W/species[j].W)))
            denom += X[j]*phi
        mu += X[k]*muk[k]/denom
    return mu


def combination_kappa(species, X, p, T, cv_molar):
    lam = [kappa_k(s, p, T, cv_molar[i]) for i, s in enumerate(species)]
    s1 = sum(X[i]*lam[i] for i in range(len(species)))
    s2 = sum(X[i]/lam[i] for i in range(len(species)))
    return 0.5*(s1 + 1.0/s2)


def main():
    gas = ct.Solution('gri30.yaml')
    names = gas.species_names

    species = []
    for sp in gas.species():
        tr = sp.transport
        species.append(Species(
            GEO[tr.geometry],
            tr.well_depth/kB,          # J -> K
            tr.diameter*1e10,          # m -> Angstrom
            tr.dipole/DEBYE,           # C*m -> Debye
            tr.polarizability*1e30,    # m^3 -> A^3
            tr.rotational_relaxation,
            sp.molecular_weight
        ))

    states = [
        ("unburnt CH4/air T=300", dict(T=300.0, X="CH4:0.095, O2:0.19, N2:0.715")),
        ("burnt T=1800", dict(T=1800.0,
             X="CO2:0.07, H2O:0.14, N2:0.72, CO:0.02, H2:0.01, O2:0.04")),
    ]

    for name, st in states:
        gas.TPX = st["T"], ct.one_atm, st["X"]
        X, Y, T, p = gas.X, gas.Y, gas.T, gas.P

        # Cantera reference (mixture-averaged)
        ct_mu = gas.viscosity
        ct_kappa = gas.thermal_conductivity
        ct_Dmix = gas.mix_diff_coeffs

        cv_molar = [sp.thermo.cp(T)/1000.0 - R for sp in gas.species()]  # J/mol/K

        our_mu = wilke_mu(species, X, p, T)
        our_kappa = combination_kappa(species, X, p, T, cv_molar)

        our_Dmix = []
        for k in range(len(species)):
            denom = 0.0
            for j in range(len(species)):
                if j != k:
                    denom += X[j]/D_jk(species[k], species[j], p, T)
            our_Dmix.append((1.0 - Y[k])/denom if denom > 0 else 0.0)

        print(f"=== {name} ===")
        print(f"  mu     : ours={our_mu:.6e}  cantera={ct_mu:.6e}  "
              f"rel err={abs(our_mu-ct_mu)/ct_mu*100:.2f}%")
        print(f"  kappa  : ours={our_kappa:.6e}  cantera={ct_kappa:.6e}  "
              f"rel err={abs(our_kappa-ct_kappa)/ct_kappa*100:.2f}%")
        for nm in ['CH4', 'O2', 'H2O', 'CO2', 'H2', 'CO', 'N2']:
            if nm in names:
                k = names.index(nm)
                if ct_Dmix[k] > 0:
                    err = abs(our_Dmix[k]-ct_Dmix[k])/ct_Dmix[k]*100
                    print(f"  D_mix({nm:4s}): ours={our_Dmix[k]:.6e}  "
                          f"cantera={ct_Dmix[k]:.6e}  rel err={err:.2f}%")


if __name__ == '__main__':
    main()
