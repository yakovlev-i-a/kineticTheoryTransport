/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 LZF / kineticTheoryTransport
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of a derivative work of OpenFOAM (GPLv3).

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kineticTheoryMulticomponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixtureType::
WilkeWeights
(
    scalar p,
    scalar T
) const
{
    forAll(mu_, i)
    {
        mu_[i] = specieThermos_[i].mu(p, T);
    }

    forAll(M_, i)
    {
        scalar sumXphi = 0;

        forAll(M_, j)
        {
            if (i != j)
            {
                const scalar phiij =
                    sqr(1 + sqrt((mu_[i]/mu_[j])*B_(i, j)))/A_(i, j);

                sumXphi += X_[j]*phiij;
            }
            else
            {
                sumXphi += X_[j];
            }
        }

        w_[i] = X_[i]/sumXphi;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::kineticTheoryMulticomponentMixture<ThermoType>::
kineticTheoryMulticomponentMixture
(
    const dictionary& dict
)
:
    multicomponentMixture<ThermoType>(dict),
    mixture_("mixture", this->specieThermos()[0]),
    transportMixture_(this->specieThermos())
{}


template<class ThermoType>
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixtureType::
transportMixtureType
(
    const PtrList<ThermoType>& specieThermos
)
:
    specieThermos_(specieThermos),
    M_(specieThermos.size()),
    A_(specieThermos.size()),
    B_(specieThermos.size()),
    X_(specieThermos.size()),
    mu_(specieThermos.size()),
    w_(specieThermos.size()),
    muCached_(false)
{
    forAll(specieThermos_, i)
    {
        M_[i] = specieThermos[i].W();
    }

    forAll(M_, i)
    {
        forAll(M_, j)
        {
            if (i != j)
            {
                A_(i, j) = ((4/sqrt(2.0))*sqrt(1 + M_[i]/M_[j]));
                B_(i, j) = sqrt(M_[j]/M_[i]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixtureType::
mu
(
    scalar p,
    scalar T
) const
{
    WilkeWeights(p, T);

    scalar mu = 0;
    forAll(w_, i)
    {
        mu += w_[i]*mu_[i];
    }

    return mu;
}


template<class ThermoType>
Foam::scalar
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixtureType::
kappa
(
    scalar p,
    scalar T
) const
{
    // Combination averaging (Mathur et al., Kee Eq. 12.119)
    scalar sum1 = 0;
    scalar sum2 = 0;

    forAll(X_, i)
    {
        const scalar kappai = specieThermos_[i].kappa(p, T);

        sum1 += X_[i]*kappai;
        sum2 += X_[i]/kappai;
    }

    return 0.5*(sum1 + 1/max(sum2, small));
}


template<class ThermoType>
const typename
Foam::kineticTheoryMulticomponentMixture<ThermoType>::thermoMixtureType&
Foam::kineticTheoryMulticomponentMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    mixture_ = Y[0]*this->specieThermos()[0];

    for (label i=1; i<Y.size(); i++)
    {
        mixture_ += Y[i]*this->specieThermos()[i];
    }

    return mixture_;
}


template<class ThermoType>
const typename
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    transportMixture_.muCached_ = false;

    scalar sumX = 0;

    forAll(Y, i)
    {
        transportMixture_.X_[i] = Y[i]/this->specieThermos()[i].W();
        sumX += transportMixture_.X_[i];
    }

    forAll(Y, i)
    {
        transportMixture_.X_[i] /= sumX;
    }

    return transportMixture_;
}


template<class ThermoType>
const typename
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixtureType&
Foam::kineticTheoryMulticomponentMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y,
    const thermoMixtureType&
) const
{
    transportMixture(Y);
    transportMixture_.muCached_ = true;
    return transportMixture_;
}


template<class ThermoType>
Foam::scalar
Foam::kineticTheoryMulticomponentMixture<ThermoType>::Djk
(
    const label k,
    const label j,
    const scalar p,
    const scalar T
) const
{
    const scalar Wk = this->specieThermo(k).W();   // [kg/kmol]
    const scalar Wj = this->specieThermo(j).W();   // [kg/kmol]

    // Защита от отрицательных/нулевых значений (резкие градиенты пламени)
    // могут привести к NaN/Inf в интеграле столкновений.
    const scalar pSafe = max(p, 1.0);      // [Pa]
    const scalar TSafe = max(T, 1.0);      // [K]

    // Reduced molecular weight, Kee Eq. (12.92)
    const scalar Wkj = Wk*Wj/(Wk + Wj);

    const scalar epsilonk = kB*this->specieThermo(k).welldepth(); // [J]
    const scalar epsilonj = kB*this->specieThermo(j).welldepth(); // [J]

    // Kee Eq. (12.36) - geometric-mean well depth
    scalar epsilonkj = sqrt(epsilonk*epsilonj);

    const scalar sigmak = 1.0e-10*this->specieThermo(k).diameter(); // [m]
    const scalar sigmaj = 1.0e-10*this->specieThermo(j).diameter(); // [m]

    // Kee Eq. (12.35) - arithmetic-mean collision diameter
    scalar sigmakj = 0.5*(sigmak + sigmaj);

    const scalar miuk =
        1.0e-21*pow(10, -7.0/2.0)*this->specieThermo(k).dipole();
    const scalar miuj =
        1.0e-21*pow(10, -7.0/2.0)*this->specieThermo(j).dipole();

    const scalar polark = 1.0e-30*this->specieThermo(k).polar();  // [m^3]
    const scalar polarj = 1.0e-30*this->specieThermo(j).polar();  // [m^3]

    scalar deltastar = 0.0;

    // Stockmayer polar/nonpolar interaction corrections (Kee Eqs. 12.37-12.42)
    // Ported 1:1 from the validated of-11 standardKinetic reference.
    if (polark == 0.0)
    {
        if (polarj == 0.0)
        {
            // both species without polarizability correction
            deltastar =
                0.5*miuk*miuj/(epsilonkj*pow3(sigmakj));
        }
        else
        {
            const scalar alphajstar = polarj/pow3(sigmaj);     // non-polar
            const scalar miukstar =
                miuk/sqrt(epsilonk*pow3(sigmak));             // polar
            const scalar ksi =
                1.0 + 0.25*alphajstar*miukstar*sqrt(epsilonk/epsilonj);

            epsilonkj *= sqr(ksi);                            // Kee Eq. (12.38)
            sigmakj *= pow(ksi, -1.0/6.0);                    // Kee Eq. (12.39)

            deltastar =
                0.5*miuk*miuj/(epsilonkj*pow3(sigmakj));
        }
    }
    else
    {
        if (polarj == 0.0)
        {
            const scalar alphakstar = polark/pow3(sigmak);     // non-polar
            const scalar miujstar =
                miuj/sqrt(epsilonj*pow3(sigmaj));             // polar
            const scalar ksi =
                1.0 + 0.25*alphakstar*miujstar*sqrt(epsilonj/epsilonk);

            epsilonkj *= sqr(ksi);                            // Kee Eq. (12.38)
            sigmakj *= pow(ksi, -1.0/6.0);                    // Kee Eq. (12.39)

            deltastar =
                0.5*miuk*miuj/(epsilonkj*pow3(sigmakj));
        }
        else
        {
            deltastar = 0.0;
        }
    }

    // Kee Eq. (12.8) - reduced temperature
    const scalar Tstar = max(TSafe*kB/epsilonkj, 1.0e-3);

    // Kee Table 12.1 - Omega^(1,1)* coefficients a_m
    const scalar a[6] = {1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5};

    // Kee Eq. (12.19) = (12.6)
    scalar omega =
        a[0]*pow(Tstar, -a[1]) + pow(Tstar + a[2], -a[3]);
    const scalar f =
        1.0
      + (exp(a[4]/Tstar) - exp(-a[5]/Tstar))*sqr(deltastar)
       /(2.0 + 2.5*deltastar);

    omega = max(omega*f, 1.0e-10);

    // Kee Eq. (12.113) - binary diffusion coefficient
    return
        max
        (
            3.0/16.0
           *sqrt(2.0*NA*pow3(kB*TSafe)/(PI*Wkj/1000.0))
           /(pSafe*sqr(sigmakj)*omega),
            1.0e-30
        );
}


// ************************************************************************* //
