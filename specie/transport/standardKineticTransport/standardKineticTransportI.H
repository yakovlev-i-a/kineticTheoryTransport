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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::standardKineticTransport<Thermo>::standardKineticTransport
(
    const Thermo& t,
    const scalar geometry,
    const scalar welldepth,
    const scalar diameter,
    const scalar dipole,
    const scalar polar,
    const scalar relaxation
)
:
    Thermo(t),
    geometry_(geometry),
    welldepth_(welldepth),
    diameter_(diameter),
    dipole_(dipole),
    polar_(polar),
    relaxation_(relaxation)
{}


template<class Thermo>
inline Foam::standardKineticTransport<Thermo>::standardKineticTransport
(
    const word& name,
    const standardKineticTransport& st
)
:
    Thermo(name, st),
    geometry_(st.geometry_),
    welldepth_(st.welldepth_),
    diameter_(st.diameter_),
    dipole_(st.dipole_),
    polar_(st.polar_),
    relaxation_(st.relaxation_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::standardKineticTransport<Thermo>>
Foam::standardKineticTransport<Thermo>::clone() const
{
    return autoPtr<standardKineticTransport<Thermo>>
    (
        new standardKineticTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    const scalar sigma = 1.0e-10*diameter_;        // collision diameter [m]
    const scalar Tstar = T/welldepth_;             // reduced temperature
    const scalar epsilon = kB*welldepth_;          // [J]
    const scalar miu = 1.0e-21*pow(10, -7.0/2.0)*dipole_; // [m^(3/2) J^(1/2)]

    // Reduced dipole moment (Kee Eq. 12.18)
    const scalar deltastar = 0.5*sqr(miu)/(epsilon*pow3(sigma));

    // Kee Eq. (12.100)
    return
        5.0/16.0
       *sqrt((this->W()/1000.0)/NA*kB*T/PI)
       /(sqr(sigma)*Omega(2, Tstar, deltastar));
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::Omega
(
    const label index,
    const scalar Tstar,
    const scalar deltastar
) const
{
    scalar omega = 0.0;

    if (index == 1)
    {
        // Kee Table 12.1, coefficients a_m
        const scalar a[6] = {1.0548, 0.15504, 0.55909, 2.1705, 0.093193, 1.5};

        // Kee Eq. (12.19) = (12.6)
        omega = a[0]*pow(Tstar, -a[1]) + pow(Tstar + a[2], -a[3]);
        const scalar f =
            1.0
          + (exp(a[4]/Tstar) - exp(-a[5]/Tstar))*sqr(deltastar)
           /(2.0 + 2.5*deltastar);

        omega *= f;
    }
    else if (index == 2)
    {
        // Kee Table 12.1, coefficients b_m
        const scalar b[6] = {1.0413, 0.11930, 0.43628, 1.6041, 0.095661, 2.0};

        // Kee Eq. (12.20) = (12.7)
        omega = b[0]*pow(Tstar, -b[1]) + pow(Tstar + b[2], -b[3]);
        const scalar f =
            1.0
          + (exp(b[4]/Tstar) - exp(-b[5]/Tstar))*sqr(deltastar)
           /(2.0 + 2.5*deltastar);

        omega *= f;
    }
    else
    {
        FatalErrorInFunction
            << "Collision integral index must be 1 or 2, got " << index << nl
            << exit(FatalError);
    }

    return omega;
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::Dkk
(
    const scalar p,
    const scalar T
) const
{
    const scalar sigma = 1.0e-10*diameter_;        // collision diameter [m]
    const scalar Tstar = T/welldepth_;
    const scalar epsilon = kB*welldepth_;          // [J]
    const scalar miu = 1.0e-21*pow(10, -7.0/2.0)*dipole_;
    const scalar deltastar = 0.5*sqr(miu)/(epsilon*pow3(sigma));

    // Kee Eq. (12.110) - self-diffusion coefficient
    return
        3.0/8.0
       *sqrt(NA*pow3(kB*T)/(PI*(this->W()/1000.0)))
       /(p*sqr(sigma)*Omega(1, Tstar, deltastar));
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::F
(
    const scalar T
) const
{
    // Kee Eq. (12.112)
    return
        1.0
      + 0.5*pow(PI, 3.0/2.0)*sqrt(welldepth_/T)
      + (sqr(PI)/4.0 + 2.0)*(welldepth_/T)
      + pow(PI, 3.0/2.0)*pow(welldepth_/T, 3.0/2.0);
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::Zrot
(
    const scalar T
) const
{
    // Kee Eq. (12.111)
    return relaxation_*F(298)/F(T);
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    // Convert the specific heat from [J/(kg K)] to [J/(mol K)]
    const scalar Cv_total = (this->W()/1000.0)*this->Cv(p, T);

    // Kee Eq. (8.124) - translational contribution
    const scalar Cv_trans = 3.0/2.0*R;

    scalar Cv_rot = 0.0;
    scalar Cv_vib = 0.0;

    if (geometry_ == 0)
    {
        // Single atom: no rotational/vibrational contribution
        Cv_rot = 0.0;
        Cv_vib = 0.0;
    }
    else if (geometry_ == 1)
    {
        Cv_rot = R;                                     // Kee Eq. (8.127)
        Cv_vib = Cv_total - Cv_trans - Cv_rot;          // Kee Eq. (12.108)
    }
    else if (geometry_ == 2)
    {
        Cv_rot = 3.0/2.0*R;                             // Kee Eq. (8.128)
        Cv_vib = Cv_total - Cv_trans - Cv_rot;          // Kee Eq. (12.108)
    }
    else
    {
        FatalErrorInFunction
            << "Molecule geometry must be 0, 1 or 2, got " << geometry_ << nl
            << exit(FatalError);
    }

    const scalar rhoDkkOverMu = this->rho(p, T)*Dkk(p, T)/mu(p, T);

    // Kee Eq. (12.106)
    const scalar A = 5.0/2.0 - rhoDkkOverMu;

    // Kee Eq. (12.107)
    const scalar B =
        Zrot(T)
      + 2.0/PI*(5.0/3.0*Cv_rot/R + rhoDkkOverMu);

    // Kee Eq. (12.103)
    const scalar ftrans = 5.0/2.0*(1.0 - 2.0/PI*Cv_rot/Cv_trans*A/B);

    // Kee Eq. (12.104)
    const scalar frot = rhoDkkOverMu*(1.0 + 2.0/PI*A/B);

    // Kee Eq. (12.105)
    const scalar fvib = rhoDkkOverMu;

    // Kee Eq. (12.102)
    return
        mu(p, T)/(this->W()/1000.0)
       *(ftrans*Cv_trans + frot*Cv_rot + fvib*Cv_vib);
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::geometry() const
{
    return geometry_;
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::welldepth() const
{
    return welldepth_;
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::diameter() const
{
    return diameter_;
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::dipole() const
{
    return dipole_;
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::polar() const
{
    return polar_;
}


template<class Thermo>
inline Foam::scalar Foam::standardKineticTransport<Thermo>::relaxation() const
{
    return relaxation_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::standardKineticTransport<Thermo>::operator=
(
    const standardKineticTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    geometry_ = st.geometry_;
    welldepth_ = st.welldepth_;
    diameter_ = st.diameter_;
    dipole_ = st.dipole_;
    polar_ = st.polar_;
    relaxation_ = st.relaxation_;
}


template<class Thermo>
inline void Foam::standardKineticTransport<Thermo>::operator+=
(
    const standardKineticTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();

        // Molecular transport parameters are not additive: keep those of the
        // most recently added species (consistent with the of-11 reference).
        geometry_ = st.geometry_;
        welldepth_ = st.welldepth_;
        diameter_ = st.diameter_;
        dipole_ = st.dipole_;
        polar_ = st.polar_;
        relaxation_ = st.relaxation_;
    }
}


template<class Thermo>
inline void Foam::standardKineticTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::standardKineticTransport<Thermo> Foam::operator+
(
    const standardKineticTransport<Thermo>& st1,
    const standardKineticTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < small)
    {
        return standardKineticTransport<Thermo>
        (
            t,
            st1.geometry_,
            st1.welldepth_,
            st1.diameter_,
            st1.dipole_,
            st1.polar_,
            st1.relaxation_
        );
    }
    else
    {
        return standardKineticTransport<Thermo>
        (
            t,
            st1.geometry_,
            st1.welldepth_,
            st1.diameter_,
            st1.dipole_,
            st1.polar_,
            st1.relaxation_
        );
    }
}


template<class Thermo>
inline Foam::standardKineticTransport<Thermo> Foam::operator*
(
    const scalar s,
    const standardKineticTransport<Thermo>& st
)
{
    return standardKineticTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.geometry_,
        st.welldepth_,
        st.diameter_,
        st.dipole_,
        st.polar_,
        st.relaxation_
    );
}


// ************************************************************************* //
