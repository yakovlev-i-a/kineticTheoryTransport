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

Class
    Foam::standardKineticTransport

Description
    Transport package for a pure species based on the standard kinetic theory
    of gases (Chapman-Enskog / Lennard-Jones + Stockmayer polar correction).

    The pure-species dynamic viscosity mu_k, thermal conductivity kappa_k and
    self-diffusion coefficient D_kk are evaluated from the six molecular
    transport parameters read from the CHEMKIN ``transport`` dictionary:

    \verbatim
    transport
    {
        geometry        1;       // 0 = atom, 1 = linear, 2 = nonlinear
        welldepth       97.53;   // Lennard-Jones potential well depth  [K]
        diameter        3.621;   // Lennard-Jones collision diameter  [Angstrom]
        dipole          0.0;     // dipole moment                      [Debye]
        polar           1.76;    // polarizability                     [A^3]
        relaxation      0.0;     // rotational relaxation number at 298 K [-]
    }
    \endverbatim

    This class is a drop-in replacement for sutherlandTransport and is the
    ThermoType-level building block of libkineticTheoryTransport.

    Reference:
        Wang, An & Chan (2025), Computer Physics Communications 309, 109472.
        Kee, Coltrin & Glarborg, "Chemically Reacting Flow" (2003).

SourceFiles
    standardKineticTransportI.H
    standardKineticTransport.C

\*---------------------------------------------------------------------------*/

#ifndef standardKineticTransport_H
#define standardKineticTransport_H

#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class standardKineticTransport;

template<class Thermo>
inline standardKineticTransport<Thermo> operator+
(
    const standardKineticTransport<Thermo>&,
    const standardKineticTransport<Thermo>&
);

template<class Thermo>
inline standardKineticTransport<Thermo> operator*
(
    const scalar,
    const standardKineticTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const standardKineticTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                 Class standardKineticTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class standardKineticTransport
:
    public Thermo
{
    // Private Data

        //- Ratio of the circumference of a circle to its diameter
        const scalar PI = constant::mathematical::pi;

        //- Avogadro constant [mol^-1]
        const scalar NA = constant::physicoChemical::NA.value();

        //- Boltzmann constant [J/K]
        const scalar kB = constant::physicoChemical::k.value();

        //- Universal gas constant [J/(mol K)]
        const scalar R = constant::physicoChemical::R.value();

        //- Molecule geometrical configuration
        //  0 = single atom, 1 = linear molecule, 2 = nonlinear molecule
        scalar geometry_;

        //- Lennard-Jones potential well depth [K]  (epsilon/kB)
        scalar welldepth_;

        //- Lennard-Jones collision diameter [Angstrom]  (10^-10 m)
        scalar diameter_;

        //- Dipole moment [Debye]
        scalar dipole_;

        //- Polarizability [cubic Angstrom]  (10^-30 m^3)
        scalar polar_;

        //- Rotational relaxation collision number at 298 K [-]
        scalar relaxation_;


    // Private Member Functions

        //- Read coefficient from the transport sub-dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        inline standardKineticTransport
        (
            const Thermo& t,
            const scalar geometry,
            const scalar welldepth,
            const scalar diameter,
            const scalar dipole,
            const scalar polar,
            const scalar relaxation
        );

        //- Construct as named copy
        inline standardKineticTransport
        (
            const word&,
            const standardKineticTransport&
        );

        //- Construct from name and dictionary
        standardKineticTransport(const word& name, const dictionary& dict);

        //- Construct from base thermo and dictionary
        standardKineticTransport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<standardKineticTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "standardKinetic<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Collision integral Omega^(l,s)* for the given reduced
        //  temperature and reduced dipole moment
        //  index: 1 -> Omega^(1,1)*, 2 -> Omega^(2,2)*
        inline scalar Omega
        (
            const label index,
            const scalar Tstar,
            const scalar deltastar
        ) const;

        //- Self-diffusion coefficient [m^2/s]
        inline scalar Dkk(const scalar p, const scalar T) const;

        //- Temperature-dependent rotational-relaxation factor F(T)
        inline scalar F(const scalar T) const;

        //- Rotational relaxation collision number Zrot(T)
        inline scalar Zrot(const scalar T) const;

        //- Molecule geometrical configuration (0/1/2)
        inline scalar geometry() const;

        //- Lennard-Jones potential well depth [K]
        inline scalar welldepth() const;

        //- Lennard-Jones collision diameter [Angstrom]
        inline scalar diameter() const;

        //- Dipole moment [Debye]
        inline scalar dipole() const;

        //- Polarizability [A^3]
        inline scalar polar() const;

        //- Rotational relaxation collision number at 298 K [-]
        inline scalar relaxation() const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator=(const standardKineticTransport&);
        inline void operator+=(const standardKineticTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend standardKineticTransport operator+ <Thermo>
        (
            const standardKineticTransport&,
            const standardKineticTransport&
        );

        friend standardKineticTransport operator* <Thermo>
        (
            const scalar,
            const standardKineticTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const standardKineticTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "standardKineticTransportI.H"

#ifdef NoRepository
    #include "standardKineticTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
