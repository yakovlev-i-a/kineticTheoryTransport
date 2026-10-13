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
    Foam::kineticTheoryMixture

Description
    Non-template runtime interface exposing the kinetic-theory binary mass
    diffusion coefficient of a multicomponent mixture.

    The templated kineticTheoryMulticomponentMixture<ThermoType> inherits this
    interface.  Because the mixture is a (non-virtual) base class of the
    concrete thermo type (through BasicThermo<MixtureType>), a thermophysical
    transport model holding only the abstract fluidMulticomponentThermo& can
    recover the mixture via a dynamic_cast to this interface.

SourceFiles
    (none - interface only)

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryMixture_H
#define kineticTheoryMixture_H

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class kineticTheoryMixture Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryMixture
{
public:

    //- Runtime type information
    TypeName("kineticTheoryMixture");


    //- Destructor
    virtual ~kineticTheoryMixture()
    {}


    // Member Functions

        //- Binary mass diffusion coefficient between species i and j [m^2/s]
        virtual scalar Djk
        (
            const label i,
            const label j,
            const scalar p,
            const scalar T
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
