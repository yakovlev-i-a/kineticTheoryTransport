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

#include "standardKineticTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::standardKineticTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").lookup<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::standardKineticTransport<Thermo>::standardKineticTransport
(
    const word& name,
    const dictionary& dict
)
:
    Thermo(name, dict),
    geometry_(readCoeff("geometry", dict)),
    welldepth_(readCoeff("welldepth", dict)),
    diameter_(readCoeff("diameter", dict)),
    dipole_(readCoeff("dipole", dict)),
    polar_(readCoeff("polar", dict)),
    relaxation_(readCoeff("relaxation", dict))
{}


template<class Thermo>
Foam::standardKineticTransport<Thermo>::standardKineticTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    geometry_(readCoeff("geometry", dict)),
    welldepth_(readCoeff("welldepth", dict)),
    diameter_(readCoeff("diameter", dict)),
    dipole_(readCoeff("dipole", dict)),
    polar_(readCoeff("polar", dict)),
    relaxation_(readCoeff("relaxation", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::standardKineticTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("geometry", geometry_);
    dict.add("welldepth", welldepth_);
    dict.add("diameter", diameter_);
    dict.add("dipole", dipole_);
    dict.add("polar", polar_);
    dict.add("relaxation", relaxation_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const standardKineticTransport<Thermo>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
