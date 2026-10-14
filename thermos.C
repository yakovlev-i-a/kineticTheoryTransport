/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 LZF / kineticTheoryTransport
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of a derivative work of OpenFOAM (GPLv3).
    See the GNU General Public License for more details.
\*---------------------------------------------------------------------------*/

#include "psiMulticomponentThermo.H"
#include "rhoFluidMulticomponentThermo.H"

#include "kineticTheoryMulticomponentMixture.H"

#include "standardKineticTransport.H"

#include "forGases.H"

#include "makeFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // ThermoPhysics = standardKineticTransport
    //                 <species::thermo<janafThermo<perfectGas<specie>>,
    //                                  sensibleEnthalpy>>
    typedefThermo
    (
        standardKineticTransport,
        sensibleEnthalpy,
        janafThermo,
        perfectGas,
        specie
    );

    // psi-based: the default route of multicomponentFluid (reactingFoam)
    makeFluidMulticomponentThermos
    (
        psiThermo,
        psiMulticomponentThermo,
        kineticTheoryMulticomponentMixture,
        standardKineticTransportsensibleEnthalpyjanafThermoperfectGasspecie
    );

    // rho-based
    makeFluidMulticomponentThermos
    (
        rhoFluidThermo,
        rhoFluidMulticomponentThermo,
        kineticTheoryMulticomponentMixture,
        standardKineticTransportsensibleEnthalpyjanafThermoperfectGasspecie
    );
}

// ************************************************************************* //
