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

#include "fluidMulticomponentThermophysicalTransportModels.H"

#include "kineticTheoryFickianFourier.H"
#include "kineticTheoryHirschfelderCurtiss.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeLaminarThermophysicalTransportModel(kineticTheoryFickianFourier);
makeLaminarThermophysicalTransportModel(kineticTheoryHirschfelderCurtiss);


// ************************************************************************* //
