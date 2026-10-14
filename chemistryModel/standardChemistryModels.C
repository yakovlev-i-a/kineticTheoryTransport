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

Description
    Static instantiation of the standard chemistryModel (and its reaction /
    chemistry-reduction machinery) templated on the kineticTheoryTransport
    thermodynamics package.

    This replaces the on-the-fly dynamicCode compilation that would otherwise
    happen in every case, i.e. avoids the `dynamicCode/` directory being
    generated at runtime.
\*---------------------------------------------------------------------------*/

#include "typedefThermo.H"

#include "specie.H"

#include "thermo.H"

// EoS
#include "perfectGas.H"

// Thermo
#include "janafThermo.H"
#include "sensibleEnthalpy.H"

// Transport
#include "standardKineticTransport.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define ThermoPhysics                                                          \
    standardKineticTransportsensibleEnthalpyjanafThermoperfectGasspecie

namespace Foam
{
    typedefThermo
    (
        standardKineticTransport,
        sensibleEnthalpy,
        janafThermo,
        perfectGas,
        specie
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeChemistryModel.H"
#include "Standard_chemistryModel.H"

namespace Foam
{
    makeChemistryModel(Standard, ThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeChemistryReductionMethod.H"

#include "noChemistryReduction.H"
#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

namespace Foam
{
    defineChemistryReductionMethod(nullArg, ThermoPhysics);

    makeChemistryReductionMethod(none, ThermoPhysics);
    makeChemistryReductionMethod(DAC, ThermoPhysics);
    makeChemistryReductionMethod(DRG, ThermoPhysics);
    makeChemistryReductionMethod(DRGEP, ThermoPhysics);
    makeChemistryReductionMethod(EFA, ThermoPhysics);
    makeChemistryReductionMethod(PFA, ThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "FallOffReactionRate.H"
#include "ChemicallyActivatedReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "MichaelisMentenReactionRate.H"
#include "LangmuirHinshelwoodReactionRate.H"
#include "fluxLimitedLangmuirHinshelwoodReactionRate.H"
#include "surfaceArrheniusReactionRate.H"

namespace Foam
{
    defineReaction(nullArg, ThermoPhysics);

    // Irreversible/reversible/non-equilibrium-reversible reactions
    makeIRNReactions(ArrheniusReactionRate, ThermoPhysics);
    makeIRNReactions(LandauTellerReactionRate, ThermoPhysics);
    makeIRNReactions(thirdBodyArrheniusReactionRate, ThermoPhysics);

    // Irreversible/reversible reactions
    makeIRReactions(JanevReactionRate, ThermoPhysics);
    makeIRReactions(powerSeriesReactionRate, ThermoPhysics);

    // Irreversible/reversible fall-off reactions
    makeIRTemplate2Reactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );

    // Irreversible/reversible chemically activated reactions
    makeIRTemplate2Reactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );
    makeIRTemplate2Reactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );

    // Michaelis-Menten Reactions
    makeIReactions(MichaelisMentenReactionRate, ThermoPhysics);

    // Langmuir-Hinshelwood Reactions
    makeIRReactions(LangmuirHinshelwoodReactionRate, ThermoPhysics);

    // Flux-limited Langmuir-Hinshelwood Reactions
    makeGeneralReaction
    (
        IrreversibleReaction,
        fluxLimitedLangmuirHinshelwoodReactionRate,
        ThermoPhysics
    );

    // Surface-Arrhenius Reactions
    makeGeneralReaction
    (
        IrreversibleReaction,
        surfaceArrheniusReactionRate,
        ThermoPhysics
    );
}


// ************************************************************************* //
