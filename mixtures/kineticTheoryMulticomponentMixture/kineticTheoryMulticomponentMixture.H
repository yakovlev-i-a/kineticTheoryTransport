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
    Foam::kineticTheoryMulticomponentMixture

Description
    Multi-component mixture based on the standard kinetic theory of gases.

    - Pure-species mu_k, kappa_k, D_kk come from standardKineticTransport.
    - Mixture viscosity: Wilke's rule.
    - Mixture thermal conductivity: combination averaging (Mathur et al.).
    - Binary mass diffusion coefficients D_jk: Chapman-Enskog with
      Stockmayer (Lennard-Jones + polar) collision integrals.

    This is the mixture-level building block of libkineticTheoryTransport,
    drop-in replacement for coefficientWilkeMulticomponentMixture.

SourceFiles
    kineticTheoryMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryMulticomponentMixture_H
#define kineticTheoryMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "kineticTheoryMixture.H"
#include "FieldListSlice.H"
#include "scalarMatrices.H"

#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class kineticTheoryMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class kineticTheoryMulticomponentMixture
:
    public multicomponentMixture<ThermoType>,
    public kineticTheoryMixture
{
    // Private Data

        //- Ratio of the circumference of a circle to its diameter
        const scalar PI = constant::mathematical::pi;

        //- Avogadro constant [mol^-1]
        const scalar NA = constant::physicoChemical::NA.value();

        //- Boltzmann constant [J/K]
        const scalar kB = constant::physicoChemical::k.value();


public:

    // Public Classes/Typedefs

        //- Mixing type for thermodynamic properties
        typedef typename ThermoType::thermoType thermoMixtureType;

        //- Mixing type for transport properties
        class transportMixtureType
        {
            // Private Data

                //- List of specie thermo
                const PtrList<ThermoType>& specieThermos_;

                //- List of molecular weights
                scalarList M_;

                //- Matrix of 4/sqrt(2)*sqrt(1 + M_[i]/M_[j])
                scalarSquareMatrix A_;

                //- Matrix of sqrt(M_[j]/M_[i])
                scalarSquareMatrix B_;

                //- List of mole fractions
                mutable scalarList X_;

                //- List of specie viscosities
                mutable scalarList mu_;

                //- List of Wilke weights
                mutable scalarList w_;

                //- Viscosity cache state (Wilke weights reused by kappa)
                mutable bool muCached_;

                //- Calculate the Wilke weights and store in w_
                void WilkeWeights(const scalar p, const scalar T) const;


        public:

            friend class kineticTheoryMulticomponentMixture;


            // Constructors

                //- Construct from list of specie thermo
                transportMixtureType(const PtrList<ThermoType>& specieThermos);


            // Transport properties

                //- Dynamic viscosity [kg/m/s] (Wilke)
                scalar mu(const scalar p, const scalar T) const;

                //- Thermal conductivity [W/m/K] (combination averaging)
                scalar kappa(const scalar p, const scalar T) const;
        };


private:

    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable thermoMixtureType mixture_;

        //- Mutable storage for the cell/face mixture transport data
        mutable transportMixtureType transportMixture_;


public:

    // Constructors

        //- Construct from a dictionary
        kineticTheoryMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        kineticTheoryMulticomponentMixture
        (
            const kineticTheoryMulticomponentMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "kineticTheoryMulticomponentMixture<"
              + ThermoType::typeName() + '>';
        }

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;

        //- Binary mass diffusion coefficient between species i and j [m^2/s]
        virtual scalar Djk
        (
            const label i,
            const label j,
            const scalar p,
            const scalar T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kineticTheoryMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
