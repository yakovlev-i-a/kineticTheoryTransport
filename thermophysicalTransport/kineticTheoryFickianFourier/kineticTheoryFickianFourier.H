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
    Foam::laminarThermophysicalTransportModels::kineticTheoryFickianFourier

Description
    Laminar mixture-averaged transport model based on the standard kinetic
    theory of gases (Chapman-Enskog / Lennard-Jones + Stockmayer).

    - Species mass diffusion coefficients w.r.t. the mixture D_{i,m} are
      assembled from binary diffusion coefficients D_{ij}, themselves computed
      from molecular transport parameters (sigma, epsilon/kB, dipole,
      polarizability) rather than read from the dictionary.
    - Mixture viscosity (Wilke) and thermal conductivity (combination
      averaging) come from the kinetic-theory mixture.
    - The energy equation receives the correct differential-diffusion
      enthalpy-flux correction (sum of j_i h_i), identical to Fickian.

    Usage
    \verbatim
    thermophysicalTransport
    {
        model    kineticTheoryFickianFourier;
        libs     ("libkineticTheoryTransport.so");
    }
    \endverbatim

    The thermo model must be instantiated with
    kineticTheoryMulticomponentMixture (see the library's thermos.C).

SourceFiles
    kineticTheoryFickianFourier.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryFickianFourier_H
#define kineticTheoryFickianFourier_H

#include "unityLewisFourier.H"
#include "kineticTheoryMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
            Class kineticTheoryFickianFourier Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class kineticTheoryFickianFourier
:
    public unityLewisFourier<laminarThermophysicalTransportModel>
{
    // Private data

        //- Mass diffusion coefficients in the mixture [m^2/s]
        mutable PtrList<volScalarField> Dm_;


    // Private Member Functions

        //- Update Dm_ from the kinetic-theory binary coefficients
        void updateDm() const;

        //- Update Dm_ if not yet constructed and return
        const PtrList<volScalarField>& Dm() const;


public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("kineticTheoryFickianFourier");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        kineticTheoryFickianFourier
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~kineticTheoryFickianFourier()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Mass diffusivity
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> D(const volScalarField& Yi) const;

        //- Mass diffusivity
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> D
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the patch heat flux [W/m^2]
        virtual tmp<scalarField> q(const label patchi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Update the diffusion coefficients
        virtual void predict();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap& map);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap& map);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kineticTheoryFickianFourier.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
