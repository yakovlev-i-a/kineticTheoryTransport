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
    Foam::laminarThermophysicalTransportModels::kineticTheoryHirschfelderCurtiss

Description
    Laminar mixture-averaged transport model using the Hirschfelder-Curtiss
    approximation of the species diffusion velocities.

    The species diffusion velocity of species i is

        V_i = - D_{i,m} grad(Y_i) - (D_{i,m} Y_i / W_m) grad(W_m)

    with the mixture-averaged mass diffusivity

        D_{i,m} = (1 - Y_i)
                  / ( W_m * sum_{j != i} Y_j / (W_j D_{ij}) )

    where W_m is the mean molecular weight and D_{ij} the binary (Chapman-
    Enskog / Lennard-Jones + Stockmayer) diffusion coefficients from
    standardKineticTransport.

    To enforce global mass conservation (sum_i Y_i V_i = 0) a correction
    velocity is added to every species equation:

        V_c = sum_k ( D_{k,m} grad(Y_k) + (D_{k,m} Y_k / W_m) grad(W_m) )
            = - sum_k V_k

    The species equation flux is then

        j_i = - rho D_{i,m} grad(Y_i)
              - rho D_{i,m} Y_i/W_m grad(W_m)
              + rho V_c Y_i

    and the enthalpy diffusion flux sum_i h_i j_i is added to the energy
    equation without the Fickian inert-specie correction (it is zero
    automatically because sum_i j_i = 0).

    Reference:
        Hirschfelder, Curtiss & Bird (1954), "Molecular Theory of Gases and
        Liquids".
        Kee, Coltrin & Glarborg (2003), "Chemically Reacting Flow".

Usage
    \verbatim
    thermophysicalTransport
    {
        model    kineticTheoryHirschfelderCurtiss;
        libs     ("libkineticTheoryTransport.so");
    }
    \endverbatim

    The thermo model must be instantiated with
    kineticTheoryMulticomponentMixture (see the library's thermos.C).

SourceFiles
    kineticTheoryHirschfelderCurtiss.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryHirschfelderCurtiss_H
#define kineticTheoryHirschfelderCurtiss_H

#include "unityLewisFourier.H"
#include "kineticTheoryMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
        Class kineticTheoryHirschfelderCurtiss Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class kineticTheoryHirschfelderCurtiss
:
    public unityLewisFourier<laminarThermophysicalTransportModel>
{
    // Private data

        //- Mass diffusion coefficients in the mixture [m^2/s]
        mutable PtrList<volScalarField> Dm_;

        //- Correction-velocity mass flux rho*V_c & S_f [kg/s]
        mutable autoPtr<surfaceScalarField> phiVcPtr_;


    // Private Member Functions

        //- Update Dm_ from the kinetic-theory binary coefficients
        void updateDm() const;

        //- Update Dm_ if not yet constructed and return
        const PtrList<volScalarField>& Dm() const;

        //- Update the correction-velocity mass flux
        void updateVc() const;

        //- Return the correction-velocity mass flux
        const surfaceScalarField& phiVc() const;


public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("kineticTheoryHirschfelderCurtiss");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        kineticTheoryHirschfelderCurtiss
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~kineticTheoryHirschfelderCurtiss()
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

        //- Return the specie flux for the given specie mass-fraction
        //  [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the specie flux
        //  for the given specie mass-fraction for patch [kg/m^2/s]
        virtual tmp<scalarField> j
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Update the diffusion coefficients and correction velocity
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
    #include "kineticTheoryHirschfelderCurtiss.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
