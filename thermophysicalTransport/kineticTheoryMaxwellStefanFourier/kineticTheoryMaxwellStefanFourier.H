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
    Foam::laminarThermophysicalTransportModels::
    kineticTheoryMaxwellStefanFourier

Description
    Multi-component Maxwell-Stefan generalised Fick's law diffusion
    coefficients with Fourier temperature-gradient heat flux, for laminar
    flow.  This is the strict (coupled) multi-component diffusion closure.

    Unlike the native MaxwellStefan/MaxwellStefanFourier models (where the
    binary mass diffusion coefficients D_ij are read as Function2<scalar> of
    (p, T) from the dictionary), this model evaluates D_ij on the fly from the
    standard kinetic theory of gases (Chapman-Enskog / Lennard-Jones +
    Stockmayer) through kineticTheoryMulticomponentMixture::Djk(i, j, p, T).

    The generalised Fick's law diffusion coefficients are obtained by the
    Taylor-Krishna transformation of the Maxwell-Stefan equation in every cell
    and patch face:

        nabla X_i = sum_{j != i} X_i X_j / D_ij (V_j - V_i)

    which is inverted (N-1 sub-system, excluding the default specie) to give

        j_i = -rho D_ii^Fick grad Y_i - sum_{j != i} rho D_ij^Fick grad Y_j

    The off-diagonal part is treated explicitly (jexp_), the diagonal part is
    implicit (unityLewisFourier base with DEff = rho D_ii^Fick).

    The heat flux includes the differential-diffusion enthalpy correction
    sum_j j_j h_j, identical to MaxwellStefan/MaxwellStefanFourier.

    Usage
    \verbatim
    thermophysicalTransport
    {
        model    kineticTheoryMaxwellStefanFourier;
        libs     ("libkineticTheoryTransport.so");
    }
    \endverbatim

    The thermo model must be instantiated with
    kineticTheoryMulticomponentMixture (see the library's thermos.C).

SourceFiles
    kineticTheoryMaxwellStefanFourier.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryMaxwellStefanFourier_H
#define kineticTheoryMaxwellStefanFourier_H

#include "unityLewisFourier.H"
#include "kineticTheoryMixture.H"
#include "MeshObjects.H"
#include "LUscalarMatrix.H"
#include "scalarMatrices.H"
#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
        Class kineticTheoryMaxwellStefanFourier Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class kineticTheoryMaxwellStefanFourier
:
    public unityLewisFourier<laminarThermophysicalTransportModel>,
    public TopoChangeableMeshObject<fvMesh>
{
    // Private data

        //- Generalised Fick's law diffusion coefficients field list.
        //  This is the diagonal of the mass diffusion coefficient matrix
        mutable PtrList<volScalarField> Dii_;

        //- List of fields of the explicit part of the mass diffusion flux
        //  of the species
        mutable PtrList<surfaceScalarField> jexp_;

        //- List of effective (mixture) mass diffusivity fields
        mutable PtrList<volScalarField> Di_;


    // Workspace for diffusion coefficient transformation

        //- Molecular weights of the species
        mutable scalarField W;

        //- List of mass-fraction field pointers
        //  for the internal mesh or a patch
        mutable List<const scalarField*> YPtrs;

        //- Matrix of mass diffusion coefficient field pointers
        //  for the internal mesh or a patch
        mutable SquareMatrix<scalarField*> DijPtrs;

        //- Mass-fractions at a cell or face
        mutable scalarField Y;

        //- Mole-fractions at a cell or face
        mutable scalarField X;

        //- Binary mass diffusion coefficient matrix at a cell or face
        mutable scalarSquareMatrix DD;

        //- Matrix form of the coefficients in the Maxwell-Stefan equation
        //  at a cell or face
        mutable LUscalarMatrix A;

        //- Matrix of composition coefficients at a cell or face
        //  used to transform the binary mass diffusion coefficients into
        //  the generalised Fick's law diffusion coefficients
        mutable scalarSquareMatrix B;

        //- Inverse of A
        mutable scalarSquareMatrix invA;

        //- Matrix of the generalised Fick's law diffusion coefficients
        //  at a cell or face
        mutable scalarSquareMatrix D_;


    // Private Member Functions

        //- Binary mass diffusion coefficient D_ij(i, j, p, T) [m^2/s] as a
        //  volScalarField, evaluated on the fly from kinetic theory
        tmp<volScalarField> binaryDjk
        (
            const label i,
            const label j,
            const volScalarField& p,
            const volScalarField& T
        ) const;

        //- Transform the Binary mass diffusion coefficient matrix DD into the
        //  matrix of the generalised Fick's law diffusion coefficients D
        void transformDiffusionCoefficient() const;

        //- Calls transformDiffusionCoefficient() for each cell and patch face
        //  and maps the resulting D into DijPtrs
        void transformDiffusionCoefficientFields() const;

        //- Calls transformDiffusionCoefficientFields() for the internal mesh
        //  and each patch and maps the resulting DijPtrs into Dii_ and Dij
        void transform(List<PtrList<volScalarField>>& Dij) const;

        //- Update Dii_
        void updateDii() const;

        //- Update Dii_ if not yet constructed and return
        const PtrList<volScalarField>& Dii() const;

        //- Update jexp_ if not yet constructed and return
        const PtrList<surfaceScalarField>& jexp() const;

        //- Update Di_ if not yet constructed and return
        const PtrList<volScalarField>& Di() const;


public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("kineticTheoryMaxwellStefanFourier");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        kineticTheoryMaxwellStefanFourier
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~kineticTheoryMaxwellStefanFourier()
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

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
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

        //- Update the diffusion coefficients and flux corrections
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
    #include "kineticTheoryMaxwellStefanFourier.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
