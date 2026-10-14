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

#include "kineticTheoryHirschfelderCurtiss.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

template<class laminarThermophysicalTransportModel>
void kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
updateDm() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    const kineticTheoryMixture& mixture =
        dynamicCast<const kineticTheoryMixture>(this->thermo());

    Dm_.setSize(Y.size());

    const volScalarField Wm(this->thermo().W());

    forAll(Dm_, i)
    {
        volScalarField sumXbyD
        (
            volScalarField::New
            (
                "sumXbyD",
                T.mesh(),
                dimless/dimKinematicViscosity/Wm.dimensions()
            )
        );

        forAll(Y, j)
        {
            if (j != i)
            {
                const scalarField& pCells = p.primitiveField();
                const scalarField& TCells = T.primitiveField();

                tmp<scalarField> tDij
                (
                    new scalarField(pCells.size())
                );

                scalarField& Dij = tDij.ref();

                forAll(pCells, celli)
                {
                    Dij[celli] =
                        mixture.Djk(i, j, pCells[celli], TCells[celli]);
                }

                sumXbyD.primitiveFieldRef() +=
                    Y[j].primitiveField()
                   /(
                        this->thermo().Wi(j).value()*Dij
                    );
            }
        }

        Dm_.set
        (
            i,
            (1 - Y[i])/Wm
           /max(sumXbyD, dimensionedScalar(sumXbyD.dimensions(), small))
        );
    }
}


template<class laminarThermophysicalTransportModel>
const PtrList<volScalarField>&
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::Dm() const
{
    if (!Dm_.size())
    {
        updateDm();
    }

    return Dm_;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
updateVc() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& rho = this->momentumTransport().rho();
    const volScalarField Wm(this->thermo().W());
    const fvMesh& mesh = rho.mesh();

    if (!phiVcPtr_.valid())
    {
        phiVcPtr_.reset
        (
            surfaceScalarField::New
            (
                "phiVc",
                mesh,
                dimensionedScalar(dimMass/dimTime, 0)
            ).ptr()
        );
    }

    volVectorField Vc
    (
        volVectorField::New
        (
            "Vc",
            mesh,
            dimensionedVector(dimVelocity, Zero)
        )
    );

    const PtrList<volScalarField>& Dmi = this->Dm();

    forAll(Y, k)
    {
        Vc += Dmi[k]*fvc::grad(Y[k]) + (Dmi[k]*Y[k]/Wm)*fvc::grad(Wm);
    }

    phiVcPtr_() = fvc::interpolate(rho*Vc) & mesh.Sf();
}


template<class laminarThermophysicalTransportModel>
const surfaceScalarField&
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
phiVc() const
{
    if (!phiVcPtr_.valid())
    {
        updateVc();
    }

    return phiVcPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
kineticTheoryHirschfelderCurtiss
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisFourier<laminarThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo
    )
{
    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
bool kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
read()
{
    return
        unityLewisFourier<laminarThermophysicalTransportModel>::read();
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::D
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "D",
        this->momentumTransport().rho()
       *Dm()[this->thermo().specieIndex(Yi)]
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::D
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Dm()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *Dm()[this->thermo().specieIndex(Yi)]
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Dm()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    const label i = this->thermo().specieIndex(Yi);
    const volScalarField& rho = this->momentumTransport().rho();
    const volScalarField Wm(this->thermo().W());
    const volScalarField& Di = Dm()[i];

    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "j(" + Yi.name() + ')',
            this->thermo().phaseName()
        ),
       -fvc::interpolate(rho*Di)*fvc::snGrad(Yi)
       -fvc::interpolate(rho*Di*Yi/Wm)*fvc::snGrad(Wm)
      + (phiVc()/rho.mesh().magSf())*fvc::interpolate(Yi)
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const label i = this->thermo().specieIndex(Yi);
    const scalarField& rhop =
        this->momentumTransport().rho().boundaryField()[patchi];
    const volScalarField Wm(this->thermo().W());
    const scalarField Wmp = Wm.boundaryField()[patchi];
    const scalarField snGradWmp = Wm.boundaryField()[patchi].snGrad();
    const scalarField& Dip = Dm()[i].boundaryField()[patchi];
    const fvMesh& mesh = Yi.mesh();

    return
      - rhop*Dip*Yi.boundaryField()[patchi].snGrad()
      - rhop*Dip*Yi.boundaryField()[patchi]/Wmp*snGradWmp
      + (
            phiVc().boundaryField()[patchi]
           /mesh.magSf().boundaryField()[patchi]
        )*Yi.boundaryField()[patchi];
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const label i = this->thermo().specieIndex(Yi);
    const volScalarField& rho = this->momentumTransport().rho();
    const volScalarField Wm(this->thermo().W());
    const volScalarField& Di = Dm()[i];

    tmp<fvScalarMatrix> tDivj
    (
        - fvm::laplacian(rho*Di, Yi)
        - fvc::laplacian(rho*Di*Yi/Wm, Wm)
    );

    tDivj.ref() += fvm::div(phiVc(), Yi, "div(phi,Yi_h)");

    return tDivj;
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->thermo().phaseName()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    if (Y.size())
    {
        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(dimEnergy/dimArea/dimTime, 0)
            )
        );

        forAll(Y, i)
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));

            sumJh += j(Y[i])*fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::q
(
    const label patchi
) const
{
    tmp<scalarField> tmpq
    (
      - (
            this->alpha().boundaryField()[patchi]
           *this->kappaEff(patchi)
           *this->thermo().T().boundaryField()[patchi].snGrad()
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const scalarField& p = this->thermo().p().boundaryField()[patchi];
    const scalarField& T = this->thermo().T().boundaryField()[patchi];

    if (Y.size())
    {
        scalarField sumJh(tmpq->size(), scalar(0));

        forAll(Y, i)
        {
            const scalarField hi(this->thermo().hsi(i, p, T));

            sumJh += j(Y[i], patchi)*hi;
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
            he
        )
    );

    tmpDivq.ref() -=
        fvm::laplacianCorrection(this->alpha()*this->alphaEff(), he);

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    surfaceScalarField sumJh
    (
        surfaceScalarField::New
        (
            "sumJh",
            he.mesh(),
            dimensionedScalar(dimEnergy/dimArea/dimTime, 0)
        )
    );

    forAll(Y, i)
    {
        const volScalarField hi(this->thermo().hsi(i, p, T));

        sumJh += j(Y[i])*fvc::interpolate(hi);
    }

    tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());

    return tmpDivq;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
predict()
{
    unityLewisFourier<laminarThermophysicalTransportModel>::predict();
    updateDm();
    updateVc();
}


template<class laminarThermophysicalTransportModel>
bool kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
movePoints()
{
    return true;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
topoChange
(
    const polyTopoChangeMap& map
)
{
    Dm_.clear();
    phiVcPtr_.clear();
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
mapMesh
(
    const polyMeshMap& map
)
{
    Dm_.clear();
    phiVcPtr_.clear();
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryHirschfelderCurtiss<laminarThermophysicalTransportModel>::
distribute
(
    const polyDistributionMap& map
)
{
    Dm_.clear();
    phiVcPtr_.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
