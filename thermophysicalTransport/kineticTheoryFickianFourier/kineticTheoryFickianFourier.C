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

#include "kineticTheoryFickianFourier.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
void kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::
updateDm() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    const kineticTheoryMixture& mixture =
        dynamicCast<const kineticTheoryMixture>(this->thermo());

    Dm_.setSize(Y.size());

    const volScalarField Wm(this->thermo().W());

    // sumXbyD 必须通过 volScalarField 运算累积（而非只改 primitiveFieldRef），
    // 否则 boundaryField 保持 0，导致 Dm 在边界除以 small -> 巨大值 -> 矩阵病态。
    volScalarField sumXbyD
    (
        volScalarField::New
        (
            "sumXbyD",
            T.mesh(),
            dimensionedScalar
            (
                dimless/dimKinematicViscosity/Wm.dimensions(),
                0
            )
        )
    );

    forAll(Dm_, i)
    {
        sumXbyD = dimensionedScalar(sumXbyD.dimensions(), 0);

        forAll(Y, j)
        {
            if (j != i)
            {
                // Dij 构造为完整 volScalarField（internal + boundary）
                volScalarField Dij
                (
                    volScalarField::New
                    (
                        "Dij",
                        T.mesh(),
                        dimensionedScalar(dimKinematicViscosity, 0)
                    )
                );

                forAll(T, celli)
                {
                    Dij[celli] =
                        mixture.Djk(i, j, p[celli], T[celli]);
                }
                forAll(T.boundaryField(), patchi)
                {
                    scalarField& bf = Dij.boundaryFieldRef()[patchi];
                    forAll(bf, facei)
                    {
                        bf[facei] = mixture.Djk
                        (
                            i, j,
                            p.boundaryField()[patchi][facei],
                            T.boundaryField()[patchi][facei]
                        );
                    }
                }

                sumXbyD += Y[j]/(this->thermo().Wi(j)*Dij);
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
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::Dm() const
{
    if (!Dm_.size())
    {
        updateDm();
    }

    return Dm_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::
kineticTheoryFickianFourier
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
bool kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::read()
{
    return
        unityLewisFourier<laminarThermophysicalTransportModel>::read();
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::D
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
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::D
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
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::DEff
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
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::DEff
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
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::q() const
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
        surfaceScalarField sumJ
        (
            surfaceScalarField::New
            (
                "sumJ",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(sumJ.dimensions()*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(this->j(Y[i]));

                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const volScalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::q
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
        scalarField sumJ(tmpq->size(), scalar(0));
        scalarField sumJh(tmpq->size(), scalar(0));

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const scalarField hi(this->thermo().hsi(i, p, T));

                const scalarField ji(this->j(Y[i], patchi));

                sumJ += ji;

                sumJh += ji*hi;
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const scalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*hi;
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::
divq(volScalarField& he) const
{
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
            he
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    tmpDivq.ref() -=
        fvm::laplacianCorrection(this->alpha()*this->alphaEff(), he);

    surfaceScalarField sumJ
    (
        surfaceScalarField::New
        (
            "sumJ",
            he.mesh(),
            dimensionedScalar(dimMass/dimArea/dimTime, 0)
        )
    );

    surfaceScalarField sumJh
    (
        surfaceScalarField::New
        (
            "sumJh",
            he.mesh(),
            dimensionedScalar(sumJ.dimensions()*he.dimensions(), 0)
        )
    );

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));

            const surfaceScalarField ji(this->j(Y[i]));

            sumJ += ji;

            sumJh += ji*fvc::interpolate(hi);
        }
    }

    {
        const label i = this->thermo().defaultSpecie();

        const volScalarField hi(this->thermo().hsi(i, p, T));

        sumJh -= sumJ*fvc::interpolate(hi);
    }

    tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());

    return tmpDivq;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::predict()
{
    unityLewisFourier<laminarThermophysicalTransportModel>::predict();
    updateDm();
}


template<class laminarThermophysicalTransportModel>
bool kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::movePoints()
{
    return true;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::topoChange
(
    const polyTopoChangeMap& map
)
{
    Dm_.clear();
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::mapMesh
(
    const polyMeshMap& map
)
{
    Dm_.clear();
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryFickianFourier<laminarThermophysicalTransportModel>::distribute
(
    const polyDistributionMap& map
)
{
    Dm_.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
