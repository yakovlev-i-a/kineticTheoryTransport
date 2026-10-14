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

#include "kineticTheoryMaxwellStefanFourier.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "surfaceInterpolate.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
binaryDjk
(
    const label i,
    const label j,
    const volScalarField& p,
    const volScalarField& T
) const
{
    const kineticTheoryMixture& mixture =
        dynamicCast<const kineticTheoryMixture>(this->thermo());

    const speciesTable& species = this->thermo().species();

    tmp<volScalarField> tD
    (
        volScalarField::New
        (
            "Djk(" + species[i] + ',' + species[j] + ')',
            p.mesh(),
            dimensionedScalar(dimKinematicViscosity, 0)
        )
    );
    volScalarField& D = tD.ref();

    // Внутренние ячейки
    forAll(p, celli)
    {
        D[celli] = mixture.Djk(i, j, p[celli], T[celli]);
    }

    // Граничные грани — ОБЯЗАТЕЛЬНО заполнять, иначе boundaryField останется
    // нулём и после преобразования Maxwell-Stefan даст вырожденную матрицу.
    forAll(p.boundaryField(), patchi)
    {
        scalarField& bf = D.boundaryFieldRef()[patchi];
        const scalarField& bp = p.boundaryField()[patchi];
        const scalarField& bT = T.boundaryField()[patchi];

        forAll(bf, facei)
        {
            bf[facei] = mixture.Djk(i, j, bp[facei], bT[facei]);
        }
    }

    return tD;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
transformDiffusionCoefficient() const
{
    const label d = this->thermo().defaultSpecie();

    // Calculate the molecular weight of the mixture and the mole fractions
    scalar Wm = 0;

    forAll(W, i)
    {
        X[i] = Y[i]/W[i];
        Wm += X[i];
    }

    Wm = 1/Wm;
    X *= Wm;

    // i counter for the specie sub-system without the default specie
    label is = 0;

    // Calculate the A and B matrices from the binary mass diffusion
    // coefficients and specie mole fractions
    forAll(X, i)
    {
        if (i != d)
        {
            A(is, is) = -X[i]*Wm/(DD(i, d)*W[d]);
            B(is, is) = -(X[i]*Wm/W[d] + (1 - X[i])*Wm/W[i]);

            // j counter for the specie sub-system without the default specie
            label js = 0;

            forAll(X, j)
            {
                if (j != i)
                {
                    A(is, is) -= X[j]*Wm/(DD(i, j)*W[i]);

                    if (j != d)
                    {
                        A(is, js) =
                            X[i]*(Wm/(DD(i, j)*W[j]) - Wm/(DD(i, d)*W[d]));

                        B(is, js) = X[i]*(Wm/W[j] - Wm/W[d]);
                    }
                }

                if (j != d)
                {
                    js++;
                }
            }

            is++;
        }
    }

    // LU decompose A and invert
    A.decompose();
    A.inv(invA);

    // Calculate the generalised Fick's law diffusion coefficients
    multiply(D_, invA, B);
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
transformDiffusionCoefficientFields() const
{
    const label d = this->thermo().defaultSpecie();

    // For each cell or patch face
    forAll(*(YPtrs[0]), pi)
    {
        forAll(W, i)
        {
            // Map YPtrs -> Y
            Y[i] = (*YPtrs[i])[pi];

            // Map DijPtrs -> DD
            forAll(W, j)
            {
                DD(i, j) = (*DijPtrs[i][j])[pi];
            }
        }

        // Transform DD -> D
        transformDiffusionCoefficient();

        // i counter for the specie sub-system without the default specie
        label is = 0;

        forAll(W, i)
        {
            if (i != d)
            {
                // j counter for the specie sub-system
                // without the default specie
                label js = 0;

                // Map D -> DijPtrs
                forAll(W, j)
                {
                    if (j != d)
                    {
                        (*DijPtrs[i][j])[pi] = D_(is, js);

                        js++;
                    }
                }

                is++;
            }
        }
    }
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
transform
(
    List<PtrList<volScalarField>>& Dij
) const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& Y0 = Y[0];

    forAll(W, i)
    {
        // Map this->thermo().Y() internal fields -> YPtrs
        YPtrs[i] = &Y[i].primitiveField();

        // Map Dii_ internal fields -> DijPtrs
        DijPtrs[i][i] = &Dii_[i].primitiveFieldRef();

        // Map Dij internal fields -> DijPtrs
        forAll(W, j)
        {
            if (j != i)
            {
                DijPtrs[i][j] = &Dij[i][j].primitiveFieldRef();
            }
        }
    }

    // Transform binary mass diffusion coefficients internal field DijPtrs ->
    // generalised Fick's law diffusion coefficients DijPtrs
    transformDiffusionCoefficientFields();

    forAll(Y0.boundaryField(), patchi)
    {
        forAll(W, i)
        {
            // Map this->thermo().Y() patch fields -> YPtrs
            YPtrs[i] = &Y[i].boundaryField()[patchi];

            // Map Dii_ patch fields -> DijPtrs
            DijPtrs[i][i] = &Dii_[i].boundaryFieldRef()[patchi];

            // Map Dij patch fields -> DijPtrs
            forAll(W, j)
            {
                if (j != i)
                {
                    DijPtrs[i][j] = &Dij[i][j].boundaryFieldRef()[patchi];
                }
            }
        }

        // Transform binary mass diffusion coefficients patch field DijPtrs ->
        // generalised Fick's law diffusion coefficients DijPtrs
        transformDiffusionCoefficientFields();
    }
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
updateDii() const
{
    const label d = this->thermo().defaultSpecie();

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField& rho = this->momentumTransport().rho();

    Dii_.setSize(Y.size());
    jexp_.setSize(Y.size());

    List<PtrList<volScalarField>> Dij(Y.size());

    // Оценить бинарные коэффициенты массовой диффузии из кинетической
    // теории и инициализировать явную часть потока массы.
    forAll(Y, i)
    {
        if (i != d)
        {
            if (jexp_.set(i))
            {
                jexp_[i] = Zero;
            }
            else
            {
                jexp_.set
                (
                    i,
                    surfaceScalarField::New
                    (
                        "jexp" + Y[i].name(),
                        Y[i].mesh(),
                        dimensionedScalar(dimensionSet(1, -2, -1, 0, 0), 0)
                    )
                );
            }
        }

        Dii_.set(i, binaryDjk(i, i, p, T).ptr());

        Dij[i].setSize(Y.size());

        forAll(Y, j)
        {
            if (j > i)
            {
                Dij[i].set(j, binaryDjk(i, j, p, T).ptr());
            }
            else if (j < i)
            {
                Dij[i].set(j, Dij[j][i].clone());
            }
        }
    }

    // Преобразовать бинарные коэффициенты в коэффициенты обобщённого
    // закона Фика.
    transform(Dij);

    // Накопить явную часть потока массы
    forAll(Y, j)
    {
        if (j != d)
        {
            const surfaceScalarField snGradYj(fvc::snGrad(Y[j]));

            forAll(Y, i)
            {
                if (i != d && i != j)
                {
                    jexp_[i] -= fvc::interpolate(rho*Dij[i][j])*snGradYj;
                }
            }
        }
    }

    // Вычислить эффективную диффузию видов, если требуется
    if (!Di_.empty())
    {
        Di_.setSize(Y.size());

        forAll(Y, i)
        {
            // Начать с диагонального коэффициента
            Di_.set(i, Dii_[i].clone().ptr());

            // Массовая доля подсистемы без этого вида
            const volScalarField sumYnoti(max(1 - Y[i], small));

            // Вычесть недиагональные коэффициенты, взвешенные их долей
            // в подсистеме
            forAll(Y, j)
            {
                if (j != i)
                {
                    Di_[i] -= Y[j]/sumYnoti*Dij[i][j];
                }
            }
        }
    }
}


template<class laminarThermophysicalTransportModel>
const PtrList<volScalarField>&
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
Dii() const
{
    if (!Dii_.size())
    {
        updateDii();
    }

    return Dii_;
}


template<class laminarThermophysicalTransportModel>
const PtrList<surfaceScalarField>&
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
jexp() const
{
    if (!jexp_.size())
    {
        updateDii();
    }

    return jexp_;
}


template<class laminarThermophysicalTransportModel>
const PtrList<volScalarField>&
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
Di() const
{
    if (!Di_.size())
    {
        const PtrList<volScalarField>& Y = this->thermo().Y();
        Di_.setSize(Y.size()); // <-- Mark as to be computed from now on
        updateDii();
    }

    return Dii_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
kineticTheoryMaxwellStefanFourier
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
    ),

    TopoChangeableMeshObject<fvMesh>(*this),

    W(this->thermo().species().size()),

    YPtrs(W.size()),
    DijPtrs(W.size()),

    Y(W.size()),
    X(W.size()),
    DD(W.size()),
    A(W.size() - 1),
    B(A.m()),
    invA(A.m()),
    D_(W.size())
{
    // Set the molecular weights of the species
    forAll(W, i)
    {
        W[i] = this->thermo().Wi(i).value();
    }

    read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
bool kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
read()
{
    // Никаких бинарных D из словаря: они вычисляются на лету из
    // кинетической теории (kineticTheoryMulticomponentMixture::Djk).
    return unityLewisFourier<laminarThermophysicalTransportModel>::read();
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
D
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "D",
        this->momentumTransport().rho()
       *Di()[this->thermo().specieIndex(Yi)]
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
D
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Di()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


template<class laminarThermophysicalTransportModel>
tmp<volScalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()*Dii()[this->thermo().specieIndex(Yi)]
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Dii()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
q() const
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

    const label d = this->thermo().defaultSpecie();

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
            if (i != d)
            {
                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(this->j(Y[i]));

                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = d;

            const volScalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
q
(
    const label patchi
) const
{
    tmp<scalarField> tmpq
    (
      - (
            this->alpha().boundaryField()[patchi]
           *this->kappaEff(patchi)
           *this->thermo().T().boundaryField()[patchi]
        )
    );

    const label d = this->thermo().defaultSpecie();

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    if (Y.size())
    {
        scalarField sumJ(tmpq->size(), scalar(0));
        scalarField sumJh(tmpq->size(), scalar(0));

        forAll(Y, i)
        {
            if (i != d)
            {
                const scalarField hi(this->thermo().hsi(i, p, T));

                const scalarField ji(this->j(Y[i], patchi));

                sumJ += ji;

                sumJh += ji*hi;
            }
        }

        {
            const label i = d;

            const scalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*hi;
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
divq
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

    const label d = this->thermo().defaultSpecie();

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
        if (i != d)
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));

            const surfaceScalarField ji(this->j(Y[i]));

            sumJ += ji;

            sumJh += ji*fvc::interpolate(hi);
        }
    }

    {
        const label i = d;

        const volScalarField hi(this->thermo().hsi(i, p, T));

        sumJh -= sumJ*fvc::interpolate(hi);
    }

    tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());

    return tmpDivq;
}


template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
j
(
    const volScalarField& Yi
) const
{
    const label d = this->thermo().defaultSpecie();

    if (this->thermo().specieIndex(Yi) == d)
    {
        const PtrList<volScalarField>& Y = this->thermo().Y();

        tmp<surfaceScalarField> tjd
        (
            surfaceScalarField::New
            (
                IOobject::groupName
                (
                    "j" + name(d),
                    this->thermo().phaseName()
                ),
                Yi.mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField& jd = tjd.ref();

        forAll(Y, i)
        {
            if (i != d)
            {
                jd -= this->j(Y[i]);
            }
        }

        return tjd;
    }
    else
    {
        return
            unityLewisFourier<laminarThermophysicalTransportModel>::j(Yi)
          + jexp()[this->thermo().specieIndex(Yi)];
    }
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const label d = this->thermo().defaultSpecie();

    if (this->thermo().specieIndex(Yi) == d)
    {
        const PtrList<volScalarField>& Y = this->thermo().Y();

        tmp<scalarField> tjd
        (
            new scalarField(Yi.boundaryField()[patchi].size(), scalar(0))
        );
        scalarField& jd = tjd.ref();

        forAll(Y, i)
        {
            if (i != d)
            {
                jd -= this->j(Y[i], patchi);
            }
        }

        return tjd;
    }
    else
    {
        return
            unityLewisFourier<laminarThermophysicalTransportModel>::
            j(Yi, patchi)
          + jexp()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
    }
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix>
kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
divj
(
    volScalarField& Yi
) const
{
    return
        unityLewisFourier<laminarThermophysicalTransportModel>::divj(Yi)
      + fvc::div(jexp()[this->thermo().specieIndex(Yi)]*Yi.mesh().magSf());
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
predict()
{
    unityLewisFourier<laminarThermophysicalTransportModel>::predict();
    updateDii();
}


template<class laminarThermophysicalTransportModel>
bool kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
movePoints()
{
    return true;
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
topoChange
(
    const polyTopoChangeMap& map
)
{
    // Delete the cached Dii and jexp, will be re-created in predict
    Dii_.clear();
    jexp_.clear();

    // Clear Di elements (if any) individually so the entire list is still
    // marked as to be computed from now on
    forAll(Di_, i)
    {
        Di_.set(i, nullptr);
    }
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
mapMesh
(
    const polyMeshMap& map
)
{
    Dii_.clear();
    jexp_.clear();

    forAll(Di_, i)
    {
        Di_.set(i, nullptr);
    }
}


template<class laminarThermophysicalTransportModel>
void kineticTheoryMaxwellStefanFourier<laminarThermophysicalTransportModel>::
distribute
(
    const polyDistributionMap& map
)
{
    Dii_.clear();
    jexp_.clear();

    forAll(Di_, i)
    {
        Di_.set(i, nullptr);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
