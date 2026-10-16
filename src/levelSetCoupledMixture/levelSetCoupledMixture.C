/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2026 Amirhossein Taran
-------------------------------------------------------------------------------
License
    This file is part of the LevelSet-VOF solver, built on OpenFOAM.

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "levelSetCoupledMixture.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(levelSetCoupledMixture, 0);
}


const Foam::Enum<Foam::levelSetCoupledMixture::couplingSource>
Foam::levelSetCoupledMixture::couplingSourceNames
({
    { couplingSource::volumeFraction, "volumeFraction" },
    { couplingSource::levelSet, "levelSet" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::levelSetCoupledMixture::minCellSpacing() const
{
    // Cell-centre spacing across internal faces. Unlike cbrt(V) this is
    // also right for 2-D and axisymmetric meshes of any thickness.
    return gMin(mag(levelSet_.mesh().delta())().primitiveField());
}


void Foam::levelSetCoupledMixture::readCoupling()
{
    const dictionary coeffs(subOrEmptyDict("levelSet"));

    cellSize_ = coeffs.getOrDefault<scalar>("cellSize", minCellSpacing());
    widthCoeff_ = coeffs.getOrDefault<scalar>("widthCoeff", 1.5);
    seedCoeff_ = coeffs.getOrDefault<scalar>("seedCoeff", 0.75);
    pseudoTimeCoeff_ = coeffs.getOrDefault<scalar>("pseudoTimeCoeff", 0.1);

    if (cellSize_ <= 0 || widthCoeff_ <= 0 || pseudoTimeCoeff_ <= 0)
    {
        FatalIOErrorInFunction(coeffs)
            << "cellSize, widthCoeff and pseudoTimeCoeff must be positive"
            << exit(FatalIOError);
    }

    // Enough pseudo-time to push the distance a full half-width away from
    // the interface, so cells outside the band always see |d| > eps
    nRedistanceIter_ = coeffs.getOrDefault<label>
    (
        "nRedistanceIter",
        label(std::ceil(widthCoeff_/pseudoTimeCoeff_ - SMALL))
    );

    curvatureFrom_ = couplingSourceNames.getOrDefault
    (
        "curvatureFrom",
        coeffs,
        couplingSource::levelSet
    );

    propertiesFrom_ = couplingSourceNames.getOrDefault
    (
        "propertiesFrom",
        coeffs,
        couplingSource::levelSet
    );

    Info<< "Level-set coupling:" << nl
        << "    cell size       " << cellSize_ << nl
        << "    half-width      " << widthCoeff_*cellSize_ << nl
        << "    redistance      " << nRedistanceIter_ << " x "
        << pseudoTimeCoeff_*cellSize_ << nl
        << "    curvature from  " << couplingSourceNames[curvatureFrom_] << nl
        << "    properties from " << couplingSourceNames[propertiesFrom_]
        << nl << endl;
}


Foam::dimensionedScalar Foam::levelSetCoupledMixture::halfWidth() const
{
    return dimensionedScalar(dimLength, widthCoeff_*cellSize_);
}


void Foam::levelSetCoupledMixture::redistance()
{
    const dimensionedScalar h(dimLength, cellSize_);
    const dimensionedScalar dTau(pseudoTimeCoeff_*h);

    const volScalarField seed((2.0*alpha1() - 1.0)*seedCoeff_*h);

    // Sign of the seed; zero on cells sitting exactly on the interface
    const volScalarField seedSign
    (
        seed/max(mag(seed), dimensionedScalar(dimLength, VSMALL))
    );

    levelSet_ == seed;
    levelSet_.correctBoundaryConditions();

    for (label iter = 0; iter < nRedistanceIter_; ++iter)
    {
        levelSet_ += seedSign*(1.0 - mag(fvc::grad(levelSet_)))*dTau;
        levelSet_.correctBoundaryConditions();
    }

    // Report how far |grad d| is from 1 inside the interface band
    const scalarField& V = levelSet_.mesh().V();
    const scalarField inBand
    (
        pos0(halfWidth().value() - mag(levelSet_.primitiveField()))
    );
    const scalarField gradError
    (
        mag(1.0 - mag(fvc::grad(levelSet_)().primitiveField()))
    );
    const scalar bandVolume = gSum(inBand*V);

    if (bandVolume > 0)
    {
        Info<< "Level set: mean |1 - |grad d|| in band = "
            << gSum(inBand*gradError*V)/bandVolume << endl;
    }
}


void Foam::levelSetCoupledMixture::updateHeaviside()
{
    using constant::mathematical::pi;

    // Clamping x to [-1, 1] makes the smooth formula exact outside the band
    const volScalarField x
    (
        max(min(levelSet_/halfWidth(), scalar(1)), scalar(-1))
    );

    heaviside_ == 0.5*(1.0 + x + sin(pi*x)/pi);
}


void Foam::levelSetCoupledMixture::updateCurvature()
{
    const surfaceVectorField gradDf(fvc::interpolate(fvc::grad(levelSet_)));

    const surfaceScalarField nf
    (
        (gradDf/(mag(gradDf) + dimensionedScalar(dimless, 1e-6)))
      & levelSet_.mesh().Sf()
    );

    curvature_ = -fvc::div(nf);
}


void Foam::levelSetCoupledMixture::blendViscosity()
{
    // Blend the dynamic viscosity, then divide by the blended density, so
    // that rho*nu is the Heaviside average of mu1 and mu2
    const volScalarField& H = heaviside_;

    nu_ =
    (
        H*rho1_*nuModel1_->nu() + (1.0 - H)*rho2_*nuModel2_->nu()
    )/(H*rho1_ + (1.0 - H)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::levelSetCoupledMixture::levelSetCoupledMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    immiscibleIncompressibleTwoPhaseMixture(U, phi),
    cellSize_(0),
    widthCoeff_(0),
    seedCoeff_(0),
    pseudoTimeCoeff_(0),
    nRedistanceIter_(0),
    curvatureFrom_(couplingSource::levelSet),
    propertiesFrom_(couplingSource::levelSet),
    sigmaModel_(surfaceTensionModel::New(*this, U.mesh())),
    levelSet_
    (
        IOobject
        (
            "levelSet",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    heaviside_
    (
        IOobject
        (
            "heaviside",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    curvature_
    (
        IOobject
        (
            "levelSetCurvature",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    )
{
    readCoupling();
    rebuildLevelSet();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::levelSetCoupledMixture::dirac() const
{
    using constant::mathematical::pi;

    const dimensionedScalar eps(halfWidth());

    const volScalarField x
    (
        max(min(levelSet_/eps, scalar(1)), scalar(-1))
    );

    return tmp<volScalarField>::New
    (
        "levelSetDirac",
        (1.0 + cos(pi*x))/(2.0*eps)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::levelSetCoupledMixture::blendedDensity() const
{
    return tmp<volScalarField>::New
    (
        "rhoBlended",
        heaviside_*rho1_ + (1.0 - heaviside_)*rho2_
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::levelSetCoupledMixture::blendedMassFlux
(
    const surfaceScalarField& phi
) const
{
    return phi*(rho2_ + (rho1_ - rho2_)*fvc::interpolate(heaviside_));
}


Foam::tmp<Foam::surfaceScalarField>
Foam::levelSetCoupledMixture::capillaryForce() const
{
    if (curvatureFrom_ == couplingSource::volumeFraction)
    {
        return surfaceTensionForce();
    }

    return
        fvc::interpolate(sigmaModel_->sigma()*curvature_)
       *fvc::snGrad(heaviside_);
}


void Foam::levelSetCoupledMixture::rebuildLevelSet()
{
    redistance();
    updateHeaviside();
    updateCurvature();

    if (propertiesFromLevelSet())
    {
        blendViscosity();
    }
}


void Foam::levelSetCoupledMixture::correct()
{
    immiscibleIncompressibleTwoPhaseMixture::correct();

    if (propertiesFromLevelSet())
    {
        blendViscosity();
    }
}


bool Foam::levelSetCoupledMixture::read()
{
    if (immiscibleIncompressibleTwoPhaseMixture::read())
    {
        sigmaModel_->readDict(*this);
        readCoupling();
        return true;
    }

    return false;
}


// ************************************************************************* //
