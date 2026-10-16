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

Class
    Foam::levelSetCoupledMixture

Description
    Immiscible, incompressible two-phase mixture whose volume fraction is
    carried by the VOF (MULES) transport, and which rebuilds a signed-distance
    level set from that volume fraction every time it is corrected.

    The level set supplies a smooth Heaviside function and an interface
    curvature. These can replace the volume-fraction based values in the
    surface tension force and in the density and viscosity blending. The
    coupling follows the simple coupled level-set/VOF idea of Albadawi et al.
    (2013), Int. J. Multiphase Flow 53:11-28.

    Rebuild, with h the reference cell size:
    \verbatim
        seed        d0  = (2 alpha1 - 1) seedCoeff h
        redistance  d  <- d + dTau sgn(d0) (1 - |grad d|),
                    dTau = pseudoTimeCoeff h,  nRedistanceIter times
        Heaviside   H   = (1 + x + sin(pi x)/pi)/2,
                    x   = clamp(d/eps, -1, 1),  eps = widthCoeff h
        curvature   kappa = -div(n_f),  n_f = (grad d)_f/|(grad d)_f|
    \endverbatim

    All settings are optional and are read from the \c levelSet
    sub-dictionary of \c constant/transportProperties:
    \verbatim
    levelSet
    {
        cellSize         0.00625;   // [m] default: smallest cell spacing
        widthCoeff       1.5;       // interface half-width / cellSize
        seedCoeff        0.75;      // |seed| / cellSize
        pseudoTimeCoeff  0.1;       // pseudo-time step / cellSize
        nRedistanceIter  15;        // default: ceil(widthCoeff/pseudoTimeCoeff)
        curvatureFrom    levelSet;  // levelSet | volumeFraction
        propertiesFrom   levelSet;  // levelSet | volumeFraction
    }
    \endverbatim

    The surface tension coefficient is the \c sigma entry that the
    interfaceProperties base class already reads.

SourceFiles
    levelSetCoupledMixture.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_levelSetCoupledMixture_H
#define Foam_levelSetCoupledMixture_H

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "surfaceTensionModel.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class levelSetCoupledMixture Declaration
\*---------------------------------------------------------------------------*/

class levelSetCoupledMixture
:
    public immiscibleIncompressibleTwoPhaseMixture
{
public:

    // Public Data Types

        //- Which interface description feeds a given coupling term
        enum class couplingSource : unsigned char
        {
            volumeFraction,
            levelSet
        };

        //- Names for couplingSource
        static const Enum<couplingSource> couplingSourceNames;


private:

    // Private Data

        //- Reference cell size h [m]
        scalar cellSize_;

        //- Interface half-width in units of h
        scalar widthCoeff_;

        //- Magnitude of the seed distance in units of h
        scalar seedCoeff_;

        //- Pseudo-time step of the redistance iterations in units of h
        scalar pseudoTimeCoeff_;

        //- Number of redistance iterations
        label nRedistanceIter_;

        //- Source of the curvature in the surface tension force
        couplingSource curvatureFrom_;

        //- Source of the density and viscosity blending
        couplingSource propertiesFrom_;

        //- Surface tension coefficient model
        autoPtr<surfaceTensionModel> sigmaModel_;

        //- Signed distance to the interface [m], positive in phase 1
        volScalarField levelSet_;

        //- Smoothed Heaviside function of the signed distance
        volScalarField heaviside_;

        //- Interface curvature from the signed distance [1/m]
        volScalarField curvature_;


    // Private Member Functions

        //- Read the coupling settings from the levelSet sub-dictionary
        void readCoupling();

        //- Smallest cell-centre spacing of the mesh
        scalar minCellSpacing() const;

        //- Interface half-width as a dimensioned length
        dimensionedScalar halfWidth() const;

        //- Seed the signed distance from alpha1 and redistance it
        void redistance();

        //- Update the Heaviside function from the signed distance
        void updateHeaviside();

        //- Update the curvature from the signed distance
        void updateCurvature();

        //- Overwrite the kinematic viscosity with the Heaviside blend
        void blendViscosity();

        //- No copy construct
        levelSetCoupledMixture(const levelSetCoupledMixture&) = delete;

        //- No copy assignment
        void operator=(const levelSetCoupledMixture&) = delete;


public:

    //- Runtime type information
    TypeName("levelSetCoupledMixture");


    // Constructors

        //- Construct from velocity and volumetric flux
        levelSetCoupledMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~levelSetCoupledMixture() = default;


    // Member Functions

        //- Signed distance field
        const volScalarField& levelSet() const noexcept
        {
            return levelSet_;
        }

        //- Smoothed Heaviside field
        const volScalarField& heaviside() const noexcept
        {
            return heaviside_;
        }

        //- Level-set curvature field
        const volScalarField& curvature() const noexcept
        {
            return curvature_;
        }

        //- True if density and viscosity are blended with the Heaviside
        bool propertiesFromLevelSet() const noexcept
        {
            return propertiesFrom_ == couplingSource::levelSet;
        }

        //- Smoothed Dirac delta of the signed distance [1/m]
        tmp<volScalarField> dirac() const;

        //- Mixture density blended with the Heaviside
        tmp<volScalarField> blendedDensity() const;

        //- Mass flux consistent with blendedDensity()
        tmp<surfaceScalarField> blendedMassFlux
        (
            const surfaceScalarField& phi
        ) const;

        //- Face surface tension force density [N/m^3], to be multiplied by
        //  the face area. Uses the curvature selected by curvatureFrom.
        tmp<surfaceScalarField> capillaryForce() const;

        //- Rebuild the level set from the current alpha1: redistance,
        //  Heaviside, curvature and, if selected, the blended viscosity.
        //  Call once per outer corrector, after the alpha equation.
        void rebuildLevelSet();

        //- Correct the VOF properties and re-apply the viscosity blend.
        //  Cheap: the alpha equation calls it after every MULES correction.
        virtual void correct();

        //- Re-read transportProperties, including the levelSet settings
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
