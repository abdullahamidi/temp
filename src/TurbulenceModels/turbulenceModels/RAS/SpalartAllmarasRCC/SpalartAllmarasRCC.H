/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::RASModels::SpalartAllmarasRCCRCC

Group
    grpRASTurbulence

Description

This version is according to the 2012 paper by Allmaras, Johnson, and Spalart: 
"Modifications and Clarifications for the Implementation of the Spalart-Allmaras Turbulence Model",
As well as the NASA Turbulence Modeling Resource Webpage on the Spalart-Allmaras Turbulence Model.
In addition, a Rotation/Curvature Correction (RCC) method is applied, according to the 2012 paper by Qiang and Yong:
"A New, Simpler Rotation/Curvature Correction Method for the Spalart-Allmaras Turbulence Model"
The model with RCC applied is referred to as SARCM.

SourceFiles
    SpalartAllmarasRCC.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasRCC_H
#define SpalartAllmarasRCC_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class SpalartAllmarasRCC Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class SpalartAllmarasRCC
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        SpalartAllmarasRCC(const SpalartAllmarasRCC&) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasRCC&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;
            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cs_;

//// S.A {29.11.2022}
            dimensionedScalar cr1_;
            dimensionedScalar cr2_;
            dimensionedScalar cr3_;

            // rotating velocity
            dimensionedVector omga_;
//// S.A {29.11.2022}


        // Fields

            //- Modified kinematic viscosity [m2/s]
            volScalarField nuTilda_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField::Internal& y_;


    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField::Internal> fv2
        (
            const volScalarField::Internal& chi,
            const volScalarField::Internal& fv1
        ) const;

        tmp<volScalarField::Internal> Stilda() const;

        tmp<volScalarField::Internal> fw
        (
            const volScalarField::Internal& Stilda
        ) const;


        //- Update nut with the latest available nuTilda
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasRCC");


    // Constructors

        //- Construct from components
        SpalartAllmarasRCC
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~SpalartAllmarasRCC() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return the (estimated) turbulent kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the (estimated) turbulent kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the (estimated) specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Solve the turbulence equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasRCC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

