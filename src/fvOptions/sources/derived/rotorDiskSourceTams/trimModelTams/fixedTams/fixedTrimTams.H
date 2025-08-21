/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fixedTrimTams

Description
    Trim model where the operating characteristics of rotor
    (e.g. blade pitch angle) are fixed, hence the possibility
    of varying thrust and torque.

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSource1
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (runtime modifiable)
        trimModel    fixedTrimTams;

        fixedTrimTamsCoeffs
        {
            // Conditional mandatory entries (runtime modifiable)

                // when trimModel=fixedTrimTams
                theta0   5.0;
                theta1c  1.0;
                theta1s  4.0;
        }
    }

See also
  - Foam::fv::rotorDiskSource
  - Foam::trimModel
  - Foam::targetCoeffTrim

SourceFiles
    fixedTrimTams.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTrimTams_H
#define fixedTrimTams_H

#include "trimModelTams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fixedTrimTams Declaration
\*---------------------------------------------------------------------------*/

class fixedTrimTams
:
    public trimModelTams
{
protected:

    // Protected Data

        //- Geometric angle of attack [rad]
        scalarField thetag_;


public:

    //- Run-time type information
    TypeName("fixedTrimTams");


    // Constructors

        //- Constructor from rotor and dictionary
        fixedTrimTams(const fv::rotorDiskSourceTams& rotor, const dictionary& dict);

        //- No copy construct
        fixedTrimTams(const fixedTrimTams&) = delete;

        //- No copy assignment
        void operator=(const fixedTrimTams&) = delete;


    //- Destructor
    virtual ~fixedTrimTams() = default;


    // Member Functions

        //- Read
        void read(const dictionary& dict);

        //- Return the geometric angle of attack [rad]
        virtual tmp<scalarField> thetag() const;

        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        );

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
