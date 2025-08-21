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
    Foam::trimModelTams

Description
    Base class for trim models for handling
    blade characteristics and thrust-torque relations.

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSourceTams1
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (runtime modifiable)
        trimModelTams    <trimModelTams>;

        <trimModelTams>Coeffs
        {
            // Conditional mandatory entries (runtime modifiable)

                // when trimModelTams=fixedTrim
                theta0   5.0;
                theta1c  1.0;
                theta1s  4.0;
        }
    }

    where the entries mean:
    \table
      Property   | Description                          | Type  | Reqd | Dflt
      trimModelTams  | Trim model specification             | word  | yes  | -
      file  | Name of file containing profile characteristics | word   | no  | -
    \endtable

    Options for the \c trimModelTams entry:
    \verbatim
      fixed       | Fixed blade/angle characteristics
      targetForce | Target thrust/torque by using time-variant blade characs.
    \endverbatim

See also
  - Foam::fv::rotorDiskSourceTams
  - Foam::fixedTrim
  - Foam::targetCoeffTrim

SourceFiles
    trimModelTams.C

\*---------------------------------------------------------------------------*/

#ifndef trimModelTams_H
#define trimModelTams_H

#include "rotorDiskSourceTams.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class trimModelTams Declaration
\*---------------------------------------------------------------------------*/

class trimModelTams
{
protected:

    // Protected Data

        //- Reference to the rotor source model
        const fv::rotorDiskSourceTams& rotor_;

        //- Name of model
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;


public:

    //- Run-time type information
    TypeName("trimModelTams");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            trimModelTams,
            dictionary,
            (
                const fv::rotorDiskSourceTams& rotor,
                const dictionary& dict
            ),
            (rotor, dict)
        );


    // Constructors

        //- Construct from components
        trimModelTams
        (
            const fv::rotorDiskSourceTams& rotor,
            const dictionary& dict,
            const word& name
        );


    // Selectors

        //- Return a reference to the selected trim model
        static autoPtr<trimModelTams> New
        (
            const fv::rotorDiskSourceTams& rotor,
            const dictionary& dict
        );


    //- Destructor
    virtual ~trimModelTams() = default;


    // Member Functions

        //- Read
        virtual void read(const dictionary& dict);

        //- Return the geometric angle of attack [rad]
        virtual tmp<scalarField> thetag() const = 0;

        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        ) = 0;

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
