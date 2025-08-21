/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::centerOfMass

Group
    grpForcesFunctionObjects

Description
    Computes the center of mass for a fluid. To be used with interFoam solver

Usage

    CoM
    {
        type            centerOfMass;
        libs            (oscfdFunctionObjects);
        writeControl   timeStep;
        writeInterval  1;
        log            true; // if you want to see the result on the terminal
        writeToFile    yes;
    }

SourceFiles
    centerOfMass.C


TAI Opensource CFD Center
Saleh Abuhanieh
Sahin Yigit
Baris Cumhur

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_centerOfMass_H
#define functionObjects_centerOfMass_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class centerOfMass Declaration
\*---------------------------------------------------------------------------*/

class centerOfMass
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    bool calculate_mmi_;

    // Protected Data

    vector com_;

    tensor mmi_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

public:

    //- Runtime type information
    TypeName("centerOfMass");


    // Constructors

        //- Construct from Time and dictionary
        centerOfMass
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        centerOfMass(const centerOfMass&) = delete;

        //- No copy assignment
        void operator=(const centerOfMass&) = delete;


    //- Destructor
    virtual ~centerOfMass() = default;


    // Member Functions

        //- Read the centerOfMass data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
