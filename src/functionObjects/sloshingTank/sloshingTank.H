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
    Foam::functionObjects::sloshingTank

Group
    grpForcesFunctionObjects

Description
    Computes the center of mass for a fluid. To be used with interFoam solver

Usage

    CoM
    {
        type            sloshingTank;
        libs            (oscfdFunctionObjects);
        writeControl   timeStep;
        writeInterval  1;
        log            true; // if you want to see the result on the terminal
        writeToFile    yes;
    }

SourceFiles
    sloshingTank.C


TAI Opensource CFD Center
Saleh Abuhanieh
Sahin Yigit
Baris Cumhur

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_centerOfMass_H
#define functionObjects_centerOfMass_H


#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Vector.H"
#include "vectorField.H"
#include "surfaceMesh.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "solidBodyMotionFunction.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sloshingTank Declaration
\*---------------------------------------------------------------------------*/

class sloshingTank
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:
    septernion fixed_ref_transformation_;

    // Protected Data

    vectorField fixed_ref_cc_;

    vector com_;
    vector com_fixed_;

    tensor mmi_cm_;
    tensor mmi_fixed_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

public:

    //- Runtime type information
    TypeName("sloshingTank");


    // Constructors

        //- Construct from Time and dictionary
        sloshingTank
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        sloshingTank(const sloshingTank&) = delete;

        //- No copy assignment
        void operator=(const sloshingTank&) = delete;


    //- Destructor
    virtual ~sloshingTank() = default;


    // Member Functions

        //- Read the sloshingTank data
    virtual bool read(const dictionary&);

        virtual bool execute();


        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
