/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    foamToTetDualMesh

Group
    grpPostProcessingUtilities

Description
    printCellID.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "Time.H"
#include "IOobjectList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "printCellID"
    );


    #include "setRootCase.H"
    #include "createTime.H"

    // Read the mesh
    #include "createNamedMesh.H"

    // Read the tetDualMesh
    Info<< "printCellID for time = "
        << runTime.timeName() << nl << endl;

    volScalarField cellID
    (
        IOobject
        (
           "cellID",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ//,
            //IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    forAll (mesh.C(), celli)
    {
        cellID[celli] = celli;
    }

    cellID.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
