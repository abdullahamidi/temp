/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    reconstructCyclic

Description
    A postprocessing utility for reconstructing the full domain from
    a cyclic sector (both mesh and fields).
    
Author

    Saleh Abuhanieh, Nov 2023.    
    

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "axisAngleRotation.H"
#include "transformGeometricField.H"
#include "IOobjectList.H"
#include "mergePolyMesh.H"
#include "polyTopoChanger.H"
#include "perfectInterface.H"
#include "zeroGradientFvPatchFields.H"


using namespace Foam;
using namespace Foam::coordinateRotations;

template<class GeoField>
void ReadAndRotateFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const dimensionedTensor& rotT
)
{
    // Objects of field type
    IOobjectList fields(objects.lookupClass<GeoField>());

    forAllConstIters(fields, fieldIter)
    {
        GeoField fld(*fieldIter(), mesh);
        Info<< "    Rotating " << fld.name() << endl;
        transform(fld, rotT, fld);
        fld.write();
    }
}


void rotateFields
(
    const fvMesh& mesh,
    const Time& runTime,
    const tensor& rotationT
)
{
    // Need dimensionedTensor for geometric fields
    const dimensionedTensor rotT(rotationT);

    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime.timeName());

    ReadAndRotateFields<volVectorField>(mesh, objects, rotT);
    ReadAndRotateFields<volSphericalTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<volSymmTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<volTensorField>(mesh, objects, rotT);

    ReadAndRotateFields<surfaceVectorField>(mesh, objects, rotT);
    ReadAndRotateFields<surfaceSphericalTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<surfaceSymmTensorField>(mesh, objects, rotT);
    ReadAndRotateFields<surfaceTensorField>(mesh, objects, rotT);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A postprocessing utility for reconstructing the full domain from a cyclic sector (both mesh and fields)"
    );
    
    timeSelector::addOptions();

    argList::addArgument
    (
        "rotate-angle",
        "Rotate <angle> degrees about <vector> - eg, '((1 0 0) 45)'"
    );
    
    argList::noParallel(); // only serial for now
    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    
    // Get a rotation specification
    List<tensor> rotList(Zero);
    const Tuple2<vector, scalar> rotAxisAngle
    (
        args.get<Tuple2<vector, scalar>>(1)
    );
    const vector& axis = rotAxisAngle.first();
    const scalar angle = rotAxisAngle.second();
    const label numberOfSectors = 360/angle;
    Info<< "Rotating points and fields " << nl
        << "    about " << axis << nl
        << "    angle " << angle << nl  
        << "    No. of sectors " << numberOfSectors << nl;
    rotList.resize(numberOfSectors);
    forAll(rotList, roti)
    {
        rotList[roti] = axisAngle::rotation(axis, angle*roti, true);
    }
    
    #include "createTime.H"
    
    if (args.found("time"))
    {
        if (args["time"] == "constant")
        {
            runTime.setTime(instant(0, "constant"), 0);
        }
        else
        {
            const scalar timeValue = args.get<scalar>("time");
            runTime.setTime(instant(timeValue), 0);
        }
    }
    
    Info<<"Reconstructing the "<<numberOfSectors<<" cyclic sectors (mesh and fields) for time: "<<runTime.timeName()<<endl;
    
    Foam::PtrList<Foam::fvMesh> meshes(numberOfSectors);
        
    forAll(meshes, meshi)
    {
        meshes.set
        (
            meshi,
            new Foam::fvMesh
            (
                Foam::IOobject
                (
                    Foam::polyMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                ),
                false  // Do not initialise
            )
        );
    }
    
    for (auto& mesh : meshes)
    {
        mesh.init(true);  // Initialise all (lower levels and current)
    }

    // start the operations
    forAll(meshes, meshi)
    {
        Info<<"meshi: "<<meshi<<endl;
        pointField& points = const_cast<pointField&>(meshes[meshi].points());
        Info<< "Rotating points of sector "<<meshi<<" by " << rotList[meshi] << endl;
        transform(points, rotList[meshi], points);

        Info<< "Rotating vector/tensor fields" << endl;        
        rotateFields(meshes[meshi], runTime, rotList[meshi]);
    }
    
    mergePolyMesh masterMesh
    (
        IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    // merge meshes   
    Foam::PtrList<Foam::polyMesh> meshesToAdd(numberOfSectors-1);
    forAll(meshesToAdd, meshii)
    {
        masterMesh.addMesh(meshes[meshii+1]);
    }
    masterMesh.merge();
    
    // stitch mesh   
    polyTopoChanger stitcher(masterMesh, IOobject::NO_READ); 
        // Patch names
    word masterPatchName, slavePatchName;
    DynamicList<word> patchNames;
    forAll(masterMesh.boundaryMesh(), patchi)
    {
        if (isA<cyclicPolyPatch>(masterMesh.boundaryMesh()[patchi]))
        {
            patchNames.append(masterMesh.boundaryMesh()[patchi].name());
        }
    }
    if(patchNames.size() == 2)
    {
        masterPatchName = patchNames[0];
        slavePatchName  = patchNames[1];        
    }
    else
    {
        FatalError
            << "Number of cyclic patches must be 2." << endl
            << exit(FatalError);
    }        
        // Zone names
    const word mergePatchName(masterPatchName + slavePatchName);
    const word cutZoneName(mergePatchName + "CutFaceZone");

        // Master/slave patches
    const polyPatch& masterPatch = masterMesh.boundaryMesh()[masterPatchName];

    masterMesh.pointZones().clearAddressing();
    masterMesh.faceZones().clearAddressing();
    masterMesh.cellZones().clearAddressing();

        // Lists of master and slave faces:
    labelList faceIds;

        // Markup master face ids
    faceIds.setSize(masterPatch.size());
    std::iota(faceIds.begin(), faceIds.end(), masterPatch.start());

    stitcher.clear();
    stitcher.setSize(1);

        // Add new (empty) zone for resulting internal faces
    masterMesh.faceZones()
    (
        cutZoneName,
        true // verbose
    ).resetAddressing(std::move(faceIds), false);


        // Add the perfect interface mesh modifier
    stitcher.set
    (
        0,
        new perfectInterface
        (
            "stitchMesh", //"couple" + Foam::name(actioni),
            0,
            stitcher,
            cutZoneName,
            masterPatchName,
            slavePatchName
        )
    );
     
        // Execute all polyMeshModifiers
    autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh(true);
    masterMesh.movePoints(morphMap->preMotionPoints());

    runTime++;    
    masterMesh.write(); // I couldn't convert the mergePolyMesh to fvMesh directly, thus, I write then read from disk

    Foam::fvMesh mesh_combined
    (
        Foam::IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    
    // fields
    Info<< nl << "creating fields for the combined mesh" <<'\n'<< endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh_combined,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_combined,
        dimensionedScalar(dimPressure, Zero)
    );
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh_combined,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_combined,
        dimensionedScalar(dimTemperature, Zero)
    );
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh_combined,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_combined,
        dimensionedScalar(dimDensity, Zero)
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh_combined,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_combined,
        dimensionedVector(dimVelocity, Zero)
    ); 
  
    label total = 0;
    word oldTimeName = Foam::name(runTime.value()-1);
    forAll(meshes, meshi)
    {
        IOobjectList objects(meshes[meshi], oldTimeName); //IOobjectList objects(meshes[meshi], runTime.timeName());
        // volScalarField
        for (const word& fieldName : objects.sortedNames<volScalarField>())
        {
            if(fieldName == "p" || fieldName == "T" || fieldName == "rho")
            {
                volScalarField field(*(objects[fieldName]), meshes[meshi]);
                volScalarField& fieldCombined(mesh_combined.lookupObjectRef<volScalarField>(fieldName));

                forAll(meshes[meshi].C(), i)
                {       
                    fieldCombined[total+i] = field[i];
                }
                forAll(meshes[meshi].boundary(), patchID)
                {
                    if (isA<cyclicPolyPatch>(meshes[meshi].boundaryMesh()[patchID]))
                    {
                        continue;
                    }
                    label patchSize = meshes[meshi].boundary()[patchID].size();
                    forAll (meshes[meshi].boundary()[patchID], facei) 
                    {
                        fieldCombined.boundaryFieldRef()[patchID][facei + meshi*patchSize] = field.boundaryField()[patchID][facei];                      
                    }
                }
            }
        }
        for (const word& fieldName : objects.sortedNames<volVectorField>())
        {
            if(fieldName == "U")
            {
                volVectorField field(*(objects[fieldName]), meshes[meshi]);
                volVectorField& fieldCombined(mesh_combined.lookupObjectRef<volVectorField>(fieldName));

                forAll(meshes[meshi].C(), i)
                {       
                    fieldCombined[total+i] = field[i];
                }
                forAll(meshes[meshi].boundary(), patchID)
                {
                    if (isA<cyclicPolyPatch>(meshes[meshi].boundaryMesh()[patchID]))
                    {
                        continue;
                    }
                    label patchSize = meshes[meshi].boundary()[patchID].size();
                    forAll (meshes[meshi].boundary()[patchID], facei) 
                    {
                        fieldCombined.boundaryFieldRef()[patchID][facei + meshi*patchSize] = field.boundaryField()[patchID][facei];                        
                    }
                }
            }
        }        
        total += meshes[meshi].nCells();   
    }
    

    p.write();
    T.write();   
    rho.write();
    U.write();               
    //runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
