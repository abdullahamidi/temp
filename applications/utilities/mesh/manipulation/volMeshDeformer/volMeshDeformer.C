/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License

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
    volMeshDeformer

Description
    This utility is a stand-alone volume mesh deformer. Its mainly
    based on the mesh deformation functionalities provided by
    the adjointOptimisationFoam library. Tested for single and
    multiple deformation boxes.

    Saleh Abuhanieh {04.2023}

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "NURBS3DVolume.H"
#include "optMeshMovement.H"
#include "volumetricBSplinesMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (runTime.value() > 0.0)
    {
        FatalErrorInFunction
            << "Make sure that you start the defomration process from time 0 to avoid  \n"
            << "accidentally deforming the mesh again!. Delete the time folders from your case."
            << exit(FatalError);
    }

    const dictionary dict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    labelHashSet patches
    (
        mesh.boundaryMesh().patchSet
        (
            dict.subDict("meshMovement").get<wordRes>("patches")
        )
    );
    if (patches.empty())
    {
        FatalErrorInFunction
            << "No valid patche(s) were selected to deform."
            << exit(FatalError);
    }
    labelList sensitivityPatchIDs = patches.toc();

    autoPtr<optMeshMovement> optMeshMovement;

    optMeshMovement.reset
    (
        optMeshMovement::New
        (
            mesh,
            dict.subDict("meshMovement"),
            sensitivityPatchIDs
        ).ptr()
    );

    const dictionary NURBSdict(dict.subDict("volumetricBSplinesMotionSolverCoeffs"));

    wordList controlBoxes(NURBSdict.toc());
    IOdictionary dict_deformation
    (
        IOobject
        (
            "deformationInputs",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    List<vectorField> noDeformationBoxes;

    for (const word& boxName : controlBoxes)
    {
        if (NURBSdict.isDict(boxName))
        {
            vectorField xyz(dict_deformation.lookup(boxName));
            noDeformationBoxes.append(xyz);
        }
    }

    List<vector> allDeformations;

    forAll(noDeformationBoxes, boxNo)
    {
        forAll(noDeformationBoxes[boxNo], i)
        {
            allDeformations.append(noDeformationBoxes[boxNo][i]);
        }
    }

    Info<<"Total number of boxes: "<<noDeformationBoxes.size()<<endl;
    Info<<"Total number of available control points: "<<allDeformations.size()<<endl;
    Info<<"Total number of available design variables: "<<(allDeformations.size())*3<<endl;

    vectorField cpMovement_temp(allDeformations.size(), Zero);  
    forAll(allDeformations, i)
    {
        cpMovement_temp[i] = allDeformations[i];
    }  

    const vectorField& cpMovement(cpMovement_temp);
    label activeVariables = 0;
    forAll(cpMovement, v)
    {
        //Info<<"cpMovement[v][0] : "<< cpMovement[v][0]<<"  cpMovement[v][1] : "<< cpMovement[v][1]<<"  cpMovement[v][2] : "<< cpMovement[v][2]<<endl;
        if(mag(cpMovement[v][0]) > SMALL)
        {
            activeVariables++;
        }
        if(mag(cpMovement[v][1]) > SMALL)
        {
            activeVariables++;
        }
        if(mag(cpMovement[v][2]) > SMALL)
        {
            activeVariables++;
        }
    }
    Info<<"Total number of active design variables: "<<activeVariables<<endl;

    autoPtr<displacementMethod>& displMethodPtr = optMeshMovement->returnDisplacementMethod();
    autoPtr<motionSolver>& motionPtr = displMethodPtr->getMotionSolver();
    refCast<volumetricBSplinesMotionSolver>(motionPtr()).setControlPointsMovement(cpMovement);

    // Move mesh
    displMethodPtr->update();


    runTime++;

    Info<< "  Writing new mesh points " << endl;

        pointIOField points
        (
            IOobject
            (
               "points",
                mesh.pointsInstance(),
                mesh.meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh.points()
        );
        points.write();

    // Check mesh quality
    mesh.checkMesh(true);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

