/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    sumoSU2ToFoam

Description
    Converts the su2 mesh file as written by SUMO (pentagrow).
    Only 3D mesh which contains tetra and prisims.
    Please make sure that the MARKER_TAG in the su2 file is
    a single string.
    S.A 

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "polyPatch.H"
#include "cellModeller.H"
#include "triFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::validArgs.append("su2 mesh file");

    #include "setRootCase.H"
    #include "createTime.H"

    if (!args.check())
    {
         FatalError.exit();
    }

    std::ifstream su2File(args.args()[1].c_str());

    // skip the first three lines
    string line;
    std::getline(su2File, line);
    std::getline(su2File, line);
    std::getline(su2File, line);

    // read the fourth line (NDIME=int)
    std::getline(su2File, line);
    int indexL4 = line.find_last_not_of("0123456789");
    string resultL4 = line.substr(indexL4 + 1);
    label NDIME = std::stoi(resultL4);
    if (NDIME != 3)
    {
        FatalErrorInFunction
            << "The SU2 Mesh is not Three Dimensional Mesh"
            << exit(FatalError);
    }

    // read the fifth line (NELEM=int)
    std::getline(su2File, line);
    int indexL5 = line.find_last_not_of("0123456789");
    string resultL5 = line.substr(indexL5 + 1);
    label NELEM = std::stoi(resultL5);
    Info<<"3D Mesh With "<<NELEM<<" Cells"<<endl;

    //
    // Read Cells.
    Info<<"Reading Cells ... "<<endl;
    //
    cellShapeList cells(NELEM);

    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& prism = *(cellModeller::lookup("prism"));

    labelList tetPoints(4);
    labelList prismPoints(6);
    label tetNo = 0; label prismNo = 0;

    for (label i=0; i<NELEM; i++)
    {
        label SU2_Element_Type_Code;
        su2File >> SU2_Element_Type_Code;
        if (SU2_Element_Type_Code == 10)
        {
            su2File >> tetPoints[0] >> tetPoints[1] >> tetPoints[2] >> tetPoints[3];
            cells[i] = cellShape(tet, tetPoints);
            tetNo += 1;
        }
        else if (SU2_Element_Type_Code == 13)
        {
            su2File >> prismPoints[0] >> prismPoints[1] >> prismPoints[2] >> prismPoints[3] >> prismPoints[4] >> prismPoints[5];
            cells[i] = cellShape(prism, prismPoints);
            prismNo += 1;
        }
        else
        {
            FatalErrorInFunction
                << "This Tool Can Handle Only Tetras and Prisms"
                << exit(FatalError);
        }
    }
    Info<<"tet Cells = "<<tetNo<<" and "<<"prism Cells = "<<prismNo<<endl;
    Info<<"Reading Cells Done "<<endl;


    // read the line (NPOIN=int)
    string lineP; su2File >> lineP;
    int indexPL = lineP.find_last_not_of("0123456789");
    string resultPL = lineP.substr(indexPL + 1);
    label NPOIN = std::stoi(resultPL);
    Info<<"The Mesh Has "<<NPOIN<<" Points"<<endl;

    //
    // Read nodes.
    Info<<"Reading Points/Vertices ... "<<endl;
    //
    pointField points(NPOIN);

    forAll(points, pointi)
    {
        scalar x,y,z,ID;

        su2File >> x >> y >> z >> ID;

        points[pointi] = point(x, y, z);
    }

    Info<<"Reading Points/Vertices Done "<<endl; //Info<<"points.size()"<<points.size()<<endl; string test; su2File >> test; Info<<test<<endl; 


    // read the line (NMARK=int)
    string lineM; su2File >> lineM;
    int indexML = lineM.find_last_not_of("0123456789");
    string resultML = lineM.substr(indexML + 1);
    label NMARK = std::stoi(resultML);
    Info<<"The Mesh Has "<<NMARK<<" Boundary Patches"<<endl;

    //
    // Read nodes.
    Info<<"Reading Boundaries ... "<<endl;
    //
    wordList patchNames(NMARK); //This is a list of patch Names
    DynamicList<label> boundaryPatch; //This is a list contains the patchID for each boundary face
    label numberOfBoundaryFaces = 0;
    DynamicList<triFace> boundaryTraingles;
    for (label i=0; i<NMARK; i++)
    {
        // read the line (MARKER_TAG=string)
        string line1; su2File >> line1;
        int index1 = line1.find("="); string MARKER_TAG = line1.substr(index1+1); 
        patchNames[i] = MARKER_TAG; //Info<<patchNames[i]<<endl;
        // read the line (MARKER_ELEMS=int)
        string line2; su2File >> line2;
        int index2 = line2.find("="); string MARKER_ELEMS_string = line2.substr(index2+1);
        label MARKER_ELEMS = std::stoi(MARKER_ELEMS_string);  //Info<<MARKER_ELEMS<<endl;
        numberOfBoundaryFaces += MARKER_ELEMS;
        // Read the faces of each boundary patch  
        for (label j=0; j<MARKER_ELEMS; j++)
        {
            label SU2_Element_Type_Code,p1,p2,p3;
            su2File >> SU2_Element_Type_Code >> p1 >> p2 >> p3;
            triFace tri(p1, p2, p3);
            boundaryTraingles.append(tri);
            boundaryPatch.append(i);
        }

    }
    Info<<"The Mesh Has "<<numberOfBoundaryFaces<<" Boundary faces"<<endl;
    Info<<"Reading Boundaries Done "<<endl;


    //
    Info<<"Adding Boundary Faces to Patches ... "<<endl;
    //
    // This is a list of patches, and each patch is a list of boundary faces
    faceListList patchFaces(NMARK);

    wordList patchTypes(NMARK, polyPatch::typeName);
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = polyPatch::typeName;
    wordList patchPhysicalTypes(NMARK, polyPatch::typeName);

    // This is a list of all the boundary faces
    faceList boundaryFaces(numberOfBoundaryFaces);
    // Boundary faces as three vertices
    HashTable<label, triFace, Hash<triFace>> vertsToBoundary(numberOfBoundaryFaces); // for changing the orientation of the triangle
    forAll(boundaryFaces, facei)
    {
        boundaryFaces[facei].setSize(3);
        boundaryFaces[facei][0] = boundaryTraingles(facei)[0];
        boundaryFaces[facei][1] = boundaryTraingles(facei)[1];
        boundaryFaces[facei][2] = boundaryTraingles(facei)[2];

        vertsToBoundary.insert(boundaryTraingles(facei), facei);
    }
/*
    // loop over the tet cells
    for (label celli=0; celli<tetNo; celli++)//forAll(cells, celli)
    {
        const cellShape& cll = cells[celli];

        // Get the four (outwards pointing) faces of the cell
        faceList tris(cll.faces());

        forAll(tris, i)
        {
            const face& f = tris[i];

            // Is there any boundary face with same vertices?
            // (uses commutative hash)
            HashTable<label, triFace, Hash<triFace>>::iterator iter =
                vertsToBoundary.find(triFace(f[0], f[1], f[2]));

            if (iter != vertsToBoundary.end())
            {
                label facei = iter();
                const triFace& tri = iter.key();

                // Determine orientation of tri v.s. cell centre.
                point cc(cll.centre(points));
                point fc(tri.centre(points));
                vector fn(tri.area(points));

                if (((fc - cc) & fn) < 0)
                {
                    // Boundary face points inwards. Flip.
                    boundaryFaces[facei].flip();
                }

                // Done this face so erase from hash
                vertsToBoundary.erase(iter);
            }
        }
    }

    if (vertsToBoundary.size())
    {
        // Didn't find cells connected to boundary faces.
        WarningInFunction
            << "There are boundary faces without attached cells."
            << "Boundary faces (as triFaces):" << vertsToBoundary.toc()
            << endl;
    }
*/
    // Sort boundaryFaces by patch.
    List<DynamicList<face>> allPatchFaces(NMARK);

    forAll(boundaryPatch, facei)
    {
        label patchi = boundaryPatch[facei];
        allPatchFaces[patchi].append(boundaryFaces[facei]);
    }

    Info<< "Patches:" << nl
        << "\tpatch ID\tPatch name\tSize" << nl
        << "\t--------\t----------\t----" << endl;

    forAll(allPatchFaces, patchi)
    {
        Info<< '\t' << patchi << "\t\t"
            << patchNames[patchi] << "\t\t"
            << allPatchFaces[patchi].size() << endl;

        patchFaces[patchi].transfer(allPatchFaces[patchi]);
    }


    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        std::move(points),
        cells,
        patchFaces,
        patchNames,
        patchTypes,
        defaultFacesName,
        defaultFacesType,
        patchPhysicalTypes
    );

    Info<< "Writing mesh to " << runTime.constant() << endl << endl;

    mesh.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //
