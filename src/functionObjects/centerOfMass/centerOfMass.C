/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "centerOfMass.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(centerOfMass, 0);
    addToRunTimeSelectionTable(functionObject, centerOfMass, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::centerOfMass::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Center of Mass");
    writeCommented(os, "Time");
    writeTabbed(os, "X");
    writeTabbed(os, "Y");
    writeTabbed(os, "Z");
    if (calculate_mmi_){
        writeTabbed(os, "IXX");
        writeTabbed(os, "IXY");
        writeTabbed(os, "IXZ");
        writeTabbed(os, "IYX");
        writeTabbed(os, "IYY");
        writeTabbed(os, "IYZ");
        writeTabbed(os, "IZX");
        writeTabbed(os, "IZY");
        writeTabbed(os, "IZZ");
    }
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::centerOfMass::centerOfMass
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict)
{
    read(dict);

    writeFileHeader(file());

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::centerOfMass::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    calculate_mmi_ = dict.getOrDefault<bool>("calculateMassMomentOfInertia", true);

    return true;
}


bool Foam::functionObjects::centerOfMass::execute()
{

    immiscibleIncompressibleTwoPhaseMixture& mix =
        mesh_.lookupObjectRef<immiscibleIncompressibleTwoPhaseMixture>("transportProperties");

    const scalarField& alpha1 = mix.alpha1().primitiveField();
    const scalarField& alpha2 = mix.alpha2().primitiveField();

    scalar rho1 = mix.rho1().value();
    scalar rho2 = mix.rho2().value();

    const vectorField& CC = mesh_.C().primitiveField();
    const scalarField& V = mesh_.V();

    vector vU = gSum((rho1 * alpha1 + rho2 * alpha2) * V * CC);
    scalar sL = gSum((rho1 * alpha1 + rho2 * alpha2) * V);
    com_ = vU / sL;

    vectorField translatedCC(CC - com_);

    if (calculate_mmi_) {
        tensor temp = gSum((rho1 * alpha1 + rho2 * alpha2) * V * translatedCC * translatedCC);
        mmi_[tensor::XX] = temp[tensor::YY] + temp[tensor::ZZ];
        mmi_[tensor::XY] = -temp[tensor::XY];
        mmi_[tensor::XZ] = -temp[tensor::XZ];

        mmi_[tensor::YX] = -temp[tensor::YX];
        mmi_[tensor::YY] = temp[tensor::XX] + temp[tensor::ZZ];
        mmi_[tensor::YZ] = -temp[tensor::YZ];

        mmi_[tensor::ZX] = -temp[tensor::ZX];
        mmi_[tensor::ZY] = -temp[tensor::ZY];
        mmi_[tensor::ZZ] = temp[tensor::XX] + temp[tensor::YY];
    }

    return true;
}

bool Foam::functionObjects::centerOfMass::write()
{
    Log
     << "Fluid Center of Mass" << endl
     << "Time: " << mesh_.time().timeName() << endl
     << "X: " << com_[vector::X] << endl
     << "Y: " << com_[vector::Y] << endl
     << "Z: " << com_[vector::Z] << endl;
    if (calculate_mmi_) {
        Log 
         << "Ixx: " << mmi_[tensor::XX] << endl
         << "Ixy: " << mmi_[tensor::XY] << endl
         << "Ixz: " << mmi_[tensor::XZ] << endl
         << "Iyx: " << mmi_[tensor::YX] << endl
         << "Iyy: " << mmi_[tensor::YY] << endl
         << "Iyz: " << mmi_[tensor::YZ] << endl
         << "Izx: " << mmi_[tensor::ZX] << endl
         << "Izy: " << mmi_[tensor::ZY] << endl
         << "Izz: " << mmi_[tensor::ZZ] << endl;
    }

    if (Pstream::master())
    {

        file()
            << mesh_.time().timeName()
            << tab << com_[vector::X]
            << tab << com_[vector::Y]
            << tab << com_[vector::Z];

        if (calculate_mmi_) {
            file()
                << tab << mmi_[tensor::XX]
                << tab << mmi_[tensor::XY]
                << tab << mmi_[tensor::XZ]
                << tab << mmi_[tensor::YX]
                << tab << mmi_[tensor::YY]
                << tab << mmi_[tensor::YZ]
                << tab << mmi_[tensor::ZX]
                << tab << mmi_[tensor::ZY]
                << tab << mmi_[tensor::ZZ];
            }
        
            file() << endl;

    }

    return true;
}


// ************************************************************************* //
