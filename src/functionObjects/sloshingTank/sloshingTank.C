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

#include "sloshingTank.H"
#include "volFields.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "dynamicMotionSolverFvMesh.H"
#include "solidBodyMotionSolver.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sloshingTank, 0);
    addToRunTimeSelectionTable(functionObject, sloshingTank, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::sloshingTank::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Center of Mass");
    writeCommented(os, "Time");
    writeTabbed(os, "CoMx");
    writeTabbed(os, "CoMy");
    writeTabbed(os, "CoMz");
    writeTabbed(os, "CoMx_fixed");
    writeTabbed(os, "CoMy_fixed");
    writeTabbed(os, "CoMz_fixed");
    writeTabbed(os, "IXX_cm");
    writeTabbed(os, "IXY_cm");
    writeTabbed(os, "IXZ_cm");
    writeTabbed(os, "IYX_cm");
    writeTabbed(os, "IYY_cm");
    writeTabbed(os, "IYZ_cm");
    writeTabbed(os, "IZX_cm");
    writeTabbed(os, "IZY_cm");
    writeTabbed(os, "IZZ_cm");
    writeTabbed(os, "IXX_body");
    writeTabbed(os, "IXY_body");
    writeTabbed(os, "IXZ_body");
    writeTabbed(os, "IYX_body");
    writeTabbed(os, "IYY_body");
    writeTabbed(os, "IYZ_body");
    writeTabbed(os, "IZX_body");
    writeTabbed(os, "IZY_body");
    writeTabbed(os, "IZZ_body");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sloshingTank::sloshingTank
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

bool Foam::functionObjects::sloshingTank::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    fixed_ref_transformation_ = septernion(
        dict.getOrDefault<vector>("refBodyFixedAxisCoord", vector {}),
        quaternion(quaternion::XYZ, dict.getOrDefault<vector>("refBodyFixedAxisRot", vector {}) * Foam::constant::mathematical::pi / 180.));

    IOdictionary dynamicMeshDict(IOobject
    (
        "dynamicMeshDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ));

    autoPtr<solidBodyMotionFunction> solid_body_motion_function_(solidBodyMotionFunction::New(dynamicMeshDict.optionalSubDict("dynamicMeshDictCoeffs"), mesh_.time()));

    vector back_translation = solid_body_motion_function_->transformation().t();
    tensor back_rotation = solid_body_motion_function_->transformation().r().R();
    vector ref_translation = fixed_ref_transformation_.t();
    tensor ref_rotation = fixed_ref_transformation_.r().R();

    fixed_ref_cc_ = mesh_.C().primitiveField();

    forAll(fixed_ref_cc_, i) {
        fixed_ref_cc_[i] = ref_rotation & (((fixed_ref_cc_[i] & back_rotation) + back_translation) - ref_translation);
    }

    return true;
}

bool Foam::functionObjects::sloshingTank::execute()
{
    immiscibleIncompressibleTwoPhaseMixture& mix =
        mesh_.lookupObjectRef<immiscibleIncompressibleTwoPhaseMixture>("transportProperties");

    const scalarField& alpha1 = mix.alpha1().primitiveField();
    const scalarField& alpha2 = mix.alpha2().primitiveField();

    scalar rho1 = mix.rho1().value();
    scalar rho2 = mix.rho2().value();

    const vectorField& cc = mesh_.C().primitiveField();
    const scalarField& V = mesh_.V();

    com_ = gSum((rho1 * alpha1 + rho2 * alpha2) * V * cc) / gSum((rho1 * alpha1 + rho2 * alpha2) * V);
    com_fixed_ = gSum((rho1 * alpha1 + rho2 * alpha2) * V * fixed_ref_cc_) / gSum((rho1 * alpha1 + rho2 * alpha2) * V);

    vectorField cc_rel_com(cc - com_);

    tensor temp = gSum((rho1 * alpha1 + rho2 * alpha2) * V * cc_rel_com * cc_rel_com);
    mmi_cm_[tensor::XX] = temp[tensor::YY] + temp[tensor::ZZ];
    mmi_cm_[tensor::XY] = -temp[tensor::XY];
    mmi_cm_[tensor::XZ] = -temp[tensor::XZ];

    mmi_cm_[tensor::YX] = -temp[tensor::YX];
    mmi_cm_[tensor::YY] = temp[tensor::XX] + temp[tensor::ZZ];
    mmi_cm_[tensor::YZ] = -temp[tensor::YZ];

    mmi_cm_[tensor::ZX] = -temp[tensor::ZX];
    mmi_cm_[tensor::ZY] = -temp[tensor::ZY];
    mmi_cm_[tensor::ZZ] = temp[tensor::XX] + temp[tensor::YY];

    temp = gSum((rho1 * alpha1 + rho2 * alpha2) * V * fixed_ref_cc_ * fixed_ref_cc_);
    mmi_fixed_[tensor::XX] = temp[tensor::YY] + temp[tensor::ZZ];
    mmi_fixed_[tensor::XY] = -temp[tensor::XY];
    mmi_fixed_[tensor::XZ] = -temp[tensor::XZ];

    mmi_fixed_[tensor::YX] = -temp[tensor::YX];
    mmi_fixed_[tensor::YY] = temp[tensor::XX] + temp[tensor::ZZ];
    mmi_fixed_[tensor::YZ] = -temp[tensor::YZ];

    mmi_fixed_[tensor::ZX] = -temp[tensor::ZX];
    mmi_fixed_[tensor::ZY] = -temp[tensor::ZY];
    mmi_fixed_[tensor::ZZ] = temp[tensor::XX] + temp[tensor::YY];

    return true;
}

bool Foam::functionObjects::sloshingTank::write()
{
    Log
     << "Fluid Center of Mass" << endl
        << "Time: " << mesh_.time().timeName() << endl
        << "CoMx: " << com_[vector::X] << endl
        << "CoMy: " << com_[vector::Y] << endl
        << "CoMz: " << com_[vector::Z] << endl
        << "CoMx_fixed: " << com_[vector::X] << endl
        << "CoMy_fixed: " << com_[vector::Y] << endl
        << "CoMz_fixed: " << com_[vector::Z] << endl
        << "Ixx_cm: " << mmi_cm_[tensor::XX] << endl
        << "Ixy_cm: " << mmi_cm_[tensor::XY] << endl
        << "Ixz_cm: " << mmi_cm_[tensor::XZ] << endl
        << "Iyx_cm: " << mmi_cm_[tensor::YX] << endl
        << "Iyy_cm: " << mmi_cm_[tensor::YY] << endl
        << "Iyz_cm: " << mmi_cm_[tensor::YZ] << endl
        << "Izx_cm: " << mmi_cm_[tensor::ZX] << endl
        << "Izy_cm: " << mmi_cm_[tensor::ZY] << endl
        << "Izz_cm: " << mmi_cm_[tensor::ZZ] << endl
        << "Ixx_fixed: " << mmi_fixed_[tensor::XX] << endl
        << "Ixy_fixed: " << mmi_fixed_[tensor::XY] << endl
        << "Ixz_fixed: " << mmi_fixed_[tensor::XZ] << endl
        << "Iyx_fixed: " << mmi_fixed_[tensor::YX] << endl
        << "Iyy_fixed: " << mmi_fixed_[tensor::YY] << endl
        << "Iyz_fixed: " << mmi_fixed_[tensor::YZ] << endl
        << "Izx_fixed: " << mmi_fixed_[tensor::ZX] << endl
        << "Izy_fixed: " << mmi_fixed_[tensor::ZY] << endl
        << "Izz_fixed: " << mmi_fixed_[tensor::ZZ] << endl;

    if (Pstream::master())
    {

        file()
            << mesh_.time().timeName()
            << tab << com_[vector::X]
            << tab << com_[vector::Y]
            << tab << com_[vector::Z]
            << tab << com_fixed_[vector::X]
            << tab << com_fixed_[vector::Y]
            << tab << com_fixed_[vector::Z]
            << tab << mmi_cm_[tensor::XX]
            << tab << mmi_cm_[tensor::XY]
            << tab << mmi_cm_[tensor::XZ]
            << tab << mmi_cm_[tensor::YX]
            << tab << mmi_cm_[tensor::YY]
            << tab << mmi_cm_[tensor::YZ]
            << tab << mmi_cm_[tensor::ZX]
            << tab << mmi_cm_[tensor::ZY]
            << tab << mmi_cm_[tensor::ZZ]
            << tab << mmi_fixed_[tensor::XX]
            << tab << mmi_fixed_[tensor::XY]
            << tab << mmi_fixed_[tensor::XZ]
            << tab << mmi_fixed_[tensor::YX]
            << tab << mmi_fixed_[tensor::YY]
            << tab << mmi_fixed_[tensor::YZ]
            << tab << mmi_fixed_[tensor::ZX]
            << tab << mmi_fixed_[tensor::ZY]
            << tab << mmi_fixed_[tensor::ZZ]
            << endl;

    }

    return true;
}


// ************************************************************************* //
