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

\*---------------------------------------------------------------------------*/

#include "movingWallMRFVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallMRFVelocityFvPatchVectorField::
movingWallMRFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::movingWallMRFVelocityFvPatchVectorField::
movingWallMRFVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::movingWallMRFVelocityFvPatchVectorField::
movingWallMRFVelocityFvPatchVectorField
(
    const movingWallMRFVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::movingWallMRFVelocityFvPatchVectorField::
movingWallMRFVelocityFvPatchVectorField
(
    const movingWallMRFVelocityFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::movingWallMRFVelocityFvPatchVectorField::
movingWallMRFVelocityFvPatchVectorField
(
    const movingWallMRFVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallMRFVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    //// S.A {Sep. 2022}
    bool MRF_NS = false;
    if(mesh.foundObject<surfaceVectorField>("omegaCrossR"))
    {
        MRF_NS = true;
    }
    ////

    if (mesh.moving())
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        const pointField& oldPoints = mesh.oldPoints();

        vectorField oldFc(pp.size());

        forAll(oldFc, i)
        {
            oldFc[i] = pp[i].centre(oldPoints);
        }

        const scalar deltaT = mesh.time().deltaTValue();

        const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());

        scalarField phip
        (
            p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
        );

        const vectorField n(p.nf());
        const scalarField& magSf = p.magSf();
        tmp<scalarField> Un = phip/(magSf + VSMALL);


        vectorField::operator=(Up + n*(Un - (n & Up)));
    }
    else if(MRF_NS)  //// S.A {Sep. 2022}
    {
        const fvPatch& p = patch();
        const vectorField n(p.nf());
        const surfaceVectorField& omegaCrossRns = mesh.lookupObject<surfaceVectorField>("omegaCrossR");

        const vectorField Up(omegaCrossRns.boundaryField()[this->patch().index()]);

        const surfaceScalarField& minusMRFFlux = mesh.lookupObject<surfaceScalarField>("minusMRFFlux");
        tmp<scalarField> tNormalVelocity = minusMRFFlux.boundaryField()[this->patch().index()]/this->patch().magSf();
        const scalarField Un = -tNormalVelocity();
        vectorField::operator=(Up + n*(Un - (n & Up)));
    }            ////

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingWallMRFVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingWallMRFVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
