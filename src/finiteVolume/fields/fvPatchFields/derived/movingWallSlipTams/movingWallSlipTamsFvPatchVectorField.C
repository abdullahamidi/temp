/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014 Johan Heyns - CSIR, South Africa
    Copyright (C) 2011 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "movingWallSlipTamsFvPatchVectorField.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingWallSlipTamsFvPatchVectorField::movingWallSlipTamsFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchField<vector>(p, iF)
{}


Foam::movingWallSlipTamsFvPatchVectorField::movingWallSlipTamsFvPatchVectorField
(
    const movingWallSlipTamsFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::movingWallSlipTamsFvPatchVectorField::movingWallSlipTamsFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<vector>(p, iF)
{
    //// S.A {Sep. 2022}
    if(this->internalField().mesh().foundObject<surfaceScalarField>("minusMRFFlux"))
    {
        //Info<<"MRF simulation is detected at movingWallSlipTams"<<endl;
        MRF = true;
    }
    ////

    evaluate();
}


Foam::movingWallSlipTamsFvPatchVectorField::movingWallSlipTamsFvPatchVectorField
(
    const movingWallSlipTamsFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingWallSlipTamsFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchField<vector>::autoMap(m);
}


void Foam::movingWallSlipTamsFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    transformFvPatchField<vector>::rmap(ptf, addr);
}


Foam::tmp<Foam::vectorField >
Foam::movingWallSlipTamsFvPatchVectorField::snGrad() const
{
    const vectorField nHat(this->patch().nf());
    const vectorField pif(this->patchInternalField());

    const fvMesh& mesh = this->internalField().mesh();

    if (mesh.changing())
    {
        tmp<scalarField> tNormalVelocity = mesh.phi().boundaryField()[this->patch().index()]/this->patch().magSf();
        const scalarField& normalVelocity = tNormalVelocity();

        return
        (
            (nHat*normalVelocity + transform(I - sqr(nHat), pif)) - pif
        )*this->patch().deltaCoeffs();
    }
    else if (MRF)  //// S.A {Sep. 2022}
    {
        const surfaceScalarField& minusMRFFlux = mesh.lookupObject<surfaceScalarField>("minusMRFFlux");
        tmp<scalarField> tNormalVelocity = minusMRFFlux.boundaryField()[this->patch().index()]/this->patch().magSf();
        const scalarField& normalVelocity = tNormalVelocity();

        return
        (
            (-nHat*normalVelocity + transform(I - sqr(nHat), pif)) - pif
        )*this->patch().deltaCoeffs();
    }             ////
    else
    {
        return
        (
            (transform(I - sqr(nHat), pif)) - pif
        )*this->patch().deltaCoeffs();
    }
}


void Foam::movingWallSlipTamsFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = this->internalField().mesh();

    const vectorField nHat(this->patch().nf());

    if (mesh.changing())
    {
        tmp<scalarField> tNormalVelocity = mesh.phi().boundaryField()[this->patch().index()]/this->patch().magSf();
        const scalarField& normalVelocity = tNormalVelocity;

        vectorField::operator=
        (
            nHat*normalVelocity
          + transform(I - sqr(nHat), this->patchInternalField())
        );
    }
    else if (MRF)  //// S.A {Sep. 2022}
    {
        const surfaceScalarField& minusMRFFlux = mesh.lookupObject<surfaceScalarField>("minusMRFFlux");
        tmp<scalarField> tNormalVelocity = minusMRFFlux.boundaryField()[this->patch().index()]/this->patch().magSf();
        const scalarField& normalVelocity = tNormalVelocity;

        vectorField::operator=
        (
            -nHat*normalVelocity
          + transform(I - sqr(nHat), this->patchInternalField())
        );
    }             ////
    else
    {
        vectorField::operator=
        (
            transform(I - sqr(nHat), this->patchInternalField())
        );

    }

    transformFvPatchField<vector>::evaluate();
}


Foam::tmp<Foam::vectorField >
Foam::movingWallSlipTamsFvPatchVectorField::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}


void Foam::movingWallSlipTamsFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchField<vector>::write(os);
#if FOUNDATION >= 7
    writeEntry(os, "value", *this);
#else
    writeEntry("value", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingWallSlipTamsFvPatchVectorField
    );
}

// ************************************************************************* //
