/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

Class
    Foam::movingWallMRFVelocityFvPatchVectorField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a velocity condition for cases with
    moving walls.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            movingWallMRFVelocity;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchVectorField
    Foam::MRFZone

SourceFiles
    movingWallMRFVelocityFvPatchVectorField.C

This boundary condition inforces the the realtive velocity = 0 condition
on the boundary. To be used with SRF/MRF simulations.
S.A {Sep. 2022}

\*---------------------------------------------------------------------------*/

#ifndef movingWallMRFVelocityFvPatchVectorField_H
#define movingWallMRFVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingWallMRFVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingWallMRFVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("movingWallMRFVelocity");


    // Constructors

        //- Construct from patch and internal field
        movingWallMRFVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallMRFVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallMRFVelocityFvPatchVectorField
        //  onto a new patch
        movingWallMRFVelocityFvPatchVectorField
        (
            const movingWallMRFVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingWallMRFVelocityFvPatchVectorField
        (
            const movingWallMRFVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingWallMRFVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingWallMRFVelocityFvPatchVectorField
        (
            const movingWallMRFVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new movingWallMRFVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
