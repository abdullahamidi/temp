/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "foamExtend_inverseDistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{
    defineTypeNameAndDebug(foamExtend_inverseDistance, 0);
    addToRunTimeSelectionTable(cellCellStencil, foamExtend_inverseDistance, mesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellCellStencils::foamExtend_inverseDistance::stencilWeights
(
    const point& sample,
    const pointList& donorCcs,
    scalarList& weights
) const
{

    // Inverse-distance weighting

    weights.setSize(donorCcs.size());
    scalar sum = 0.0;
    forAll(donorCcs, i)
    {
        scalar d = mag(sample-donorCcs[i]);

        if (d > ROOTVSMALL)
        {
            weights[i] = 1.0/d;
            sum += weights[i];
        }
        else
        {
            // Short circuit
            weights = 0.0;
            weights[i] = 1.0;
            return;
        }
    }
    forAll(weights, i)
    {
        weights[i] /= sum;
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCellStencils::foamExtend_inverseDistance::foamExtend_inverseDistance
(
    const fvMesh& mesh,
    const dictionary& dict,
    const bool doUpdate
)
:
    foamExtendStencil(mesh, dict, false)
{
    if (doUpdate)
    {
        update();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellCellStencils::foamExtend_inverseDistance::~foamExtend_inverseDistance()
{}


// ************************************************************************* //
