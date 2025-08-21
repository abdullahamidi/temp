/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

Class
    Foam::cellCellStencils::tiogaStencil

Description
    A wrapper class for using the TIOGA overset library as a cellCellStencil
    within the OpenFOAM (ESI) framework. Tested with OpenFOAM v2206
    
Author
    Saleh Abuhanieh {Sep. 2023}    

SourceFiles
    tiogaStencil.C

\*---------------------------------------------------------------------------*/

#ifndef tiogaStencil_H
#define tiogaStencil_H

#include "cellCellStencil.H"
#include "volFields.H"
#include "labelVector.H"
#include "pointList.H"
#include "tioga.h"
#include "tiogaMeshAtrributes.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshSubset;

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                     Class tiogaStencil Declaration
\*---------------------------------------------------------------------------*/

class tiogaStencil
:
    public cellCellStencil
{
protected:

    // Protected data

        //- Dictionary of motion control parameters
        const dictionary dict_;

        //- Per cell the cell type
        labelList cellTypes_;

        //- Indices of interpolated cells
        labelList interpolationCells_;

        //- Fetch interpolated cells
        autoPtr<mapDistribute> cellInterpolationMap_;

        //- Interpolation stencil
        labelListList cellStencil_;

        //- Interpolation weights
        scalarListList cellInterpolationWeights_;

        //- Amount of interpolation
        volScalarField cellInterpolationWeight_;

        //- tioga
        TIOGA::tioga* tg_; 
        //

        //- Restart condition flag
        bool restartCondition_;


private:

    // Private Member Functions

        //- No copy construct
        tiogaStencil(const tiogaStencil&) = delete;

        //- No copy assignment
        void operator=(const tiogaStencil&) = delete;

public:

    //- Runtime type information
    TypeName("tiogaStencil");


    // Constructors

        //- Construct from fvMesh
        tiogaStencil(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~tiogaStencil();


    // Member Functions

        //- Update stencils. Return false if nothing changed.
        virtual bool update();

        //- Return the cell type list
        virtual const labelUList& cellTypes() const
        {
            return cellTypes_;
        }

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const
        {
            return interpolationCells_;
        }

        //- Return a communication schedule
        virtual const mapDistribute& cellInterpolationMap() const
        {
            if (!cellInterpolationMap_.valid())
            {
                const_cast<tiogaStencil&>(*this).update();
            }
            return cellInterpolationMap_();
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            return cellStencil_;
        }

        //- Weights for cellStencil
        virtual const scalarListList& cellInterpolationWeights() const
        {
            return cellInterpolationWeights_;
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            return cellInterpolationWeight_;
        }

        //- Calculate inverse distance weights for a single acceptor
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
