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
    

Description

SourceFiles
    foamExtend_averageValue.C

\*---------------------------------------------------------------------------*/

#ifndef foamExtend_averageValue_H
#define foamExtend_averageValue_H

#include "foamExtendStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                        Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class foamExtend_averageValue
:
    public foamExtendStencil
{
    // Private Member Functions

        //- No copy construct
        foamExtend_averageValue(const foamExtend_averageValue&) = delete;

        //- No copy assignment
        void operator=(const foamExtend_averageValue&) = delete;


public:

    //- Runtime type information
    TypeName("foamExtend_averageValue");


    // Constructors

        //- Construct from fvMesh
        foamExtend_averageValue(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~foamExtend_averageValue();


    // Member Functions

        //- Calculate lsq weights for single acceptor
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
