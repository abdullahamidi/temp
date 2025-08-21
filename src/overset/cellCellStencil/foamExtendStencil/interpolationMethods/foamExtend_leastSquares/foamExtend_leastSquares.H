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
    foamExtend_leastSquares.C

\*---------------------------------------------------------------------------*/

#ifndef foamExtend_leastSquares_H
#define foamExtend_leastSquares_H

#include "foamExtendStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                        Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class foamExtend_leastSquares
:
    public foamExtendStencil
{
    // Private Member Functions

        //- No copy construct
        foamExtend_leastSquares(const foamExtend_leastSquares&) = delete;

        //- No copy assignment
        void operator=(const foamExtend_leastSquares&) = delete;


public:

    //- Runtime type information
    TypeName("foamExtend_leastSquares");


    // Constructors

        //- Construct from fvMesh
        foamExtend_leastSquares(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~foamExtend_leastSquares();


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
