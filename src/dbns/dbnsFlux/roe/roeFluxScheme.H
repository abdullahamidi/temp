/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

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

Class
    Foam::roeFluxScheme

Description

    Roe-Pike flux splitting scheme according to:
    
    [1] P. Roe, J. Pike, Efficient Construction and Utilisation of Approximate
        Riemann Solutions, 1985.
    
    For the implementation:
    
    [2] E. F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics,
        3rd Edition, 2008. (section 11.3)
        
    [3] J. Blazek, Computational Fluid Dynamics Principles and Applications,
        3rd Edition, 2015. (section 4.3.3)
        
    [4] S. Chun et al., Chinese Journal of Aeronautics,
        Analysis on capabilities of density-based solvers within
        OpenFOAM to distinguish aerothermal variables in diffusion
        boundary layer. 26(6): 1370-1379.
        
    [5] foam-extend rhoFlux implementation.   

SourceFiles
    roeFluxScheme.C

Author
    Saleh Abuhanieh

\*---------------------------------------------------------------------------*/

#ifndef roeFluxScheme_H
#define roeFluxScheme_H

#include "fluxScheme.H"
#include "surfaceMesh.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class roeFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class roeFluxScheme
:
    public fluxScheme
{

protected:
    // Protected data
        const fvMesh& mesh_;
        const fluidThermo& thermo_;
        const volScalarField& rho_;
        const volVectorField& U_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;

        surfaceScalarField pos_;
        surfaceScalarField neg_;

public:
    //- Runtime type information
    TypeName("Roe");


	// Constructors

        //- Construct from components
        roeFluxScheme
        (
            const dictionary& dict,
            const fluidThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& rhoU,
            const volScalarField& rhoE
        );


    // Destructor

        ~roeFluxScheme();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
