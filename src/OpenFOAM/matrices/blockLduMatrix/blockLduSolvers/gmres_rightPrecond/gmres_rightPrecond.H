/*---------------------------------------------------------------------------*\

    TAMS-AERO: Is part of the Turkish Aerospace Multiphysics Solver
               which is dedicated for the aerodynamic applications.

-------------------------------------------------------------------------------
License
    This file is part of TAMS-AERO.

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
    gmres_rightPrecond

Description
    Righ-preconditioned GMRES solver as described by Y. Saad. 2003
    ("Iterative methods for sparse linear systems"), secion 9.3.2.

Author
    Saleh Abuhanieh {July 2024}
        Turkish Aerospace

\*---------------------------------------------------------------------------*/

#ifndef gmres_rightPrecond_H
#define gmres_rightPrecond_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "jacobianMatrix.H"
#include "residualIO.H"
#include "preconditioner.H"
#include "hisaSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class gmres_rightPrecondSolver Declaration
\*---------------------------------------------------------------------------*/

template <int nScalar, int nVector>
class gmres_rightPrecond
:
    public hisaSolver<nScalar,nVector>
{
    protected:

        inline void givensRotation(const scalar& h, const scalar& beta, scalar& c, scalar& s) const;

    public:
        //- Runtime type information
        TypeName("GMRES_rightPrecond");

        gmres_rightPrecond
        (
            const dictionary& dict,
            const jacobianMatrix<nScalar,nVector>& jacobian,
            const preconditioner<nScalar,nVector>* preconditioner,
            residualIO& defaultTol // Default residual tolerance and ordering of residual I/O
        );

        // Main solve functions:

        // Solve J dU = Res (J is LHS Jacobian):

        //- Returns number of iterations
        virtual label solve
        (
            PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW,  //Initial/returned solution variables
            const PtrList<volScalarField>& sR, const PtrList<volVectorField>& vR,  //Residuals
            autoPtr< residualIO >& pInitRes
        ) const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gmres_rightPrecondTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

