/*---------------------------------------------------------------------------*\

    TAMS-AERO: Is part of the Turkish Aerospace Multiphysics Solver
               which is dedicated for the aerodynamic applications.
               
-------------------------------------------------------------------------------
License
    This file is part of TAMS-AERO.

Class
    lusgsTams

Description
    This is an optimized version of the LU-SGS preconditioning
    as per Luo, Baum, Lohner, JCP 146 664-690 (1998).    

Author
    Saleh Abuhanieh {Oct. 2024}
        Turkish Aerospace

\*---------------------------------------------------------------------------*/

#ifndef lusgsTams_H
#define lusgsTams_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "preconditioner.H"
#include "jacobianMatrix.H"
#include "fvjMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lusgsTams Declaration
\*---------------------------------------------------------------------------*/

template <int nScalar, int nVector>
class lusgsTams
:
    public preconditioner<nScalar, nVector>
{

    protected:

        //- Protected data
            autoPtr<scalarField> rDiagCoeff_; 
            autoPtr<scalarField> rDiagCoeffRemappedUp_;
            autoPtr<scalarField> rDiagCoeffRemappedLo_;                                             

    public:
        //- Runtime type information
        TypeName("LUSGS_Tams");

        lusgsTams
        (
            const dictionary& dict,
            const jacobianMatrix<nScalar,nVector>& jacobian,
            const preconditioner<nScalar,nVector>* prePreconditioner
        );

        //- Precondition a vector
        virtual void precondition
        (
            PtrList<volScalarField>& sVec,
            PtrList<volVectorField>& vVec
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "lusgsTamsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
