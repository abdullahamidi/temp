/*---------------------------------------------------------------------------*\

    TAMS-AERO: Is part of the Turkish Aerospace Multiphysics Solver
               which is dedicated for the aerodynamic applications.
               
-------------------------------------------------------------------------------
License
    This file is part of TAMS-AERO.

\*---------------------------------------------------------------------------*/

#include "lusgsTams.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "fvjOperators.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

template <int nScalar, int nVector>
lusgsTams<nScalar, nVector>::lusgsTams
(
    const dictionary& dict,
    const jacobianMatrix<nScalar, nVector>& jacobian,
    const preconditioner<nScalar,nVector>* prePreconditioner
)
:
    preconditioner<nScalar, nVector>
    (
        typeName,
        dict,
        jacobian,
        prePreconditioner
    )
{
    // Generate a scalar diagonal coefficient based on the max of the diagonal
    // of the Jacobian

    rDiagCoeff_.reset(new scalarField(this->jacobian_.mesh().nCells(), GREAT));

    forN(jacobian.mesh().nCells(), celli)
    {
        forN(nScalar, i)
        {
            if (this->jacobian_.dSBySExists(i,i))
            {
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(this->jacobian_.dSByS(i,i).diag()[celli])
                    );
            }
            else
            {
                FatalErrorInFunction
                    << "Diagonal S" << i << " of Jacobian not populated."
                    << exit(FatalError);
            }
        }
        forN(nVector, i)
        {
            if (this->jacobian_.dVByVExists(i,i))
            {
                const tensor& diag = this->jacobian_.dVByV(i,i).diag()[celli];
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(diag.xx())
                    );
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(diag.yy())
                    );
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(diag.zz())
                    );
            }
            else
            {
                FatalErrorInFunction
                    << "Diagonal V" << i << " of Jacobian not populated."
                    << exit(FatalError);
            }
        }
        if (rDiagCoeff_()[celli] < VSMALL)
        {
            FatalErrorInFunction << "All diagonals of Jacobian are zero." << endl 
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template <int nScalar, int nVector>
void lusgsTams<nScalar, nVector>::precondition
(
    PtrList<volScalarField>& sVec,
    PtrList<volVectorField>& vVec
) const
{

    // Call base class to apply any pre-preconditioner
    preconditioner<nScalar,nVector>::precondition(sVec, vVec);
    
    const fvMesh& msh(this->mesh_);
    scalarField V(msh.V());

    // Residual is still in strong form
    forN(nScalar, i) sVec[i].primitiveFieldRef() *= V;
    forN(nVector, i) vVec[i].primitiveFieldRef() *= V;
        
    // Step I: divide by the diagonal
    forN(nScalar, i) sVec[i].primitiveFieldRef() *= rDiagCoeff_();
    forN(nVector, i) vVec[i].primitiveFieldRef() *= rDiagCoeff_();
    
    
    label nIntFaces(msh.nInternalFaces());
    
    // Step II: forwardSweep
        forN(nScalar, i)
        {
            forN(nScalar, j)
            {
                const fvjMatrix<scalar>& matrix(this->jacobian_.dSByS(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());
                if (this->jacobian_.dSBySExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ losortPtr =
                        matrix.lduAddr().losortAddr().begin();
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const scalar* const __restrict__ lowerPtr =
                        matrix.lower().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();  
                    volScalarField& sVec_i = sVec[i];
                    volScalarField& sVec_j = sVec[j];                                               
                    for (label face=0; face<nIntFaces; face++)
                    {
                        const label sface = losortPtr[face];
                        sVec_i[uPtr[sface]] -=
                            rDiagCoeffPtr[uPtr[sface]]*dot(lowerPtr[sface], sVec_j[lPtr[sface]]);  
                    }
                }
            }
        }

        forN(nScalar, i)
        {
            forN(nVector, j)
            {
                const fvjMatrix<vector>& matrix(this->jacobian_.dSByV(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());
                if (this->jacobian_.dSByVExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ losortPtr =
                        matrix.lduAddr().losortAddr().begin();
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const vector* const __restrict__ lowerPtr =
                        matrix.lower().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();
                    volScalarField& sVec_i = sVec[i];
                    volVectorField& vVec_j = vVec[j];                                                                          
                    for (label face=0; face<nIntFaces; face++)
                    {
                        const label sface = losortPtr[face];
                        sVec_i[uPtr[sface]] -=
                            rDiagCoeffPtr[uPtr[sface]]*dot(lowerPtr[sface], vVec_j[lPtr[sface]]);   
                    }                                    
                }
            }
        }
        
        forN(nVector, i)
        {
            forN(nScalar, j)
            {
                const fvjMatrix<vector>& matrix(this->jacobian_.dVByS(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());            
                if (this->jacobian_.dVBySExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ losortPtr =
                        matrix.lduAddr().losortAddr().begin();
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const vector* const __restrict__ lowerPtr =
                        matrix.lower().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();
                    volVectorField& vVec_i = vVec[i];
                    volScalarField& sVec_j = sVec[j];                                                                           
                    for (label face=0; face<nIntFaces; face++)
                    {
                        const label sface = losortPtr[face];
                        vVec_i[uPtr[sface]] -=
                            rDiagCoeffPtr[uPtr[sface]]*dot(lowerPtr[sface], sVec_j[lPtr[sface]]);   
                    }                                    
                } 
            }
        }

        forN(nVector, i)
        {
            forN(nVector, j)
            {
                const fvjMatrix<tensor>& matrix(this->jacobian_.dVByV(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());            
                if (this->jacobian_.dVByVExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ losortPtr =
                        matrix.lduAddr().losortAddr().begin();
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const tensor* const __restrict__ lowerPtr =
                        matrix.lower().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();
                    volVectorField& vVec_i = vVec[i];
                    volVectorField& vVec_j = vVec[j];                                                                          
                    for (label face=0; face<nIntFaces; face++)
                    {
                        const label sface = losortPtr[face];
                        vVec_i[uPtr[sface]] -=
                            rDiagCoeffPtr[uPtr[sface]]*dot(lowerPtr[sface], vVec_j[lPtr[sface]]);   
                    }                                     
                }
            }
        }

    // Step III: reverseSweep
        forN(nScalar, i)
        {
            forN(nScalar, j)
            {
                const fvjMatrix<scalar>& matrix(this->jacobian_.dSByS(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());
                if (this->jacobian_.dSBySExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const scalar* const __restrict__ upperPtr =
                        matrix.upper().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();
                    volScalarField& sVec_i = sVec[i];
                    volScalarField& sVec_j = sVec[j];                                                                          
                    for (label face=nIntFaces - 1; face>=0; face--)
                    {
                        sVec_i[lPtr[face]] -=
                            rDiagCoeffPtr[lPtr[face]]*dot(upperPtr[face], sVec_j[uPtr[face]]);      
                    }
                }
            }
        }

        forN(nScalar, i)
        {
            forN(nVector, j)
            {
                const fvjMatrix<vector>& matrix(this->jacobian_.dSByV(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());
                if (this->jacobian_.dSByVExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const vector* const __restrict__ upperPtr =
                        matrix.upper().begin(); 
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin(); 
                    volScalarField& sVec_i = sVec[i];
                    volVectorField& vVec_j = vVec[j];                                                                        
                    for (label face=nIntFaces - 1; face>=0; face--)
                    {
                        sVec_i[lPtr[face]] -=
                            rDiagCoeffPtr[lPtr[face]]*dot(upperPtr[face], vVec_j[uPtr[face]]);      
                    }                                
                }
            }
        }

        forN(nVector, i)
        {
            forN(nScalar, j)
            {
                const fvjMatrix<vector>& matrix(this->jacobian_.dVByS(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());            
                if (this->jacobian_.dVBySExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const vector* const __restrict__ upperPtr =
                        matrix.upper().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();
                    volVectorField& vVec_i = vVec[i];
                    volScalarField& sVec_j = sVec[j];                                                                          
                    for (label face=nIntFaces - 1; face>=0; face--)
                    {
                        vVec_i[lPtr[face]] -=
                            rDiagCoeffPtr[lPtr[face]]*dot(upperPtr[face], sVec_j[uPtr[face]]);      
                    }                                  
                }
            }
        }

        forN(nVector, i)
        {
            forN(nVector, j)
            {
                const fvjMatrix<tensor>& matrix(this->jacobian_.dVByV(i,j));
                bool hasLorU = (matrix.hasLower() || matrix.hasUpper());            
                if (this->jacobian_.dVByVExists(i,j) && hasLorU)
                {
                    const label* const __restrict__ uPtr =
                        matrix.lduAddr().upperAddr().begin();
                    const label* const __restrict__ lPtr =
                        matrix.lduAddr().lowerAddr().begin();   
                    const tensor* const __restrict__ upperPtr =
                        matrix.upper().begin();
                    const scalar* const __restrict__ rDiagCoeffPtr =
                        rDiagCoeff_().begin();
                    volVectorField& vVec_i = vVec[i];
                    volVectorField& vVec_j = vVec[j];                                                                         
                    for (label face=nIntFaces - 1; face>=0; face--)
                    {
                        vVec_i[lPtr[face]] -=
                            rDiagCoeffPtr[lPtr[face]]*dot(upperPtr[face], vVec_j[uPtr[face]]);      
                    }                                   
                }
            }
        }                         
                                                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
