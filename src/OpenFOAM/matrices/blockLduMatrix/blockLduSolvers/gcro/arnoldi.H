/*---------------------------------------------------------------------------*\

    OSCFD Group

-------------------------------------------------------------------------------
Class
    gcro

Description
    Solves linear system J dU = Res where J is approximate LHS Jacobian

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

\*---------------------------------------------------------------------------*/

#ifndef Arnoldi_H
#define Arnoldi_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "jacobianMatrix.H"
#include "preconditioner.H"
#include "algebra.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

template<typename Htype>
inline void gram_schmidt_process
(
    PtrList<volScalarField>& s, PtrList<volVectorField>& v,
    const List<PtrList<volScalarField>>& sPtr, const List<PtrList<volVectorField>>& vPtr,
    label nVecs,
    Htype&& H,
    label column
) {

    labelList reduceRequest(nVecs, -1);

    for (label i = 0; i < nVecs; i++) {
        forAll(s, n)                         // h_ij = <v_i, v>  
        {
            H(i, column) += sumProd(s[n].primitiveField(), sPtr[i][n].primitiveField());
        }
        forAll(v, n) {
            H(i, column) += sumProd(v[n].primitiveField(), vPtr[i][n].primitiveField());
        }
        reduce(H(i, column), sumOp<scalar>(), Pstream::msgType(), UPstream::worldComm, reduceRequest[i]);
    }
    for (label i = 0; i < nVecs; i++) {
        if (Pstream::parRun() && reduceRequest[i] != -1) {
            Pstream::waitRequest(reduceRequest[i]);
            reduceRequest[i] = -1;
        }
        forAll(s, n)                          // v = v - h_ij v_i
        {
            s[n].primitiveFieldRef() -= H(i, column) * sPtr[i][n].primitiveField();
        }
        forAll(v, n) {
            v[n].primitiveFieldRef() -= H(i, column) * vPtr[i][n].primitiveField();
        }
    }
}

template <int nScalar, int nVector, typename HMatrixT>
inline void arnoldi_iteration
(
    const jacobianMatrix<nScalar, nVector>& J,
    const preconditioner<nScalar, nVector>* preconditioner,
    const PtrList<volScalarField>& sResidual, const PtrList<volVectorField>& vResidual, // residuals and temporary storage
    label nKrylov, // nKrylov
    label nRitz, // nRitz
    List<PtrList<volScalarField>>& sVPtr, List<PtrList<volVectorField>>& vVPtr, // Krylov vectors
    HMatrixT&& H,
    scalarField& c
) {
    forAll(sVPtr[nRitz], n) {
        sVPtr[nRitz][n].primitiveFieldRef() = sResidual[n].primitiveField();
    }
    forAll(vVPtr[nRitz], n) {
        vVPtr[nRitz][n].primitiveFieldRef() = vResidual[n].primitiveField();
    }

    // Execute preconditioning
    preconditioner->precondition(sVPtr[nRitz], vVPtr[nRitz]);
    c[nRitz] = normalize(sVPtr[nRitz], vVPtr[nRitz]);

    for (label j = nRitz; j < nKrylov; j++)	// Search directions
    {
        // Set search direction - delay scaling vector to allow parallel comms overlap
        PtrList<volScalarField>& sV = sVPtr[j + 1];
        PtrList<volVectorField>& vV = vVPtr[j + 1];

        // Matrix vector product
        J.matrixMul(sVPtr[j], vVPtr[j], sV, vV);
        preconditioner->precondition(sV, vV);       // v_j+1 = P^-1 A v_j

        // Gram-Schmidt step - Orthogonalization on previous Krylov vectors
        gram_schmidt_process(sV, vV, sVPtr, vVPtr, j + 1, H, j - nRitz);

        H(j + 1, j - nRitz) = normalize(sV, vV);
    }
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif      

// ************************************************************************* //

