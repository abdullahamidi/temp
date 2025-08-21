/*---------------------------------------------------------------------------*\

    OSCFD Group

-------------------------------------------------------------------------------
Class
    gcro

Description
    Solves linear system J dU = Res as an alternative to GMRES solver under HiSA
    project where J is approximate LHS Jacobian.
    Reference articles listed below;
        - M. L. Parks, E. Sturler, D. D. Johnson, S. Maiti, "Recycling Krylov
        subspaces for sequences of linear systems", SIAM Journal on Scientific
        Computing, 2004
        - S. Xu, S.Timme, and K. J. Badcock, "Krylov subspace recycling for
        linearised aerodynamics analysis using DLR-TAU", University of
        Liverpool, 2015

Authors
    Barış Cumhur
        Turkish Aerospace Industries, Turkey

\*---------------------------------------------------------------------------*/

#ifndef gcro_H
#define gcro_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "QRMatrix.H"
#include "fvCFD.H"
#include "hisaSolver.H"
#include "jacobianMatrix.H"
#include "preconditioner.H"
#include "residualIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                  Class gcroSolver Declaration
\*---------------------------------------------------------------------------*/

template <int nScalar, int nVector>
class gcro : public hisaSolver<nScalar, nVector> {
protected:
public:
  //- Runtime type information
  TypeName("GCRO");

  gcro(const dictionary &dict, const jacobianMatrix<nScalar, nVector> &jacobian,
       const preconditioner<nScalar, nVector> *preconditioner,
       residualIO &
           defaultTol // Default residual tolerance and ordering of residual I/O
  );

  // Main solve functions:

  // Solve J dU = Res (J is LHS Jacobian):

  //- Returns number of iterations
  virtual label
  solve(PtrList<volScalarField> &sW,
        PtrList<volVectorField> &vW, // Initial/returned solution variables
        const PtrList<volScalarField> &sR,
        const PtrList<volVectorField> &vR, // Residuals
        autoPtr<residualIO> &pInitRes) const;

private:
  void update_solution(const List<PtrList<volScalarField>> &sVPtr,
                       const List<PtrList<volVectorField>> &vVPtr,
                       const scalarRectangularMatrix &H, const scalarField &c,
                       label nKrylov, const PtrList<volScalarField> &sR,
                       const PtrList<volVectorField> &vR,
                       PtrList<volScalarField> &dsW,
                       PtrList<volVectorField> &dvW,
                       PtrList<volScalarField> &sResidual,
                       PtrList<volVectorField> &vResidual) const;

  void update_solution(const List<PtrList<volScalarField>> &sVPtr,
                       const List<PtrList<volVectorField>> &vVPtr,
                       const List<PtrList<volScalarField>> &sUk,
                       const List<PtrList<volVectorField>> &vUk,
                       const List<PtrList<volScalarField>> &sCk,
                       const List<PtrList<volVectorField>> &vCk,
                       const scalarRectangularMatrix &H,
                       const scalarRectangularMatrix &B, label nKrylov,
                       label nRitz, const PtrList<volScalarField> &sR,
                       const PtrList<volVectorField> &vR,
                       PtrList<volScalarField> &dsW,
                       PtrList<volVectorField> &dvW,
                       PtrList<volScalarField> &sResidual,
                       PtrList<volVectorField> &vResidual) const;

  void update_solution(PtrList<volScalarField> &sW, PtrList<volVectorField> &vW,
                       const PtrList<volScalarField> &dsW,
                       const PtrList<volVectorField> &dvW) const;

  bool check_residual(residualIO &finalRes, const residualIO &initRes,
                      const residualIO &tol, const residualIO &tolRel,
                      const PtrList<volScalarField> &sResidual,
                      const PtrList<volVectorField> &vResidual,
                      const FixedList<scalar, nScalar> &sNormFactor,
                      const FixedList<scalar, nVector> &vNormFactor) const;

  static scalarRectangularMatrix Yk;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "gcroI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
