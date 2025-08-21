/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa

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

\*---------------------------------------------------------------------------*/

#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "algebra.H"
#include "arnoldi.H"
#include "emptyFvPatch.H"
#include "fvjOperators.H"
#include "gcro.H"
#include "scalarMatrices.H"
#include "transformFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

template <int nScalar, int nVector>
gcro<nScalar, nVector>::gcro(
    const dictionary &dict, const jacobianMatrix<nScalar, nVector> &jacobian,
    const preconditioner<nScalar, nVector> *preconditioner,
    residualIO
        &defaultTol // Default residual tolerance and ordering of residual I/O
    )
    : hisaSolver<nScalar, nVector>(typeName, dict, jacobian, preconditioner,
                                   defaultTol) {}

template <int nScalar, int nVector>
label gcro<nScalar, nVector>::solve(
    PtrList<volScalarField> &sW,
    PtrList<volVectorField> &vW, // Initial/returned solution variables
    const PtrList<volScalarField> &sR,
    const PtrList<volVectorField> &vR, // Residuals
    autoPtr<residualIO> &pInitRes) const {
  // Read settings from dictionary
  const dictionary &dict = this->dict_;
  // Number of Krylov-space vectors
  label nKrylov = dict.lookupOrDefault<label>("nKrylov", 8);
  label nRitz = dict.lookupOrDefault<label>("nRitz", 4);
  bool EnableRecycling = dict.lookupOrDefault<bool>("EnableRecycling", false);
  // Number of GCRO iterations (restarts)
  label maxIter = dict.lookupOrDefault<label>("maxIter", 20);
  // Solver absolute tolerance
  residualIO tol(this->defaultTol_);
  // Solver relative tolerance
  residualIO tolRel(this->defaultTol_, 0.01);
  if (dict.found("solverTol")) {
    dict.lookup("solverTol") >> tol;
  }
  if (dict.found("solverTolRel")) {
    dict.lookup("solverTolRel") >> tolRel;
  }
  // Make sure the relative tolerance is ignored in non-solved directions
  vector::labelType validComponents =
      this->mesh_.solutionD(); //-1 for empty directions
  forAll(validComponents, cmpt) {
    if (validComponents[cmpt] == -1) {
      forN(nVector, vectorI) {
        vector v = tolRel.getVector(vectorI);
        v[cmpt] = GREAT;
        tolRel.setVector(vectorI, v);
      }
    }
  }

  Info << "Solving for (";
  const labelList &residualOrdering = tol.residualOrdering();
  for (label i = 0; i < residualOrdering.size(); i++) {
    if (i) {
      Info << " ";
    }
    forAll(sW, j) {
      if (residualOrdering[j] == i) {
        Info << sW[j].name();
      }
    }
    forAll(vW, j) {
      if (residualOrdering[nScalar + j] == i) {
        Info << vW[j].name();
      }
    }
  }
  Info << ")" << endl;

  // Allocate temp storage
  PtrList<volScalarField> dsW(nScalar);
  PtrList<volVectorField> dvW(nVector);
  // Allocate variables to hold solved increment
  PtrList<volScalarField> sResidual(nScalar);
  PtrList<volVectorField> vResidual(nVector);
  forN(nScalar, k) {
    dsW.set(k, new volScalarField(
                   IOobject("d" + sW[k].name(), this->mesh_.time().timeName(),
                            this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
                   sW[k], zeroGradientFvPatchScalarField::typeName));
    sResidual.set(
        k, new volScalarField(
               IOobject("tempVector", this->mesh_.time().timeName(),
                        this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
               this->mesh_, dimless, zeroGradientFvPatchScalarField::typeName));
  }
  forN(nVector, k) {
    dvW.set(k, new volVectorField(
                   IOobject("d" + vW[k].name(), this->mesh_.time().timeName(),
                            this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
                   vW[k], zeroGradientFvPatchVectorField::typeName));
    vResidual.set(
        k, new volVectorField(
               IOobject("tempVector", this->mesh_.time().timeName(),
                        this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
               this->mesh_, dimless, zeroGradientFvPatchVectorField::typeName));
  }

  // Allocate residual change storage
  List<PtrList<volScalarField>> sCk(nRitz);
  List<PtrList<volVectorField>> vCk(nRitz);
  List<PtrList<volScalarField>> sUk(nRitz);
  List<PtrList<volVectorField>> vUk(nRitz);
  for (label i = 0; i < nRitz; i++) {
    sCk[i].setSize(nScalar);
    sUk[i].setSize(nScalar);
    forN(nScalar, k) {
      sCk[i].set(
          k,
          new volScalarField(
              IOobject("ResidualOptVector", this->mesh_.time().timeName(),
                       this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
              this->mesh_, dimless, zeroGradientFvPatchScalarField::typeName));
      sUk[i].set(
          k,
          new volScalarField(
              IOobject("ResidualOptVector", this->mesh_.time().timeName(),
                       this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
              this->mesh_, dimless, zeroGradientFvPatchScalarField::typeName));
    }
    vCk[i].setSize(nVector);
    vUk[i].setSize(nVector);
    forN(nVector, k) {
      vCk[i].set(
          k,
          new volVectorField(
              IOobject("ResidualOptVector", this->mesh_.time().timeName(),
                       this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
              this->mesh_, dimless, zeroGradientFvPatchVectorField::typeName));
      vUk[i].set(
          k,
          new volVectorField(
              IOobject("ResidualOptVector", this->mesh_.time().timeName(),
                       this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
              this->mesh_, dimless, zeroGradientFvPatchVectorField::typeName));
    }
  }

  // Allocate Krylov vectors
  List<PtrList<volScalarField>> sVPtr(nKrylov + 1);
  List<PtrList<volVectorField>> vVPtr(nKrylov + 1);
  for (label i = 0; i < nKrylov + 1; i++) {
    sVPtr[i].setSize(nScalar);
    forN(nScalar, k) {
      sVPtr[i].set(
          k,
          new volScalarField(
              IOobject("krylovVector", this->mesh_.time().timeName(),
                       this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
              this->mesh_, dimless, zeroGradientFvPatchScalarField::typeName));
    }
    vVPtr[i].setSize(nVector);
    forN(nVector, k) {
      vVPtr[i].set(
          k,
          new volVectorField(
              IOobject("krylovVector", this->mesh_.time().timeName(),
                       this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
              this->mesh_, dimless, zeroGradientFvPatchVectorField::typeName));
    }
  }

  // Norm factor: more invariant representation of the residual
  // OF and Hrv suggest the following
  // Type xRef = gAverage(x);
  // pA = A xRef;
  // scalar NormFactor = gSum(cmptMag(A x - pA) + cmptMag(source - pA)) +
  // matrix.small; NOTE: pA ensures initial residual normalises to 1 when
  // initialised with const fields

  // Calc A(W-WAve)

  forN(nScalar, k) {
    scalar avg = gAverage(sW[k].primitiveField());
    dsW[k].primitiveFieldRef() -= avg;
    dsW[k].boundaryFieldRef() -= avg;
  }
  forN(nVector, k) {
    vector avg = gAverage(vW[k].primitiveField());
    dvW[k].primitiveFieldRef() -= avg;
    dvW[k].boundaryFieldRef() -= avg;
  }

  // Matrix multiplication
  this->jacobian_.matrixMul(dsW, dvW, sResidual, vResidual);

  forN(nScalar, k) {
    dsW[k].primitiveFieldRef() = 0.0;
    dsW[k].boundaryFieldRef() = 0.0;
  }
  forN(nVector, k) {
    dvW[k].primitiveFieldRef() = vector::zero;
    dvW[k].boundaryFieldRef() = vector::zero;
  }

  pInitRes.set(new residualIO(tol));
  residualIO &initRes(pInitRes());
  // Use vector magnitude for normalisation
  FixedList<scalar, nScalar> sNormFactor;
  forAll(sNormFactor, i) {
    sNormFactor[i] =
        gSum(mag(sResidual[i].primitiveField()) + mag(sR[i].primitiveField())) +
        VSMALL;
    initRes.setScalar(i, gSumMag(sR[i].primitiveField()) / sNormFactor[i]);
  }
  FixedList<scalar, nVector> vNormFactor;
  forAll(vNormFactor, i) {
    vNormFactor[i] =
        gSum(mag(vResidual[i].primitiveField()) + mag(vR[i].primitiveField())) +
        VSMALL;
    initRes.setVector(i, gSumCmptMag(vR[i].primitiveField()) / vNormFactor[i]);
  }

  residualIO finalRes(initRes);

  // Solver admin
  bool solverStop = false;
  label solverIter = 0;

  // Initialisation
  using QRMatrix = QRMatrix<scalarRectangularMatrix>;

  // Initial dW is zero so r_0 is equal to R
  forN(nScalar, k) {
    sResidual[k].primitiveFieldRef() = sR[k].primitiveField();
  }
  forN(nVector, k) {
    vResidual[k].primitiveFieldRef() = vR[k].primitiveField();
  }

  Info << "  Initial residual: ";
  Info << finalRes;
  Info << endl;

  scalarRectangularMatrix H(nKrylov + 1, nKrylov, Zero);
  scalarField c(nKrylov + 1, Zero);
  scalarRectangularMatrix Pk;
  scalarDiagonalMatrix lambda;
  scalarRectangularMatrix Q;
  label eigen_count = 0;
  label prev_eigen_count = 0;
  labelList reduceRequests(nKrylov + 1);

  if (EnableRecycling && Yk.size() > 0) {
    prev_eigen_count = Yk.nCols();

    for (label i = 0; i < prev_eigen_count; ++i) {
      collect(Yk.subMatrix(0, i, Yk.nRows(), 1), sCk[i], vCk[i]);
      this->jacobian_.matrixMul(sCk[i], vCk[i], sUk[i], vUk[i]);
      expand(sUk[i], vUk[i], Yk.subMatrix(0, i, Yk.nRows(), 1));
    }

    QRMatrix decomposition(Yk, QRMatrix::modes::FULL,
                           QRMatrix::outputs::BOTH_QR,
                           QRMatrix::pivoting::FALSE);

    scalarSquareMatrix R(
        decomposition.R().subMatrix(0, 0, Yk.nCols(), Yk.nCols()));

    Yk = Yk * inverse(R);

    for (label i = 0; i < prev_eigen_count; ++i) {
      collect(decomposition.Q().subMatrix(0, i, Yk.nRows(), 1), sCk[i], vCk[i]);
      collect(Yk.subColumn(i), sUk[i], vUk[i]);
    }

    scalarField Ckr0(prev_eigen_count);
    for (label i = 0; i < prev_eigen_count; i++) {
      forAll(sCk[i], n) {
        Ckr0[i] +=
            sumProd(sCk[i][n].primitiveField(), sResidual[n].primitiveField());
      }
      forAll(vCk[i], n) {
        Ckr0[i] +=
            sumProd(vCk[i][n].primitiveField(), vResidual[n].primitiveField());
      }
      reduce(Ckr0[i], sumOp<scalar>(), Pstream::msgType(), UPstream::worldComm,
             reduceRequests[i]);
    }

    for (label i = 0; i < prev_eigen_count; i++) {
      if (Pstream::parRun() && reduceRequests[i] != -1)
        Pstream::waitRequest(reduceRequests[i]);

      forN(nScalar, n) {
        dsW[n].primitiveFieldRef() += sUk[i][n].primitiveField() * Ckr0[i];
        sResidual[n].primitiveFieldRef() -=
            sCk[i][n].primitiveField() * Ckr0[i];
      }
      forN(nVector, n) {
        dvW[n].primitiveFieldRef() += vUk[i][n].primitiveField() * Ckr0[i];
        vResidual[n].primitiveFieldRef() -=
            vCk[i][n].primitiveField() * Ckr0[i];
      }
    }

    prev_eigen_count = min(prev_eigen_count, nRitz);
  }

  for (solverIter = 0; solverIter < maxIter && !solverStop;
       ++solverIter) // Outer loop
  {
    // VPtr stands for W_m+1
    for (label j = 0; j < prev_eigen_count; j++) {
      forN(nScalar, n) {
        sVPtr[j][n].primitiveFieldRef() = sCk[j][n].primitiveField();
      }
      forN(nVector, n) {
        vVPtr[j][n].primitiveFieldRef() = vCk[j][n].primitiveField();
      }
    }

    // Ck = W_m+1 Q
    for (label j = 0; j < eigen_count; j++) {
      forN(nScalar, n) { sCk[j][n].primitiveFieldRef() = Zero; }
      forN(nVector, n) { vCk[j][n].primitiveFieldRef() = Zero; }

      for (label i = 0; i <= nKrylov; i++) {
        const PtrList<volScalarField> &sVi = sVPtr[i];
        const PtrList<volVectorField> &vVi = vVPtr[i];

        forN(nScalar, n) {
          sCk[j][n].primitiveFieldRef() += sVi[n].primitiveField() * Q[i][j];
        }
        forN(nVector, n) {
          vCk[j][n].primitiveFieldRef() += vVi[n].primitiveField() * Q[i][j];
        }
      }
    }

    // Update W_m+1 from Ck for this iteration
    // Note : new Krylov vectors not computed yet.
    for (label j = 0; j < eigen_count; j++) {
      forN(nScalar, n) {
        sVPtr[j][n].primitiveFieldRef() = sCk[j][n].primitiveField();
      }
      forN(nVector, n) {
        vVPtr[j][n].primitiveFieldRef() = vCk[j][n].primitiveField();
      }
    }

    if (eigen_count == 0)
      eigen_count = prev_eigen_count; // for deflated restart

    H = scalarRectangularMatrix(nKrylov + 1, nKrylov,
                                Zero); // G = [Dk, Bk; 0 H_m - k]
    c = scalarField(nKrylov + 1, Zero);
    auto Hmp1 = H.subMatrix(0, eigen_count, nKrylov + 1, nKrylov);

    // Begin reducing the Uk vector to calculate Dk
    for (label n = 0; n < eigen_count; ++n) {
      reduceRequests[n] = magsqr(sUk[n], vUk[n], H[n][n]);
    }

    arnoldi_iteration<nScalar, nVector>(this->jacobian_, this->preconditioner_,
                                        sResidual, vResidual, nKrylov,
                                        eigen_count, sVPtr, vVPtr, Hmp1, c);

    for (label n = 0; n < eigen_count; ++n) {
      if (Pstream::parRun() && reduceRequests[n] != -1) {
        Pstream::waitRequest(reduceRequests[n]);
      }
      H[n][n] = scale(sUk[n], vUk[n], 1 / sqrt(H[n][n])); // Dk
    }

    // VPtr now stands for V^_m
    for (label j = 0; j < eigen_count; j++) {
      forN(nScalar, n) {
        sVPtr[j][n].primitiveFieldRef() = sUk[j][n].primitiveField();
      }
      forN(nVector, n) {
        vVPtr[j][n].primitiveFieldRef() = vUk[j][n].primitiveField();
      }
    }

    // x = x + V^_m * d_m
    update_solution(sVPtr, vVPtr, H, c, nKrylov, sR, vR, dsW, dvW, sResidual,
                    vResidual);

    solverStop = check_residual(finalRes, initRes, tol, tolRel, sResidual,
                                vResidual, sNormFactor, vNormFactor);

    Info << "  GCRO-DR iteration: " << solverIter + 1 << "   Residual: ";
    Info << finalRes;
    Info << endl;

    forN(nScalar, i) { Info << " dS" << i << "=" << gSum(dsW[i]) << nl; }
    forN(nVector, i) { Info << " dV" << i << "=" << gSum(dvW[i]) << nl; }

    // Calculate new Uk
    // We dont need Ck yet so its going to be computed from a new iteration
    prev_eigen_count = eigen_count;
    if (eigen_count > 0) {
      scalarSquareMatrix eigenProblemA = scalarSquareMatrix(H.T() * H);
      List<labelList> reduceRequests2d(nKrylov);
      scalarRectangularMatrix wv(nKrylov + 1, nKrylov, Zero);

      // VPtr now stands for V^_m. And W_m+1 consists from 2 parts Ck and Vm+1
      // W_m+1 V^_m
      for (label j = 0; j < nKrylov; j++) {
        reduceRequests2d[j] = labelList(nKrylov + 1, -1);
        for (label i = 0; i < eigen_count; i++) {
          forAll(sCk[i], n) {
            wv(i, j) += sumProd(sCk[i][n].primitiveField(),
                                sVPtr[j][n].primitiveField());
          }
          forAll(vCk[i], n) {
            wv(i, j) += sumProd(vCk[i][n].primitiveField(),
                                vVPtr[j][n].primitiveField());
          }
          reduce(wv(i, j), sumOp<scalar>(), Pstream::msgType(),
                 UPstream::worldComm, reduceRequests2d[j][i]);
        }
        for (label i = eigen_count; i < nKrylov + 1; i++) {
          forAll(sVPtr[i], n) {
            wv(i, j) += sumProd(sVPtr[i][n].primitiveField(),
                                sVPtr[j][n].primitiveField());
          }
          forAll(vVPtr[i], n) {
            wv(i, j) += sumProd(vVPtr[i][n].primitiveField(),
                                vVPtr[j][n].primitiveField());
          }
          reduce(wv(i, j), sumOp<scalar>(), Pstream::msgType(),
                 UPstream::worldComm, reduceRequests2d[j][i]);
        }
      }

      if (Pstream::parRun()) {
        for (label j = 0; j < nKrylov; j++) {
          for (label i = 0; i < nKrylov + 1; i++) {
            if (reduceRequests2d[j][i] != -1)
              Pstream::waitRequest(reduceRequests2d[j][i]);
          }
        }
      }

      scalarSquareMatrix eigenProblemB = scalarSquareMatrix(H.T() * wv);

      eigen_count =
          eig_smallest(eigenProblemA, eigenProblemB, Pk, lambda, nRitz);

    } else {
      scalarSquareMatrix Hm = H.subMatrix(0, 0, nKrylov, nKrylov);
      scalarField e(nKrylov, Zero);
      e[nKrylov - 1] = 1;
      scalarSquareMatrix eigenProblemA = Hm + sqr(H(nKrylov, nKrylov - 1)) *
                                                  inverse(Hm.T()) *
                                                  transposeProduct(e);
      eigen_count = eig_smallest(eigenProblemA, Pk, lambda, nRitz);
    }

    if (eigen_count > 0) {
      QRMatrix decomposition(H * Pk, QRMatrix::modes::FULL,
                             QRMatrix::outputs::BOTH_QR,
                             QRMatrix::pivoting::FALSE);
      Q = decomposition.Q().subMatrix(0, 0, nKrylov + 1, eigen_count);
      auto Pk_invR = Pk * inverse(decomposition.R().subMatrix(0, 0, eigen_count,
                                                              eigen_count));

      // Uk = Vm Pk R^-1
      for (label j = 0; j < eigen_count; j++) {
        forN(nScalar, n) { sUk[j][n].primitiveFieldRef() = Zero; }
        forN(nVector, n) { vUk[j][n].primitiveFieldRef() = Zero; }

        for (label i = 0; i < nKrylov; i++) {
          const PtrList<volScalarField> &sVi = sVPtr[i];
          const PtrList<volVectorField> &vVi = vVPtr[i];

          forN(nScalar, n) {
            sUk[j][n].primitiveFieldRef() +=
                sVi[n].primitiveField() * Pk_invR[i][j];
          }
          forN(nVector, n) {
            vUk[j][n].primitiveFieldRef() +=
                vVi[n].primitiveField() * Pk_invR[i][j];
          }
        }
      }
    }
  }

  update_solution(sW, vW, dsW, dvW);

  forN(nScalar, i) { Info << " S" << i << "=" << gSum(sW[i]) << nl; }
  forN(nVector, i) { Info << " V" << i << "=" << gSum(vW[i]) << nl; }

  if (EnableRecycling) {
    label count = 0;
    forN(nScalar, n) { count += sResidual[n].size(); }
    forN(nVector, n) { count += vResidual[n].size() * 3; }
    Yk.resize(count, sUk.size());

    for (label i = 0; i < sUk.size(); ++i) {
      expand(sUk[i], vUk[i], Yk.subColumn(i));
    }
  }

  return solverIter;
}

template <int nScalar, int nVector>
void gcro<nScalar, nVector>::update_solution(
    const List<PtrList<volScalarField>> &sVPtr,
    const List<PtrList<volVectorField>> &vVPtr,
    const scalarRectangularMatrix &H, const scalarField &c, label nKrylov,
    const PtrList<volScalarField> &sR, const PtrList<volVectorField> &vR,
    PtrList<volScalarField> &dsW, PtrList<volVectorField> &dvW,
    PtrList<volScalarField> &sResidual,
    PtrList<volVectorField> &vResidual) const {
  scalarField y = lsq(H, c);
  // scalarField Hy = H * y;

  forN(nScalar, n) { sResidual[n].primitiveFieldRef() = Zero; }
  forN(nVector, n) { vResidual[n].primitiveFieldRef() = Zero; }

  for (label i = 0; i < nKrylov; i++) {
    const PtrList<volScalarField> &sVi = sVPtr[i];
    const PtrList<volVectorField> &vVi = vVPtr[i];
    const scalar &yi = y[i];
    // const scalar& Hyi = Hy[i];

    forN(nScalar, n) {
      dsW[n].primitiveFieldRef() += sVi[n].primitiveField() * yi;
      // sResidual[n].primitiveFieldRef() -= sVi[n].primitiveField() * Hyi;
    }
    forN(nVector, n) {
      dvW[n].primitiveFieldRef() += vVi[n].primitiveField() * yi;
      // vResidual[n].primitiveFieldRef() -= vVi[n].primitiveField() * Hyi;
    }
  }

  /*forN(nScalar, n) {
      sResidual[n].primitiveFieldRef() -= sVPtr[nKrylov][n].primitiveField() *
  Hy[nKrylov];
  }
  forN(nVector, n) {
      vResidual[n].primitiveFieldRef() -= vVPtr[nKrylov][n].primitiveField() *
  Hy[nKrylov];
  }*/

  this->jacobian_.matrixMul(dsW, dvW, sResidual, vResidual);

  forN(nScalar, n) {
    sResidual[n].primitiveFieldRef() =
        sR[n].primitiveField() - sResidual[n].primitiveField();
  }
  forN(nVector, n) {
    vResidual[n].primitiveFieldRef() =
        vR[n].primitiveField() - vResidual[n].primitiveField();
  }
}

template <int nScalar, int nVector>
void gcro<nScalar, nVector>::update_solution(
    PtrList<volScalarField> &sW, PtrList<volVectorField> &vW,
    const PtrList<volScalarField> &dsW,
    const PtrList<volVectorField> &dvW) const {

  forN(nScalar, k) {
    sW[k].primitiveFieldRef() += dsW[k].primitiveField();
    sW[k].correctBoundaryConditions();
  }
  forN(nVector, k) {
    vW[k].primitiveFieldRef() += dvW[k].primitiveField();

    // Zero the variable in non-solved directions
    vector::labelType validComponents =
        this->mesh_.solutionD(); //-1 for empty directions
    forAll(validComponents, cmpt) {
      if (validComponents[cmpt] == -1) {
        vW[k].replace(cmpt, dimensionedScalar("0", vW[k].dimensions(), 0.0));
      }
    }

    vW[k].correctBoundaryConditions();
  }
}

template <int nScalar, int nVector>
bool gcro<nScalar, nVector>::check_residual(
    residualIO &finalRes, const residualIO &initRes, const residualIO &tol,
    const residualIO &tolRel, const PtrList<volScalarField> &sResidual,
    const PtrList<volVectorField> &vResidual,
    const FixedList<scalar, nScalar> &sNormFactor,
    const FixedList<scalar, nVector> &vNormFactor) const {
  forN(nScalar, k) {
    finalRes.setScalar(k,
                       gSumMag(sResidual[k].primitiveField()) / sNormFactor[k]);
  }
  forN(nVector, k) {
    vector res = gSumCmptMag(vResidual[k].primitiveField()) / vNormFactor[k];
    // Zero the residual in non-solved directions
    vector::labelType validComponents =
        this->mesh_.solutionD(); //-1 for empty directions
    forAll(validComponents, cmpt) {
      if (validComponents[cmpt] == -1) {
        res[cmpt] = 0.0;
      }
    }
    finalRes.setVector(k, res);
  }

  return (max(finalRes / tol) < 1.0) ||
         (max(finalRes / (initRes + ROOTVSMALL) / (tolRel + ROOTVSMALL)) < 1.0);
}

} // namespace Foam

// ************************************************************************* //
