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

#ifndef Eigen_H
#define Eigen_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "jacobianMatrix.H"
#include "preconditioner.H"
#include "EigenMatrix.H"
#include "QRMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

inline scalarSquareMatrix inverse(scalarSquareMatrix A) {
    scalarSquareMatrix invA(A.n());
    labelList pivotIndices(A.n());

    LUDecompose(A, pivotIndices);

    for (label i = 0; i < A.n(); ++i) {
        scalarField source(A.n(), 0);
        source[i] = 1;
        LUBacksubstitute(A, pivotIndices, source);
        invA.subColumn(i) = source;
    }

    return invA;
}

inline label __get_smallest_real_eigs(scalarRectangularMatrix& v, scalarDiagonalMatrix& lambdaRe, const scalarDiagonalMatrix& lambdaIm, label k) {
    auto unary_op_reel = [](scalar val) {
        return val == 0;
    };

    label reel_eig_count = std::count_if(lambdaIm.begin(), lambdaIm.end(), unary_op_reel);

    if (reel_eig_count < k) {
        Info << "There isn't enough real eigen pairs (" << k << " requested " << reel_eig_count << " available)" << nl;
        k = reel_eig_count;
    }

    if (k == 0) return 0;

    auto binary_op_reel = [&](scalar lhs, scalar rhs) {
        return lhs == 0 && rhs != 0;
    };

    labelList order = lambdaIm.sortPermutation(binary_op_reel);

    lambdaRe.applyPermutation(order);
    List<bool> pass(order.size(), false);
    for (label i = 0; i < order.size(); ++i) {
        if (pass[i]) continue;
        pass[i] = true;
        label prev = i;
        label j = order[i];
        while (i != j) {
            for (label k = 0; k < v.nRows(); ++k)
                Swap(v[k][prev], v[k][j]);
            pass[j] = true;
            prev = j;
            j = order[j];
        }
    }

    lambdaRe.resize(reel_eig_count);

    if (k < reel_eig_count) {
        auto op_smallest = [](scalar lhs, scalar rhs) {
            return lhs > rhs;
        };
        order = lambdaRe.sortPermutation(op_smallest);
        lambdaRe.applyPermutation(order);
        pass = List<bool>(order.size(), false);
        for (label i = 0; i < order.size(); ++i) {
            if (pass[i]) continue;
            pass[i] = true;
            label prev = i;
            label j = order[i];
            while (i != j) {
                for (label k = 0; k < v.nRows(); ++k)
                    Swap(v[k][prev], v[k][j]);
                pass[j] = true;
                prev = j;
                j = order[j];
            }
        }
    }
    lambdaRe.resize(k);
    v = v.subMatrix(0, 0, v.n(), k);

    return k;
}

inline label eig_smallest(const scalarSquareMatrix& A, scalarRectangularMatrix& v, scalarDiagonalMatrix& lambda, label k) {
    using EigenMatrix = EigenMatrix < scalar>;

    EigenMatrix decomposition(A);

    lambda = decomposition.EValsRe();
    v = decomposition.EVecs();

    return __get_smallest_real_eigs(v, lambda, decomposition.EValsIm(), k);
}

inline label eig_smallest(const scalarSquareMatrix& A, const scalarSquareMatrix& B, scalarRectangularMatrix& v, scalarDiagonalMatrix& lambda, label k) {
    using EigenMatrix = EigenMatrix < scalar>;

    EigenMatrix decomposition(inverse(B) * A);

    lambda = decomposition.EValsRe();
    v = decomposition.EVecs();

    return __get_smallest_real_eigs(v, lambda, decomposition.EValsIm(), k);
}

inline scalarField lsq(const scalarRectangularMatrix& H, const scalarField& c) {
    using QRMatrix = QRMatrix < scalarRectangularMatrix>;

    QRMatrix decomposition(H, QRMatrix::modes::FULL, QRMatrix::outputs::BOTH_QR, QRMatrix::pivoting::FALSE);
    scalarField y( decomposition.solve(c));
    y.resize(H.n());
    return y;
}

inline scalarSquareMatrix transposeProduct(const scalarField& v) {
    scalarSquareMatrix result(v.size());

    for (label i = 0; i < v.size(); ++i) {
        for (label j = 0; j < v.size(); ++j) {
            result(i, j) = v[i] * v[j];
        }
    }

    return result;
}

inline label magsqr(const PtrList<volScalarField>& s, const PtrList<volVectorField>& v, scalar& magnitude) {
    label reduceRequest = -1;

    magnitude = 0.0;
    forAll(s, n) {
        magnitude += sumSqr(s[n].primitiveField());
    }
    forAll(v, n) {
        magnitude += sum(magSqr(v[n].primitiveField()));
    }
    reduce(magnitude, sumOp<scalar>(), Pstream::msgType(), UPstream::worldComm, reduceRequest);

    return reduceRequest;
}

inline scalar scale(PtrList<volScalarField>& s, PtrList<volVectorField>& v, scalar scaling) {
    forAll(s, n) {
        s[n].primitiveFieldRef() *= scaling;
    }
    forAll(v, n) {
        v[n].primitiveFieldRef() *= scaling;
    }

    return scaling;
}

inline scalar normalize(PtrList<volScalarField>& s, PtrList<volVectorField>& v, scalar norm_factor) {
    scale(s, v, 1 / norm_factor);

    return norm_factor;
}

inline scalar normalize(PtrList<volScalarField>& s, PtrList<volVectorField>& v) {
    scalar beta;

    label reduceRequest = magsqr(s, v, beta);

    if (Pstream::parRun() && reduceRequest != -1) {
        Pstream::waitRequest(reduceRequest);
    }
    beta = Foam::sqrt(beta);

    return normalize(s, v, beta);
}

template<typename T>
inline void expand(const PtrList<volScalarField>& s, const PtrList<volVectorField>& v, T&& vector) {
    label size = 0;
    label count = s.size() + v.size() * 3;

    forAll(s, i) {
        size += s[i].size();
    }

    forAll(v, i) {
        size += v[i].size() * 3;
    }

    for (label i = 0; i < size / count; ++i) {
        forAll(s, j) {
            vector(i * count + j, 0) = s[j][i];
        }
        forAll(v, j) {
            vector(i * count + s.size() + j * v.size() + 0, 0) = v[j][i].x();
            vector(i * count + s.size() + j * v.size() + 1, 0) = v[j][i].y();
            vector(i * count + s.size() + j * v.size() + 2, 0) = v[j][i].z();
        }
    }
}

template<typename T>
inline void collect(const T&& vector, PtrList<volScalarField>& s, PtrList<volVectorField>& v) {
    label size = 0;
    label count = s.size() + v.size() * 3;

    forAll(s, i) {
        size += s[i].size();
    }

    forAll(v, i) {
        size += v[i].size() * 3;
    }

    for (label i = 0; i < size / count; ++i) {
        forAll(s, j) {
            s[j][i] = vector(i * count + j, 0);
        }
        forAll(v, j) {
            v[j][i].x() = vector(i * count + s.size() + j * v.size() + 0, 0);
            v[j][i].y() = vector(i * count + s.size() + j * v.size() + 1, 0);
            v[j][i].z() = vector(i * count + s.size() + j * v.size() + 2, 0);
        }
    }
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

