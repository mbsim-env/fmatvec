/* Copyright (C) 2003-2005  Martin Förg

 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact:
 *   martin.o.foerg@googlemail.com
 *
 */

#include "config.h"
#include "linear_algebra.h"
#include "linear_algebra_fixed.h"
#include "linear_algebra_var.h"
#include "linear_algebra_fixed_var.h"
#include "linear_algebra_var_fixed.h"
#include "linear_algebra_double.h"
#include "wrapper.h"
#include <stdexcept>
#include <sstream>

#define CVT_TRANSPOSE(c) \
   (((c) == CblasNoTrans) ? 'N' : \
    ((c) == CblasTrans) ? 'T' : \
    ((c) == CblasConjTrans) ? 'C' : \
    -1)

#define CVT_UPLO(c) \
   (((c) == CblasUpper) ? 'U' : \
    ((c) == CblasLower) ? 'L' : \
    -1)

#define CVT_SIDE(c) \
   (((c) == CblasLeft) ? 'L' : \
    ((c) == CblasRight) ? 'R' : \
    -1)

template <class T>
std::string toStr(const T& val) {
  std::ostringstream oss;
  oss << val;
  return oss.str();
}

//-------------------------------------
// Matrix operations
//-------------------------------------
namespace fmatvec {

  //-------------------------------------
  // operations
  //-------------------------------------

  Matrix<General, Ref, Ref, double> slvLU(const SquareMatrix<Ref, double> &A, const Matrix<General, Ref, Ref, double> &X) {

    assert(A.size() == X.rows());

    Matrix<General, Ref, Ref, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    SquareMatrix<Ref, double> B = A;

    auto *ipiv = new int[A.size()];

    int info = dgesv(B.blasOrder(), B.size(), Y.cols(), B(), B.ldim(), ipiv, Y(), Y.ldim());

    delete[] ipiv;

    if(info != 0)
      throw std::runtime_error("Exception in slvLU: dgesv exited with info="+toStr(info));

    return Y;
  }

  Matrix<General, Var, Var, double> slvLU(const SquareMatrix<Var, double> &A, const Matrix<General, Var, Var, double> &X, int &info) {

    assert(A.size() == X.rows());

    Matrix<General, Var, Var, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    SquareMatrix<Var, double> B = A;

    auto *ipiv = new int[A.size()];

    info = dgesv(B.blasOrder(), B.size(), Y.cols(), B(), B.ldim(), ipiv, Y(), Y.ldim());

    delete[] ipiv;

    return Y;
  }

  Matrix<General, Ref, Ref, double> slvLUFac(const SquareMatrix<Ref, double> &A, const Matrix<General, Ref, Ref, double> &X, const Vector<Ref, int> &ipiv) {

    assert(A.size() == X.rows());

    Matrix<General, Ref, Ref, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    SquareMatrix<Ref, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dgetrs(B.blasOrder(), B.blasTrans(), B.size(), Y.cols(), B(), B.ldim(), ipiv(), Y(), Y.ldim());
#else
    int info = dgetrs(B.blasOrder(), CVT_TRANSPOSE(B.blasTrans()), B.size(), Y.cols(), B(), B.ldim(), ipiv(), Y(), Y.ldim());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLUFac: dgetrs exited with info="+toStr(info));

    return Y;
  }

  Matrix<General, Var, Var, double> slvLUFac(const SquareMatrix<Var, double> &A, const Matrix<General, Var, Var, double> &X, const Vector<Var, int> &ipiv) {

    assert(A.size() == X.rows());

    Matrix<General, Var, Var, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    SquareMatrix<Var, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dgetrs(B.blasOrder(), B.blasTrans(), B.size(), Y.cols(), B(), B.ldim(), ipiv(), Y(), Y.ldim());
#else
    int info = dgetrs(B.blasOrder(), CVT_TRANSPOSE(B.blasTrans()), B.size(), Y.cols(), B(), B.ldim(), ipiv(), Y(), Y.ldim());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLUFac: dgetrs exited with info="+toStr(info));

    return Y;
  }

  Matrix<General, Ref, Ref, double> slvQR(const SquareMatrix<Ref, double> &A, const Matrix<General, Ref, Ref, double> &X) {

    assert(A.size() == X.rows());

    Matrix<General, Ref, Ref, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    SquareMatrix<Ref, double> B = A;

    int info = dgels(B.blasTrans(), B.rows(), B.cols(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());

    if(info != 0)
      throw std::runtime_error("Exception in slvQR: dgels exited with info="+toStr(info));

    return Y;
  }

  int Doolittle_LU_with_Pivoting_Solve(double *A, double B[], const int pivot[], double x[], int n) {
    int i, k;
    double *p_k;
    double dum;

    //         Solve the linear equation Lx = B for x, where L is a lower
    //         triangular matrix with an implied 1 along the diagonal.

    for (k = 0, p_k = A; k < n; p_k += n, k++) {
      if (pivot[k] != k) {
        dum = B[k];
        B[k] = B[pivot[k]];
        B[pivot[k]] = dum;
      }
      x[k] = B[k];
      for (i = 0; i < k; i++)
        x[k] -= x[i] * *(p_k + i);
    }

    //         Solve the linear equation Ux = y, where y is the solution
    //         obtained above of Lx = B and U is an upper triangular matrix.

    for (k = n - 1, p_k = A + n * (n - 1); k >= 0; k--, p_k -= n) {
      if (pivot[k] != k) {
        dum = B[k];
        B[k] = B[pivot[k]];
        B[pivot[k]] = dum;
      }
      for (i = k + 1; i < n; i++)
        x[k] -= x[i] * *(p_k + i);
      if (fabs(*(p_k + k)) < 1e-10) //TODO: use eps here
        return -1;
      x[k] /= *(p_k + k);
    }

    return 0;
  }

  Vector<Ref, double> slvLU(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x) {

    assert(A.size() == x.size());

    Vector<Ref, double> y = x;

    if (x.size() == 0)
      return y;

    SquareMatrix<Ref, double> B = A;

    auto *ipiv = new int[A.size()];

    int info = dgesv(B.blasOrder(), B.size(), 1, B(), B.ldim(), ipiv, y(), y.size());

    delete[] ipiv;

    if(info != 0)
      throw std::runtime_error("Exception in slvLU: dgesv exited with info="+toStr(info));

    return y;
  }

  Vector<Ref, double> slvLU(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x, int & info) {

    assert(A.size() == x.size());

    Vector<Ref, double> y = x;

    if (x.size() == 0)
      return y;

    SquareMatrix<Ref, double> B = A;

    auto *ipiv = new int[A.size()];

    info = dgesv(B.blasOrder(), B.size(), 1, B(), B.ldim(), ipiv, y(), y.size());

    delete[] ipiv;

    return y;
  }

  Vector<Var, double> slvLU(const SquareMatrix<Var, double> &A, const Vector<Var, double> &x, int & info) {

    assert(A.size() == x.size());

    Vector<Var, double> y = x;

    if (x.size() == 0)
      return y;

    SquareMatrix<Var, double> B = A;

    auto *ipiv = new int[A.size()];

    info = dgesv(B.blasOrder(), B.size(), 1, B(), B.ldim(), ipiv, y(), y.size());

    delete[] ipiv;

    return y;
  }

  Vector<Ref, double> slvLUFac(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x, const Vector<Ref, int> &ipiv) {

    assert(A.size() == x.size());

    Vector<Ref, double> y = x;

    if (x.size() == 0)
      return y;

    SquareMatrix<Ref, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dgetrs(B.blasOrder(), B.blasTrans(), B.size(), 1, B(), B.ldim(), ipiv(), y(), y.size());
#else
    int info = dgetrs(B.blasOrder(), CVT_TRANSPOSE(B.blasTrans()), B.size(), 1, B(), B.ldim(), ipiv(), y(), y.size());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLUFac: dgetrs exited with info="+toStr(info));

    return y;
  }

  Vector<Var, double> slvLUFac(const SquareMatrix<Var, double> &A, const Vector<Var, double> &x, const Vector<Var, int> &ipiv) {

    assert(A.size() == x.size());

    Vector<Var, double> y = x;

    if (x.size() == 0)
      return y;

    SquareMatrix<Var, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dgetrs(B.blasOrder(), B.blasTrans(), B.size(), 1, B(), B.ldim(), ipiv(), y(), y.size());
#else
    int info = dgetrs(B.blasOrder(), CVT_TRANSPOSE(B.blasTrans()), B.size(), 1, B(), B.ldim(), ipiv(), y(), y.size());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLUFac: dgetrs exited with info="+toStr(info));

    return y;
  }

  Vector<Ref, double> slvQR(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x) {

    assert(A.size() == x.size());

    Vector<Ref, double> y = x;

    if (x.size() == 0)
      return y;

    SquareMatrix<Ref, double> B = A;

    int info = dgels(B.blasTrans(), B.rows(), B.cols(), y.cols(), B(), B.ldim(), y(), y.size());

    if(info != 0)
      throw std::runtime_error("Exception in slvQR: dgels exited with info="+toStr(info));

    return y;
  }

  SquareMatrix<Ref, double> inv(const SquareMatrix<Ref, double> &A) {

    SquareMatrix<Ref, double> B = A;

    if (A.size() == 0)
      return B;

    auto *ipiv = new int[A.size()];

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv);

    if(info != 0)
      throw std::runtime_error("Exception in inv: dgetrf exited with info="+toStr(info));

    dgetri(B.blasOrder(), B.size(), B(), B.ldim(), ipiv);

    delete[] ipiv;

    return B;
  }

  Matrix<Symmetric, Ref, Ref, double> inv(const Matrix<Symmetric, Ref, Ref, double> &A) {

    Matrix<Symmetric, Ref, Ref, double> B = A;

    if (A.size() == 0)
      return B;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dpotrf(B.blasOrder(), B.blasUplo(), B.size(), B(), B.ldim());
#else
    int info = dpotrf(B.blasOrder(),CVT_UPLO(B.blasUplo()) , B.size(), B(), B.ldim());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in inv: dpotrf exited with info="+toStr(info));

#ifndef HAVE_LIBMKL_INTEL_LP64
    dpotri(B.blasOrder(), B.blasUplo(), B.rows(), B(), B.ldim());
#else
    dpotri( B.blasOrder(),CVT_UPLO(B.blasUplo()) , B.rows(), B(), B.ldim());
#endif

    return B;
  }

  Matrix<Diagonal, Ref, Ref, double> inv(const Matrix<Diagonal, Ref, Ref, double> &A) {

    Matrix<Diagonal, Ref, Ref, double> B(A.size());

    if (A.size() == 0)
      return B;

    for (int i = 0; i < A.size(); i++)
      B(i) = 1.0 / A(i);

    return B;
  }

  Matrix<General, Ref, Ref, double> facLU(const Matrix<General, Ref, Ref, double> &A, Vector<Ref, int> &ipiv) {

    Matrix<General, Ref, Ref, double> B = A;

    if (A.rows() == 0 || A.cols() == 0)
      return B;

    int is = A.rows() < A.cols() ? A.rows() : A.cols();
    if (ipiv.size() != is) {
      ipiv.resize(is);
    }

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv());

    if(info != 0)
      throw std::runtime_error("Exception in facLU: dgetrf exited with info="+toStr(info));

    return B;
  }

  SquareMatrix<Ref, double> facLU(const SquareMatrix<Ref, double> &A, Vector<Ref, int> &ipiv) {

    SquareMatrix<Ref, double> B = A;

    if (A.size() == 0)
      return B;

    int is = A.size();
    if (ipiv.size() != is) {
      ipiv.resize(is);
    }

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv());

    if(info != 0)
      throw std::runtime_error("Exception in facLU: dgetrf exited with info="+toStr(info));

    return B;
  }

  SquareMatrix<Var, double> facLU(const SquareMatrix<Var, double> &A, Vector<Var, int> &ipiv) {

    SquareMatrix<Var, double> B = A;

    if (A.size() == 0)
      return B;

    int is = A.size();
    if (ipiv.size() != is) {
      ipiv.resize(is);
    }

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv());

    if(info != 0)
      throw std::runtime_error("Exception in facLU: dgetrf exited with info="+toStr(info));

    return B;
  }

  int facLU(double *A, int pivot[], int n) {
      int i, j, k;
      double *p_k, *p_row, *p_col = nullptr;
      double max;

      //         For each row and column, k = 0, ..., n-1,

      for (k = 0, p_k = A; k < n; p_k += n, k++) {

        //            find the pivot row

        pivot[k] = k;
        max = fabs(*(p_k + k));
        for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
          if (max < fabs(*(p_row + k))) {
            max = fabs(*(p_row + k));
            pivot[k] = j;
            p_col = p_row;
          }
        }

        //     and if the pivot row differs from the current row, then
        //     interchange the two rows.

        if (pivot[k] != k)
          for (j = 0; j < n; j++) {
            max = *(p_k + j);
            *(p_k + j) = *(p_col + j);
            *(p_col + j) = max;
          }

        //                and if the matrix is singular, return error

        if (fabs(*(p_k + k)) < 1.e-10) //TODO: use eps
          return -1;

        //      otherwise find the lower triangular matrix elements for column k.

        for (i = k + 1, p_row = p_k + n; i < n; p_row += n, i++) {
          *(p_row + k) /= *(p_k + k);
        }

        //            update remaining matrix

        for (i = k + 1, p_row = p_k + n; i < n; p_row += n, i++)
          for (j = k + 1; j < n; j++)
            *(p_row + j) -= *(p_row + k) * *(p_k + j);

      }

      return 0;
    }


  Matrix<Symmetric, Ref, Ref, double> facLL(const Matrix<Symmetric, Ref, Ref, double> &A) {

    Matrix<Symmetric, Ref, Ref, double> B = A;

    if (A.size() == 0)
      return B;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dpotrf(B.blasOrder(), B.blasUplo(), B.size(), B(), B.ldim());
#else
    int info = dpotrf(B.blasOrder(), CVT_UPLO(B.blasUplo()), B.size(), B(), B.ldim());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in facLL: dpotrf exited with info="+toStr(info));

    return B;
  }

  double nrm1(const Vector<Ref, double> &x) {

    if (x.size() == 0)
      return 0.0;

    return dasum(x.size(), x(), x.inc());
  }

  double nrmInf(const Vector<Ref, double> &x) {

    if (x.size() == 0)
      return 0.0;

    int id = idamax(x.size(), x(), x.inc());
    return fabs(x(id));
  }

  double nrm2(const Vector<Ref, double> &x) {

    if (x.size() == 0)
      return 0.0;

    return dnrm2(x.size(), x(), x.inc());
  }

  Vector<Ref, double> slvLL(const Matrix<Symmetric, Ref, Ref, double> &A, const Vector<Ref, double> &x) {

    assert(A.size() == x.size());

    Vector<Ref, double> y = x;

    if (x.size() == 0)
      return y;

    Matrix<Symmetric, Ref, Ref, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dposv(B.blasOrder(), B.blasUplo(), B.size(), 1, B(), B.ldim(), y(), y.size());
#else
    int info = dposv(B.blasOrder(), CVT_UPLO(B.blasUplo()), B.size(), 1, B(), B.ldim(), y(), y.size());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLL: dposv exited with info="+toStr(info));

    return y;
  }

  Matrix<General, Ref, Ref, double> slvLL(const Matrix<Symmetric, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &X) {

    assert(A.size() == X.rows());

    Matrix<General, Ref, Ref, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    Matrix<Symmetric, Ref, Ref, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dposv(B.blasOrder(), B.blasUplo(), B.size(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());
#else
    int info = dposv(B.blasOrder(), CVT_UPLO(B.blasUplo()), B.size(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLL: dposv exited with info="+toStr(info));

    return Y;
  }

  Vector<Ref, double> slvLLFac(const Matrix<Symmetric, Ref, Ref, double> &A, const Vector<Ref, double> &x) {

    assert(A.size() == x.size());

    Vector<Ref, double> y = x;

    if (x.size() == 0)
      return y;

    Matrix<Symmetric, Ref, Ref, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dpotrs(B.blasOrder(), B.blasUplo(), B.size(), 1, B(), B.ldim(), y(), y.size());
#else
    int info = dpotrs(B.blasOrder(),CVT_UPLO(B.blasUplo()) , B.size(), 1, B(), B.ldim(), y(), y.size());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLLFac: dpotrs exited with info="+toStr(info));

    return y;
  }

  Matrix<General, Ref, Ref, double> slvLLFac(const Matrix<Symmetric, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &X) {

    assert(A.size() == X.rows());

    Matrix<General, Ref, Ref, double> Y = X;

    if (X.rows() == 0 || X.cols() == 0)
      return Y;

    Matrix<Symmetric, Ref, Ref, double> B = A;

#ifndef HAVE_LIBMKL_INTEL_LP64
    int info = dpotrs(B.blasOrder(), B.blasUplo(), B.size(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());
#else
    int info = dpotrs(B.blasOrder(),CVT_UPLO(B.blasUplo()) , B.size(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());
#endif

    if(info != 0)
      throw std::runtime_error("Exception in slvLLFac: dpotrs exited with info="+toStr(info));

    return Y;
  }

  Vector<Ref, std::complex<double>> eigval(const SquareMatrix<Ref, double> &A) {

    double *vl = nullptr, *vr = nullptr;
    auto *wr = new double[A.size()];
    auto *wi = new double[A.size()];

    SquareMatrix<Ref, double> B = A;

    dgeev('N', 'N', A.size(), B(), B.ldim(), wr, wi, vl, B.size(), vr, B.size());

    Vector<Ref, std::complex<double>> w(A.size());
    for (int i = 0; i < A.size(); i++)
      w(i) = std::complex<double>(wr[i], wi[i]);

    delete[] wr;
    delete[] wi;

    return w;
  }

  int eigvec(const SquareMatrix<Ref, double> &A, SquareMatrix<Ref, std::complex<double>> &V, Vector<Ref, std::complex<double>> &w) {

    SquareMatrix<Ref, double> B = A;
    w.resize(A.size(),NONINIT);
    V.resize(A.size(),NONINIT);
    SquareMatrix<Ref, double> Vreal(A.size(),NONINIT);

    double *vl = nullptr;
    auto *wr = new double[A.size()];
    auto *wi = new double[A.size()];

    int info = dgeev('N', 'V', A.size(), B(), B.ldim(), wr, wi, vl, B.size(), Vreal(), B.size());

    for (int i = 0; i < A.size(); i++) 
      w(i) = std::complex<double>(wr[i], wi[i]);
    for (int i = 0; i < A.size(); i++) {
      if((i < A.size()-1) and (w(i+1)==conj(w(i)))) {
        for(int j = 0; j < A.size(); j++) {
          V(j,i) = std::complex<double>(Vreal(j,i), Vreal(j,i+1));
          V(j,i+1) = std::complex<double>(Vreal(j,i), -Vreal(j,i+1));
        }
        i++;
      }
      else {
        for(int j = 0; j < A.size(); j++) {
          V(j,i) = std::complex<double>(Vreal(j,i), 0);
        }
      }
    }

    delete[] wr;
    delete[] wi;

    return info;
  }
  
  int eigvec(const Matrix<Symmetric, Ref, Ref, double> &A, const Matrix<Symmetric, Ref, Ref, double> &B, SquareMatrix<Ref, double> &eigenvectors, Vector<Ref, double> &eigenvalues) {
    const int dim = A.size();
    auto *w = new double[dim];
    SquareMatrix<Ref, double> B_(dim);
    eigenvectors.resize(dim);
    eigenvalues.resize(dim);
    for (int z = 0; z < dim; z++)
      for (int s = 0; s <= z; s++) {
        eigenvectors(z, s) = A(z, s);
        eigenvectors(s, z) = A(z, s);
        B_(z, s) = B(z, s);
        B_(s, z) = B(z, s);
      }

    int info = dsygv(1, 'V', 'L', dim, eigenvectors(), eigenvectors.ldim(), B_(), B_.ldim(), w);

    for (int i = 0; i < dim; i++) {
      Vector<Ref, double> evTmp = eigenvectors.col(i);
      if (nrm2(evTmp) > 0)
        eigenvectors.set(i, evTmp / nrm2(evTmp));
      eigenvalues(i) = w[i];
    }

    delete[] w;

    return info;
  }

  int svd(Matrix<General, Ref, Ref, double> &A,    Matrix<General, Ref, Ref, double> &S,  SquareMatrix<Ref, double> &U, SquareMatrix<Ref, double> &VT, int Rueckgabe) {
	  const int dim_m = A.rows();
	  const int dim_n = A.cols();

	  U.resize(dim_m);
	  S.resize(dim_m,dim_n);
	  VT.resize(dim_n);

	  const int lda = A.ldim();
	  const int ldu = U.ldim();
	  const int ldvt = VT.ldim();

	  int info = dgesvd('A','A',dim_m, dim_n, A(), lda, S(), U(), ldu, VT(), ldvt);

	  if(Rueckgabe==1) {
		  int i=0;
		  int minimum;

		  if(dim_m < dim_n) minimum = dim_m;
		  else minimum = dim_n;

		  for (i=1; i<minimum; i++) {
			  S(i,i)=S(i,0);
			  S(i,0) = 0.;
		  }
	  }

	  return info;
  }

  Vector<Ref, double> eigval(const Matrix<Symmetric, Ref, Ref, double> &A) {

    Vector<Ref, double> w(A.size(), NONINIT);
    Matrix<Symmetric, Ref, Ref, double> B = A;

    dsyev('N', 'L', B.size(), B(), B.ldim(), w());

    return w;

  }

  Vector<Ref, double> eigvalSel(const Matrix<Symmetric, Ref, Ref, double> &A, int il, int iu, double abstol) {

    assert(il >= 1);
    assert(iu >= il);
    assert(iu <= A.size());

    Matrix<Symmetric, Ref, Ref, double> B = A;
    B.ldim();
    double vl = 0, vu = 0;
    int m;
    Vector<Ref, double> w(A.size());
    double* z = nullptr;
    int ldz = 1;

    dsyevx('N', 'I', B.blasUplo(), B.size(), B(), B.ldim(), vl, vu, il, iu, abstol, &m, w(), z, ldz);

    return w(Range<Var,Var>(0, m - 1));

  }

  double rho(const SquareMatrix<Ref, double> &A) {

    if (A.rows() == 0 || A.cols() == 0)
      return 0.0;
    Vector<Ref, std::complex<double>> v = eigval(A);
    double buf = abs(v(0));
    for (int i = 0; i < v.size(); i++) {
      double absi = abs(v(i));
      if (absi >= buf)
        buf = absi;
    }
    return buf;
  }

  double rho(const Matrix<Symmetric, Ref, Ref, double> &A) {

    if (A.rows() == 0 || A.cols() == 0)
      return 0.0;
    Vector<Ref, double> v = eigval(A);
    double buf = fabs(v(0));
    for (int i = 0; i < v.size(); i++) {
      double absi = fabs(v(i));
      if (absi >= buf)
        buf = absi;
    }
    return buf;
  }

  double nrmInf(const Matrix<General, Ref, Ref, double> &A) {

    if (A.rows() == 0 || A.cols() == 0)
      return 0.0;

    return dlange('I', A.rows(), A.cols(), A(), A.ldim());
  }

  double nrmInf(const Matrix<Symmetric, Ref, Ref, double> &A) {

      if (A.rows() == 0 || A.cols() == 0)
        return 0.0;
    
      return dlansy('I', 'L', A.cols(), A(), A.ldim());
    }

  double nrm1(const Matrix<General, Ref, Ref, double> &A) {

    if (A.rows() == 0 || A.cols() == 0)
      return 0.0;

    return dlange('1', A.rows(), A.cols(), A(), A.ldim());
  }

  double nrmFro(const Matrix<General, Ref, Ref, double> &A) {

    if (A.rows() == 0 || A.cols() == 0)
      return 0.0;

    return dlange('F', A.rows(), A.cols(), A(), A.ldim());
  }

  double nrm2(const Matrix<General, Ref, Ref, double> &A) {

    if (A.rows() == 0 || A.cols() == 0)
      return 0.0;

    return sqrt(rho(JTJ(A)));
  }

  Matrix<General, Ref, Ref, double> slvLS(const Matrix<General, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &B, double rcond) {

    assert(A.rows() == B.rows());

    //    if(X.rows() == 0 || X.cols() == 0)
    //      return Y;

    Matrix<General, Ref, Ref, double> A_ = A;

    Matrix<General, Ref, Ref, double> B_(A.rows() > A.cols() ? A.rows() : A.cols(), B.cols(), NONINIT);
    B_.set(Range<Var,Var>(0, B.rows() - 1), Range<Var,Var>(0, B.cols() - 1), B);

    int info = dgelss(A.rows(), A.cols(), B_.cols(), A_(), A_.ldim(), B_(), B_.ldim(), rcond);

    if(info != 0)
      throw std::runtime_error("Exception in slvLS: dgelss exited with info="+toStr(info));

    return B_(Range<Var,Var>(0, A.cols() - 1), Range<Var,Var>(0, B.cols() - 1));
  }

  Vector<Ref, double> slvLS(const Matrix<General, Ref, Ref, double> &A, const Vector<Ref, double> &b, double rcond) {

    assert(A.rows() == b.size());

    //if(b.size() == 0)
    //return y;

    Matrix<General, Ref, Ref, double> A_ = A;

    Vector<Ref, double> b_(A.rows() > A.cols() ? A.rows() : A.cols(), NONINIT);
    b_.set(Range<Var,Var>(0, b.size() - 1), b);

    int info = dgelss(A.rows(), A.cols(), 1, A_(), A_.ldim(), b_(), b_.size(), rcond);

    if(info != 0)
      throw std::runtime_error("Exception in slvLS: dgelss exited with info="+toStr(info));

    return b_(Range<Var,Var>(0, A.cols() - 1));
  }

//  Matrix<General,Ref,Ref,double> swap(const Matrix<General,Ref,Ref,double> &X, const Vector<Ref,int> &ipiv ) {
//
//    //    assert(A.size() == X.rows());
//
//    Matrix<General,Ref,Ref,double> Y = X;
//
//    ATL_dlaswp( Y.cols(), Y(), Y.ldim(), 0, Y.rows(), ipiv(), 1 );
//
//    return Y;
//  }
//  Matrix<General,Ref,Ref,double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_DIAG unit, const SquareMatrix<Ref,double> &A, const Matrix<General,Ref,Ref,double> &X, const Vector<Ref,int> &ipiv ) {
//
//    //    assert(A.size() == X.rows());
//
//    Matrix<General,Ref,Ref,double> Y = X;
//
//    //    if(X.rows() == 0 || X.cols() == 0)
//    //      return Y;
//
//
//    ATL_dlaswp( Y.cols(), Y(), Y.ldim(), 0, Y.rows(), ipiv(), 1 );
//
//    SquareMatrix<Ref,double> B = A;
//
//    dtrsm(CblasColMajor,side,uplo,B.blasTrans(),unit,Y.rows(), Y.cols(), 1, A(),A.ldim(), Y(), Y.ldim());
//
//    return Y;  
//  }
}
