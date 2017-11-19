/* Copyright (C) 2007  Martin Förg, Jan Clauberg

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

#ifndef wrapper_h
#define wrapper_h

#ifndef HAVE_LIBMKL_INTEL_LP64
#ifndef CBLAS_ENUM_DEFINED_H
   #define CBLAS_ENUM_DEFINED_H
   enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
   enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
   enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
   enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
   enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};
#endif
#define CBLAS_INDEX int

#ifndef ATLAS_ENUM_H
   #define ATLAS_ENUM_H
   #define ATLAS_ORDER CBLAS_ORDER
      #define AtlasRowMajor CblasRowMajor
      #define AtlasColMajor CblasColMajor
   #define ATLAS_TRANS CBLAS_TRANSPOSE
      #define AtlasNoTrans CblasNoTrans
      #define AtlasTrans CblasTrans
      #define AtlasConjTrans CblasConjTrans
   #define ATLAS_UPLO CBLAS_UPLO
      #define AtlasUpper CblasUpper
      #define AtlasLower CblasLower
   #define ATLAS_DIAG CBLAS_DIAG
      #define AtlasNonUnit CblasNonUnit
      #define AtlasUnit CblasUnit
   #define ATLAS_SIDE CBLAS_SIDE
      #define AtlasLeft  CblasLeft
      #define AtlasRight CblasRight
#endif
#else
#include "mkl_cblas.h"
#endif

#include <complex>

typedef std::complex<double> doublecomplex;

#if defined(HAVE_LIBBLAS)
namespace fmatvec {
  void dscal(int N, double alpha, double *X, int incX) ;
  void dgemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, int M, int N,
                 int K, double alpha, const double *A,
                 int lda, const double *B, int ldb,
                 double beta, double *C, int ldc);
  void dsymm(CBLAS_ORDER Order, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, int M, int N,
                 double alpha, const double *A, int lda,
                 const double *B, int ldb, double beta,
                 double *C, int ldc);
  void daxpy(int N, double alpha, const double *X,
                 int incX, double *Y, int incY);
  void dcopy(int N, const double *X, int incX,
                 double *Y, int incY);
  void dgemv(CBLAS_ORDER Order,
                 CBLAS_TRANSPOSE TransA, int M, int N,
                 double alpha, const double *A, int lda,
                 const double *X, int incX, double beta,
                 double *Y, int incY);
  void dsymv(CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                 int N, double alpha, const double *A,
                 int lda, const double *X, int incX,
                 double beta, double *Y, int incY);
  double ddot(int N, const double *X, int incX,
                  const double *Y, int incY);
  double dasum(int N, const double *X, int incX);
  CBLAS_INDEX idamax(int N, const double *X, int incX);
  double dnrm2(int N, const double *X, int incX);
  int dgesv(CBLAS_ORDER Order, int N, int NRHS,
                  double *A, int lda, int *ipiv,
                  double *B, int ldb);
  int zgesv(CBLAS_ORDER Order, int N, int NRHS,
                  doublecomplex *A, int lda, int *ipiv,
                  doublecomplex *B, int ldb);
  int dgetrs(CBLAS_ORDER Order, CBLAS_TRANSPOSE Trans,
    int N, int NRHS, const double *A, int lda,
    const int *ipiv, double *B, int ldb);
  int dgetrf(CBLAS_ORDER Order, int M, int N,
                   double *A, int lda, int *ipiv);
  int dgetri(CBLAS_ORDER Order, int N, double *A,
                   int lda, const int *ipiv);
  int dpotrf(ATLAS_ORDER Order, ATLAS_UPLO Uplo, int N, double *A, int lda);
  int dpotri(ATLAS_ORDER Order, ATLAS_UPLO Uplo,
                   int N, double *A, int lda);
  int dposv(ATLAS_ORDER Order, ATLAS_UPLO Uplo,
                  int N, int NRHS, double *A, int lda,
                  double *B, int ldb);
  int dpotrs(CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                   int N, int NRHS, const double *A, int lda,
                   double *B, int ldb);
}
#endif

#if defined(HAVE_LIBATLAS)
namespace fmatvec {

extern "C" {
#include "cblas.h"
#include "clapack.h"
}

#define dscal  cblas_dscal
#define dcopy  cblas_dcopy
#define ddot   cblas_ddot
#define dnrm2  cblas_dnrm2
#define dasum  cblas_dasum
#define idamax cblas_idamax
#define daxpy  cblas_daxpy
#define dgemv  cblas_dgemv
#define dsymv  cblas_dsymv
#define dgemm  cblas_dgemm
#define dsymm  cblas_dsymm

#define dgesv  clapack_dgesv
#define zgesv  clapack_zgesv
#define dgetrs clapack_dgetrs
#define dgetrf clapack_dgetrf
#define dgetri clapack_dgetri
#define dpotrf clapack_dpotrf
#define dpotri clapack_dpotri
#define dpotrs clapack_dpotrs
#define dposv  clapack_dposv

}
#endif

#if defined(HAVE_LIBMKL_INTEL_LP64)
extern "C" {
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#include "mkl_lapack.h"
}

#define dscal  cblas_dscal
#define dcopy  cblas_dcopy
#define ddot   cblas_ddot
#define dnrm2  cblas_dnrm2
#define dasum  cblas_dasum
#define idamax cblas_idamax
#define daxpy  cblas_daxpy
#define dgemv  cblas_dgemv
#define dsymv  cblas_dsymv
#define dgemm  cblas_dgemm
#define dsymm  cblas_dsymm

#define dgesv  LAPACKE_dgesv
#define zgesv  LAPACKE_zgesv
#define dgetrs LAPACKE_dgetrs
#define dgetrf LAPACKE_dgetrf
#define dgetri LAPACKE_dgetri
#define dpotrf LAPACKE_dpotrf
#define dpotri LAPACKE_dpotri
#define dpotrs LAPACKE_dpotrs
#define dposv  LAPACKE_dposv

#endif

namespace fmatvec {

#ifndef HAVE_LIBMKL_INTEL_LP64
  int dsyevx(char jobz, char range, CBLAS_UPLO cuplo, int n, double *a, int lda, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz);
#else
   int dsyevx(const char jobz, const char range, const CBLAS_UPLO cuplo, const int n, double *a, const int lda, const double vl, const double vu, const int il, const int iu, const double abstol, int *m, double *w, double *z, const int ldz); 
#endif

  int dgels( CBLAS_TRANSPOSE ctr, int m, int n, int nrhs, double* a, int lda, double* b, int ldb);

  int dgelss(int m, int n, int nrhs, double *a, int lda, double *b, int ldb, double rcond);

  int dgeev(char jobvl, char jobvr, int n, double *a, int lda, double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr);

  int dsygv(int itype, char jobz, char uplo, int n, double *a, int lda, double *b, int ldb, double *w);
  
  int dgesvd(char jobu, char jobvt, int m, int n, double *a, int lda, double *s, double *u, int ldu, double *vt, int ldvt);

  int dsyev(char jobz, char ul, int n, double *a, int lda, double *w);

  double dlange(char norm, int m, int n, const double* a, int lda);
  
  double dlansy(char norm, char uplo, int n, const double* a, int lda);
}
#endif
