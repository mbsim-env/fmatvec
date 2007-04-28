/* Copyright (C) 2007  Martin Förg

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
 *   mfoerg@users.berlios.de
 *
 */

#ifndef wrapper_h
#define wrapper_h

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

namespace fmatvec {

#define HAVE_ATLAS

#ifndef HAVE_ATLAS
  void dscal(const int N, const double alpha, double *X, const int incX) ;
  void dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
  void dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc);
  void daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);
  void dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);
  void dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);
  void dsymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
  double ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
  double dasum(const int N, const double *X, const int incX);
  CBLAS_INDEX idamax(const int N, const double *X, const int incX);
  double dnrm2(const int N, const double *X, const int incX);
  int dgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  double *A, const int lda, int *ipiv,
                  double *B, const int ldb);
  int dgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const double *A, const int lda,
    const int *ipiv, double *B, const int ldb);
  int dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   double *A, const int lda, int *ipiv);
  int dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
                   const int lda, const int *ipiv);
  int dpotrf(const enum ATLAS_UPLO Uplo, const int N, double *A, const int lda);
  int dpotri(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, double *A, const int lda);
  int dposv(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, double *A, const int lda,
                  double *B, const int ldb);
  int dpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const double *A, const int lda,
                   double *B, const int ldb);
#else

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
#define dgetrs clapack_dgetrs
#define dgetrf clapack_dgetrf
#define dgetri clapack_dgetri
#define dpotrf clapack_dpotrf
#define dpotri clapack_dpotri
#define dpotrs clapack_dpotrs
#define dposv  clapack_dposv

#endif

  int dgels( const enum CBLAS_TRANSPOSE ctr, const int m, const int n, const int nhrs, double* a, const int lda, double* b, const int ldb);

  int dgelss(const int m, const int n, const int nrhs, double *a, const int lda, double *b, const int ldb, const double rcond);

  int dgeev(const char jobvl, const char jobvr, const int n, double *a, const int lda, double *wr, double *wi, double *vl, const int ldvl, double *vr, const int ldvr);

  int dsyev(const char jobz, const char ul, const int n, double *a, const int lda, double *w);

  int dsyevx(const char jobz, const char range, const enum CBLAS_UPLO cuplo, const int n, double *a, const int lda, const double vl, const double vu, const int il, const int iu, const double abstol, const int *m, double *w, double *z, const int ldz);

  double dlange(const char norm , const int m, const int n, const double* a, const int lda);
}

#endif
