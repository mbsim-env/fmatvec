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

#include "config.h"
#include "wrapper.h"

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




extern "C" {

#ifndef HAVE_LIBATLAS

  void dscal_(const int *n, const double *alpha, double *X, const int* incX);
  void dgemm_(const char *TransA, const char* TransB, const int *M, const int *N,
                 const int *K, const double *alpha, const double *A,
                 const int *lda, const double *B, const int *ldb,
                 const double *beta, double *C, const int *ldc);
  void dsymm_(const char *Side,
                 const char *Uplo, const int *M, const int *N,
                 const double *alpha, const double *A, const int *lda,
                 const double *B, const int *ldb, const double *beta,
                 double *C, const int *ldc);
  void daxpy_(const int *N, const double *alpha, const double *X,
                 const int *incX, double *Y, const int *incY);
  void dcopy_(const int *N, const double *X, const int *incX,
                 double *Y, const int *incY);
  void dgemv_(const char *TransA, const int *M, const int *N,
                 const double *alpha, const double *A, const int *lda,
                 const double *X, const int *incX, const double *beta,
                 double *Y, const int *incY);
  void dsymv_(const char* Uplo, const int *N, const double *alpha, const double *A,
                 const int *lda, const double *X, const int *incX,
                 const double *beta, double *Y, const int *incY);

  double ddot_(const int *N, const double *X, const int *incX,
                  const double *Y, const int *incY);
  double dasum_(const int *N, const double *X, const int *incX);
  int idamax_(const int *N, const double *X, const int *incX);
  double dnrm2_(const int *N, const double *X, const int *incX);


  void dgesv_(const int *N, const int *NRHS,
                  double *A, const int *lda, int *ipiv,
                  double *B, const int *ldb, int *info);

  void dgetrs_(const char *Trans, const int *N, const int *NRHS, const double
      *A, const int *lda, const int *ipiv, double *B, const int *ldb, int *info);

  void dgetrf_(const int *M, const int *N, double *A, const int *lda, int *ipiv,
      int *info);

  void dgetri_(const int *N, double *A, const int *lda, const int *ipiv, double
      *work, const int *lwork, int *info);

  void dpotrf_(const char *Uplo, const int *N, double *A, const int *lda, int *info);

  void dpotri_(const char *Uplo, const int *N, double *A, const int *lda, int *info);

  void dposv_(const char *Uplo, const int *N, const int *NRHS, double *A, const
      int *lda, double *B, const int *ldb, int *info);

  void dpotrs_(const char *Uplo, const int *N, const int *NRHS, const double *A,
      const int *lda, double *B, const int *ldb, int *info);

#endif

  void dgeev_(const char* jobvl, const char* jobvr, const int *n, double *a, const int *lda, double *wr, double *wi, double *vl, const int *ldvl, double* vr, const int *ldvr, double* work, const int *lwork, int *info);
  void dsyev_( const char* jobz, const char*, const int *n, double *a, const int *lda, double *w, double *work, const int *lwork, int *info);
  void dgels_( const char*, const int*, const int*, const int*, double* a, const int*, double* b, const int*, double*, const int*, int*);
  double dlange_(const char *norm , const int *m, const int *n, const double* A, const int* lda , double* work);
  void dsyevx_(const char *jobz, const char *range, const char *uplo, const int *n, double *a, const int *lda, const double *vl, const double *vu, const int *il, const int *iu, const double *abstol, const int* m, double *w, double *z, const int *ldz, double *work, const int *lwork, int *iwork, int *ifail, int *info);
  void dgelss_(const int *m, const int *n, const int *nrhs, double *A, const int *lda, double *b, const int *ldb, const double *s, const double *rcond, int *rank, double *work, const int *lwork, int *info);
}

namespace fmatvec {

#ifndef HAVE_LIBATLAS

  void dscal(const int N, const double alpha, double *X, const int incX) {
    dscal_(&N, &alpha, X, &incX);
  }
  
  void dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {

    const char trA=CVT_TRANSPOSE(TransA);
    const char trB=CVT_TRANSPOSE(TransB);
    dgemm_(&trA, &trB, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  }

  void dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc) {


    const char side = CVT_SIDE(Side);
    const char uplo = CVT_UPLO(Uplo);

    dsymm_(&side, &uplo, &M, &N, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
  }

  void daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY) {
    daxpy_(&N, &alpha, X, &incX, Y, &incY);
  }

  void dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY) {
    dcopy_(&N, X, &incX, Y, &incY);
  }

  void dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY) {

    const char trA=CVT_TRANSPOSE(TransA);

    dgemv_(&trA, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
  }

  void dsymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY) {

    const char uplo = CVT_UPLO(Uplo);

    dsymv_(&uplo, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);


  }

  double ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY) {
    return ddot_(&N, X, &incX, Y, &incY);
  }

  double dasum(const int N, const double *X, const int incX) {
    return dasum_(&N, X, &incX);
  }

  CBLAS_INDEX idamax(const int N, const double *X, const int incX) {
    return idamax_(&N, X, &incX)-1;
  }

  double dnrm2(const int N, const double *X, const int incX) {
    return dnrm2_(&N, X, &incX);
  }



  int dgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                  double *A, const int lda, int *ipiv,
                  double *B, const int ldb) {

    int info;
    dgesv_(&N, &NRHS, A, &lda, ipiv, B, &ldb, &info);

    return info;
  }

  int dgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const double *A, const int lda,
    const int *ipiv, double *B, const int ldb) {

    const char trA=CVT_TRANSPOSE(Trans);
    int info;
    dgetrs_(&trA, &N, &NRHS, A, &lda, ipiv, B, &ldb, &info);
    return info;
  }

  int dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   double *A, const int lda, int *ipiv) {
    int info;
    dgetrf_(&M, &N, A, &lda, ipiv, &info);
    return info;
  }

  int dgetri(const enum CBLAS_ORDER Order, const int N, double *A,
                   const int lda, const int *ipiv) {
    int info;
    int lwork = 2*N;
    double *work = new double[lwork];
    dgetri_(&N, A, &lda, ipiv, work, &lwork, &info);
    delete [] work;
    return info;
  }

  int dpotrf(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo, const int N, double *A, const int lda) {

    int info;
    const char uplo = CVT_UPLO(Uplo);
    dpotrf_(&uplo, &N, A, &lda, &info);
    return info;
  }

  int dpotri(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, double *A, const int lda) {
    int info;
    const char uplo = CVT_UPLO(Uplo);
    dpotri_(&uplo, &N, A, &lda, &info);
    return info;
  }

  int dposv(const enum ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, double *A, const int lda,
                  double *B, const int ldb) {
    int info;
    const char uplo = CVT_UPLO(Uplo);
    dposv_(&uplo, &N, &NRHS, A, &lda, B, &ldb, &info);
    return info;
  }

  int dpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const double *A, const int lda,
                   double *B, const int ldb) {
    int info;
    const char uplo = CVT_UPLO(Uplo);
    dpotrs_(&uplo, &N, &NRHS, A, &lda, B, &ldb, &info);
    return info;
  }

#endif


  int dgels( const enum CBLAS_TRANSPOSE ctr, const int m, const int n, const int nrhs, double* a, const int lda, double* b, const int ldb) {

    const char tr=CVT_TRANSPOSE(ctr);
    const int lwork =  m*10;
    double *work = new double[lwork];
    int info;
    dgels_( &tr, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info );

    delete [] work;

    return info;
  }

  int dgelss(const int m, const int n, const int nrhs, double *a, const int lda, double *b, const int ldb, const double rcond) {

    int minmn = m<n?m:n;
    double *s = new double[minmn]; 
    int lwork = 2*(3*minmn + m+n);
    double *work = new double[lwork];
    int rank;
    int info;

    dgelss_( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info );
    delete [] s;
    delete [] work;
    return info;
  }

  int dgeev(const char jobvl, const char jobvr, const int n, double *a, const int lda, double *wr, double *wi, double *vl, const int ldvl, double *vr, const int ldvr) {

    int info;
    const int lwork = 10*n;
    double *work = new double[lwork];

    dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &n, vr, &n, work, &lwork, &info);

    delete [] work;

    return info;

  }

  int dsyev(const char jobz, const char ul, const int n, double *a, const int lda, double *w) {

    int lwork = 10*n;
    double *work = new double[lwork];
    int info;

    dsyev_(&jobz, &ul, &n, a, &lda, w, work, &lwork, &info);

    delete [] work;

    return info;

  }

  int dsyevx(const char jobz, const char range, const enum CBLAS_UPLO cuplo, const int n, double *a, const int lda, const double vl, const double vu, const int il, const int iu, const double abstol, const int *m, double *w, double *z, const int ldz) {

    int lwork = 8*n;  // opt(NB+3)*N
    double *work = new double[lwork];
    int *iwork = new int[5*n];
    int info;
    int *ifail = new int[n];
    char uplo = CVT_UPLO(cuplo);

    dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, iwork, ifail, &info );

    delete [] work;
    delete [] iwork;
    delete [] ifail;

    return info; 

  }
  
  double dlange(const char norm , const int m, const int n, const double* a, const int lda) {

    double *work = new double[2*m];
    double res = dlange_(&norm , &m, &n, a, &lda , work );
    delete [] work;
    return res;

  }


}
