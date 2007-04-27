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

extern "C" {
#include "cblas.h"
  void dgeev_( char* jobvl, char* jobvr, int *n, double *a, int *lda,double * wr,double * wi,double * vl, int *ldvl,double* vr,int * ldvr,double* work,int * lwork, int *info );
  void dsyev_( char* jobz, char*, int *n, double *a, int *lda,double * w, double* work,int * lwork, int *info );
  void dgels_( char*, int*, int*, int*,double* a, int*,double* b, int*, double*, int*, int*);
  double dlange_( char *norm , int *m, int *n, const double* A, int* lda , double* work );
  void dsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int* m, double *w, double *z, int *ldz, double *work, int *lwork, int *iwork, int *ifail, int *info );
  void dgelss_( int *m, int *n, int *nrhs, double *A, int *lda, double *b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *info );
}

namespace fmatvec {

#define CVT_TRANSPOSE(c) \
   (((c) == CblasNoTrans) ? 'N' : \
    ((c) == CblasTrans) ? 'T' : \
    ((c) == CblasConjTrans) ? 'C' : \
    -1)

#define CVT_UPLO(c) \
   (((c) == CblasUpper) ? 'U' : \
    ((c) == CblasLower) ? 'L' : \
    -1)

  int clapack_dgels( CBLAS_TRANSPOSE ctr, int m, int n, int nrhs, double* a, int lda, double* b, int ldb) {

    char tr=CVT_TRANSPOSE(ctr);
    int lwork =  m*10;
    double *work = new double[lwork];
    int info;
    dgels_( &tr, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info );

    delete [] work;

    return info;
  }

  int clapack_dgelss(int m, int n, int nrhs, double *a, int lda, double *b, int ldb, double rcond) {

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

  int clapack_dgeev(char jobvl, char jobvr, int n, double *a, int lda, double *wr, double *wi, double *vl, int ldvl, double* vr, int ldvr) {

    int info;
    int lwork = 10*n;
    double *work = new double[lwork];

    dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &n, vr, &n, work, &lwork, &info);

    delete [] work;

    return info;

  }

  int clapack_dsyev( char jobz, char ul, int n, double *a, int lda, double *w) {

    int lwork = 10*n;
    double *work = new double[lwork];
    int info;

    dsyev_(&jobz, &ul, &n, a, &lda, w, work, &lwork, &info);

    delete [] work;

    return info;

  }

  int clapack_dsyevx(char jobz, char range, CBLAS_UPLO cuplo, int n, double *a, int lda, double vl, double vu, int il, int iu, double abstol, int &m, double *w, double *z, int ldz) {

    int lwork = 8*n;  // opt(NB+3)*N
    double *work = new double[lwork];
    int *iwork = new int[5*n];
    int info;
    int *ifail = new int[n];
    char uplo = CVT_UPLO(cuplo);

    dsyevx_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info );

    delete [] work;
    delete [] iwork;
    delete [] ifail;

    return info; 

  }
  
  double clapack_dlange(char norm , int m, int n, const double* a, int lda) {

    double *work = new double[2*m];
    double res = dlange_(&norm , &m, &n, a, &lda , work );
    delete [] work;
    return res;

  }


}
