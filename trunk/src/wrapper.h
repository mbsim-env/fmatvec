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

namespace fmatvec {

  int clapack_dgels( CBLAS_TRANSPOSE ctr, int m, int n, int nhrs, double* a, int lda, double* b, int ldb);

  int clapack_dgelss(int m, int n, int nrhs, double *a, int lda, double *b, int ldb, double rcond);

  int clapack_dgeev(char jobvl, char jobvr, int n, double *a, int lda, double *wr, double *wi, double *vl, int ldvl, double* vr, int ldvr);

  int clapack_dsyev( char jobz, char ul, int n, double *a, int lda, double *w);

  int clapack_dsyevx(char jobz, char range, CBLAS_UPLO cuplo, int n, double *a, int lda, double vl, double vu, int il, int iu, double abstol, int &m, double *w, double *z, int ldz);

  double clapack_dlange(char norm , int m, int n, const double* a, int lda);
}

#endif
