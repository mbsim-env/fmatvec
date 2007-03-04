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
 *   martin.foerg@web.de
 *
 */

#ifndef blas_extensions_double_h
#define blas_extensions_double_h

#ifndef CBLAS_ENUM_DEFINED_H
   #define CBLAS_ENUM_DEFINED_H
   enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
   enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
   enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
   enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
   enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};
#endif

void daxpy(const enum CBLAS_TRANSPOSE  transa, const enum CBLAS_TRANSPOSE transb, int m, int n, double alpha, const double *a, int lda1, double *b, int ldb1);

void dcopy(const enum CBLAS_TRANSPOSE  transa, const enum CBLAS_TRANSPOSE transb, int m, int n, const double *a, int lda1, double *b, int ldb1);

void dscal(int n, double da, const double *dx, int incx, double *c);

void dscal(const enum CBLAS_TRANSPOSE  transa, int m, int n, double alpha, 
	   const double *a, int lda1, double *c, int ldc);

void dscal(char transa, int m, int n, double alpha, double *a, int lda1);

void daxpy(int n, double da, const double *dx, int incx, const double *dy, int incy, double *c);

void daxpy(const enum CBLAS_TRANSPOSE transa, const enum CBLAS_TRANSPOSE transb, int m, int n, double alpha, const double *a, int lda1, const double *b, int ldb1, double *c, int ldc);

#endif
