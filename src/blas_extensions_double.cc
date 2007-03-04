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

#include "config.h"
#include "blas_extensions_double.h"

extern "C" {
#include "cblas.h"
}

void daxpy(const enum CBLAS_TRANSPOSE  transa, const enum CBLAS_TRANSPOSE transb, int m, int n, double alpha, const double *a, int lda1, double *b, int ldb1) {

  int i, nn, one = 1;
  bool nota = (transa == CblasNoTrans) ? true : false;
  bool notb = (transb == CblasNoTrans) ? true : false;

  if(nota && notb && lda1==m && ldb1==m) {
    nn = (m)*(n);
    cblas_daxpy(nn, alpha, a, one, b, one);
  } else {
    if(nota) {
      if(notb) {
	for(i = 0; i < n; i++)
	  cblas_daxpy(m, alpha, a+i*(lda1), one, b+i*(ldb1), one);
      } else {
	for(i = 0; i < n; i++)
	  cblas_daxpy(m, alpha, a+i*(lda1), one, b+i, ldb1);
      }
    } else {
      if(notb) {
	for(i = 0; i < n; i++)
	  cblas_daxpy(m, alpha, a+i, lda1, b+i*(ldb1), one);
      } else {
	for(i = 0; i < n; i++)
	  cblas_daxpy(m, alpha, a+i, lda1, b+i, ldb1);
      }
    }
  }
}

void dcopy(const enum CBLAS_TRANSPOSE  transa, const enum CBLAS_TRANSPOSE transb, int m, int n, const double *a, int lda1, double *b, int ldb1) {

  bool nota = (transa == CblasNoTrans) ? true : false;
  bool notb = (transb == CblasNoTrans) ? true : false;
  int i, one = 1, nn;
  
  if(nota && notb && lda1==m && ldb1==m) {
    nn = m*n;
    cblas_dcopy(nn,  a, one, b, one);
    return;
  }
  if(nota) {
      if(notb) {
	for(i = 0; i < n; i++)
	  cblas_dcopy(m,  a+i*lda1, one,  b+i*ldb1, one);
	return;
      } else {
	for(i = 0; i < n; i++)
	  cblas_dcopy(m,  a+i*lda1, one, b+i, ldb1);
	return;
      }
  } else {
    if(notb) {
      for(i = 0; i < n; i++)
	cblas_dcopy(m,  a+i, lda1, b+i*ldb1, one);
      return;
    } else {
      for(i = 0; i < n; i++)
	cblas_dcopy(m,  a+i, lda1, b+i, ldb1);
      return;
    }
  }
}

void dscal(int n, double da, const double *dx, int incx, double *c) {

    int i, m;

    if (incx == 1) {
      m = n % 4;
      if (m != 0) {
	for (i = 0; i < m; ++i) 
	  c[i]= da * dx[i];
	
	if (n < 4) 
	  return;
      }
      for (i = m; i < n; i +=4)	{
	c[i]= da * dx[i];
	c[i+1]= da * dx[i+1];
	c[i+2]= da * dx[i+2];
	c[i+3]= da * dx[i+3];
      }
    } else
      for (i = 0; i < n; i++) 
	c[i] =da * dx[i*incx];
} 

void dscal(const enum CBLAS_TRANSPOSE  transa, int m, int n, double alpha, const double *a, int lda1, double *c, int ldc) {

  int i, nn, one = 1;
  bool nota = (transa == CblasNoTrans) ? true : false;

  if(nota && lda1==m && ldc==m) {
    nn = m*n;
    dscal(nn, alpha, a, one, c);
    return;
  }
  if(nota) {
    for(i = 0; i < n; i++)
      dscal(m, alpha, a+i*lda1, one, c+i*ldc);
  } else {
    for(i = 0; i < n; i++)
      dscal(m, alpha, a+i, lda1, c+i*ldc);
  }
}

void dscal(char transa, int m, int n, double alpha, double *a, int lda1) {
    int i, nn, one = 1;
    bool nota = (transa == CblasNoTrans) ? true : false;

    if(nota && lda1==m) {
      nn = m*n;
      cblas_dscal(nn, alpha, a, one);
      return;
    }
    if(nota) {
      for(i = 0; i < n; i++)
	cblas_dscal(m, alpha, a+i*lda1, one);
    }
}

void daxpy(int n, double da, const double *dx, int incx, const double *dy, int incy, double *c)
{
  int i, m;
/* for(i=0; i<n; i++)
  *c++=*dx++ + *dy++;
 return;
for(i=0; i<n; i++)
  c[i]=dx[i] + dy[i];
return;*/
 if (incx == 1)
  {
    if(incy == 1)
    {
      m = n % 4;
      if (m != 0)
      {
	for (i = 0; i < m; ++i)
	  c[i]=dy[i] + da * dx[i];

	if (n < 4)
	  return;
      }
      for (i = m; i < n; i +=4)
      {
	c[i]=dy[i] + da * dx[i];
	c[i+1]=dy[i+1] + da * dx[i+1];
	c[i+2]=dy[i+2] + da * dx[i+2];
	c[i+3]=dy[i+3] + da * dx[i+3];
      }
    }
    else
    {
      for (i = 0; i < n; i++)
	c[i]=dy[i*incy] + da * dx[i];
    }
  }
  else
  {
    if(incy==1)
    {
      for (i = 0; i < n; i++)
	c[i]=dy[i] + da * dx[i*incx];
    }
    else
    {
      for (i = 0; i < n; i++)
	c[i]=dy[i*incy] + da * dx[i*incx];
    }
  }
}
void daxpy(const enum CBLAS_TRANSPOSE transa, const enum CBLAS_TRANSPOSE transb, int m, int n, double alpha, const double *a, int lda1, const double *b, int ldb1, double *c, int ldc)
{
  int i, nn, one = 1;
  bool nota = (transa == CblasNoTrans) ? true : false;
  bool notb = (transb == CblasNoTrans) ? true : false;


  if(nota && notb && lda1==m && ldb1==m && ldc==m)
  {
    nn = (m)*(n);
    daxpy(nn, alpha, a, one, b, one, c);
  }
  else
  {
    if(nota)
    {
      if(notb)
      {
	for(i = 0; i < n; i++)
	  daxpy(m, alpha, a+i*(lda1), one, b+i*(ldb1), one, c+i* ldc);
      }
      else
      {
	for(i = 0; i < n; i++)
	  daxpy(m, alpha, a+i*(lda1), one, b+i, ldb1, c+i* ldc);
      }
    }
    else
    {
      if(notb)
      {
	for(i = 0; i < n; i++)
	  daxpy(m, alpha, a+i, lda1, b+i*(ldb1), one, c+i* ldc);
      }
      else
      {
	for(i = 0; i < n; i++)
	  daxpy(m, alpha, a+i, lda1, b+i, ldb1,c+i* ldc);
      }
    }
  }
}


