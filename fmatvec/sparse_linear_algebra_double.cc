/* Copyright (C) 2003-2022  Martin Förg

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
#include "sparse_linear_algebra_double.h"
#include "linear_algebra.h"
#include "spooles.h"
#include "arpack.h"
//#include <stdexcept>
//#include <sstream>

namespace fmatvec {

  Vector<Ref, double> slvLU(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Vector<Ref, double> &x) {
    assert(A.size() == x.size());
    Spooles spooles(A);
    spooles.factorize();
    return spooles.solve(x);
  }

  Matrix<General, Ref, Ref, double> slvLU(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &X) {
    assert(A.size() == X.rows());
    Spooles spooles(A);
    spooles.factorize();
    return spooles.solve(X);
  }

  int eigvec(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Matrix<SymmetricSparse, Ref, Ref, double> &M, int nev, double sigma, Matrix<General, Ref, Ref, double> &eigenvectors, Vector<Ref, double> &eigenvalues) {
    const char *which = "LM";
    int ido = 0;
    char bmat = 'G';
    int n = A.rows();
    int ldv=n;
    double tol = 0;
    double* resid = new double[n];
    for(int i=0; i<n; i++)
      resid[i] = 0;
    int ncv = 2*nev;
    if (ncv > n)
      ncv = n;
    double *v = new double[ldv*ncv];
    int *iparam = new int[11];
    for(int i=0; i<11; i++)
      iparam[i] = 0;
    iparam[0] = 1;
    iparam[2] = 300;
    iparam[3] = 1;
    iparam[6] = 3;
    int *ipntr = new int[11];
    double *workd = new double[3*n];
    int lworkl = (ncv+8)*ncv;
    double *workl = new double[lworkl];
    int info = 0;

    Vector<Ref,int> ipiv(A.rows(),NONINIT);
    Spooles spooles(A,M,-sigma);
    spooles.factorize();

    while(true) {
      dsaupd_c(&ido, &bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
      if (ido==-1) {
	Vector<Ref,double> x(n,&workd[ipntr[0]-1]);
	Vector<Ref,double> y(n,&workd[ipntr[1]-1]);
	y = spooles.solve(M*x);
      }
      else if (ido==1) {
	Vector<Ref,double> t(n,&workd[ipntr[0]-1]);
	Vector<Ref,double> x(n,&workd[ipntr[1]-1]);
	Vector<Ref,double> y(n,&workd[ipntr[2]-1]);
	x = spooles.solve(y);
      }
      else if (ido==2) {
	Vector<Ref,double> x(n,&workd[ipntr[0]-1]);
	Vector<Ref,double> y(n,&workd[ipntr[1]-1]);
	y = M*x;
      }
      else
	break;
    } 

    eigenvectors.resize(n,nev,NONINIT);
    eigenvalues.resize(nev,NONINIT);
    int rvec = true;
    char howmny = 'A';
    int *select = new int[ncv];
    dseupd_c( rvec, &howmny, select, eigenvalues(), eigenvectors(), eigenvectors.ldim(), sigma, &bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);

    delete [] resid;
    delete [] v;
    delete [] iparam;
    delete [] ipntr;
    delete [] workd;
    delete [] workl;
    delete [] select;

    return info;
  }

}
