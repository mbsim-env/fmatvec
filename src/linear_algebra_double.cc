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
 *   mfoerg@users.berlios.de
 *
 */

#include "config.h"
#include "linear_algebra.h"
#include "linear_algebra_double.h"
#include "blas_extensions_double.h"
#include "wrapper.h"

#define FMATVEC_NO_INITIALIZATION
#define FMATVEC_NO_BOUNDS_CHECK

//-------------------------------------
// Matrix operations
//-------------------------------------
namespace fmatvec {

  Matrix<General, double> operator*(const Matrix<General, double> &A, double alpha) {

    Matrix<General, double> B(A.rows(),A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),alpha,A(),A.ldim(),B(),
	B.ldim());

    return B;
  }

  Matrix<General, double> operator/(const Matrix<General, double> &A, double alpha) {

    Matrix<General, double> B(A.rows(),A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),1./alpha,A(),A.ldim(),
	B(),B.ldim());

    return B;
  }

  Matrix<General, double> operator*(double alpha, const Matrix<General, double> &A) {

    Matrix<General, double> B(A.rows(),A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),alpha,A(),A.ldim(),B(),
	B.ldim());

    return B;
  }

  SquareMatrix<double> operator*(const SquareMatrix<double> &A, double alpha) {

    SquareMatrix<double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),alpha,A(),A.ldim(),B(),
	B.ldim());

    return B;
  }

  SquareMatrix<double> operator/(const SquareMatrix<double> &A, double alpha) {

    SquareMatrix<double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),1./alpha,A(),A.ldim(),
	B(), B.ldim());

    return B;
  }

  SquareMatrix<double> operator*(double alpha, const SquareMatrix<double> &A) {

    SquareMatrix<double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),alpha,A(),A.ldim(),B(),
	B.ldim());

    return B;
  }

  Matrix<Symmetric, double> operator*(const Matrix<Symmetric, double> &A, double alpha) {

    Matrix<Symmetric, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    for(int i=0; i<A.size(); i++) 
      for(int j=0; j<=i; j++) 
	B(i,j) = A(i,j)*alpha;

    return B;
  }

  Matrix<Symmetric, double> operator/(const Matrix<Symmetric, double> &A, double alpha) {

    Matrix<Symmetric, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    for(int i=0; i<A.size(); i++) 
      for(int j=0; j<=i; j++) 
	B(i,j) = A(i,j)/alpha;

    return B;
  }

  Matrix<Symmetric, double> operator*(double alpha, const Matrix<Symmetric, double> &A) {

    Matrix<Symmetric, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    for(int i=0; i<A.size(); i++) 
      for(int j=0; j<=i; j++) 
	B(i,j) = A(i,j)*alpha;

    return B;
  }

  Matrix<Diagonal, double> operator*(double alpha, const Matrix<Diagonal, double> &A) {

    Matrix<Diagonal, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.size(),alpha,A(),1,B());

    return B;

  }

  Matrix<Diagonal, double> operator*(const Matrix<Diagonal, double> &A, double alpha) {

    Matrix<Diagonal, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.size(),alpha,A(),1,B());

    return B;

  }

  Matrix<Diagonal, double> operator/(const Matrix<Diagonal, double> &A, double alpha) {

    Matrix<Diagonal, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return B;
    }
#endif

    myblas_dscal(A.size(),1./alpha,A(),1,B());

    return B;

  }

  Matrix<Diagonal, double>& operator/=(const Matrix<Diagonal, double >& A_, double alpha) {

    Matrix<Diagonal, double> &A = const_cast<Matrix<Diagonal, double> &>(A_);

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return A;
    }
#endif

    dscal(A.size(), 1./alpha, A(),1);

    return A;
  }

  Matrix<Diagonal, double>& operator*=(const Matrix<Diagonal, double >& A_, double alpha) {

    Matrix<Diagonal, double> &A = const_cast<Matrix<Diagonal, double> &>(A_);

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return A;
    }
#endif

    dscal(A.size(), alpha, A(),1);

    return A;
  }

  SquareMatrix<double>& operator*=(const SquareMatrix<double>& A_, double alpha) {

    SquareMatrix<double> &A = const_cast<SquareMatrix<double> &>(A_);

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return A;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),alpha,A(),A.ldim());

    return A;
  }

  SquareMatrix<double>& operator/=(const SquareMatrix<double>& A_, double alpha) {

    SquareMatrix<double> &A = const_cast<SquareMatrix<double> &>(A_);

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return A;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),1./alpha,A(),A.ldim());

    return A;
  }

  Matrix<General, double>& operator*=(const Matrix<General, double > &A_, const double &alpha) {

    Matrix<General, double> &A = const_cast<Matrix<General, double> &>(A_);

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return A;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),alpha,A(),A.ldim());

    return A;
  }

  Matrix<General, double>& operator/=(const Matrix<General, double > &A_, const double &alpha) {

    Matrix<General, double> &A = const_cast<Matrix<General, double> &>(A_);

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return A;
    }
#endif

    myblas_dscal(A.blasTrans(),A.rows(),A.cols(),1./alpha,A(),A.ldim());

    return A;
  }

  //-------------------------------------
  // Matrix/Matrix operations
  //-------------------------------------

  Matrix<Symmetric, double > operator+(const Matrix<Symmetric, double > &A, const Matrix<Symmetric, double > &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == B.size());
#endif

    Matrix<Symmetric, double > C(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return C;
    }
#endif

    for(int i=0; i<A.size(); i++)
      for(int j=0; j<=i; j++)
	C(i,j) = A(i,j) + B(i,j);
    return C;
  }

  SquareMatrix<double> operator+(const SquareMatrix<double> &A, const SquareMatrix<double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == B.size());
#endif

    SquareMatrix<double> C(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return C;
    }
#endif

    myblas_dcopy(A.blasTrans(), C.blasTrans(), A.rows(), A.cols(), A(), A.ldim(),
	C(),C.ldim());
    myblas_daxpy(B.blasTrans(), C.blasTrans(), B.rows(), B.cols(),1., B(),
	B.ldim(), C(),C.ldim());

    return C;
  }

  SquareMatrix<double> operator-(const SquareMatrix<double> &A, const SquareMatrix<double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == B.size());
#endif

    SquareMatrix<double> C(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return C;
    }
#endif

    myblas_dcopy(A.blasTrans(), C.blasTrans(), A.rows(), A.cols(), A(), A.ldim(),
	C(),C.ldim());
    myblas_daxpy(B.blasTrans(), C.blasTrans(), B.rows(), B.cols(),-1., B(),
	B.ldim(), C(),C.ldim());

    return C;
  }

   Matrix<Symmetric, double > operator-(const Matrix<Symmetric, double > &A, const Matrix<Symmetric, double > &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == B.size());
#endif

    Matrix<Symmetric, double > C(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return C;
    }
#endif

    for(int i=0; i<A.size(); i++)
      for(int j=0; j<=i; j++)
	C(i,j) = A(i,j) - B(i,j);
    return C;
  }

  SquareMatrix<double> operator*(const SquareMatrix<double> &A, const SquareMatrix<double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == B.size());
#endif

    SquareMatrix<double> C(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0) {
      return C;
    }
#endif

    dgemm(A.blasOrder(),A.blasTrans(),B.blasTrans(),A.rows(),B.cols(),
	A.cols(),1.,A(),A.ldim(),B(),B.ldim(),0.,
	C(),C.ldim());

    return C;
  }

  Matrix<General, double> operator*(const Matrix<General, double> &A,const Matrix<Diagonal, double> &B) { 

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.cols() == B.rows());
#endif

    Matrix<General, double> C(A.rows(),B.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return C;
    }
#endif

    for(int m=0; m<C.rows(); m++)
      for(int n=0; n<C.cols(); n++)
	C(m,n)=A(m,n)*B(n);

    return C;
  }

  Matrix<General, double> operator*(const Matrix<Diagonal, double> &A, const Matrix<General, double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.cols() == B.rows());
#endif

    Matrix<General, double> C(A.rows(),B.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return C;
    }
#endif

    for(int m=0; m<C.rows(); m++)
      for(int n=0; n<C.cols(); n++)
	C(m,n)=A(m)*B(m,n);

    return C;
  }


  Matrix<General, double>& operator+=(const Matrix<General, double > &A_, const Matrix<General, double> &B) {

    Matrix<General, double> &A = const_cast<Matrix<General, double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return A;
    }
#endif

    myblas_daxpy(B.blasTrans(), A.blasTrans(), A.rows(), A.cols(),1., B(),
	B.ldim(), A(), A.ldim());

    return A;
  }

  Matrix<General, double>& operator-=(const Matrix<General, double > &A_, const Matrix<General, double> &B) {

    Matrix<General, double> &A = const_cast<Matrix<General, double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return A;
    }
#endif

    myblas_daxpy(B.blasTrans(), A.blasTrans(), A.rows(), A.cols(),-1., B(),
	B.ldim(), A(), A.ldim());

    return A;
  }

  Matrix<General, double > operator+(const Matrix<General, double > &A, const Matrix<General, double > &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols()); 
#endif

    Matrix<General, double > C(A.rows(),A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return C;
    }
#endif

    myblas_daxpy(A.blasTrans(), B.blasTrans(), A.rows(), A.cols(),1., A(), A.ldim(), B(),B.ldim(),C(),C.ldim());
    return C;
  }

  Matrix<General, double> operator-(const Matrix<General, double> &A, const Matrix<General, double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());
#endif

    Matrix<General, double> C(A.rows(),A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0) {
      return C;
    }
#endif

    myblas_dcopy(A.blasTrans(), C.blasTrans(), A.rows(), A.cols(), A(), A.ldim(),
	C(),C.ldim());
    myblas_daxpy(B.blasTrans(), C.blasTrans(), B.rows(), B.cols(),-1., B(),
	B.ldim(), C(),C.ldim());

    return C;
  }

  Matrix<General, double> operator*(const Matrix<General, double> &A, 
      const Matrix<General, double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.cols() == B.rows());
#endif

    Matrix<General, double> C(A.rows(),B.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || B.cols() == 0)
      return C;
    else if(A.cols() == 0) {
      C.init(0);
      return C;
    }
#endif

    dgemm(A.blasOrder(), A.blasTrans(), B.blasTrans(), A.rows() ,B.cols(), A.cols(), 1., A(), A.ldim(), B(), B.ldim(), 0., C(),C.ldim());

    return C;
  }

  Matrix<General, double> operator*(const Matrix<General, double> &A, const Matrix<Symmetric, double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.cols() == B.rows());
#endif

    int m, n;
    CBLAS_SIDE side;
    if(A.transposed()) {
      n = A.rows();
      m = B.cols();
      side = CblasLeft;
    } else {
      m = A.rows();
      n = B.cols();
      side = CblasRight;
    }

    Matrix<General, double> C(m,n);

#ifdef FMATVEC_VOID_CHECK
    if(m == 0 || n == 0)
      return A.transposed() ? trans(C) : C;
#endif

    dsymm(A.blasOrder(), side, B.blasUplo(), C.rows() , C.cols(), 1., B(), B.ldim(), A(), A.ldim(), 0., C(), C.ldim());
    return A.transposed() ? trans(C) : C;
  }

  Matrix<General, double> operator*(const Matrix<Symmetric, double> &A, const Matrix<General, double> &B) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.cols() == B.rows());
#endif

    int m, n;
    CBLAS_SIDE side;
    if(B.transposed()) {
      n = A.rows();
      m = B.cols();
      side = CblasRight;
    } else {
      m = A.rows();
      n = B.cols();
      side = CblasLeft;
    }
    Matrix<General, double> C(m,n);

#ifdef FMATVEC_VOID_CHECK
    if(m == 0 || n == 0)
      return B.transposed() ? trans(C) : C;
#endif

    dsymm(A.blasOrder(), side, A.blasUplo(), C.rows() , C.cols(), 1., A(), A.ldim(), B(), B.ldim(), 0., C(), C.ldim());
    return B.transposed() ? trans(C) : C;
  }

  Matrix<Symmetric, double> JTJ(const Matrix<General, double> &J) {

    Matrix<Symmetric, double> R(J.cols());

#ifdef FMATVEC_VOID_CHECK
    if(J.cols() == 0)
      return R;
    else if(J.rows() == 0) {
      R.init(0);
      return R;
    }
#endif

    CBLAS_TRANSPOSE jtt = J.transposed() ? CblasNoTrans : CblasTrans;

    dgemm(J.blasOrder(), jtt, J.blasTrans(), J.cols() ,J.cols(), J.rows(), 1., J(), J.ldim(), J(), J.ldim(), 0., R(),R.ldim());

    return R;
  }

  Matrix<Symmetric, double> JTMJ(const Matrix<Symmetric, double> &M, const Matrix<General, double> &J) {

#ifdef FMMTVEC_SIZE_CHECK 
    assert(M.cols() == J.rows());
#endif

    int m, n;
    CBLAS_SIDE side;
    CBLAS_TRANSPOSE ct,jtt;
    if(J.transposed()) {
      n = M.rows();
      m = J.cols();
      ct = CblasTrans;
      jtt = CblasNoTrans;
      side = CblasRight;
    } else {
      m = M.rows();
      n = J.cols();
      ct = CblasNoTrans;
      jtt = CblasTrans;
      side = CblasLeft;
    }

    //  Matrix<General, double> C(m,n);
    //

    Matrix<Symmetric, double> R(J.cols());

#ifdef FMATVEC_VOID_CHECK
    if(J.cols() == 0)
      return R;
    else if(J.rows() == 0) {
      R.init(0);
      return R;
    }
#endif

    double *C = new double[m*n];

    dsymm(M.blasOrder(), side, M.blasUplo(), m , n, 1., M(), M.ldim(), J(), J.ldim(), 0., C, m);

    dgemm(J.blasOrder(), jtt , ct, J.cols() ,J.cols(), J.rows(), 1., J(), J.ldim(), C, m, 0., R(),R.ldim());

    delete [] C;

    return R;
  }

  Matrix<Symmetric, double> JTMJ(const Matrix<Diagonal, double> &M, const Matrix<General, double> &J) {

#ifdef FMMTVEC_SIZE_CHECK 
    assert(M.size() == J.rows());
#endif

    Matrix<Symmetric, double> R(J.cols());

#ifdef FMATVEC_VOID_CHECK
    if(J.cols() == 0)
      return R;
    else if(J.rows() == 0) {
      R.init(0);
      return R;
    }
#endif

    CBLAS_TRANSPOSE jtt = J.transposed() ? CblasNoTrans : CblasTrans;
    Matrix<General, double> MJ = M * J;

    dgemm(J.blasOrder(), jtt, J.blasTrans(), J.cols() ,MJ.cols(), J.rows(), 1., J(), J.ldim(), MJ(), MJ.ldim(), 0., R(),R.ldim());

    return R;

  }

  Matrix<Diagonal, double>& operator+=(const Matrix<Diagonal, double >& A_, const Matrix<Diagonal, double> &B) {

    Matrix<Diagonal, double> &A = const_cast<Matrix<Diagonal, double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.size() == B.size()); 
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return A;
#endif

    daxpy(B.size(), 1., B(), 1, A(), 1);

    return A;
  }

  Matrix<Diagonal, double>& operator-=(const Matrix<Diagonal, double >& A_, const Matrix<Diagonal, double> &B) {

    Matrix<Diagonal, double> &A = const_cast<Matrix<Diagonal, double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.size() == B.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return A;
#endif

    daxpy(B.size(), -1., B(), 1, A(),1);

    return A;
  }


  SquareMatrix<double>& operator+=(const SquareMatrix<double>& A_, const SquareMatrix<double> &B) {

    SquareMatrix<double> &A = const_cast<SquareMatrix<double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.size() == B.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return A;
#endif

    myblas_daxpy(B.blasTrans(), A.blasTrans(), A.rows(), A.cols(), 1., B(), B.ldim(), A(), A.ldim());

    return A;
  }

  Matrix<Symmetric, double>& operator+=(const Matrix<Symmetric, double >& A_, const Matrix<Symmetric, double> &B) {

    Matrix<Symmetric, double> &A = const_cast<Matrix<Symmetric, double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.size() == B.size()); 
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return A;
#endif

    // TODO Blas FUnktion suchen!
    for(int i=0; i<A.size(); i++)
      for(int j=0; j<=i; j++)
	A(i,j) += B(i,j);

    return A;
  }

  Matrix<Symmetric, double>& operator-=(const Matrix<Symmetric, double >& A_, const Matrix<Symmetric, double> &B) {

    Matrix<Symmetric, double> &A = const_cast<Matrix<Symmetric, double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.size() == B.size()); 
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return A;
#endif

    // TODO Blas FUnktion suchen!
    for(int i=0; i<A.size(); i++)
      for(int j=0; j<=i; j++)
	A(i,j) -= B(i,j);

    return A;
  }

  SquareMatrix<double>& operator-=(const SquareMatrix<double>& A_, const SquareMatrix<double> &B) {

    SquareMatrix<double> &A = const_cast<SquareMatrix<double> &>(A_);

#ifdef FMATVEC_SIZE_CHECK
    assert(A.size() == B.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return A;
#endif

    myblas_daxpy(B.blasTrans(), A.blasTrans(), A.rows(), A.cols(), -1., B(), B.ldim(), A(), A.ldim());

    return A;
  }

  //-------------------------------------
  // Vector operations
  //-------------------------------------

  Vector<double>& operator/=(const Vector<double> &x_, const double &alpha) {

    Vector<double> &x = const_cast<Vector<double> &>(x_);

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    dscal(x.size(), 1./alpha, x(), x.inc());

    return x;
  }

  Vector<double>& operator*=(const Vector<double> &x_, const double &alpha) {

    Vector<double> &x = const_cast<Vector<double> &>(x_);

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    dscal(x.size(), alpha, x(), x.inc());

    return x;
  }

  Vector<double> operator*(const Vector<double> &x, double alpha) {

    Vector<double> y(x.size());//=x.copy();
    //Vector<double> y=x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    myblas_dscal(x.size(),alpha,x(),x.inc(),y());

    return y;
  }

  Vector<double> operator*(double alpha, const Vector<double> &x) {

    Vector<double> y(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    myblas_dscal(x.size(),alpha,x(),x.inc(),y());

    return y;
  }

  Vector<double> operator/(const Vector<double> &x, double alpha) {

    Vector<double> y(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    myblas_dscal(x.size(),1./alpha,x(),x.inc(),y());

    return y;
  }

  //-------------------------------------
  // Vector/Vector operations
  //-------------------------------------


  Vector<double>& operator+=(const Vector<double> &x_, const Vector<double> &y) {

    Vector<double> &x = const_cast<Vector<double> &>(x_);

#ifdef FMATVEC_SIZE_CHECK
    assert(x.size() == y.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    daxpy(y.size(), 1., y(), y.inc(), x(), x.inc());

    return x;
  }

  Vector<double>& operator-=(const Vector<double> &x_, const Vector<double> &y) {

    Vector<double> &x = const_cast<Vector<double> &>(x_);

#ifdef FMATVEC_SIZE_CHECK
    assert(x.size() == y.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    daxpy(y.size(), -1., y(), y.inc(), x(), x.inc());

    return x;
  }

  Vector<double> operator+(const Vector<double> &x, const Vector<double> &y) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(x.size() == y.size());
#endif

    Vector<double> z(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return z;
#endif

    myblas_daxpy(x.size(),1, x(), x.inc(), y(), y.inc(), z());

    /* dcopy(x.size(), x(), x.inc(), z(),z.inc());
       daxpy(y.size(), 1., y(), y.inc(), z(),z.inc());
       */
    return z;
  }

  Vector<double> operator-(const Vector<double> &x, const Vector<double> &y) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(x.size() == y.size());
#endif

    Vector<double> z(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return z;
#endif

    dcopy(x.size(), x(), x.inc(), z(),z.inc());
    daxpy(y.size(), -1., y(), y.inc(), z(),z.inc());

    return z;
  }

  //-------------------------------------
  // RowVector operations
  //-------------------------------------

  RowVector<double>& operator/=(const RowVector<double> &x_, const double &alpha) {

    RowVector<double> &x = const_cast<RowVector<double> &>(x_);

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    dscal(x.size(), 1./alpha, x(), x.inc());

    return x;
  }

  RowVector<double>& operator*=(const RowVector<double> &x_, const double &alpha) {

    RowVector<double> &x = const_cast<RowVector<double> &>(x_);

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    dscal(x.size(), alpha, x(), x.inc());

    return x;
  }

  RowVector<double> operator*(const RowVector<double> &x, double alpha) {

    RowVector<double> y(x.size());//=x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    myblas_dscal(x.size(),alpha,x(),x.inc(),y());

    return y;
  }

  RowVector<double> operator*(double alpha, const RowVector<double> &x) {

    RowVector<double> y(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    myblas_dscal(x.size(),alpha,x(),x.inc(),y());

    return y;
  }

  RowVector<double> operator/(const RowVector<double> &x, double alpha) {

    RowVector<double> y(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    myblas_dscal(x.size(),1./alpha,x(),x.inc(),y());

    return y;
  }

  //-------------------------------------
  // RowVector/RowVector operations
  //-------------------------------------


  RowVector<double>& operator+=(const RowVector<double> &x_, const RowVector<double> &y) {

    RowVector<double> &x = const_cast<RowVector<double> &>(x_);

#ifdef FMATVEC_SIZE_CHECK
    assert(x.size() == y.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    daxpy(y.size(), 1., y(), y.inc(), x(), x.inc());

    return x;
  }

  RowVector<double>& operator-=(const RowVector<double> &x_, const RowVector<double> &y) {

    RowVector<double> &x = const_cast<RowVector<double> &>(x_);

#ifdef FMATVEC_SIZE_CHECK
    assert(x.size() == y.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return x;
#endif

    daxpy(y.size(), -1., y(), y.inc(), x(), x.inc());

    return x;
  }

  RowVector<double> operator+(const RowVector<double> &x, const RowVector<double> &y) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(x.size() == y.size());
#endif 

    RowVector<double> z(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return z;
#endif

    myblas_daxpy(x.size(),1, x(), x.inc(), y(), y.inc(), z());

    /* dcopy(x.size(), x(), x.inc(), z(),z.inc());
       daxpy(y.size(), 1., y(), y.inc(), z(),z.inc());
       */
    return z;
  }

  RowVector<double> operator-(const RowVector<double> &x, const RowVector<double> &y) {


#ifdef FMATVEC_SIZE_CHECK 
    assert(x.size() == y.size());
#endif

    RowVector<double> z(x.size());

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return z;
#endif

    dcopy(x.size(), x(), x.inc(), z(),z.inc());
    daxpy(y.size(), -1., y(), y.inc(), z(),z.inc());

    return z;
  }


  //-------------------------------------
  // Matrix/Vector operations
  //-------------------------------------

  Vector<double> operator*(const Matrix<General, double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.cols() == x.size());
#endif

    Vector<double> y(A.rows());

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0)
      return y;
    else if(A.cols() == 0) {
      y.init(0);
      return y;
    }
#endif

    int m, n;

    if(A.transposed()) {
      m=A.cols();
      n=A.rows();
    } else { 
      m=A.rows();
      n=A.cols();
    }

    dgemv(A.blasOrder(),A.blasTrans(),m,n, 1 ,A(), A.ldim(),x(),x.inc(),0, y(), y.inc());

    return y;
  }

  Vector<double> operator*(const Matrix<Symmetric, double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return y;
#endif

    dsymv(A.blasOrder(),A.blasUplo(),A.size(), 1.0 ,A(), A.ldim(),x(),x.inc(),0, y(), y.inc());

    return y;
  }

  Vector<double> operator*(const Matrix<Diagonal, double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return y;
#endif

    for(int m=0; m<y.size(); m++)
      y(m)=A(m)*x(m);

    return y;
  }

  //-------------------------------------
  // RowVector/Matrix operations
  //-------------------------------------

  RowVector<double> operator*(const RowVector<double> &x, const Matrix<General, double> &A) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.rows() == x.size());
#endif

    RowVector<double> y(A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.cols() == 0)
      return y;
    else if(A.rows() == 0) {
      y.init(0);
      return y;
    }
#endif

    int m, n;
    int incx = x.transposed()?x.ldim():1;
    CBLAS_TRANSPOSE transa = A.transposed()?CblasNoTrans:CblasTrans;

    if(A.transposed()) {
      m=A.cols();
      n=A.rows();
    } else { 
      m=A.rows();
      n=A.cols();
    }

    dgemv(A.blasOrder(),transa,m,n, 1 ,A(), A.ldim(),x(),x.inc(),0, y(), y.inc());

    return y;
  }

  RowVector<double> operator*(const RowVector<double> &x, const Matrix<Symmetric, double> &A) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    RowVector<double> y(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return y;
#endif

    dsymv(A.blasOrder(),A.blasUplo(),A.size(), 1.0 ,A(), A.ldim(),x(),x.inc(),0, y(), y.inc());

    return y;
  }

  RowVector<double> operator*(const RowVector<double> &x, const Matrix<Diagonal, double> &A) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.rows() == x.size());
#endif

    RowVector<double> y(A.cols());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return y;
#endif

    for(int m=0; m<y.size(); m++)
      y(m)=A(m)*x(m);
    return y;
  }

  //-------------------------------------
  // RowVector/Vector operations
  //-------------------------------------

  double operator*(const RowVector<double> &x, const Vector<double> &y) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(x.size() == y.size());
#endif

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return 0.0;
#endif

    return ddot(x.size(), x(), x.inc(), y(),y.inc());
  }

  //-------------------------------------
  // operations
  //-------------------------------------

  Matrix<General, double> slvLU(const SquareMatrix<double> &A, const Matrix<General, double> &X) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == X.rows());
#endif

    Matrix<General, double> Y = X.copy();

#ifdef FMATVEC_VOID_CHECK
    if(X.rows() == 0 || X.cols() == 0)
      return Y;
#endif

    SquareMatrix<double> B = A.copy();

    int *ipiv = new int[A.size()];

    int info = dgesv(B.blasOrder(), B.size(), Y.cols(), B(), B.ldim(), ipiv, Y(), Y.ldim());

    delete [] ipiv;

    assert(info==0);

    return Y;  
  }

  Matrix<General, double> slvLUFac(const SquareMatrix<double> &A, const Matrix<General, double> &X, const Vector<int> &ipiv) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == X.rows());
#endif

    Matrix<General, double> Y = X.copy();

#ifdef FMATVEC_VOID_CHECK
    if(X.rows() == 0 || X.cols() == 0)
      return Y;
#endif

    SquareMatrix<double> B = A.copy();

    int info = dgetrs(B.blasOrder(), B.blasTrans(), B.size(), Y.cols(), B(), B.ldim(), ipiv(), Y(), Y.ldim());

    assert(info==0);

    return Y;  
  }

  Matrix<General, double> slvQR(const SquareMatrix<double> &A, const Matrix<General, double> &X) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == X.rows());
#endif

    Matrix<General, double> Y = X.copy();

#ifdef FMATVEC_VOID_CHECK
    if(X.rows() == 0 || X.cols() == 0)
      return Y;
#endif

    SquareMatrix<double> B = A.copy();

    int info = dgels( B.blasTrans(), B.rows(), B.cols(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());

    assert(info==0);

    return Y;
  }

  Vector<double> slvLU(const SquareMatrix<double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y = x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    SquareMatrix<double> B = A.copy();

    int *ipiv = new int[A.size()];

    int info = dgesv(B.blasOrder(), B.size(), 1, B(), B.ldim(), ipiv, y(), y.size());

    delete [] ipiv;

    assert(info==0);

    return y;  
  }

  Vector<double> slvLUFac(const SquareMatrix<double> &A, const Vector<double> &x, const Vector<int> &ipiv) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y = x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    SquareMatrix<double> B = A.copy();

    int info = dgetrs(B.blasOrder(), B.blasTrans(), B.size(), 1, B(), B.ldim(), ipiv(), y(), y.size());

    assert(info==0);

    return y;  
  }

  Vector<double> slvQR(const SquareMatrix<double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y = x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    SquareMatrix<double> B = A.copy();

    int info = dgels( B.blasTrans(), B.rows(), B.cols(), y.cols(), B(), B.ldim(), y(), y.size());

    assert(info==0);

    return y;  
  }

  SquareMatrix<double> inv(const SquareMatrix<double> &A) {

    SquareMatrix<double> B=A.copy();

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return B;
#endif

    int *ipiv = new int[A.size()];

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv);

    assert(info==0);

    dgetri( B.blasOrder(), B.size(), B(), B.ldim(), ipiv);

    delete [] ipiv;

    return B;
  }

  Matrix<Symmetric, double> inv(const Matrix<Symmetric, double> &A) {

    Matrix<Symmetric, double> B=A.copy();

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return B;
#endif

    int info = dpotrf(B.blasOrder(), B.blasUplo(), B.size(), B(), B.ldim());

    assert(info==0);

    dpotri( B.blasOrder(), B.blasUplo(), B.rows(), B(), B.ldim());

    return B;
  }

  Matrix<Diagonal, double> inv(const Matrix<Diagonal, double> &A) {

    Matrix<Diagonal, double> B(A.size());

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return B;
#endif

    for(int i=0; i<A.size(); i++)
      B(i) = 1.0/A(i);

    return B;
  }

  Matrix<General, double> facLU(const Matrix<General, double> &A, Vector<int> &ipiv) {

    Matrix<General, double> B=A.copy();

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return B;
#endif

    int is = A.rows() < A.cols() ? A.rows() : A.cols();
    if(ipiv.size() != is) {
      ipiv.resize(is);
    }

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv());

    assert(info==0);

    return B;
  }

  SquareMatrix<double> facLU(const SquareMatrix<double> &A, Vector<int> &ipiv) {

    SquareMatrix<double> B=A.copy();

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return B;
#endif

    int is = A.size();
    if(ipiv.size() != is) {
      ipiv.resize(is);
    }

    int info = dgetrf(B.blasOrder(), B.rows(), B.cols(), B(), B.ldim(), ipiv());

    assert(info==0);

    return B;
  }

  Matrix<Symmetric, double> facLL(const Matrix<Symmetric, double> &A) {

    Matrix<Symmetric, double> B=A.copy();

#ifdef FMATVEC_VOID_CHECK
    if(A.size() == 0)
      return B;
#endif

    int info = dpotrf(B.blasOrder(), B.blasUplo(), B.size(), B(), B.ldim());

    assert(info==0);

    return B;
  }

  double nrm1(const Vector<double> &x) {

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return 0.0;
#endif

    return dasum(x.size(), x(), x.inc());
  }

  double nrmInf(const Vector<double> &x) {

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return 0.0;
#endif

    int id = idamax(x.size(), x(), x.inc());
    return fabs(x(id));
  }

  double nrm2(const Vector<double> &x) {

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return 0.0;
#endif

    return dnrm2(x.size(), x(), x.inc());
  }

  Vector<double> slvLL(const Matrix<Symmetric, double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y = x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    Matrix<Symmetric, double> B = A.copy();

    int info = dposv(B.blasOrder(), B.blasUplo(), B.size(), 1, B(), B.ldim(), y(), y.size());

    assert(info==0);

    return y;  
  }

  Matrix<General, double> slvLL(const Matrix<Symmetric, double> &A, const Matrix<General, double> &X) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == X.rows());
#endif

    Matrix<General, double> Y = X.copy();

#ifdef FMATVEC_VOID_CHECK
    if(X.rows() == 0 || X.cols() == 0)
      return Y;
#endif

    Matrix<Symmetric, double> B = A.copy();

    int info = dposv(B.blasOrder(), B.blasUplo(), B.size(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());

    assert(info==0);

    return Y;  
  }

  Vector<double> slvLLFac(const Matrix<Symmetric, double> &A, const Vector<double> &x) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == x.size());
#endif

    Vector<double> y = x.copy();

#ifdef FMATVEC_VOID_CHECK
    if(x.size() == 0)
      return y;
#endif

    Matrix<Symmetric, double> B = A.copy();

    int info = dpotrs(B.blasOrder(), B.blasUplo(), B.size(), 1, B(), B.ldim(), y(), y.size());

    assert(info==0);

    return y;  
  }

  Matrix<General, double> slvLLFac(const Matrix<Symmetric, double> &A, const Matrix<General, double> &X) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.size() == X.rows());
#endif

    Matrix<General, double> Y = X.copy();

#ifdef FMATVEC_VOID_CHECK
    if(X.rows() == 0 || X.cols() == 0)
      return Y;
#endif

    Matrix<Symmetric, double> B = A.copy();

    int info = dpotrs(B.blasOrder(), B.blasUplo(), B.size(), Y.cols(), B(), B.ldim(), Y(), Y.ldim());

    assert(info==0);

    return Y;  
  }

  Vector<complex<double> > eigval(const SquareMatrix<double> &A) {

    double *vl, *vr;
    double *wr = new double[A.size()];
    double *wi = new double[A.size()];

    SquareMatrix<double> B = A.copy();

    int info = dgeev('N','N', A.size(), B(), B.ldim(), wr, wi, vl, B.size(), vr, B.size());

    Vector<complex<double> > w(A.size());
    for(int i=0; i<A.size(); i++)
      w(i)=complex<double>(wr[i],wi[i]);

    delete [] wr;
    delete [] wi;

    return w;

  }

  Vector<double> eigval(const Matrix<Symmetric, double> &A) {

    Vector<double> w(A.size(),NONINIT);
    Matrix<Symmetric, double> B = A.copy();

    int info = dsyev('N', 'L', B.size(), B(), B.ldim(), w());

    return w;

  }

  Vector<double> eigvalSel(const Matrix<Symmetric, double> &A, int il, int iu, double abstol) {

    assert(il>=1);
    assert(iu>=il);
    assert(iu<=A.size());

    Matrix<Symmetric, double> B = A.copy();
    int ldb = B.ldim();
    double vl,vu;
    int m;
    Vector<double> w(A.size());
    double* z;
    int ldz = 1;

    int info = dsyevx('N', 'I', B.blasUplo(), B.size(), B(), B.ldim(), vl, vu, il, iu, abstol, &m, w(), z, ldz);

    return w(0,m-1);

  }


  double rho(const SquareMatrix<double> &A) {

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return 0.0;
#endif
    Vector<complex<double> > v = eigval(A);
    double buf = abs(v(0));
    for(int i=0; i<v.size(); i++) {
      double absi = abs(v(i));
      if(absi >= buf)
	buf = absi;
    }
    return buf;
  }

  double rho(const Matrix<Symmetric, double> &A) {

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return 0.0;
#endif
    Vector<double> v = eigval(A);
    double buf = fabs(v(0));
    for(int i=0; i<v.size(); i++) {
      double absi = fabs(v(i));
      if(absi >= buf)
	buf = absi;
    }
    return buf;
  }

  double nrmInf(const Matrix<General,double> &A) {

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return 0.0;
#endif

    return dlange('I' , A.rows(), A.cols(), A(), A.ldim());
  }

  double nrm1(const Matrix<General,double> &A) {

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return 0.0;
#endif

    return dlange('1' , A.rows(), A.cols(), A(), A.ldim());
  }

  double nrmFro(const Matrix<General,double> &A) {

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return 0.0;
#endif

    return dlange('F' , A.rows(), A.cols(), A(), A.ldim());
  }

  double nrm2(const Matrix<General,double> &A) {

#ifdef FMATVEC_VOID_CHECK
    if(A.rows() == 0 || A.cols() == 0)
      return 0.0;
#endif

    return sqrt(rho(JTJ(A)));
  }


  Vector<double> slvLS(const Matrix<General,double> &A, const Vector<double> &b, double rcond) {

#ifdef FMATVEC_SIZE_CHECK 
    assert(A.rows() == b.size());
#endif

    Matrix<General,double> A_ = A.copy();

    Vector<double> b_(A.rows()>A.cols()?A.rows():A.cols() ,NONINIT);
    b_(0,b.size()-1) = b;

    int info = dgelss( A.rows(), A.cols(), 1, A_(), A_.ldim(), b_(), b_.size(), rcond);

    assert(info == 0);

    return b_(0,A.cols()-1);
  }

//  Matrix<General, double> swap(const Matrix<General, double> &X, const Vector<int> &ipiv ) {
//
//    //#ifdef FMATVEC_SIZE_CHECK 
//    //    assert(A.size() == X.rows());
//    //#endif
//
//    Matrix<General, double> Y = X.copy();
//
//    ATL_dlaswp( Y.cols(), Y(), Y.ldim(), 0, Y.rows(), ipiv(), 1 );
//
//    return Y;
//  }
//  Matrix<General, double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_DIAG unit, const SquareMatrix<double> &A, const Matrix<General, double> &X, const Vector<int> &ipiv ) {
//
//    //#ifdef FMATVEC_SIZE_CHECK 
//    //    assert(A.size() == X.rows());
//    //#endif
//
//    Matrix<General, double> Y = X.copy();
//
//    //#ifdef FMATVEC_VOID_CHECK
//    //    if(X.rows() == 0 || X.cols() == 0)
//    //      return Y;
//    //#endif
//
//
//    ATL_dlaswp( Y.cols(), Y(), Y.ldim(), 0, Y.rows(), ipiv(), 1 );
//
//    SquareMatrix<double> B = A.copy();
//
//    dtrsm(CblasColMajor,side,uplo,B.blasTrans(),unit,Y.rows(), Y.cols(), 1, A(),A.ldim(), Y(), Y.ldim());
//
//    return Y;  
//  }
}
