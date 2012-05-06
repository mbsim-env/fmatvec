/* Copyright (C) 2003-2005  Martin Förg, Rober Huber

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

#ifndef linear_algebra_h
#define linear_algebra_h

#include "square_matrix.h"
#include "fixed_vector.h"

namespace fmatvec {

  /*! \brief Transpose of a vector.
   *
   * This function computes the transpose of a vector. The result is a rowvector.
   * \return The transpose.
   * */
  template <class AT>
    RowVector<General, AT> trans(const Vector<General, AT> &x) {

      return RowVector<General, AT>(x.m,x.lda,x.tp?false:true,x.memory,x.ele).copy();
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class AT>
    Vector<General, AT> operator-(const Vector<General, AT> &x) {

      Vector<General, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
	y(i)=-x(i);

      return y;
    }

  /*! \brief Vector-scalar multiplication.
   *
   * This function computes the product of a vector 
   * and a scalar.
   * \return The product.
   * */
  template <class AT>
    Vector<General, AT> operator*(const Vector<General, AT> &x, const AT& alpha) {

      Vector<General, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
	y(i) = x(i)*alpha;

      return y;
    }

  /*! \brief Scalar-vector multiplication.
   *
   * \see operator*(const Vector<General, AT>&x,const AT&).
   * */
  template <class AT>
    Vector<General, AT> operator*(const AT& alpha, const Vector<General, AT> &x) {

      Vector<General, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
	y(i) = x(i)*alpha;

      return y;
    }

  /*! \brief Transpose of a rowvector.
   *
   * This function computes the transpose of a rowvector.
   * The result is a vector.
   * \param x A vector.
   * \return The transpose.
   * */
  template <class AT>
    Vector<General, AT> trans(const RowVector<General, AT> &x) {

      return Vector<General, AT>(x.n,x.lda,x.tp?false:true,x.memory,x.ele).copy();
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a rowvector.
   * \return The negation.
   * */
  template <class AT>
    RowVector<General, AT> operator-(const RowVector<General, AT> &a) {

      RowVector<General, AT> c(a.size(),NONINIT);

      for(int i=0; i<a.size(); i++)
	c(i)=-a(i);

      return c;
    }

  /*! \brief Rowvector-scalar multiplication.
   *
   * This function computes the product of a rowvector 
   * and a scalar.
   * \return The product.
   * */
  template <class AT>
    RowVector<General, AT> operator*(const RowVector<General, AT> &x, const AT& alpha) {

      RowVector<General, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
	y(i) = x(i)*alpha;

      return y;
    }

  /*! \brief Scalar-rowvector multiplication.
   *
   * \see operator*(const RowVector<General, AT>&, const AT&).
   * */
  template <class AT>
    RowVector<General, AT> operator*(const AT &alpha, const RowVector<General, AT> &x) {

      RowVector<General, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
	y(i) = x(i)*alpha;

      return y;
    }

  /*! \brief Vector-vector addition.
   *
   * This function computes the sum of two vectors
   * \return The sum.
   * */
  template <class AT>
    inline Vector<General, AT> operator+(const Vector<General, AT> &x, const Vector<General, AT> &y) {

      assert(y.size() == x.size());

      Vector<General, AT> z(x.size(),NONINIT);

      if(y.transposed()) {
	if(x.transposed())
	  for(int i=0; i<x.size(); i++)
	    z.er(i) = x.et(i) + y.et(i);
	else
	  for(int i=0; i<x.size(); i++)
	    z.er(i) = x.er(i) + y.et(i);
      }
      else {
	if(x.transposed())
	  for(int i=0; i<x.size(); i++)
	    z.er(i) = x.et(i) + y.er(i);
	else 
	  for(int i=0; i<x.size(); i++)
	    z.er(i) = x.er(i) + y.er(i);
      } 

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT>
    Vector<General, AT> operator-(const Vector<General, AT> &x, const Vector<General, AT> &y) {

      assert(y.size() == x.size());

      Vector<General, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
	z(i) = x(i) - y(i);

      return z;
    }

  /*! \brief Scalar product.
   *
   * This function computes the scalar product of two vectors.
   * \return The scalar product.
   * */
  template <class AT>
    AT scalarProduct(const Vector<General, AT> &x, const Vector<General, AT> &y) {

      assert(x.size()==y.size());

      AT res = 0;

      for(int i=0; i<x.size(); i++)
	res += x(i)*y(i);

      return res;
    }

  /*! \brief Cross product.
   *
   * This function computes the cross product of two vectors.
   * \return The cross product.
   * */
  template <class AT>
    Vector<General, AT> crossProduct(const Vector<General, AT> &x, const Vector<General, AT> &y) {

      assert(x.size()==3);
      assert(y.size()==3);

      Vector<General, AT> z(3,NONINIT);

      z(0) = x(1)*y(2) - x(2)*y(1);
      z(1) = x(2)*y(0) - x(0)*y(2);
      z(2) = x(0)*y(1) - x(1)*y(0);

      return z;
    }

  /*! \brief Triple product.
   *
   * This function computes the triple product of three vectors.
   * \return The triple product.
   * */
  template <class AT>
    double tripleProduct(const Vector<General, AT> &a, const Vector<General, AT> &x, const Vector<General, AT> &y) {

      assert(a.size()==3);
      assert(x.size()==3);
      assert(y.size()==3);

      return a(0)*(x(1)*y(2)-x(2)*y(1)) + a(1)*(x(2)*y(0)-x(0)*y(2)) + a(2)*(x(0)*y(1)-x(1)*y(0));
    }

  /*! \brief Rowvector-rowvector subtraction.
   *
   * This function computes the difference of two rowvectors
   * \return The difference.
   * */
  template <class AT>
    RowVector<General, AT> operator-(const RowVector<General, AT> &x, const RowVector<General, AT> &y) {

      assert(y.size() == x.size());

      RowVector<General, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
	z(i) = x(i) - y(i);

      return z;
    }

  /*! \brief Scalar product.
   *
   * This function computes the product of a vector and a rowvector. The result is
   * a scalar.
   * \return The product.
   * */
  template <class AT> 
    AT operator*(const RowVector<General, AT> &x, const Vector<General, AT> &y) {

      assert(x.size() == y.size());

      AT res=0;

      for(int i=0; i<x.size(); i++)
	res+=x(i)*y(i);

      return res;
    }

  /*! \brief Transpose of a matrix.
   *
   * This function computes the transpose of a matrix.
   * \f[ \boldsymbol{A} \rightarrow \boldsymbol{A}^T  \f]
   * The result is a matrix.
   * \param A A general matrix.
   * \return The transpose.
   * */
  template <class AT>
    Matrix<General, AT> trans(const Matrix<General, AT> &A) {

      return Matrix<General, AT>(A.n,A.m,A.lda,A.tp?false:true,A.memory,A.ele).copy();
    }

  /*! \brief Transpose of a matrix.
   *
   * This function computes the transpose of a square matrix.
   * \f[ \boldsymbol{A} \rightarrow \boldsymbol{A}^T  \f]
   * The result is a square matrix.
   * \param A A square matrix.
   * \return The transpose.
   * */
  template <class AT>
    SquareMatrix<General, AT> trans(const SquareMatrix<General, AT> &A) {

      return SquareMatrix<General, AT>(A.n, A.lda, A.tp?false:true, A.memory, A.ele).copy();
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class AT>
    SquareMatrix<General, AT> operator-(const SquareMatrix<General, AT> &A) {

      SquareMatrix<General, AT> B(A.size(),NONINIT);

      for(int i=0; i<A.rows(); i++)
	for(int j=0; j<A.cols(); j++)
	  B(i,j)=-A(i,j);

      return B;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a matrix
   * \return The negation. 
   * */
  template <class AT, class Type>
    Matrix<Type, AT> operator-(const Matrix<Type, AT> &A) {

      Matrix<Type, AT> C(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++)
	for(int j=0; j<A.cols(); j++)
	  C(i,j)=-A(i,j);

      return C;
    }

  /*! \brief Matrix-scalar multiplication.
   *
   * This function computes the product of a matrix 
   * and a scalar.
   * \param A The matrix. 
   * \param alpha The scalar. 
   * \return The product.
   * */
  template <class AT, class Type>
    Matrix<Type, AT> operator*(const Matrix<Type, AT > &A, const AT &alpha) {

      Matrix<Type, AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++) 
	for(int j=0; j<A.cols(); j++) 
	  B(i,j) = A(i,j)*alpha;

      return B;
    }

  /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<Type, AT >&, const AT&).
   * */
  template <class AT, class Type>
    Matrix<Type, AT> operator*(const AT &alpha, const Matrix<Type, AT > &A) {

      Matrix<Type, AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++) 
	for(int j=0; j<A.cols(); j++) 
	  B(i,j) = A(i,j)*alpha;

      return B;
    }

  /*! \brief Tilde operator
   *
   *  Example:
   * \code 
   * tilde(x);
   * \endcode 
   *  \f[
   *  x=\begin{pmatrix}x_1 \\ x_2 \\ x_3\end{pmatrix} \quad \Rightarrow \quad
   *  \tilde x = \begin{pmatrix} 0 & -x_3 & x_2 \\ x_3 &  0 & -x_1 \\ -x_2 & x_1 & 0\end{pmatrix}
   *  \f]
   * */
  template <class AT>
    SquareMatrix<General, AT> tilde(const Vector<General, AT> &x) {

      SquareMatrix<General, AT> B(x.size(),NONINIT);

      B(0,0) =  0;
      B(1,1) =  0;
      B(2,2) =  0;
      B(0,1) = -x(2);
      B(0,2) =  x(1);
      B(1,0) =  x(2);
      B(1,2) = -x(0);
      B(2,0) = -x(1);
      B(2,1) =  x(0);

      return B;
    }

  /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a matrix 
   * and a vector.
   * \return The product.
   * */
  template <class AT>
    inline Vector<General, AT> operator*(const Matrix<General, AT> &A, const Vector<General, AT> &x) {

      assert(A.cols() == x.size());

      Vector<General, AT> z(A.rows(),NONINIT);

      if(x.transposed()) {
	if(A.transposed()) {
	  for(int i=0; i<z.size(); i++) {
	    z.er(i)=0;
	    for(int j=0; j<x.size(); j++)
	      z.er(i) += A.et(i,j) * x.et(j);
	  }
	} 
	else {
	  for(int i=0; i<z.size(); i++) {
	    z.er(i)=0;
	    for(int j=0; j<x.size(); j++)
	      z.er(i) += A.er(i,j) * x.et(j);
	  }
	}
      }
      else {
	if(A.transposed()) {
	  for(int i=0; i<z.size(); i++) {
	    z.er(i)=0;
	    for(int j=0; j<x.size(); j++)
	      z.er(i) += A.et(i,j) * x.er(j);
	  }
	} 
	else {
	  for(int i=0; i<z.size(); i++) {
	    z.er(i)=0;
	    for(int j=0; j<x.size(); j++)
	      z.er(i) += A.er(i,j) * x.er(j);
	  }
	}
      } 

      return z;
    }

  /*! \brief Vector-matrix multiplication.
   *
   * This function computes the product of a vector 
   * and a matrix.
   * \return The product.
   * */
  template <class AT, class Type> 
    RowVector<General, AT> operator*(const RowVector<General, AT> &x, const Matrix<Type, AT> &A) {

      assert(x.size() == A.rows());

      RowVector<General, AT> z(A.cols(),NONINIT);

      for(int i=0; i<z.size(); i++) {
	z(i)=0;
	for(int j=0; j<x.size(); j++)
	  z(i) += x(j) * A(j,i);
      }

      return z;
    }

   /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
 template <class AT, class Type1, class Type2> 
    Matrix<General, AT> operator+(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {

      assert(A.rows() == B.rows());
      assert(A.cols() == B.cols());

      Matrix<General, AT> C(A.rows(),A.cols(),NONINIT);
      for(int i=0; i<A.rows(); i++)
	for(int j=0; j<A.cols(); j++)
	  C(i,j) = A(i,j) + B(i,j);

      return C;
    }

   /*! \brief Matrix-matrix subtraction.
   *
   * This function computes the difference of two matrices
   * \return The difference.
   * */
  template <class AT, class Type1, class Type2>
    Matrix<General, AT> operator-(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {

      assert(A.rows() == B.rows());
      assert(A.cols() == B.cols());

      Matrix<General, AT> C(A.rows(),A.cols(),NONINIT);
      for(int i=0; i<A.rows(); i++)
	for(int j=0; j<A.cols(); j++)
	  C(i,j) = A(i,j) - B(i,j);

      return C;
    }

   /*! \brief Matrix-matrix multiplication.
   *
   * This function computes the product of two matrices
   * \return The product.
   * */
  template <class AT, class Type1, class Type2> 
    Matrix<General, AT> operator*(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {

      assert(A.cols() == B.rows());

      Matrix<General, AT> C(A.rows(),B.cols(),NONINIT);

      for(int i=0; i<C.rows(); i++) {
	for(int j=0; j<C.cols(); j++) {
	  C(i,j)=0;
	  for(int k=0; k<A.cols(); k++)
	    C(i,j) += A(i,k) * B(k,j);
	}
      }

      return C;
    }

   /*! \brief Matrix-matrix comparison.
   *
   * This function compares two matrices
   * \return True, if the matrices are identical, false otherwise.
   * */
  template <class AT, class Type1, class Type2> 
    bool operator==(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {

      if(&A == &B)
	return true;

      if(A.rows() != B.rows() || A.cols() != B.cols())
	return false;

      for(int i=0; i<A.rows(); i++) 
	for(int j=0; j<A.cols(); j++)
	  if(A(i,j)!=B(i,j))
	    return false;

      return true;
    }

   /*! \brief Matrix-matrix comparison.
   *
   * This function compares two matrices.
   * \return True, if the matrices are different, false otherwise.
   * */
  template <class AT, class Type1, class Type2> 
    bool operator!=(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {

      if(A.rows() != B.rows() || A.cols() != B.cols())
	return true;

      for(int i=0; i<A.rows(); i++) 
	for(int j=0; j<A.cols(); j++)
	  if(A(i,j)!=B(i,j))
	    return true;

      return false;
    }

   /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a sparse matrix
   * and a vector.
   * \return The product.
   * */
  template <class AT> 
    Vector<General, AT> operator*(const Matrix<Sparse, AT> &A, const Vector<General, AT> &x) {
      assert(A.cols() == x.size());

      Vector<General, AT> z(A.rows(),NONINIT);

      const AT *a = A();
      const int *ia = A.Ip();
      const int *ja = A.Jp();
      for(int i=0; i<A.rows(); i++) {
	z(i) = 0;
	for(int j=ia[i]; j<ia[i+1]; j++)
	  z(i) += a[j]*x(ja[j]);
      }

      return z;
    }

   /*! \brief Maximum value.
   *
   * This function computes the maximum value of a vector.
   * \return The maximum value.
   * */
  template <class AT>
    AT max(const Vector<General, AT> &x) {
      assert(x.size() > 0);
      AT maximum = x(0);
      for(int i=1; i<x.size(); i++) {
	if(x(i) > maximum)
	  maximum = x(i);
      }
      return maximum;
    }

   /*! \brief Index of maximum value.
   *
   * This function computes the index of the maximum value of a vector
   * \return The index of the maximum value.
   * */
  template <class AT>
    int maxIndex(const Vector<General, AT> &x) {
      assert(x.size() > 0);
      AT maximum = x(0);
      int index = 0;
      for(int i=1; i<x.size(); i++) {
	if(x(i) > maximum) {
	  maximum = x(i);
	  index = i;
	}
      }
      return index;
    }

   /*! \brief Minimum value.
   *
   * This function computes the minimum value of a vector.
   * \return The minimum value.
   * */
  template <class AT>
    AT min(const Vector<General, AT> &x) {
      assert(x.size() > 0);
      AT minimum = x(0);
      for(int i=1; i<x.size(); i++) {
	if(x(i) < minimum)
	  minimum = x(i);
      }
      return minimum;
    }

   /*! \brief Index of minimum value.
   *
   * This function computes the index of the minimum value of a vector
   * \return The index of the minimum value.
   * */
  template <class AT>
    int minIndex(const Vector<General, AT> &x) {
      assert(x.size() > 0);
      AT minimum = x(0);
      int index = 0;
      for(int i=1; i<x.size(); i++) {
	if(x(i) < minimum) {
	  minimum = x(i);
	  index = i;
	}
      }
      return index;
    }

  // HR 28.09.2006
  /*! \brief Bubble Sort Algorithm (stable sorting Algorithm )
   *
   * Values of rowvectors of Matrix A are sorted in ascending order 
   * \param A_ Matrix to be sorted
   * \param PivotCol Column of A used as sorting index 
   */
  template <class AT>
    Matrix<General, AT> bubbleSort(const Matrix<General, AT> &A_, int PivotCol) {
      assert(A_.rows()>0);
      assert(A_.cols()>0);
      assert(A_.cols()> PivotCol);
      Matrix<General, AT> A = A_.copy();
      int i,j,N;
      RowVector<General, AT> tmp;
      N = A.rows();
      for (i=1; i<= N-1; i++) {
        for (j=0; j< N-1; j++) {
          if (A(j,PivotCol) >  A(j+1,PivotCol)) {
            tmp = A.row(j);  
            A.row(j) = A.row(j+1);
            A.row(j+1) = tmp;
          }
        }
      }
      return A;
    }

  // HR 28.09.2006
  /*! internal function of QuickSortMedian */
  template <class AT>
    void quicksortmedian_intern(Matrix<General, AT> &A, int PivotCol, RowVector<General, AT> &tmp, int l, int r){

      if(r>l){
	int i=l-1, j=r;
	//RowVec tmp;     
	if(r-l > 3){ //Median of three
	  int m=l+(r-l)/2;
	  if(A(l,PivotCol)>A(m,PivotCol))
	  { tmp= A.row(l);
	    A.row(l) = A.row(m);
	    A.row(m) =tmp; }
	    if(A(l,PivotCol)> A(r,PivotCol))
	    { tmp= A.row(l);
	      A.row(l) = A.row(r);
	      A.row(r) =tmp; }
	    else if(A(r,PivotCol) > A(m,PivotCol) )
	    { tmp= A.row(r);
	      A.row(r) = A.row(m);
	      A.row(m) =tmp; }
	}

	for(;;){
	  while(A(++i,PivotCol) < A(r,PivotCol));
	  while(A(--j,PivotCol) > A(r,PivotCol) && j>i);
	  if(i>=j) break;
	  tmp= A.row(i);
	  A.row(i) = A.row(j);
	  A.row(j) =tmp;
	}
	tmp= A.row(i);
	A.row(i) = A.row(r);
	A.row(r) =tmp;
	quicksortmedian_intern(A, PivotCol,tmp, l, i-1);
	quicksortmedian_intern(A, PivotCol,tmp, i+1, r);
      }
    }

  // HR 28.09.2006
  /*! \brief Modified QuickSort Algorithm (unstable sorting Algorithm )
   *
   * Values of rowvectors of Matrix A are sorted in ascending order 
   * unstabel but very quick sorting Algorithm
   * Pivot Elements as 'Median of three'
   * \param A_ Matrix to be sorted
   * \param PivotCol Column of A used as sorting index 
   */
  template <class AT>
    Matrix<General, AT> quickSortMedian(const Matrix<General, AT> &A_, int PivotCol) {
      Matrix<General, AT> A = A_.copy();
      int N = A.rows();
      RowVector<General, AT> tmp;
      quicksortmedian_intern(A, PivotCol,tmp, 0, N-1);
      return A;
    }

  /*! \brief Count nonzero elements.
   *
   * This function counts the nonzero elements of a matrix. ALL diagonal
   * elements are treated as NONZERO!!! (See the storage format for sparse matrix)
   * \return The number of nonzero or diagonal elemets.
   * */
  template <class AT> int countElements(const SquareMatrix<General, AT> &A) { 
    int k=0;
    for(int i=0; i<A.size(); i++) {
      for(int j=0; j<A.size(); j++) {
	if(fabs(A(i,j))>1e-16 || i==j) {
	  k++;
	}
      }
    }
    return k;
  }

  template <int M, class AT>
    inline Vector<General, AT>& operator+=(const Vector<General, AT> &x_, const Vector<GeneralFixed<M,1>, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x_.size() == M);
#endif
      Vector<General, AT> &x = const_cast<Vector<General, AT> &>(x_);
      for(int i=0; i<M; i++) 
	x.e(i) += y.e(i);
      return x;
    }

  template <class AT>
    inline Vector<General, AT>& operator+=(const Vector<General, AT> &x_, const Vector<GeneralVar, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x_.size() == y.size());
#endif
      Vector<General, AT> &x = const_cast<Vector<General, AT> &>(x_);
      for(int i=0; i<x.size(); i++) 
	x.e(i) += y.e(i);
      return x;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const Matrix<GeneralFixedVar<M>, AT> &A, const Vector<General, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<A.cols(); j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

}

#endif
