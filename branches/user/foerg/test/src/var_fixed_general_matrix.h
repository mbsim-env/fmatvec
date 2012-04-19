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
 *   martin.o.foerg@googlemail.com
 *
 */

#ifndef var_fixed_general_matrix_h
#define var_fixed_general_matrix_h

#include "types.h"
#include <stdlib.h>

namespace fmatvec {

  //template <int N, class AT> class FixedVector;

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type FixedGeneral and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <int N, class AT> class Matrix<VarFixedGeneral<N>, AT> {

    public:

 /// @cond NO_SHOW

    protected:

      int M;

      AT *ele;

      template <class Type> inline void deepCopy(const Matrix<Type, AT> &A); 
      inline void deepCopy(const Matrix<VarFixedGeneral<N>, AT> &A); 

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * */
      Matrix(int m=0) : M(m), ele(new AT[M*N]) {
#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      Matrix(int m, Initialization ini, const AT &a=0) : M(m), ele(new AT[M*N]) {  

	if(ini == INIT) {
	  for(int i=0; i<M*N; i++) 
	    e(i) = a;
	} else if(ini == EYE ) {
	  for(int i=0; i<M; i++) {
	    for(int j=0; j<N; j++) {
	      if (i==j) e(i,j) = 1; // operator()(i,i) = 1;
	      else e(i,j) = 0; // operator()(i,j) = 0;
	    }
	  }
	}
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      Matrix(const Matrix<VarFixedGeneral<N>, AT> &A) : M(A.M), ele(new AT[M*N]) {
	deepCopy(A);
      }

      template<class Type>
      explicit Matrix(const Matrix<Type, AT> &A) : M(A.rows()), ele(new AT[M*N]) {

#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == M); 
#endif

	deepCopy(A);
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<FixedGeneral,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const char *str);

      /*! \brief Destructor. 
       * */
      ~Matrix() {
	delete[] ele;
      }

      /*! \brief Matrix resizing. 
       *
       * Resizes the matrix to size m x n.    
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \return A reference to the calling matrix.
       * \remark The matrix will be initialised to
       * zero by default. To change this behavior, define
       * FMATVEC_NO_INITIALIZATION.
       * */
      Matrix<VarFixedGeneral<N>, AT>& resize(int n) {
	delete[] ele;
	N=n;
	ele = new AT[M*N];

#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      Matrix<VarFixedGeneral<N>, AT>& operator=(const Matrix<VarFixedGeneral<N>, AT> &A) {
	return operator<<(A);
      }

      template <class Type>
      Matrix<VarFixedGeneral<N>, AT>& operator=(const Matrix<Type, AT> &A) {
	return operator<<(A);
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> 
	inline Matrix<VarFixedGeneral<N>, AT>& operator<<(const Matrix<T, AT> &A);

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<VarFixedGeneral<N>, AT>& operator<<(const Matrix<VarFixedGeneral<N>, AT> &A);

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \remark The bounds are checked by default. 
       * To change this behavior, define
       * FMATVEC_NO_BOUNDS_CHECK.
       * \sa operator()(int,int) const
       * */
      AT& operator()(int i, int j) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<M);
	assert(j<N);
#endif

	return e(i,j);
      };

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<M);
	assert(j<N);
#endif

	return e(i,j);
      };

      AT& e(int i, int j) {
	return ele[i+j*M];
      };

      /*! \brief Element operator
       *
       * See e(int,int) 
       * */
      const AT& e(int i, int j) const {
	return ele[i+j*M];
      };

      AT& e(int i) {
	return ele[i];
      };

      /*! \brief Element operator
       *
       * See e(int,int) 
       * */
      const AT& e(int i) const {
	return ele[i];
      };

      /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return ele;};

      /*! \brief Pointer operator
       *
       * See operator()() 
       * */
      const AT* operator()() const {return ele;};

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix.
       * */
      int rows() const {return M;};

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix.
       * */
      int cols() const {return N;};

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix
       * */
      int ldim() const {return M;};

      /*! \brief Transposed status.
       *
       * Returns the blas-conform transposed status.
       * \return CblasTrans if the matrix is in transposed form, CblasNoTrans
       * otherwise. 
       * */
      const CBLAS_TRANSPOSE blasTrans() const {
	return CblasNoTrans;
      };

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasColMajor.
       * */
      const CBLAS_ORDER blasOrder() const {
	return  CblasColMajor;
      };

      /*! \brief Column operator.
       *
       * Returns a vector containing the i-th column of the calling matrix. 
       * \attention The vector and the calling matrix will share the same physical memory.
       * \param i The column, that will be returned.  
       * \return A vector containing the i-th column of the calling matrix.
       * */
      inline VarVector<AT> col(int j);

      /*! \brief Column operator.
       *
       * see col(int)
       * */
      inline const VarVector<AT> col(int j) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<VarFixedGeneral<N>, AT>& init(const AT &a);

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >();

      /*! \brief std::vector<std::vector<AT> > Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT> > object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT> > the matrix will be initialized with. 
       * */
      inline Matrix(std::vector<std::vector<AT> > m);

      inline Matrix<FixedVarGeneral<N>, AT> T();

      inline const Matrix<FixedVarGeneral<N>, AT> T() const;

      template<int K> 
	inline void set(const Index &I, const Index &J, const Matrix<VarFixedGeneral<K>, AT> &A);

      template<int K> 
	inline void add(const Index &I, const Index &J, const Matrix<VarFixedGeneral<K>, AT> &A);


  };

  template <int N, class AT> 
    Matrix<VarFixedGeneral<N>, AT>::Matrix(const char *strs) {
      // if 'strs' is a single scalar, surround it first with '[' and ']'.
      // This is more Matlab-like, because e.g. '5' and '[5]' is just the same.
      // (This functionallitiy is needed e.g. by MBXMLUtils (OpenMBV,MBSim))
      std::istringstream iss(strs);
      char c;
      iss>>c;
      if(c=='[') iss.str(strs);
      else iss.str(std::string("[")+strs+"]");

      int M = 0, n=0;
      int buf=0;
      iss >> c;
      AT x;
      do {
        iss >> x;
        iss >> c;
        if(c==';') {
          if(buf)
            assert(buf == n);

          buf=n;
          n=0;
          M++;
        }
        else if(c==',')
          n++;
        c='0';
      } while(iss);

      n++; M++;
      assert(n==N);
      ele = new AT[M*N];
      iss.clear();
      iss.seekg(0);
      iss >> c;
      for(int i=0; i<M; i++)
        for(int j=0; j<N; j++) {
          iss >> e(i,j);
          iss >> c;
        }
    }

  template <int N, class AT> template< class Type>
    inline Matrix<VarFixedGeneral<N>, AT>& Matrix<VarFixedGeneral<N>, AT>::operator<<(const Matrix<Type, AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.rows() == M); 
      assert(A.cols() == N);
#endif

      deepCopy(A);

      return *this;
    }

  template <int N, class AT>
    inline Matrix<VarFixedGeneral<N>, AT>& Matrix<VarFixedGeneral<N>, AT>::operator<<(const Matrix<VarFixedGeneral<N>, AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.rows() == M);
#endif

      deepCopy(A);

      return *this;
    }

  template <int N, class AT>
    inline  Matrix<VarFixedGeneral<N>, AT>& Matrix<VarFixedGeneral<N>, AT>::init(const AT& val) {

      for(int i=0; i<M*N; i++) 
        e(i) = val;

      return *this;
    }

  template <int N, class AT>
    inline VarVector<AT> Matrix<VarFixedGeneral<N>, AT>::col(int j) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

      VarVector<AT> x(M,NONINIT);

      for(int i=0; i<M; i++)
        x.e(i) = e(i,j);

      return x;
    }

  template <int N, class AT>
    inline const VarVector<AT> Matrix<VarFixedGeneral<N>, AT>::col(int j) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

      VarVector<AT> x(M,NONINIT);

      for(int i=0; i<M; i++)
        x.e(i) = e(i,j);

      return x;

    }

  template <int N, class AT>
    inline Matrix<FixedVarGeneral<N>, AT> Matrix<VarFixedGeneral<N>, AT>::T() {
      Matrix<FixedVarGeneral<N>, AT> A(rows(),NONINIT);
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int N, class AT>
    inline const Matrix<FixedVarGeneral<N>, AT> Matrix<VarFixedGeneral<N>, AT>::T() const {
      Matrix<FixedVarGeneral<N>, AT> A(rows(),NONINIT);
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int N, class AT> template<int K>
    inline void Matrix<VarFixedGeneral<N>, AT>::set(const Index &I, const Index &J, const Matrix<VarFixedGeneral<K>, AT> &A) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(J.end()<N);
      assert(I.size()==A.rows());
      assert(J.size()==K);
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <int N, class AT> template<int K>
    inline void Matrix<VarFixedGeneral<N>, AT>::add(const Index &I, const Index &J, const Matrix<VarFixedGeneral<K>, AT> &A) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(J.end()<N);
      assert(I.size()==A.rows());
      assert(J.size()==K);
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <int N, class AT>
    inline Matrix<VarFixedGeneral<N>, AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <int N, class AT>
    inline Matrix<VarFixedGeneral<N>, AT>::Matrix(std::vector<std::vector<AT> > m) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(m.size() == M);
      assert(m[0].size() == N);
#endif
      for(int r=0; r<rows(); r++) {
        assert(m[r].size()==cols());
        for(int c=0; c<cols(); c++)
          e(r,c)=m[r][c];
      }
    }

  /// @cond NO_SHOW

  template <int N, class AT> template <class Type>
    inline void Matrix<VarFixedGeneral<N>, AT>::deepCopy(const Matrix<Type, AT> &A) { 
      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++)
          e(i,j) = A.e(i,j);
    }

  template<int N, class AT>
    inline void Matrix<VarFixedGeneral<N>,AT>::deepCopy(const Matrix<VarFixedGeneral<N>,AT> &A) {
      for(int i=0; i<M*N; i++) 
        e(i) = A.e(i);
    }

  /// @endcond

}

#endif
