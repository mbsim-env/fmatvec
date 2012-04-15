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

#ifndef fixed_var_general_matrix_h
#define fixed_var_general_matrix_h

#include "types.h"
#include <stdlib.h>

#define FMATVEC_SIZE_CHECK
#define FMATVEC_VOID_CHECK

namespace fmatvec {

  //template <int M, class AT> class FixedVector;

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type FixedGeneral and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <int M, class AT> class Matrix<FixedVarGeneral<M>, AT> {

    public:

 /// @cond NO_SHOW

    protected:

      int N;

      AT *ele;

      template <class Type> inline void deepCopy(const Matrix<Type, AT> &A); 
      inline void deepCopy(const Matrix<FixedVarGeneral<M>, AT> &A); 

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * */
      Matrix(int n) : N(n), ele(new AT[M*N]) {
#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      Matrix(int n, Initialization ini, const AT &a=0) : N(n), ele(new AT[M*N]) {  

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
      Matrix(const Matrix<FixedVarGeneral<M>, AT> &A) : N(A.N), ele(new AT[M*N]) {
	deepCopy(A);
      }

      template<class Type>
      explicit Matrix(const Matrix<Type, AT> &A) : N(A.cols()), ele(new AT[M*N]) {

#ifdef FMATVEC_SIZE_CHECK
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
      Matrix<FixedVarGeneral, AT>& resize(int n) {
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
      Matrix<FixedVarGeneral<M>, AT>& operator=(const Matrix<FixedVarGeneral<M>, AT> &A) {
	return operator<<(A);
      }

      template <class Type>
      Matrix<FixedVarGeneral<M>, AT>& operator=(const Matrix<Type, AT> &A) {
	return operator<<(A);
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> Matrix<FixedVarGeneral<M>, AT>& operator<<(const Matrix<T, AT> &A);

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      Matrix<FixedVarGeneral<M>, AT>& operator<<(const Matrix<FixedVarGeneral<M>, AT> &A);

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
      FixedVector<M,AT> col(int j);

      /*! \brief Column operator.
       *
       * see col(int)
       * */
      const FixedVector<M,AT> col(int j) const;

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      Matrix<FixedVarGeneral<M>, AT> copy() const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      Matrix<FixedVarGeneral<M>, AT>& init(const AT &a);

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      operator std::vector<std::vector<AT> >();

      /*! \brief std::vector<std::vector<AT> > Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT> > object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT> > the matrix will be initialized with. 
       * */
      Matrix(std::vector<std::vector<AT> > m);

      Matrix<FixedVarGeneral<M>, AT> T() {
	Matrix<FixedVarGeneral<M>, AT> A(NONINIT);
	for(int i=0; i<N; i++)
	  for(int j=0; j<M; j++)
	    A.e(i,j) = e(j,i);
	return A;
      }

      const Matrix<FixedVarGeneral<M>, AT> T() const {
	Matrix<FixedVarGeneral<M>, AT> A(NONINIT);
	for(int i=0; i<N; i++)
	  for(int j=0; j<M; j++)
	    A.e(i,j) = e(j,i);
	return A;
      }

      void setCol(int j, FixedVector<M,AT> &x) {
	for(int i=0; i<M; i++)
	  e(i,j) = x.e(i);
      }

  };
  // ------------------------- Constructors -------------------------------------
  template <int M, class AT> 
    Matrix<FixedVarGeneral<M>, AT>::Matrix(const char *strs) {
    // if 'strs' is a single scalar, surround it first with '[' and ']'.
    // This is more Matlab-like, because e.g. '5' and '[5]' is just the same.
    // (This functionallitiy is needed e.g. by MBXMLUtils (OpenMBV,MBSim))
    std::istringstream iss(strs);
    char c;
    iss>>c;
    if(c=='[') iss.str(strs);
    else iss.str(std::string("[")+strs+"]");

    int m = 0, N=0;
    int buf=0;
    iss >> c;
    AT x;
    do {
      iss >> x;
      iss >> c;
      if(c==';') {
	if(buf)
	  assert(buf == N);

	buf=N;
	N=0;
	m++;
      }
      else if(c==',')
	N++;
      c='0';
    } while(iss);

    N++; m++;
    assert(m==M);
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
  // ----------------------------------------------------------------------------

   template <int M, class AT> template< class Type>
    Matrix<FixedVarGeneral<M>, AT>& Matrix<FixedVarGeneral<M>, AT>::operator<<(const Matrix<Type, AT> &A) { 

#ifdef FMATVEC_SIZE_CHECK
      assert(A.rows() == M); 
      assert(A.cols() == N);
#endif

      deepCopy(A);

      return *this;
    }

   template <int M, class AT>
    Matrix<FixedVarGeneral<M>, AT>& Matrix<FixedVarGeneral<M>, AT>::operator<<(const Matrix<FixedVarGeneral<M>, AT> &A) { 

#ifdef FMATVEC_SIZE_CHECK
      assert(A.cols() == N);
#endif

      deepCopy(A);

      return *this;
    }

  template <int M, class AT>
    Matrix<FixedVarGeneral<M>, AT>&  Matrix<FixedVarGeneral<M>, AT>::init(const AT& val) {

      for(int i=0; i<M*N; i++) 
	e(i) = val;

      return *this;
    }

  template <int M, class AT>
    inline FixedVector<M,AT>  Matrix<FixedVarGeneral<M>, AT>::col(int j) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

      FixedVector<M,AT> x(NONINIT);

      for(int i=0; i<M; i++)
	x.e(i) = e(i,j);

      return x;
    }

  template <int M, class AT>
    inline const FixedVector<M,AT>  Matrix<FixedVarGeneral<M>, AT>::col(int j) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

      FixedVector<M,AT> x(NONINIT);

      for(int i=0; i<M; i++)
	x.e(i) = e(i,j);

      return x;

    }

  template <int M, class AT>
    Matrix<FixedVarGeneral<M>, AT> Matrix<FixedVarGeneral<M>, AT>::copy() const {

      Matrix<FixedVarGeneral<M>, AT> A(N,NONINIT);
      A.deepCopy(*this);

      return A;
    }


  template <int M, class AT>
    Matrix<FixedVarGeneral<M>, AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <int M, class AT>
    Matrix<FixedVarGeneral<M>, AT>::Matrix(std::vector<std::vector<AT> > m) {
#ifdef FMATVEC_SIZE_CHECK
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
   
  template <int M, class AT> template <class Type>
  inline void Matrix<FixedVarGeneral<M>, AT>::deepCopy(const Matrix<Type, AT> &A) { 
      for(int i=0; i<M; i++) 
	for(int j=0; j<N; j++)
	  e(i,j) = A.e(i,j);
    }

  template<int M, class AT>
  inline void Matrix<FixedVarGeneral<M>,AT>::deepCopy(const Matrix<FixedVarGeneral<M>,AT> &A) {
      for(int i=0; i<M*N; i++) 
	e(i) = A.e(i);
    }

   /// @endcond

}

#endif
