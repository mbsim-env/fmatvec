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

#ifndef general_fixed_matrix_h
#define general_fixed_matrix_h

#include "types.h"
#include <stdlib.h>

#define FMATVEC_SIZE_CHECK
#define FMATVEC_VOID_CHECK

namespace fmatvec {

//  template <class AT> class Vector;
//  template <class AT> class RowVector;
//  template <class AT> class SquareMatrix;
  //template <class AT> class Matrix<Symmetric, AT>;

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type FixedSize and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <int M, int N, class AT> class Matrix<FixedSize<M,N>, AT> {

    public:

 /// @cond NO_SHOW

      friend class Matrix<Symmetric, AT>;
      
//      template <class T> friend Matrix<FixedSize, T>  trans(const Matrix<FixedSize, T> &A);

    protected:

      AT ele[M*N];

      template <class Type> void deepCopy(const Matrix<Type, AT> &A); 
      void deepCopy(const Matrix<FixedSize<M,N>, AT> &A); 

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      Matrix() {
#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      Matrix(Initialization ini, const AT &a=0) {  

	if(ini == INIT) {
	  for(int i=0; i<M*N; i++) 
	    ele[i] = a;
	} else if(ini == EYE ) {
	  for(int i=0; i<M; i++) {
	    for(int j=0; j<N; j++) {
	      if (i==j) ele[i+j*M] = 1; // operator()(i,i) = 1;
	      else ele[i+j*M] = 0; // operator()(i,j) = 0;
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
      Matrix(const Matrix<FixedSize<M,N>, AT> &A) {
	  for(int i=0; i<M*N; i++) 
	    ele[i] = A.ele[i];
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<FixedSize,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const char *str);

      /*! \brief Destructor. 
       * */
      ~Matrix() {
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      Matrix<FixedSize<M,N>, AT>& operator=(const Matrix<FixedSize<M,N>, AT> &A) {
	return operator<<(A);
      }

      template <class Type>
      Matrix<FixedSize<M,N>, AT>& operator=(const Matrix<Type, AT> &A) {
	return operator<<(A);
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> Matrix<FixedSize<M,N>, AT>& operator<<(const Matrix<T, AT> &A);

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      Matrix<FixedSize<M,N>, AT>& operator<<(const Matrix<FixedSize<M,N>, AT> &A);

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

	return ele[i+j*M];
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

	return ele[i+j*M];//  return ele[i*lda+j*ldb];
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

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      Matrix<FixedSize<M,N>, AT> copy() const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      Matrix<FixedSize<M,N>, AT>& init(const AT &a);

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

      Matrix<FixedSize<M,N>, AT> T() {
	Matrix<FixedSize<N,M>, AT> A(NONINIT);
	for(int i=0; i<N; i++)
	  for(int j=0; j<M; j++)
	    A.ele[i+j*M] = ele[j+i*M];
	return A;
      }

      const Matrix<FixedSize<M,N>, AT> T() const {
	Matrix<FixedSize<N,M>, AT> A(NONINIT);
	for(int i=0; i<N; i++)
	  for(int j=0; j<M; j++)
	    A.ele[i+j*M] = ele[j+i*M];
	return A;
      }
  };
  // ------------------------- Constructors -------------------------------------
  template <int M, int N, class AT> 
    Matrix<FixedSize<M,N>, AT>::Matrix(const char *strs) {
    // if 'strs' is a single scalar, surround it first with '[' and ']'.
    // This is more Matlab-like, because e.g. '5' and '[5]' is just the same.
    // (This functionallitiy is needed e.g. by MBXMLUtils (OpenMBV,MBSim))
    std::istringstream iss(strs);
    char c;
    iss>>c;
    if(c=='[') iss.str(strs);
    else iss.str(std::string("[")+strs+"]");

    int m = 0,n=0;
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
	m++;
      }
      else if(c==',')
	n++;
      c='0';
    } while(iss);

    n++; m++;
    assert(m==M);
    assert(n==N);
    iss.clear();
    iss.seekg(0);
    iss >> c;
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++) {
	iss >> ele[i+j*M];
	iss >> c;
      }
  }
  // ----------------------------------------------------------------------------

   template <int M, int N, class AT> template< class Type>
    Matrix<FixedSize<M,N>, AT>& Matrix<FixedSize<M,N>, AT>::operator<<(const Matrix<Type, AT> &A) { 

      if(A.rows() == 0 || A.cols() == 0)
	return *this;

#ifdef FMATVEC_SIZE_CHECK
      if(A.rows() != M || A.cols() != N)
	throw;
      //assert(M == A.rows()); Funktioniert nicht
      //assert(N == A.cols());
#endif

      deepCopy(A);

      return *this;
    }

   template <int M, int N, class AT>
    Matrix<FixedSize<M,N>, AT>& Matrix<FixedSize<M,N>, AT>::operator<<(const Matrix<FixedSize<M,N>, AT> &A) { 

      deepCopy(A);

      return *this;
    }

  template <int M, int N, class AT>
    Matrix<FixedSize<M,N>, AT>&  Matrix<FixedSize<M,N>, AT>::init(const AT& val) {

      for(int i=0; i<M*N; i++) 
	ele[i] = val;

      return *this;
    }

  template <int M, int N, class AT>
    Matrix<FixedSize<M,N>, AT> Matrix<FixedSize<M,N>, AT>::copy() const {

      Matrix<FixedSize<M,N>, AT> A(NONINIT);
      A.deepCopy(*this);

      return A;
    }

  template <int M, int N, class AT> template <class Type>
    void Matrix<FixedSize<M,N>, AT>::deepCopy(const Matrix<Type, AT> &A) { 
      for(int i=0; i<M; i++) 
	for(int j=0; j<N; j++)
          ele[i+j*M] = A.operator()(i,j);
    }

//  template <class AT>
//    Matrix<FixedSize<M,N>, AT>::operator std::vector<std::vector<AT> >() {
//      std::vector<std::vector<AT> > ret(rows());
//      for(int r=0; r<rows(); r++) {
//        ret[r].resize(cols());
//        for(int c=0; c<cols(); c++)
//          ret[r][c]=operator()(r,c);
//      }
//      return ret;
//    }
//
//  template <class AT>
//    Matrix<FixedSize<M,N>, AT>::Matrix(std::vector<std::vector<AT> > m) : memory(m.size()*m[0].size()), ele((AT*)memory.get()), m(m.size()), n(m[0].size()), lda(m.size()), tp(false) {
//#ifndef FMATVEC_NO_INITIALIZATION 
//      init(0);
//#endif
//      for(int r=0; r<rows(); r++) {
//        assert(m[r].size()==cols());
//        for(int c=0; c<cols(); c++)
//          operator()(r,c)=m[r][c];
//      }
//    }

   /// @cond NO_SHOW
   
  template<int M, int N, class AT>
    void Matrix<FixedSize<M,N>,AT>::deepCopy(const Matrix<FixedSize<M,N>,AT> &A) {
      for(int i=0; i<M*N; i++) 
	ele[i] = A.ele[i];
    }

   /// @endcond

  template <int M, int N, class AT>
  inline Matrix<FixedSize<M,N>, AT> operator+(const Matrix<FixedSize<M,N>, AT> &A1, const Matrix<FixedSize<M,N>, AT> &A2) {
    Matrix<FixedSize<M,N>, AT> A3;
    for(int i=0; i<M*N; i++) 
      A3()[i] = A1()[i] + A2()[i];
    return A3;
  }

  template <int M, int N, class AT>
  inline Matrix<FixedSize<M,N>, AT> trans(const Matrix<FixedSize<M,N>, AT> &A) {
    Matrix<FixedSize<M,N>, AT> B(NONINIT);
    for(int i=0; i<N; i++)
      for(int j=0; j<M; j++)
	B()[i+j*M] = A()[j+i*M];
    return B;
  }

}

#endif
