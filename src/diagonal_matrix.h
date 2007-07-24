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

#ifndef diagonal_matrix_h
#define diagonal_matrix_h

#include "matrix.h"
#include "types.h"
#include "memory.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for diagonal matrices.
   *
   * Template class Matrix of shape type Diagonal.
   * The template parameter AT defines the
   * atomic type of the matrix. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class Matrix<Diagonal, AT> {

      protected:

    /// @cond NO_SHOW

	Memory<AT> memory;
	AT *ele;
	int n;

	void deepCopy(const Matrix<Diagonal, AT> &x);

	Matrix(int n_, Memory<AT> memory_, const AT* ele_) : memory(memory_), ele((AT*)ele_), n(n_) {
	}

	const AT* elePtr(int i) const {
	  return ele+i;
	}

	AT* elePtr(int i) {
	  return ele+i;
	}

    /// @endcond

      public:

	/*! \brief Standard constructor
	 *
	 * Construckts a matrix with no size.
	 * */
	Matrix() : memory(), ele(0), n(0) {
	}

	/*! \brief Regular Constructor
	 *
	 * Constructs a diagonal matrix of size n x n.
	 * \param n_ The number of rows and columns.
	 * \remark The matrix will be initialised to
	 * zero by default. This default behavior can be changed by defining 
	 * FMATVEC_NO_INITIALIZATION.
	 * */
	Matrix(int n_) : memory(n_), ele((AT*)memory.get()), n(n_) {  

#ifndef FMATVEC_NO_INITIALIZATION 
	  init(0);
#endif
	}

	/*! \brief Matrix resizing. 
	 *
	 * Resizes the matrix to size n x n.    
	 * \param n_ The number of rows and columns.
	 * \return A reference to the calling matrix.
	 * \remark The matrix will be initialised to
	 * zero by default. To change this behavior, define
	 * FMATVEC_NO_INITIALIZATION.
	 * */
	Matrix<Diagonal, AT>& resize(int n_) {
	  n = n_;

	  memory.resize(n);
	  ele = (AT*)memory.get();

#ifndef FMATVEC_NO_INITIALIZATION 
	  init(0);
#endif
	  return *this;
	}

	/*! \brief Matrix resizing. 
	 *
	 * Resizes the matrix to size n x n. The matrix will be initialized to 
	 * the value given by \em a
	 * (default 0), if ini is set to INIT. If init is set to NONINIT, the
	 * matrix will not be initialized.
	 * \param n_ The number of rows and columns.
	 * \param ini INIT means initialization, NONINIT means no initialization.
	 * \param a The value, the matrix will be initialized with (default 0)
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Diagonal, AT>& resize(int n_, Initialization ini, const AT &a=0) {
	  n = n_;

	  memory.resize(n);
	  ele = (AT*)memory.get();

	  if(ini == INIT)
	    init(a);
	  else if(ini == EYE ) 
	    init(1);

	  return *this;
	}

	/*! \brief Regular Constructor
	 *
	 * Constructs a diagonal matrix of size n x n. The matrix will be initialized to
	 * the value given by \em a
	 * (default 0), if ini is set to INIT. If init is set to NONINIT, the
	 * matrix will not be initialized.
	 * \param n_ The number of rows and columns.
	 * \param ini INIT means initialization, NONINIT means no initialization.
	 * \param a The value, the matrix will be initialized with (default 0)
	 * */
	Matrix(int n_, Initialization ini, const AT& a=0) : memory(n_), ele((AT*)memory.get()), n(n_) {  

	  if(ini == INIT) {
	    AT *el=ele;
	    for(int i=0; i<n; i++)
	      *el++=a;
	  } else if(ini == EYE ) {
	    AT *el=ele;
	    for(int i=0; i<n; i++) 
	      *el++= 1;
	  }  
	}

	/*! \brief Copy Constructor
	 *
	 * Constructs a reference to the diagonal matrix \em A.
	 * \attention The physical memory of the matrix \em A 
	 * will not be copied, only referenced.
	 * \param A The matrix that will be referenced.
	 * */
	Matrix(const Matrix<Diagonal, AT> &A) : memory(A.memory), ele(A.ele) ,n(A.n) {
	}

	/*! \brief Destructor. 
	 * */
	~Matrix() {
	}

	/*! \brief Copy operator
	 *
	 * Copies the diagonal matrix given by \em A.
	 * \param A The matrix to be copied. 
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Diagonal, AT>& operator<<(const Matrix<Diagonal, AT> &A);

	/*! \brief Reference operator
	 *
	 * References the diagoanl matrix given by \em A.
	 * \param A The matrix to be referenced. 
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Diagonal, AT>& operator>>(const Matrix<Diagonal, AT> &A);

	/*! \brief Assignment operator
	 *
	 * Copies the square matrix given by \em A by calling operator<<().
	 * \param A The matrix to be assigned. 
	 * \return A reference to the calling matrix.
	 * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
	 * \sa operator<<(), operator>>()
	 * */
	Matrix<Diagonal, AT>& operator=(const Matrix<Diagonal, AT> &A) {
#ifndef FMATVEC_NO_DEEP_ASSIGNMENT 
	  return operator<<(A);
#else
	  return operator>>(A);
#endif
	}

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
	const AT& operator()(int i, int j) const {

	  static AT zero=0;

#ifndef FMATVEC_NO_BOUNDS_CHECK
	  assert(i>=0);
	  assert(j>=0);
	  assert(i<n);
	  assert(j<n);
#endif

	  return i==j ? ele[i] : zero;
	};

	/*! \brief Element operator
	 *
	 * See operator()(int) 
	 * */
	const AT& operator()(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
	  assert(i>=0);
	  assert(i<n);
#endif

	  return ele[i];
	};

	/*! \brief Element operator
	 *
	 * Returns a reference to the element in the i-th row and the i-th column. 
	 * \param i The i-th row and column of the matrix
	 * \return A reference to the element A(i,i).
	 * \remark The bounds are checked by default. 
	 * To change this behavior, define
	 * FMATVEC_NO_BOUNDS_CHECK.
	 * \sa operator()(int) const
	 * */
	AT& operator()(int i) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
	  assert(i>=0);
	  assert(i<n);
#endif

	  return ele[i];
	};

	/*! \brief Pointer operator
	 *
	 * See operator()() 
	 * */
	const AT* operator()() const {
	  return ele;
	};

	/*! \brief Pointer operator.
	 *
	 * Returns the pointer to the first element.
	 * \return The pointer to the first element.
	 * */
	AT* operator()() {
	  return ele;
	};

	/*! \brief Number of rows.
	 *
	 * \return The number of rows of the matrix
	 * */
	int rows() const {return n;};

	/*! \brief Number of columns.
	 *
	 * \return The number of columns of the matrix
	 * */
	int cols() const {return n;};

	/*! \brief Number of rows and columns.
	 *
	 * \return The number of rows and columns of the matrix
	 * */
	int size() const {return n;};

	/*! \brief Storage convention.
	 *
	 * Returns the blas-conform storage convention. 
	 * The elements are stored * in columnmajor form, 
	 * i.e. the elements are stored columnwise. 
	 * \return CblasColMajor.
	 * */
	const enum CBLAS_ORDER blasOrder() const {
	  return  CblasColMajor;
	};

	/*! \brief Matrix duplicating.
	 *
	 * The calling matrix returns a \em deep copy of itself.  
	 * \return The duplicate.
	 * */
	Matrix<Diagonal, AT> copy() const;

	/*! \brief Initialization.
	 *
	 * Initializes all elements of the calling matrix with 
	 * the value given by \em a.
	 * \param a Value all elements will be initialized with.
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Diagonal, AT>& init(const AT &a);

    };

  template <class AT>
    Matrix<Diagonal, AT>& Matrix<Diagonal, AT>::operator>>(const Matrix<Diagonal, AT> &A) { 
      if(n==0) {
	n=A.n;
      } else 
	assert(n == A.n);

      memory = A.memory;
      ele = A.ele;

      return *this;
    }

  template <class AT>
    Matrix<Diagonal, AT>& Matrix<Diagonal, AT>::operator<<(const Matrix<Diagonal, AT> &A) { 

      if(A.size() == 0)
	return *this;
	
      if(n==0) {
	n = A.n;
	memory.resize(n);
	ele = (AT*)memory.get();
      } else 
	assert(n == A.n);

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<Diagonal, AT>&  Matrix<Diagonal, AT>::init(const AT &val) {

      for(int i=0; i<rows(); i++) 
	operator()(i) = val;

      return *this;
    }

 
  template <class AT>
    Matrix<Diagonal, AT> Matrix<Diagonal, AT>::copy() const {

      Matrix<Diagonal, AT> A(n);

      A.deepCopy(*this);

      return A;
    }

  template <class AT>
    void Matrix<Diagonal, AT>::deepCopy(const Matrix<Diagonal, AT> &A) { 
      for(int i=0; i<n; i++) 
	operator()(i) = A(i);
    }

  template <> extern void Matrix<Diagonal, double>::deepCopy(const Matrix<Diagonal, double> &A);

}

#endif
