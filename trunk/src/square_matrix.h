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

#ifndef square_matrix_h
#define square_matrix_h

#include "general_matrix.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class of general quadratic
   *  matrices.
   *
   * Template class SquareMatrix with shape type General and
   * atomic type AT. The storage form is dense. The template
   * parameter AT defines the atomic type of the vector. Valid
   * types are int, float, double, complex<float> and
   * complex<double> */
  template <class AT> class SquareMatrix : public Matrix<General, AT> {

    using Matrix<General, AT>::m;
    using Matrix<General, AT>::n;
    using Matrix<General, AT>::lda;
    using Matrix<General, AT>::ele;
    using Matrix<General, AT>::tp;
    using Matrix<General, AT>::memory;
    using Matrix<General, AT>::elePtr;

    public:

    /// @cond NO_SHOW

    template <class T> friend const SquareMatrix<T> trans(const SquareMatrix<T> &A);
    template <class T> friend SquareMatrix<T> trans(SquareMatrix<T> &A);

    friend const SquareMatrix<AT> Matrix<General, AT>::operator()(const Index &I) const;
    friend SquareMatrix<AT> Matrix<General, AT>::operator()(const Index &I);

    protected:

    SquareMatrix(int n, int lda, int tp, Memory<AT> memory, const AT* ele) : Matrix<General, AT>(n, n, lda, tp, memory, ele) {
    }

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a squarematrix with no size. 
       * */
      SquareMatrix() : Matrix<General, AT>() {
      }

     /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x m.
       * \param m The number of rows and columns.
       * \remark The matrix will be initialised to
       * zero by default. This default behavior can be changed by defining 
       * FMATVEC_NO_INITIALIZATION.
       * */
      SquareMatrix(int m) : Matrix<General, AT>(m,m) {
      }

     /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x n, where m must be equal to n.
       * \param m The number of rows and columns.
       * \param n The number of rows and columns.
       * \remark The matrix will be initialised to
       * zero by default. This default behavior can be changed by defining 
       * FMATVEC_NO_INITIALIZATION.
       * */
      SquareMatrix(int m, int n) : Matrix<General, AT>(m,m) {
	assert(m == n);
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x m with the pyhsical
       * memory given by \em ele_;
       * \param m The number of rows and columns.
       * \param ele The physical memory the matrix will point to.
       * */
      SquareMatrix(int m, AT* ele) : Matrix<General, AT>(m,m,ele) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x m. The matrix will be 
       * initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param m The number of rows and columns.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * */
      SquareMatrix(int m, Initialization ini, const AT &a=0) : Matrix<General, AT>(m,m,ini,a) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x n, where m must be equal to n.
       * The matrix will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param m The number of rows and columns.
       * \param n The number of rows and columns.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * */
      SquareMatrix(int m, int n, Initialization ini, const AT &a=0) : Matrix<General, AT>(m,m,ini,a) {
	assert(m == n);
      }

      /*! \brief Matrix resizing. 
       *
       * Resizes the matrix to size n x n.    
       * \param n The number of rows and columns.
       * \return A reference to the calling matrix.
       * \remark The matrix will be initialised to
       * zero by default. To change this behavior, define
       * FMATVEC_NO_INITIALIZATION.
       * */
      SquareMatrix<AT>& resize(int n) {
	Matrix<General, AT>::resize(n,n);
	return *this;
      }

      /*! \brief Matrix resizing. 
       *
       * Resizes the matrix to size n x n. The matrix will be initialized
       * to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param n The number of rows and columns.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * \return A reference to the calling matrix.
       * */
      SquareMatrix<AT>& resize(int n,Initialization ini, const AT &a=0) {
	Matrix<General, AT>::resize(n,n,ini,a);
	return *this;
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      SquareMatrix(const SquareMatrix<AT>&  A) : Matrix<General, AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      explicit SquareMatrix(const Matrix<General, AT>&  A) : Matrix<General, AT>(A) {
#ifdef FMATVEC_SIZE_CHECK
	assert(A.rows() == A.cols());
#endif
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      SquareMatrix<AT>& operator=(const SquareMatrix<AT>&  A) {
#ifndef FMATVEC_NO_DEEP_ASSIGNMENT 
	return operator<<(A);
#else
	return operator>>(A);
#endif
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      SquareMatrix<AT>& operator<<(const SquareMatrix<AT>&  A) {
	Matrix<General,AT>::operator<<(A);
	return *this;
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> SquareMatrix<AT>& operator<<(const Matrix<T, AT> &A) {
	Matrix<General,AT>::operator<<(A);
	return *this;
      }

      /*! \brief Reference operator
       *
       * References the matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      SquareMatrix<AT>& operator>>(const SquareMatrix<AT>&  A) {
	Matrix<General,AT>::operator>>(A);
	return *this;
      }

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix.
       * */
      int size() const {return m;};

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      SquareMatrix<AT> copy() const;

      using Matrix<General, AT>::operator();
      using Matrix<General, AT>::resize;

  };

  template <class AT>
    SquareMatrix<AT> SquareMatrix<AT>::copy() const {

      SquareMatrix<AT> A(m);
      A.deepCopy(*this);

      return A;
    }


}

#endif
