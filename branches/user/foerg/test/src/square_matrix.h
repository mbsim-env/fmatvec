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
  template <class AT> class SquareMatrix<Ref,AT> : public Matrix<General,Ref,Ref,AT> {

    using Matrix<General,Ref,Ref,AT>::m;
    using Matrix<General,Ref,Ref,AT>::n;
    using Matrix<General,Ref,Ref,AT>::lda;
    using Matrix<General,Ref,Ref,AT>::ele;
    using Matrix<General,Ref,Ref,AT>::tp;
    using Matrix<General,Ref,Ref,AT>::memory;
    using Matrix<General,Ref,Ref,AT>::elePtr;

    public:

    /// @cond NO_SHOW

    template <class T> friend SquareMatrix<Ref,T> trans(const SquareMatrix<Ref,T> &A);

    friend const SquareMatrix<Ref,AT> Matrix<General,Ref,Ref,AT>::operator()(const Index &I) const;
    friend SquareMatrix<Ref,AT> Matrix<General,Ref,Ref,AT>::operator()(const Index &I);

    protected:

    SquareMatrix(int n, int lda, int tp, Memory<AT> memory, const AT* ele) : Matrix<General,Ref,Ref,AT>(n, n, lda, tp, memory, ele) {
    }

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a squarematrix with no size. 
       * */
      SquareMatrix() : Matrix<General,Ref,Ref,AT>() {
      }

     /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x m.
       * \param m The number of rows and columns.
       * \remark The matrix will be initialised to
       * zero by default. This default behavior can be changed by defining 
       * FMATVEC_NO_INITIALIZATION.
       * */
      SquareMatrix(int m) : Matrix<General,Ref,Ref,AT>(m,m) {
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
      SquareMatrix(int m, int n) : Matrix<General,Ref,Ref,AT>(m,m) {
	assert(m == n);
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x m with the pyhsical
       * memory given by \em ele_;
       * \param m The number of rows and columns.
       * \param ele The physical memory the matrix will point to.
       * */
      SquareMatrix(int m, AT* ele) : Matrix<General,Ref,Ref,AT>(m,m,ele) {
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
      SquareMatrix(int m, Initialization ini, const AT &a=0) : Matrix<General,Ref,Ref,AT>(m,m,ini,a) {
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
      SquareMatrix(int m, int n, Initialization ini, const AT &a=0) : Matrix<General,Ref,Ref,AT>(m,m,ini,a) {
	assert(m == n);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      SquareMatrix(const SquareMatrix<Ref,AT>&  A) : Matrix<General,Ref,Ref,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<General,Ref,Ref,AT>&) 
       * */
      explicit SquareMatrix(const Matrix<General,Ref,Ref,AT>&  A) : Matrix<General,Ref,Ref,AT>(A) {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == A.cols());
#endif
      }

     template<class Type, class Row, class Col>
      explicit SquareMatrix(const Matrix<Type,Row,Col,AT> &x) : Matrix<General,Ref,Ref,AT>(x)  {
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
      SquareMatrix<Ref,AT>& resize(int n=0,Initialization ini=INIT, const AT &a=0) {
	Matrix<General,Ref,Ref,AT>::resize(n,n,ini,a);
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
      SquareMatrix<Ref,AT>& operator=(const SquareMatrix<Ref,AT>&  A) {
	Matrix<General,Ref,Ref,AT>::operator=(A);
	return *this;
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col> SquareMatrix<Ref,AT>& operator<<(const Matrix<Type,Row,Col,AT> &A) {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(A.rows() == A.cols());
#endif
	Matrix<General,Ref,Ref,AT>::operator<<(A);
	return *this;
      }

      /*! \brief Reference operator
       *
       * References the matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      SquareMatrix<Ref,AT>& operator>>(const SquareMatrix<Ref,AT>&  A) {
	Matrix<General,Ref,Ref,AT>::operator>>(A);
	return *this;
      }

      SquareMatrix<Ref,AT>& operator>>(const Matrix<General,Ref,Ref,AT>&  A) {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(A.rows() == A.cols());
#endif
	Matrix<General,Ref,Ref,AT>::operator>>(A);
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
      inline SquareMatrix<Ref,AT> copy() const;

      using Matrix<General,Ref,Ref,AT>::operator();

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >();

      SquareMatrix<Ref,AT> T() {
	return SquareMatrix<Ref,AT>(n, lda, tp?false:true, memory, ele);
      }

      const SquareMatrix<Ref,AT> T() const {
	return SquareMatrix<Ref,AT>(n, lda, tp?false:true, memory, ele);
      }
  };

  template <class AT>
    inline SquareMatrix<Ref,AT> SquareMatrix<Ref,AT>::copy() const {

      SquareMatrix<Ref,AT> A(m,NONINIT);
      A.deepCopy(*this);

      return A;
    }

  template <class AT>
    inline SquareMatrix<Ref,AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(size());
      if(tp) {
	for(int r=0; r<size(); r++) {
	  ret[r].resize(size());
	  for(int c=0; c<size(); c++)
	    ret[r][c]= this->et(r,c); 
	}
      }
      else {
	for(int r=0; r<size(); r++) {
	  ret[r].resize(size());
	  for(int c=0; c<size(); c++)
	    ret[r][c]= this->er(r,c); 
	}
      }
      return ret;
    }

}

#endif
