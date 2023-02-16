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
#include <stdexcept>

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
    using Matrix<General,Ref,Ref,AT>::memory;
    using Matrix<General,Ref,Ref,AT>::elePtr;

    public:

    static constexpr bool isVector {false};
    using value_type = AT;
    using shape_type = Square;

    public:

      /*! \brief Standard constructor
       *
       * Constructs a squarematrix with no size. 
       * */
      explicit SquareMatrix() : Matrix<General,Ref,Ref,AT>() { }

//      template<class Ini=All<AT>>
//      SquareMatrix(int m, Ini ini=All<AT>()) : Matrix<General,Ref,Ref,AT>(m,m,ini) { } 

      explicit SquareMatrix(int m, Noinit ini) : Matrix<General,Ref,Ref,AT>(m,m,ini) { } 
      explicit SquareMatrix(int m, Init ini=INIT, const AT &a=AT()) : Matrix<General,Ref,Ref,AT>(m,m,ini,a) { } 
      explicit SquareMatrix(int m, Eye ini, const AT &a=1) : Matrix<General,Ref,Ref,AT>(m,m,ini,a) { } 

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will copied.
       * */
      SquareMatrix(const SquareMatrix<Ref,AT>&  A) : Matrix<General,Ref,Ref,AT>(A) { }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will copied.
       * */
      template<class Row>
      SquareMatrix(const SquareMatrix<Row,AT>&  A) : Matrix<General,Ref,Ref,AT>(A) { }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will copied.
       * */
      template<class Type, class Row, class Col>
      explicit SquareMatrix(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Ref,Ref,AT>(A)  {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
      }
     
      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x m with the pyhsical
       * memory given by \em ele_;
       * \param m The number of rows and columns.
       * \param ele The physical memory the matrix will point to.
       * */
      explicit SquareMatrix(int m, AT* ele) : Matrix<General,Ref,Ref,AT>(m,m,ele) { }

      SquareMatrix(const char *str) : Matrix<General,Ref,Ref,AT>(str) { }

      SquareMatrix<Ref,AT>& resize(int m, Noinit) {
        Matrix<General,Ref,Ref,AT>::resize(m,m,Noinit());
        return *this;
      }

      SquareMatrix<Ref,AT>& resize(int m, Init ini=INIT, const AT &a=AT()) {
        Matrix<General,Ref,Ref,AT>::resize(m,m,ini,a);
        return *this;
      }

      SquareMatrix<Ref,AT>& resize(int m, Eye ini, const AT &a=1) {
        Matrix<General,Ref,Ref,AT>::resize(m,m,ini,a);
        return *this;
      }

      //! Resize a square matrix.
      //! Throw if the row and col dimension does not match or resize to this dimension.
      void resize(int m, int n) {
        if(n!=m)
          throw std::runtime_error("Cannot resize a square matrix with different dimensions for rows and columns.");
        resize(m);
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      inline SquareMatrix<Ref,AT>& operator=(const SquareMatrix<Ref,AT> &A) {
	Matrix<General,Ref,Ref,AT>::operator=(A);
	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col> inline SquareMatrix<Ref,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
	Matrix<General,Ref,Ref,AT>::operator=(A);
	return *this;
      }

      /*! \brief Reference operator
       *
       * References the matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      SquareMatrix<Ref,AT>& operator&=(SquareMatrix<Ref,AT>&  A) {
	Matrix<General,Ref,Ref,AT>::operator&=(A);
	return *this;
      }

      SquareMatrix<Ref,AT>& operator&=(Matrix<General,Ref,Ref,AT>&  A) {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
	Matrix<General,Ref,Ref,AT>::operator&=(A);
	return *this;
      }

      /*! \brief Matrix assignment
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col> SquareMatrix<Ref,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
	Matrix<General,Ref,Ref,AT>::operator<<=(A);
	return *this;
      }

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix.
       * */
      int size() const {return m;}
      int rows() const {return m;}
      int cols() const {return m;}

      using Matrix<General,Ref,Ref,AT>::operator();
      using Matrix<General,Ref,Ref,AT>::e;

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit inline operator std::vector<std::vector<AT>>() const;

      /*! \brief std::vector<std::vector<AT>> Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT>> object.
       * An FMATVEC_ASSERT checks for constant length of each ro, UseExceptions<AT>::EXw.
       * \param m The std::vector<std::vector<AT>> the matrix will be initialized with.
       * */
      explicit SquareMatrix(const std::vector<std::vector<AT>> &m);

      const SquareMatrix<Ref,AT> T() const;
  };

  template <class AT>
    inline SquareMatrix<Ref,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows(),std::vector<AT>(cols()));
      for(int r=0; r<size(); r++) {
        for(int c=0; c<size(); c++)
          ret[r][c]= e(r,c);
      }
      return ret;
    }

  template <class AT>
    inline SquareMatrix<Ref,AT>::SquareMatrix(const std::vector<std::vector<AT>> &m) : Matrix<General,Ref,Ref,AT>(static_cast<int>(m.size()),static_cast<int>(m[0].size()),NONINIT) {
      for(int r=0; r<rows(); r++) {
        if(static_cast<int>(m[r].size())!=cols())
          throw std::runtime_error("The rows of the input have different length.");
        for(int c=0; c<cols(); c++)
          e(r,c)=m[r][c];
      }
    }

  template <class AT>
    inline const SquareMatrix<Ref,AT> SquareMatrix<Ref,AT>::T() const {
      SquareMatrix<Ref,AT> A(size(),NONINIT);
      for(int i=0; i<m; i++)
        for(int j=0; j<m; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

}

#endif
