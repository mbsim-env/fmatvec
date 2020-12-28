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

#ifndef var_square_matrix_h
#define var_square_matrix_h

#include "fixed_general_matrix.h"

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
  template <class AT> class SquareMatrix<Var,AT> : public Matrix<General,Var,Var,AT> {

    using Matrix<General,Var,Var,AT>::M;

    public:
      static constexpr bool isVector {false};

      typedef AT value_type;
      typedef Square shape_type;

      /*! \brief Standard constructor
       *
       * Constructs a squarematrix with no size. 
       * */
      explicit SquareMatrix() : Matrix<General,Var,Var,AT>() { }

//      template<class Ini=All<AT>>
//      SquareMatrix(int m, Ini ini=All<AT>()) : Matrix<General,Var,Var,AT>(m,m,ini) { } 

      explicit SquareMatrix(int m, Noinit ini) : Matrix<General,Var,Var,AT>(m,m,ini) { } 
      explicit SquareMatrix(int m, Init ini=INIT, const AT &a=AT()) : Matrix<General,Var,Var,AT>(m,m,ini,a) { } 
      explicit SquareMatrix(int m, Eye ini, const AT &a=1) : Matrix<General,Var,Var,AT>(m,m,ini,a) { } 

      // move
      SquareMatrix(SquareMatrix<Var,AT> &&src) : Matrix<General,Var,Var,AT>(std::move(src)) {}
      SquareMatrix<Var,AT>& operator=(SquareMatrix<Var,AT> &&src) {
        M=src.M;
        src.M=0;
        this->N=src.N;
        src.N=0;
        delete[]this->ele;
        this->ele=src.ele;
        src.ele=nullptr;
        return *this;
      }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      SquareMatrix(const SquareMatrix<Var,AT> &A) : Matrix<General,Var,Var,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      template<class Row>
      SquareMatrix(const SquareMatrix<Row,AT> &A) : Matrix<General,Var,Var,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      template<class Type, class Row, class Col>
      explicit SquareMatrix(const Matrix<Type,Row,Col,AT>&  A) : Matrix<General,Var,Var,AT>(A) {
	FMATVEC_ASSERT(A.rows() == A.cols(), AT); 
      }

      SquareMatrix(const char *str) : Matrix<General,Var,Var,AT>(str) { }

      SquareMatrix<Var,AT>& resize(int m, Noinit) {
        Matrix<General,Var,Var,AT>::resize(m,m,Noinit());
        return *this;
      }

      SquareMatrix<Var,AT>& resize(int m, Init ini=INIT, const AT &a=AT()) {
        Matrix<General,Var,Var,AT>::resize(m,m,ini,a);
        return *this;
      }

      SquareMatrix<Var,AT>& resize(int m, Eye ini, const AT &a=1) {
        Matrix<General,Var,Var,AT>::resize(m,m,ini,a);
        return *this;
      }

      //! Resize a var square matrix.
      //! Throw if the dimensions does not match or resize to this dimension.
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
      SquareMatrix<Var,AT>& operator=(const SquareMatrix<Var,AT> &A) {
	Matrix<General,Var,Var,AT>::operator=(A);
	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      template <class Type, class Row, class Col> SquareMatrix<Var,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
	Matrix<General,Var,Var,AT>::operator=(A);
	return *this;
      }

      /*! \brief Matrix assignment
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col> SquareMatrix<Var,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
	Matrix<General,Var,Var,AT>::operator<<=(A);
	return *this;
      }
      // move
      SquareMatrix<Var,AT>& operator<<=(SquareMatrix<Var,AT> &&src) {
        FMATVEC_ASSERT(src.rows() == src.cols(), AT);
	Matrix<General,Var,Var,AT>::operator<<=(std::move(src));
	return *this;
      }

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix.
       * */
      int size() const {return M;}
      int rows() const {return M;}
      int cols() const {return M;}

      using Matrix<General,Var,Var,AT>::operator();
      using Matrix<General,Var,Var,AT>::e;

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit inline operator std::vector<std::vector<AT>>() const;

      inline const SquareMatrix<Var,AT> T() const;
  };

  template <class AT>
    inline const SquareMatrix<Var,AT> SquareMatrix<Var,AT>::T() const {
      SquareMatrix<Var,AT> A(size(),NONINIT);
      for(int i=0; i<M; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <class AT>
    inline SquareMatrix<Var,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(size());
      for(int r=0; r<size(); r++) {
        ret[r].resize(size());
        for(int c=0; c<size(); c++)
          ret[r][c]= e(r,c);
      }
      return ret;
    }

}

#endif
