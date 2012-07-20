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

      /*! \brief Standard constructor
       *
       * Constructs a squarematrix with no size. 
       * */
      SquareMatrix() : Matrix<General,Var,Var,AT>() { }

//      template<class Ini=All<AT> >
//      SquareMatrix(int m, Ini ini=All<AT>()) : Matrix<General,Var,Var,AT>(m,m,ini) { } 

      SquareMatrix(int m) : Matrix<General,Var,Var,AT>(m,m) { } 
      SquareMatrix(int m, const Noinit &ini) : Matrix<General,Var,Var,AT>(m,m,ini) { } 
      SquareMatrix(int m, const All<AT> &ini) : Matrix<General,Var,Var,AT>(m,m,ini) { } 
      SquareMatrix(int m, const Eye<AT> &ini) : Matrix<General,Var,Var,AT>(m,m,ini) { } 

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      template<class Type, class Row, class Col>
      explicit SquareMatrix(const Matrix<Type,Row,Col,AT>&  A) : Matrix<General,Var,Var,AT>(A) {
      }

      template<class Row>
      SquareMatrix(const SquareMatrix<Row,AT> &A) : Matrix<General,Var,Var,AT>(A) {
      }

//      template<class Ini=All<AT> >
//        SquareMatrix<Var,AT>& resize(int m=0, Ini ini=All<AT>()) {
//          Matrix<General,Var,Var,AT>::resize(m,m,ini);
//          return *this;
//        }

      SquareMatrix<Var,AT>& resize(int m=0) { return resize(m,All<AT>()); }

      template<class Ini>
        SquareMatrix<Var,AT>& resize(int m, const Ini &ini) {
          Matrix<General,Var,Var,AT>::resize(m,m,ini);
          return *this;
        }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template <class Type, class Row, class Col> SquareMatrix<Var,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
	Matrix<General,Var,Var,AT>::operator=(A);
	return *this;
      }

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix.
       * */
      int size() const {return M;};

      using Matrix<General,Var,Var,AT>::operator();
      using Matrix<General,Var,Var,AT>::e;

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >();

      inline const SquareMatrix<Var,AT> T() const;

  };

  template <class AT>
    inline const SquareMatrix<Var,AT> SquareMatrix<Var,AT>::T() const {
      SquareMatrix<Var,AT> A(NONINIT);
      for(int i=0; i<M; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <class AT>
    inline SquareMatrix<Var,AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(size());
      for(int r=0; r<size(); r++) {
        ret[r].resize(size());
        for(int c=0; c<size(); c++)
          ret[r][c]= e(r,c);
      }
      return ret;
    }

}

#endif
