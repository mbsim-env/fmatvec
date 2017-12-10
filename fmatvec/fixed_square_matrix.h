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

#ifndef fixed_square_matrix_h
#define fixed_square_matrix_h

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
  template <int M, class AT> class SquareMatrix<Fixed<M>,AT> : public Matrix<General,Fixed<M>,Fixed<M>,AT> {

    public:
      typedef AT value_type;

//      template<class Ini=All<AT> >
//      SquareMatrix(Ini ini=All<AT>()) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini) { }
//      template<class Ini=All<AT> >
//      SquareMatrix(int m, Ini ini=All<AT>()) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini) { } 

      SquareMatrix(Noinit ini) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini) { }
      SquareMatrix(Init ini=INIT, const AT &a=0) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini,a) { }
      SquareMatrix(Eye ini, const AT &a=1) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini,a) { }
      SquareMatrix(int m, Noinit ini) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini) { }
      SquareMatrix(int m, Init ini=INIT, const AT &a=0) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini,a) { }
      SquareMatrix(int m, Eye ini, const AT &a=1) : Matrix<General,Fixed<M>,Fixed<M>,AT>(ini,a) { }

      SquareMatrix(const char *str) : Matrix<General,Fixed<M>,Fixed<M>,AT>(str) { }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      template<class Type, class Row, class Col>
      explicit SquareMatrix(const Matrix<Type,Row,Col,AT>&  A) : Matrix<General,Fixed<M>,Fixed<M>,AT>(A) {
      }

      template<class Row>
      SquareMatrix(const SquareMatrix<Row,AT> &A) : Matrix<General,Fixed<M>,Fixed<M>,AT>(A) {
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template <class Type, class Row, class Col> SquareMatrix<Fixed<M>,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
	Matrix<General,Fixed<M>,Fixed<M>,AT>::operator=(A);
	return *this;
      }

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix.
       * */
      int size() const {return M;};

      using Matrix<General,Fixed<M>,Fixed<M>,AT>::operator();
      using Matrix<General,Fixed<M>,Fixed<M>,AT>::e;

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >() const;

      inline const SquareMatrix<Fixed<M>,AT> T() const;

  };

  template <int M, class AT>
    inline const SquareMatrix<Fixed<M>,AT> SquareMatrix<Fixed<M>,AT>::T() const {
      SquareMatrix<Fixed<M>,AT> A(NONINIT);
      for(int i=0; i<M; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int M, class AT>
    inline SquareMatrix<Fixed<M>,AT>::operator std::vector<std::vector<AT> >()  const{
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
