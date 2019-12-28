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

#ifndef matrix_h
#define matrix_h

#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <vector>
#include <boost/scope_exit.hpp>
#include "range.h"
#include "toString.h"
#include "types.h"

/*! 
 * \brief Namespace fmatvec.
 *
 * */
namespace fmatvec {

  /*! Enumerate for initialization of matrices
  */
  //enum Initialization{INIT,NONINIT,EYE};

  /*! 
   *  \brief This is the basic matrix class for arbitrary matrices.
   *
   * Template class Matrix with shape type ST and atomic type AT. The first template parameter defines the shape type, the second parameter the
   * atomic type of the matrix. Valid shape types are General, Symmetric, GeneralBand, Diagonal and Sparse. Valid atomic types are int, float, double, complex<float> and complex<double> 
   * */
  template <class Type, class Row, class Col, class AT> class Matrix {

    protected:

      /// @cond NO_SHOW

      //AT* ele;

      inline Matrix<Type,Row,Col,AT>& copy(const Matrix<Type,Row,Col,AT> &A);

      /// @endcond

    public:
      typedef AT value_type;
//      /*! \brief Standard constructor
//       * 
//       * The standard constructor.
//       * \param m The number of rows
//       * \param n The number of columns
//       * */
//      Matrix(int m, int n) {};

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \sa operator()(int,int) const
       * */
      AT& operator()(int i, int j) {
	assert(i>=0);
	assert(j>=0);
	assert(i<rows());
	assert(j<cols());

	return e(i,j);
      }

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
	assert(i>=0);
	assert(j>=0);
	assert(i<rows());
	assert(j<cols());

	return e(i,j);
      }

      AT& e(int i, int j);

      const AT& e(int i, int j) const;

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix.
       * */
      int rows() const;

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix.
       * */
      int cols() const;

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      operator std::vector<std::vector<AT>>() const;
  };

  template <class Type, class Row, class Col, class AT> inline Matrix<Type,Row,Col,AT>& Matrix<Type,Row,Col,AT>::copy(const Matrix<Type,Row,Col,AT> &A) {
    for(int i=0; i<rows(); i++) 
      for(int j=0; j<cols(); j++) 
	e(i,j) = A.e(i,j);
    return *this;
  }

  /*! \brief Matrix output 
   *
   * This function writes a matrix into a stream.
   * \param os An output stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the output stream.
   * */
  template <class Type, class Row, class Col, class AT> std::ostream& operator<<(std::ostream &os, const Matrix<Type,Row,Col,AT> &A) {
    os << toString(static_cast<std::vector<std::vector<AT>>>(A), os.precision());
    return os;
  }

  /*! \brief Matrix input
   *
   * This function loads a matrix from a stream.
   * \param is An input stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the input stream.
   * */
  template <class Type, class Row, class Col, class AT> std::istream& operator>>(std::istream &is, Matrix<Type,Row,Col,AT> &A) {
    std::ios_base::iostate oldEx=is.exceptions();
    std::ios_base::fmtflags oldFlags=is.flags();
    BOOST_SCOPE_EXIT_TPL(&is, oldEx, oldFlags) {
      is.exceptions(oldEx);
      is.flags(oldFlags);
    } BOOST_SCOPE_EXIT_END
    is.exceptions(std::ios::failbit | std::ios::badbit);
    is.flags(std::ios_base::skipws);
    char c;
    AT e;
    while((c=is.peek())==' ' || c=='\n') is.get();
    if(c!='[') {
      // its a scalar (written without the [ ] )
      is>>e;
      A<<=Matrix<Type,Row,Col,AT>(std::vector<std::vector<AT>>(1, std::vector<AT>(1, e)));
      return is;
    }
    // its vector or matrix
    is.get();
    std::vector<std::vector<AT>> m(1);
    int r=0;
    while(true) {
      is>>e;
      m[r].push_back(e);
      while((c=is.peek())==' ') is.get();
      if     (c==','           ) { is.get(); }
      else if(c==';' || c=='\n') { is.get(); m.resize(++r + 1); }
      else if(c==']'           ) { is.get(); break; }
    }
    A<<=Matrix<Type,Row,Col,AT>(m);
    return is;
  }

  /*! \brief Matrix dump
   *
   * This function dumps a matrix to a file in ASCII-format.
   * \param is The filename.
   * \param A A matrix of any shape and type.
   * */
  template <class Type, class Row, class Col, class AT> void dump(const char* str, const Matrix<Type,Row,Col,AT> &A) {
    std::ofstream os(str);
    for (int i=0; i < A.rows(); ++i) {
      for (int j=0; j < A.cols(); ++j) 
	os << std::setw(14) << A.e(i,j);

      if (i != A.rows() - 1)
	os << std::endl;
    }
    os.close();
  }

  template <class Type, class Row, class Col, class AT>
    Matrix<Type,Row,Col,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows());
      for(int r=0; r<rows(); r++) {
	ret[r].resize(cols());
	for(int c=0; c<cols(); c++)
	  ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <class Row, class AT> class SquareMatrix {
  };

  template <class Row, class AT> class Vector {
  };

  template <class Col, class AT> class RowVector {
  };
}

#endif

