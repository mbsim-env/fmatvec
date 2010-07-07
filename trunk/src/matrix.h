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

#ifndef matrix_h
#define matrix_h

#define FMATVEC_SIZE_CHECK
#define FMATVEC_VOID_CHECK

#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

/*! 
 * \brief Namespace fmatvec.
 *
 * */
namespace fmatvec {

  /*! Enumerate for initialization of matrices
  */
  enum Initialization{INIT,NONINIT,EYE};

  /*! 
   *  \brief This is the basic matrix class for arbitrary matrices.
   *
   * Template class Matrix with shape type ST and atomic type AT. The first template parameter defines the shape type, the second parameter the
   * atomic type of the matrix. Valid shape types are General, Symmetric, GeneralBand, Diagonal and Sparse. Valid atomic types are int, float, double, complex<float> and complex<double> 
   * */
  template <class ST, class AT> class Matrix {

    protected:

      /// @cond NO_SHOW

      AT* ele;
      int m, n;

      void deepCopy(const Matrix<ST, AT> &A);

      /// @endcond

    public:
      /*! \brief Standard constructor
       * 
       * The standard constructor.
       * \param m The number of rows
       * \param n The number of columns
       * */
      Matrix(int m, int n) {};

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \sa operator()(int,int) const
       * */
      AT& operator()(int i, int j);

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const;

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix.
       * */
      int rows() const {return m;};

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix.
       * */
      int cols() const {return n;};

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      operator std::vector<std::vector<AT> >();
  };



  template <class ST, class AT> void Matrix<ST, AT>::deepCopy(const Matrix<ST, AT> &A) { 
    for(int i=0; i<m; i++) 
      for(int j=0; j<n; j++) 
	operator()(i,j) = A(i,j);
  }

  /*! \brief Matrix output 
   *
   * This function writes a matrix into a stream.
   * \param os An output stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the output stream.
   * */
  template <class ST, class AT> std::ostream& operator<<(std::ostream &os, const Matrix<ST, AT> &A) {
    os << A.rows() << " x " << A.cols() << std::endl;
    os << "[ ";
    for (int i=0; i < A.rows(); ++i) {
      for (int j=0; j < A.cols(); ++j) 
	os << std::setw(14) << A(i,j);

      if (i != A.rows() - 1)
	os << std::endl  << "  ";
    }
    os << " ]";
    return os;
  }

  /*! \brief Matrix input
   *
   * This function loads a matrix from a stream.
   * \param is An input stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the input stream.
   * */
  template <class ST, class AT> std::istream& operator>>(std::istream &is, Matrix<ST, AT> &A) {
    int m, n;
    char c;
    is >> m >> c >> n >> c;
    A.resize(m,n,NONINIT);
    for (int i=0; i < A.rows(); ++i) 
      for (int j=0; j < A.cols(); ++j) 
	is >> A(i,j);
    is >> c;
    return is;
  }

  /*! \brief Matrix dump
   *
   * This function dumps a matrix to a file in ASCII-format.
   * \param is The filename.
   * \param A A matrix of any shape and type.
   * */
  template <class ST, class AT> void dump(const char* str, const Matrix<ST, AT> &A) {
    std::ofstream os(str);
    for (int i=0; i < A.rows(); ++i) {
      for (int j=0; j < A.cols(); ++j) 
	os << std::setw(14) << A(i,j);

      if (i != A.rows() - 1)
	os << std::endl;
    }
    os.close();
  }

  template <class ST, class AT>
    Matrix<ST, AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
	ret[r].resize(cols());
	for(int c=0; c<cols(); c++)
	  ret[r][c]=operator()(r,c);
      }
      return ret;
    }

}

#endif

