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
#include <iostream>
#include <iomanip>

using namespace std;

/*! 
 * \brief Namespace fmatvec.
 *
 * */
namespace fmatvec {

  /*! Enumerate for initialization of matrices
   */
  enum Initialization{INIT,NONINIT};

  /*! 
   *  \brief This is the basic matrix class for arbitrary matrices.
   *
   * Template class Matrix with shape type ST and atomic type AT. The first template parameter defines the shape type, the second parameter the
   * atomic type of the matrix. Valid shape types are General, Symmetric, GeneralBand, Diagonal and Sparse. Valid atomic types are int, float, double, complex<float> and complex<double> 
   * */
  template <class ST, class AT>
    class Matrix {

      protected:
	AT* ele;
	int m, n;

      public:
	/*! \brief Standard constructor
	 *
	 * */
	Matrix(int m, int n) {};
	void deepCopy(const Matrix<ST, AT> &A);
	const AT& operator()(int i, int j) const;
	AT& operator()(int i, int j);

	int rows() const {return m;};
	int cols() const {return n;};
    };

  template <class ST, class AT>
    void Matrix<ST, AT>::deepCopy(const Matrix<ST, AT> &A) { 
      for(int i=0; i<m; i++) 
	for(int j=0; j<n; j++) 
	  operator()(i,j) = A(i,j);
    }

  template <class ST, class AT>
    ostream& operator<<(ostream &os,const Matrix<ST, AT> &A) {
      os << A.rows() << " x " << A.cols() << endl;
      os << "[ ";
      for (int i=0; i < A.rows(); ++i) {
	for (int j=0; j < A.cols(); ++j) 
	  os << setw(14) << A(i,j);

	if (i != A.rows() - 1)
	  os << endl  << "  ";
      }
      os << " ]";
      return os;
    }

  template <class ST, class AT>
    istream& operator>>(istream &is, Matrix<ST, AT> &A) {

      int m, n;
      char c;
      is >> m >> c >> n >> c;
      A.resize(m,n,NONINIT);
      for (int i=0; i < A.rows(); ++i) 
	for (int j=0; j < A.cols(); ++j) 
	  is >> A(i,j);
      return is;
    }

}

#endif
