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

#include "general_fixed_matrix.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class of general quadratic
   *  matrices.
   *
   * Template class FixedSquareMatrix with shape type General and
   * atomic type AT. The storage form is dense. The template
   * parameter AT defines the atomic type of the vector. Valid
   * types are int, float, double, complex<float> and
   * complex<double> */
  template <int M, class AT> class FixedSquareMatrix : public Matrix<GeneralFixed<M,M>, AT> {

    public:

    public:

      /*! \brief Standard constructor
       *
       * Constructs a squarematrix with no size. 
       * */
      FixedSquareMatrix() : Matrix<GeneralFixed<M,M>, AT>() {
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
      FixedSquareMatrix(Initialization ini, const AT &a=0) : Matrix<GeneralFixed<M,M>, AT>(ini,a) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      FixedSquareMatrix(const FixedSquareMatrix<M,AT>&  A) : Matrix<GeneralFixed<M,M>, AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See SquareMatrix(const SquareMatrix<AT>&) 
       * */
      explicit FixedSquareMatrix(const Matrix<GeneralFixed<M,M>, AT>&  A) : Matrix<GeneralFixed<M,M>, AT>(A) {
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      FixedSquareMatrix<M,AT>& operator=(const FixedSquareMatrix<M,AT>&  A) {
	return operator<<(A);
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      FixedSquareMatrix<M,AT>& operator<<(const FixedSquareMatrix<M,AT>&  A) {
	Matrix<GeneralFixed<M,M>,AT>::operator<<(A);
	return *this;
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> FixedSquareMatrix<M,AT>& operator<<(const Matrix<T, AT> &A) {
	Matrix<GeneralFixed<M,M>,AT>::operator<<(A);
	return *this;
      }

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix.
       * */
      int size() const {return M;};

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      FixedSquareMatrix<M,AT> copy() const;

      using Matrix<GeneralFixed<M,M>, AT>::operator();

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      operator std::vector<std::vector<AT> >();

  //    FixedSquareMatrix<AT> T() {
  //      return FixedSquareMatrix<AT>(n, lda, tp?false:true, memory, ele);
  //    }

  //    const FixedSquareMatrix<AT> T() const {
  //      return FixedSquareMatrix<AT>(n, lda, tp?false:true, memory, ele);
  //    }
  };

  template <int M, class AT>
    FixedSquareMatrix<M,AT> FixedSquareMatrix<M,AT>::copy() const {

      FixedSquareMatrix<M,AT> A(NONINIT);
      A.deepCopy(*this);

      return A;
    }

  template <int M, class AT>
    FixedSquareMatrix<M,AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(size());
      for(int r=0; r<size(); r++) {
	ret[r].resize(size());
	for(int c=0; c<size(); c++)
	  ret[r][c]= this->ele[r+c*M]; // ret[r][c]=operator()(r,c);
      }
      return ret;
    }

}

#endif
