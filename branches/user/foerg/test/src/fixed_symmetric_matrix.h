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

#ifndef fixed_symmetric_matrix_h
#define fixed_symmetric_matrix_h

#include "types.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for symmetric matrices.
   *
   * Template class Matrix of shape type SymmetricFixed. 
   * The template parameter AT defines the
   * atomic type of the matrix. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <int M, class AT> class Matrix<Symmetric,Fixed<M>,Fixed<M>,AT> {

    protected:

    /// @cond NO_SHOW

      AT ele[M][M];

      template <class Type, class Row, class Col> inline void deepCopy(const Matrix<Type,Row,Col,AT> &A); 
      template <class Row> inline void deepCopy(const Matrix<Symmetric,Row,Row,AT> &A); 

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      Matrix() {
#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      /*! \brief Regular Constructor
       *
       * Constructs a symmetric matrix of size n x n. The matrix will be 
       * initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param n_ The number of rows and columns.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * */
      Matrix(Initialization ini, const AT &a=0) {  

	if(ini == INIT)
	  init(a);
	else if(ini == EYE) {	 
	  init(0);
	  for(int i=0; i<M; i++) 
	    ej(i,i) = 1;
	}
      }

      Matrix(NOINIT) { }
      Matrix(int n, NOINIT) { }
      Matrix(int m, int n, NOINIT) { }
      Matrix(SCALAR, const AT &a=0) { init(a); }
      Matrix(int n, SCALAR, const AT &a=0) { init(a); }
      Matrix(int m, int n, SCALAR, const AT &a=0) { init(a); }

      Matrix(int m, int n, Initialization ini, const AT &a=0) {  

	if(ini == INIT)
	  init(a);
	else if(ini == EYE) {	 
	  init(0);
	  for(int i=0; i<M; i++) 
	    ej(i,i) = 1;
	}
      }

      explicit Matrix(const Matrix<General,Fixed<M>,Fixed<M>,AT>&  A) {
	deepCopy(A);
      }

      template<class Row>
      Matrix(const Matrix<Symmetric,Row,Row,AT> &A) {

#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.size() == M); 
#endif

	deepCopy(A);
      }

      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) {

#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == M); 
	assert(A.cols() == M);
#endif

	deepCopy(A);
      }

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \remark The bounds are checked by default. 
       * To change this behavior, define
       * FMATVEC_NO_BOUNDS_CHECK.
       * \sa operator()(int,int) const
       * */
      AT& operator()(int i, int j) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<M);
	assert(j<M);
#endif

	return e(i,j);
      };

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<M);
	assert(j<M);
#endif

	return e(i,j);
      };

      AT& ei(int i, int j) {
	return ele[i][j];
      };

      const AT& ei(int i, int j) const {
	return ele[i][j];
      };

      AT& ej(int i, int j) {
	return ele[j][i];
      };

      const AT& ej(int i, int j) const {
	return ele[j][i];
      };

       AT& e(int i, int j) {
	return j > i ? ej(i,j) : ei(i,j);
      };

      const AT& e(int i, int j) const {
	return j > i ? ej(i,j) : ei(i,j);
      };

     /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return ele;};

      /*! \brief Pointer operator
       *
       * See operator()() 
       * */
      const AT* operator()() const {return ele;};

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix
       * */
      int size() const {return M;};

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix
       * */
      int rows() const {return M;};

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix
       * */
      int cols() const {return M;};

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix
       * */
      int ldim() const {return M;};

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasColMajor.
       * */
      const CBLAS_ORDER blasOrder() const {
	return  CblasColMajor;
      };

      /*! \brief Symmetry convention.
       *
       * Returns the blas-conform symmetry convention. 
       * The elements are stored in the lower triangular
       * part of the array,
       * \return CblasLower.
       * */
      const CBLAS_UPLO blasUplo() const {
	return  CblasLower;
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>& init(const AT &a);

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >();
  };

  template <int M, class AT>
    inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>&  Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::init(const AT& val) {

      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++) 
          ej(i,j) = val; 

      return *this;
    }

  template <int M, class AT>
    inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=operator()(r,c);
      }
      return ret;
    }

  /// @cond NO_SHOW

  template <int M, class AT> template <class Type, class Row, class Col>
    inline void Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::deepCopy(const Matrix<Type,Row,Col,AT> &A) { 
      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++)
          ej(i,j) = A.e(i,j);
    }

  template <int M, class AT> template <class Row>
    inline void Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::deepCopy(const Matrix<Symmetric,Row,Row,AT> &A) { 
      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++)
          ej(i,j) = A.ej(i,j);
    }

  /// @endcond

}

#endif
