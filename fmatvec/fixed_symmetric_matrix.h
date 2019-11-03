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
      typedef AT value_type;

//      template<class Ini=All<AT> >
//      Matrix(Ini ini=All<AT>()) {
//        init(ini);
//      }
//      template<class Ini=All<AT> >
//      Matrix(int m_, int n_, Ini ini=All<AT>()) {
//        init(ini);
//      }

      explicit Matrix(Noinit ini) { }
      explicit Matrix(Init ini=INIT, const AT &a=AT()) { init(a); }
      explicit Matrix(Eye ini, const AT &a=1) { init(ini,a); }
      explicit Matrix(int m, int n, Noinit ini) { assert(m==M && n==M); }
      explicit Matrix(int m, int n, Init ini, const AT &a=AT()) { assert(m==M && n==M); init(a); }
      explicit Matrix(int m, int n, Eye ini, const AT &a=1) { assert(m==M && n==M); init(ini,a); }

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

      //! Resize a fixed matrix
      //! Do nothing for the fixed dimension and throw for any other dimension.
      void resize(int n, int m) {
        if(n!=M || m!=M)
          throw std::runtime_error("A fixed symmetric matrix cannot be resized.");
      }

      // A fixed matrix is stored in c-storage order -> transposed is always true
      bool transposed() {
        return true;
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

      Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>(const std::vector<std::vector<AT>> &A) {
        for(int r=0; r<rows(); r++) {
          if(static_cast<int>(A[r].size())!=cols())
            throw std::runtime_error("The rows of the input have different length.");
          for(int c=0; c<cols(); c++) {
            e(r,c)=A[r][c];
            if(r>c && fabs(A[r][c]-A[c][r])>fabs(A[r][c])*1e-13+1e-13)
              throw std::runtime_error("The input is not symmetric.");
          }
        }
      }

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
      inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>& init(const AT &val=AT()); 
      inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>& init(Eye, const AT &val=1);
      inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>& init(Noinit, const AT &a=AT()) { return *this; }

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >() const;
  };

  template <int M, class AT>
    inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>&  Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::init(const AT &val) {

      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++) 
          ej(i,j) = val;

      return *this;
    }

  template <int M, class AT>
    inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>&  Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::init(Eye, const AT &val) {

      for(int i=0; i<size(); i++) {
        ej(i,i) = val;
        for(int j=i+1; j<size(); j++) {
          ej(i,j) = 0; 
        }
      }
      return *this;
    }

  template <int M, class AT>
    inline Matrix<Symmetric,Fixed<M>,Fixed<M>,AT>::operator std::vector<std::vector<AT> >()  const{
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
