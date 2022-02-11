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

#ifndef var_symmetric_matrix_h
#define var_symmetric_matrix_h

#include "types.h"
#include "indices.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for symmetric matrices.
   *
   * Template class Matrix of shape type SymmetricVar. 
   * The template parameter AT defines the
   * atomic type of the matrix. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class Matrix<Symmetric,Var,Var,AT> {

    protected:

    /// @cond NO_SHOW

      int M{0};

      AT *ele;

      template <class Type, class Row, class Col> inline Matrix<Symmetric,Var,Var,AT>& copy(const Matrix<Type,Row,Col,AT> &A);
      template <class Row> inline Matrix<Symmetric,Var,Var,AT>& copy(const Matrix<Symmetric,Row,Row,AT> &A);

    /// @endcond

    public:
      static constexpr bool isVector {false};
      typedef AT value_type;
      typedef Symmetric shape_type;

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      explicit Matrix() :  ele(nullptr) {
      }

      explicit Matrix(int m, Noinit) : M(m), ele(new AT[M*M]) { }
      explicit Matrix(int m, Init ini=INIT, const AT &a=AT()) : M(m), ele(new AT[M*M]) { init(a); }
      explicit Matrix(int m, Eye ini, const AT &a=1) : M(m), ele(new AT[M*M]) { init(ini,a); }
      explicit Matrix(int m, int n, Noinit) : M(m), ele(new AT[M*M]) { }
      explicit Matrix(int m, int n, Init ini=INIT, const AT &a=AT()) : M(m), ele(new AT[M*M]) { init(a); }
      explicit Matrix(int m, int n, Eye ini, const AT &a=1) : M(m), ele(new AT[M*M]) { init(ini,a); }

      // move
      Matrix(Matrix<Symmetric,Var,Var,AT> &&src) {
        M=src.M;
        src.M=0;
        ele=src.ele;
        src.ele=nullptr;
      }
      Matrix<Symmetric,Var,Var,AT>& operator=(Matrix<Symmetric,Var,Var,AT> &&src) {
        FMATVEC_ASSERT(M == src.M, AT);
        src.M=0;
        delete[]ele;
        ele=src.ele;
        src.ele=nullptr;
        return *this;
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      Matrix(const Matrix<Symmetric,Var,Var,AT> &A) : M(A.M), ele(new AT[M*M])  {
	copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Row>
      Matrix(const Matrix<Symmetric,Row,Row,AT> &A) : M(A.size()), ele(new AT[M*M])  {
	copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : M(A.rows()), ele(new AT[M*M])  {
	FMATVEC_ASSERT(A.rows() == A.cols(), AT); 
	copy(A);
      }

      /*! \brief Destructor. 
       * */
      ~Matrix() {
	delete[] ele;
      }

      Matrix<Symmetric,Var,Var,AT>& resize(int m, Noinit) { 
        delete[] ele;
        M = m;
        ele = new AT[M*M];
        return *this;
      }

      Matrix<Symmetric,Var,Var,AT>& resize(int m, Init ini=INIT, const AT &a=AT()) { return resize(m,Noinit()).init(a); }

      Matrix<Symmetric,Var,Var,AT>& resize(int m, Eye ini, const AT &a=1) { return resize(m,Noinit()).init(ini,a); }

      /*! \brief Assignment operator
       *
       * Copies the symmetric matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Symmetric,Var,Var,AT>& operator=(const Matrix<Symmetric,Var,Var,AT> &A) {
        FMATVEC_ASSERT(M == A.M, AT);
        return copy(A);
      }

      /*! \brief Assignment operator
       *
       * Copies the symmetric matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
      inline Matrix<Symmetric,Var,Var,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
        FMATVEC_ASSERT(M == A.rows(), AT);
        return copy(A);
      }

      /*! \brief Matrix assignment
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
      inline Matrix<Symmetric,Var,Var,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
        FMATVEC_ASSERT(A.rows() == A.cols(), AT);
        if(M!=A.rows()) resize(A.rows(),NONINIT);
        return copy(A);
      }
      // move
      template<class Type, class Row, class Col>
      inline Matrix<Symmetric,Var,Var,AT>& operator<<=(Matrix<Type,Row,Col,AT> &&src) {
        FMATVEC_ASSERT(src.rows() == src.cols(), AT);
        M=src.M;
        src.M=0;
        delete[]ele;
        ele=src.ele;
        src.ele=nullptr;
        return *this;
      }

      //! Resize a var symmetric matrix.
      //! Throw if the dimension does not match or resize to this dimension
      void resize(int n, int m) {
        if(n!=m)
          throw std::runtime_error("A symmetric matrix cannot have different dimensions for rows and columns.");
        resize(n);
      }

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \remark The bounds are checked in debug mode.
       * \sa operator()(int,int) const
       * */
      AT& operator()(int i, int j) {
	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(j>=0, AT);
	FMATVEC_ASSERT(i<M, AT);
	FMATVEC_ASSERT(j<M, AT);

	return e(i,j);
      }

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(j>=0, AT);
	FMATVEC_ASSERT(i<M, AT);
	FMATVEC_ASSERT(j<M, AT);

	return e(i,j);
      }

      AT& ei(int i, int j) {
	return ele[i*M+j];
      }

      const AT& ei(int i, int j) const {
	return ele[i*M+j];
      }

      AT& ej(int i, int j) {
	return ele[i+j*M];
      }

      const AT& ej(int i, int j) const {
	return ele[i+j*M];
      }

       AT& e(int i, int j) {
	return j > i ? ej(i,j) : ei(i,j);
      }

      const AT& e(int i, int j) const {
	return j > i ? ej(i,j) : ei(i,j);
      }

     /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return ele;}

      /*! \brief Pointer operator
       *
       * See operator()() 
       * */
      const AT* operator()() const {return ele;}

      /*! \brief Size.
       *
       * \return The number of rows and columns of the matrix
       * */
      constexpr int size() const {return M;}

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix
       * */
      constexpr int rows() const {return M;}

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix
       * */
      constexpr int cols() const {return M;}

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix
       * */
      int ldim() const {return M;}

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasRowMajor.
       * */
      const CBLAS_ORDER blasOrder() const {
	return CblasRowMajor;
      }

      /*! \brief Symmetry convention.
       *
       * Returns the blas-conform symmetry convention. 
       * The elements are stored in the lower triangular
       * part of the array,
       * \return CblasLower.
       * */
      const CBLAS_UPLO blasUplo() const {
	return  CblasLower;
      }

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Symmetric,Var,Var,AT>& init(const AT &val=AT()); 
      inline Matrix<Symmetric,Var,Var,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<Symmetric,Var,Var,AT>& init(Eye eye, const AT &val=1);
      inline Matrix<Symmetric,Var,Var,AT>& init(Noinit, const AT &a=AT()) { return *this; }

      /*! \brief Submatrix operator.
       * */
      inline const Matrix<General,Var,Var,AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const;

      /*! \brief Submatrix operator.
       * */
      inline const Matrix<Symmetric,Var,Var,AT> operator()(const Range<Var,Var> &I) const;

      template<class Type, class Row, class Col> inline void set(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);
      template<class Type, class Row, class Col> inline void add(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      template<class Row> inline void set(const fmatvec::Range<Var,Var> &I, const Matrix<Symmetric,Row,Row,AT> &A);
      template<class Row> inline void add(const fmatvec::Range<Var,Var> &I, const Matrix<Symmetric,Row,Row,AT> &A);

      inline const Matrix<General,Var,Var,AT> operator()(const Indices &I, const Indices &J) const;

      inline const Matrix<Symmetric,Var,Var,AT> operator()(const Indices &I) const;

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit inline operator std::vector<std::vector<AT>>() const;

      inline int countElements() const;
  };

  template <class AT>
    inline Matrix<Symmetric,Var,Var,AT>&  Matrix<Symmetric,Var,Var,AT>::init(const AT &val) {

      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++) 
          ej(i,j) = val;

      return *this;
    }

  template <class AT>
    inline Matrix<Symmetric,Var,Var,AT>&  Matrix<Symmetric,Var,Var,AT>::init(Eye, const AT &val) {

      for(int i=0; i<size(); i++) {
        ej(i,i) = val;
        for(int j=i+1; j<size(); j++) {
          ej(i,j) = 0; 
        }
      }
      return *this;
    }

  template <class AT>
    inline const Matrix<General,Var,Var,AT> Matrix<Symmetric,Var,Var,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
      FMATVEC_ASSERT(I.end()<M, AT);
      FMATVEC_ASSERT(J.end()<M, AT);
      Matrix<General,Var,Var,AT> A(I.end()-I.start()+1,J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

  template <class AT>
    inline const Matrix<Symmetric,Var,Var,AT> Matrix<Symmetric,Var,Var,AT>::operator()(const Range<Var,Var> &I) const {
      FMATVEC_ASSERT(I.end()<M, AT);

      Matrix<Symmetric,Var,Var,AT> A(I.end()-I.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=i; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,I.start()+j);

      return A;
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<Symmetric,Var,Var,AT>::set(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      if(I.start()>=J.start()) FMATVEC_ASSERT(I.start()>=J.end(), AT)
      else FMATVEC_ASSERT(J.start()>=I.end(), AT);
      FMATVEC_ASSERT(I.end()<rows(), AT);
      FMATVEC_ASSERT(J.end()<cols(), AT);
      FMATVEC_ASSERT(I.size()==A.rows(), AT);
      FMATVEC_ASSERT(J.size()==A.cols(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<Symmetric,Var,Var,AT>::add(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      if(I.start()>=J.start()) FMATVEC_ASSERT(I.start()>=J.end(), AT)
      else FMATVEC_ASSERT(J.start()>=I.end(), AT);
      FMATVEC_ASSERT(I.end()<rows(), AT);
      FMATVEC_ASSERT(J.end()<cols(), AT);
      FMATVEC_ASSERT(I.size()==A.rows(), AT);
      FMATVEC_ASSERT(J.size()==A.cols(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <class AT> template<class Row>
    inline void Matrix<Symmetric,Var,Var,AT>::set(const fmatvec::Range<Var,Var> &I, const Matrix<Symmetric,Row,Row,AT> &A) {

      FMATVEC_ASSERT(I.end()<size(), AT);
      FMATVEC_ASSERT(I.size()==A.size(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=i, jj=ii; j<=I.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <class AT> template<class Row>
    inline void Matrix<Symmetric,Var,Var,AT>::add(const fmatvec::Range<Var,Var> &I, const Matrix<Symmetric,Row,Row,AT> &A) {

      FMATVEC_ASSERT(I.end()<size(), AT);
      FMATVEC_ASSERT(I.size()==A.size(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=i, jj=ii; j<=I.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <class AT>
    inline const Matrix<General,Var,Var,AT> Matrix<Symmetric,Var,Var,AT>::operator()(const Indices &I, const Indices &J) const {
      FMATVEC_ASSERT(I.max()<size(), AT);
      FMATVEC_ASSERT(J.max()<size(), AT);

      Matrix<General,Var,Var,AT> A(I.size(),J.size(),NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
	A.e(i,j) = e(I[i],J[j]);

      return A;
    }

  template <class AT>
    inline const Matrix<Symmetric,Var,Var,AT> Matrix<Symmetric,Var,Var,AT>::operator()(const Indices &I) const {
      FMATVEC_ASSERT(I.max()<size(), AT);

      Matrix<Symmetric,Var,Var,AT> A(I.size(),NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
	A.e(i,j) = e(I[i],I[j]);

      return A;
    }

  template <class AT>
    inline Matrix<Symmetric,Var,Var,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=operator()(r,c);
      }
      return ret;
    }

  /// @cond NO_SHOW

  template <class AT> template <class Type, class Row, class Col>
    inline Matrix<Symmetric,Var,Var,AT>& Matrix<Symmetric,Var,Var,AT>::copy(const Matrix<Type,Row,Col,AT> &A) {
      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++) 
          ej(i,j) = A.e(i,j);
      return *this;
    }

  template <class AT> template <class Row>
    inline Matrix<Symmetric,Var,Var,AT>& Matrix<Symmetric,Var,Var,AT>::copy(const Matrix<Symmetric,Row,Row,AT> &A) {
      for(int i=0; i<M; i++) 
        for(int j=i; j<M; j++) 
          ej(i,j) = A.ej(i,j);
      return *this;
    }

  template <class AT>
    inline int Matrix<Symmetric,Var,Var,AT>::countElements() const {
      int k = M;
      for (int j = 0; j < M; j++) {
        for (int i = j+1; i < M; i++) {
          if (fabs(e(i, j)) > 1e-16) {
            k+=2;
          }
        }
      }
      return k;
    }

  /// @endcond

}

#endif
