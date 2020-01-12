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

#ifndef symmetric_matrix_h
#define symmetric_matrix_h

#include "types.h"
#include "matrix.h"
#include "_memory.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for symmetric matrices.
   *
   * Template class Matrix of shape type Symmetric. 
   * The template parameter AT defines the
   * atomic type of the matrix. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class Matrix<Symmetric,Ref,Ref,AT> {

    protected:

    /// @cond NO_SHOW

      Memory<AT> memory;
      AT *ele;
      int n{0};
      int lda{0};

      friend class Matrix<General,Ref,Ref,AT>;

      template <class Type, class Row, class Col> inline Matrix<Symmetric,Ref,Ref,AT>& copy(const Matrix<Type,Row,Col,AT> &A);
      inline Matrix<Symmetric,Ref,Ref,AT>& copy(const Matrix<Symmetric,Ref,Ref,AT> &A);

      const AT* elePtr(int i, int j) const {
	return  j > i ? ele+i*lda+j : ele+i+j*lda; 
      }

      AT* elePtr(int i, int j) {
	return  j > i ? ele+i*lda+j : ele+i+j*lda; 
      }

    /// @endcond

    public:

      typedef AT value_type;

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      explicit Matrix() : memory(), ele(nullptr) {
      }

      explicit Matrix(int n_, Noinit) : memory(n_*n_), ele((AT*)memory.get()), n(n_), lda(n_) { }
      explicit Matrix(int n_, Init ini=INIT, const AT &a=AT()) : memory(n_*n_), ele((AT*)memory.get()), n(n_), lda(n_) { init(a); }
      explicit Matrix(int n_, Eye ini, const AT &a=1) : memory(n_*n_), ele((AT*)memory.get()), n(n_), lda(n_) { init(ini,a); }
      explicit Matrix(int m_, int n_, Noinit) : memory(n_*n_), ele((AT*)memory.get()), n(n_), lda(n_) { }
      explicit Matrix(int m_, int n_, Init ini=INIT, const AT &a=AT()) : memory(n_*n_), ele((AT*)memory.get()), n(n_), lda(n_) { init(a); }
      explicit Matrix(int m_, int n_, Eye ini, const AT &a=1) : memory(n_*n_), ele((AT*)memory.get()), n(n_), lda(n_) { init(ini,a); }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      Matrix(const Matrix<Symmetric,Ref,Ref,AT> &A) : memory(A.n*A.n), ele((AT*)memory.get()) , n(A.n), lda(n) {
        copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Row>
      Matrix(const Matrix<Symmetric,Row,Row,AT> &A) : memory(A.size()*A.size()), ele((AT*)memory.get()), n(A.size()), lda(n) {
        copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : memory(A.rows()*A.cols()), ele((AT*)memory.get()), n(A.cols()), lda(n) {
	assert(A.rows() == A.cols());
	copy(A);
      }

      /*! \brief Regular Constructor
       *
       * Constructs a symmetric matrix of size n x n with 
       * the pyhsical memory given by \em ele_.
       * \param n_ The number of rows and columns.
       * \param ele_ The physical memory the matrix will point to.
       * */
      explicit Matrix(int n_, AT* ele_) : memory(), ele(ele_), n(n_), lda(n_) { }

      /*! \brief Destructor. 
       * */
      ~Matrix() = default;

      Matrix<Symmetric,Ref,Ref,AT>& resize(int n_, Noinit) { 
	n = n_; lda = n;
	memory.resize(n*n);
	ele = (AT*)memory.get();
        return *this;
      }

      Matrix<Symmetric,Ref,Ref,AT>& resize(int n, Init ini=INIT, const AT &a=AT()) { return resize(n,Noinit()).init(a); }

      Matrix<Symmetric,Ref,Ref,AT>& resize(int n, Eye ini, const AT &a=1) { return resize(n,Noinit()).init(ini,a); }

      //! Resize a symmetric matrix
      //! Throw if the row and col dimension does not match or resize to this dimension.
      void resize(int n, int m) {
        if(n!=m)
          throw std::runtime_error("A symmetric matrix cannot have different dimensions for rows and columns.");
        resize(n);
      }

      /*! \brief Assignment operator
       *
       * Copies the symmetric matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Symmetric,Ref,Ref,AT>& operator=(const Matrix<Symmetric,Ref,Ref,AT> &A) {
        assert(n == A.rows());
        return copy(A);
      }

      /*! \brief Assignment operator
       *
       * Copies the symmetric matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
      inline Matrix<Symmetric,Ref,Ref,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
        assert(A.rows() == A.cols());
        assert(n == A.rows());
        return copy(A);
      }

      /*! \brief Reference operator
       *
       * References the symmetric matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Symmetric,Ref,Ref,AT>& operator&=(const Matrix<Symmetric,Ref,Ref,AT> &A) {
        n=A.n;
        memory = A.memory;
        ele = A.ele;
        lda = A.lda;
        return *this;
      }

      /*! \brief Matrix assignment
       *
       * Copies the symmetric matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
        inline Matrix<Symmetric,Ref,Ref,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
          assert(A.rows() == A.cols());
          if(n!=A.rows()) resize(A.rows(),NONINIT);
          return copy(A);
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
	assert(i>=0);
	assert(j>=0);
	assert(i<n);
	assert(j<n);
	return e(i,j);
      }

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
	assert(i>=0);
	assert(j>=0);
	assert(i<n);
	assert(j<n);

	return e(i,j);
      }

      AT& ei(int i, int j) {
	return ele[i+j*lda];
      }

      const AT& ei(int i, int j) const {
	return ele[i+j*lda];
      }

      AT& ej(int i, int j) {
	return ele[i*lda+j];
      }

      const AT& ej(int i, int j) const {
	return ele[i*lda+j];
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
      int size() const {return n;}

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix
       * */
      int rows() const {return n;}

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix
       * */
      int cols() const {return n;}

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix
       * */
      int ldim() const {return lda;}

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasColMajor.
       * */
      const CBLAS_ORDER blasOrder() const {
	return CblasColMajor;
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
      inline Matrix<Symmetric,Ref,Ref,AT>& init(const AT &val); 
      inline Matrix<Symmetric,Ref,Ref,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<Symmetric,Ref,Ref,AT>& init(Eye, const AT &val=1);
      inline Matrix<Symmetric,Ref,Ref,AT>& init(Noinit, const AT &a=AT()) { return *this; }

      /*! \brief Submatrix operator.
       * */
      inline const Matrix<General,Ref,Ref,AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const;

      /*! \brief Submatrix operator.
       * */
      inline const Matrix<Symmetric,Ref,Ref,AT> operator()(const Range<Var,Var> &I) const;

      template<class Type, class Row, class Col> inline void set(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);
      template<class Type, class Row, class Col> inline void add(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);
      inline void ref(Matrix<Symmetric,Ref,Ref,AT> &A, const fmatvec::Range<Var,Var> &I);

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit inline operator std::vector<std::vector<AT>>() const;

      /*! \brief std::vector<std::vector<AT>> Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT>> object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT>> the matrix will be initialized with.
       * */
      explicit Matrix(const std::vector<std::vector<AT>> &m);

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the matrix
//       * */
//      explicit operator AT() const {
//        assert(n==1);
//        return e(0);
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a matrix with a AT object.
//       * \param x The AT the matrix will be initialized with.
//       * */
//      explicit Matrix(const AT &x) : memory(1), ele((AT*)memory.get()), n(1), lda(1) { ele[0] = x; }

      inline int countElements() const;
  };

  template <class AT>
    inline Matrix<Symmetric,Ref,Ref,AT>&  Matrix<Symmetric,Ref,Ref,AT>::init(const AT &val) {

      for(int i=0; i<rows(); i++) 
        for(int j=i; j<cols(); j++) 
          ej(i,j) = val;

      return *this;
    }

  template <class AT>
    inline Matrix<Symmetric,Ref,Ref,AT>&  Matrix<Symmetric,Ref,Ref,AT>::init(Eye, const AT &val) {

      for(int i=0; i<size(); i++) {
        ej(i,i) = val;
        for(int j=i+1; j<size(); j++) {
          ej(i,j) = 0; 
        }
      }
      return *this;
    }

  template <class AT>
    inline const Matrix<General,Ref,Ref,AT> Matrix<Symmetric,Ref,Ref,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
      assert(I.end()<n);
      assert(J.end()<n);
      Matrix<General,Ref,Ref,AT> A(I.end()-I.start()+1,J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

  template <class AT> 
    inline const Matrix<Symmetric,Ref,Ref,AT> Matrix<Symmetric,Ref,Ref,AT>::operator()(const Range<Var,Var> &I) const {
      assert(I.end()<n);

      Matrix<Symmetric,Ref,Ref,AT> A(I.end()-I.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=i; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,I.start()+j);

      return A;
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<Symmetric,Ref,Ref,AT>::set(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      assert(I.end()<rows());
      assert(J.end()<cols());
      assert(I.size()==A.rows());
      assert(J.size()==A.cols());

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<Symmetric,Ref,Ref,AT>::add(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      assert(I.end()<rows());
      assert(J.end()<cols());
      assert(I.size()==A.rows());
      assert(J.size()==A.cols());

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <class AT>
    inline void Matrix<Symmetric,Ref,Ref,AT>::ref(Matrix<Symmetric,Ref,Ref,AT> &A, const fmatvec::Range<Var,Var> &I) {
      assert(I.end()<A.size());
      n=I.size();
      memory = A.memory;
      ele = A.elePtr(I.start(),I.start());
      lda = A.lda;
    }

  template <class AT>
    inline Matrix<Symmetric,Ref,Ref,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=operator()(r,c);
      }
      return ret;
    }

  template <class AT>
    Matrix<Symmetric,Ref,Ref,AT>::Matrix(const std::vector<std::vector<AT>> &m) : Matrix<Symmetric,Ref,Ref,AT>(m.size()) {
      for(int r=0; r<rows(); r++) {
        if(static_cast<int>(m[r].size())!=cols())
          throw std::runtime_error("The rows of the input have different length.");
        for(int c=0; c<cols(); c++) {
          e(r,c)=m[r][c];
          if(r>c && fabs(m[r][c]-m[c][r])>fabs(m[r][c])*1e-13+1e-13)
            throw std::runtime_error("The input is not symmetric.");
        }
      }
    }

  template <class AT>
    inline int Matrix<Symmetric,Ref,Ref,AT>::countElements() const {
      int k = n;
      for (int j = 0; j < n; j++) {
        for (int i = j+1; i < n; i++) {
          if (fabs(e(i, j)) > 1e-16) {
            k+=2;
          }
        }
      }
      return k;
    }

  /// @cond NO_SHOW

  template <class AT>
    inline Matrix<Symmetric,Ref,Ref,AT>& Matrix<Symmetric,Ref,Ref,AT>::copy(const Matrix<Symmetric,Ref,Ref,AT> &A) {
      for(int i=0; i<n; i++) 
        for(int j=i; j<n; j++) 
          ej(i,j) = A.ej(i,j);
      return *this;
    }

  template <class AT> template <class Type, class Row, class Col>
    inline Matrix<Symmetric,Ref,Ref,AT>& Matrix<Symmetric,Ref,Ref,AT>::copy(const Matrix<Type,Row,Col,AT> &A) {
      for(int i=0; i<n; i++) 
        for(int j=i; j<n; j++) 
          ej(i,j) = A.e(i,j);
      return *this;
    }

  /// @endcond

}

#endif
