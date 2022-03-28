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

#ifndef general_matrix_h
#define general_matrix_h

#include "_memory.h"
#include "types.h"
#include "indices.h"
#include "matrix.h"
#include <cstdlib>
#include <stdexcept>

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type General and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <class AT> class Matrix<General,Ref,Ref,AT> {

    public:
      static constexpr bool isVector {false};

      typedef AT value_type;
      typedef General shape_type;

 /// @cond NO_SHOW

      friend class Vector<Ref,AT>;
      friend class RowVector<Ref,AT>;
      friend class SquareMatrix<Ref,AT>;
      friend class Matrix<Symmetric,Ref,Ref,AT>;
      
    protected:

      template<bool Const>
      class Iterator {
        typedef Matrix<General,Ref,Ref,AT> Mat;
        public:
          Iterator(std::conditional_t<Const, const Mat, Mat> &mat_, int row_, int col_) : mat(mat_), row(row_), col(col_) {}
          bool operator!=(const Iterator& x) const { return row!=x.row || col!=x.col; }
          Iterator& operator++() { ++row; if(row>=mat.m) { row=0; ++col; } return *this; }
          std::conditional_t<Const, const AT, AT>& operator*() { return mat.e(row,col); }
          std::conditional_t<Const, const AT, AT>& operator->() { return mat.e(row,col); }
          const AT& operator*() const { return mat.e(row,col); }
          const AT& operator->() const { return mat.e(row,col); }
        private:
          std::conditional_t<Const, const Mat, Mat> &mat;
          int row;
          int col;
      };

      Memory<AT> memory;
      AT *ele;
      int m{0};
      int n{0};
      int lda{0};

      template <class Type, class Row, class Col> inline Matrix<General,Ref,Ref,AT>& copy(const Matrix<Type,Row,Col,AT> &A);

      const AT* elePtr(int i, int j) const {
	return ele+i+j*lda;
      }

      AT* elePtr(int i, int j) {
	return ele+i+j*lda;
      }

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      explicit Matrix() : memory(), ele(nullptr) { }

// Works with -std=gnu++0x only
//      template<class Ini=All<AT>>
//      Matrix(int m_, int n_, Ini ini=All<AT>()) : memory(m_*n_), ele((AT*)memory.get()), m(m_), n(n_), lda(m_) {
//        init(ini);
//      }

      explicit Matrix(int m_, int n_, Noinit) : memory(m_*n_), ele((AT*)memory.get()), m(m_), n(n_), lda(m_) { }
      explicit Matrix(int m_, int n_, Init ini=INIT, const AT &a=AT()) : memory(m_*n_), ele((AT*)memory.get()), m(m_), n(n_), lda(m_) { init0(a); }
      explicit Matrix(int m_, int n_, Eye ini, const AT &a=1) : memory(m_*n_), ele((AT*)memory.get()), m(m_), n(n_), lda(m_) { init(ini,a); }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      Matrix(const Matrix<General,Ref,Ref,AT> &A) : memory(A.m*A.n), ele((AT*)memory.get()) ,m(A.m), n(A.n), lda(m) {
        copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Row, class Col>
      Matrix(const Matrix<General,Row,Col,AT> &A) : memory(A.rows()*A.cols()), ele((AT*)memory.get()), m(A.rows()), n(A.cols()), lda(m) {
        copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : memory(A.rows()*A.cols()), ele((AT*)memory.get()), m(A.rows()), n(A.cols()), lda(m) {
	copy(A);
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x n with the pyhsical memory given by \em ele_.
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \param ele_ The physical memory the matrix will point to.
       * */
      explicit Matrix(int m_, int n_, AT* ele_) : memory(), ele(ele_), m(m_), n(n_), lda(m_) { 
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<General<Ref,Ref>double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const std::string &strs);
      Matrix(const char *strs);

      /*! \brief Destructor. 
       * */
      ~Matrix() = default;

      typedef Iterator<false> iterator;
      typedef Iterator<true> const_iterator;
      iterator begin() { return iterator(*this,0,0); }
      iterator end() { return iterator(*this,0,n); }
      const_iterator begin() const { return const_iterator(*this,0,0); }
      const_iterator end() const { return const_iterator(*this,0,n); }

      Matrix<General,Ref,Ref,AT>& resize(int m_, int n_, Noinit) {
	m = m_; n = n_; lda = m;
	memory.resize(m*n);
	ele = (AT*)memory.get();
        return *this;
      }

      Matrix<General,Ref,Ref,AT>& resize(int m, int n, Init ini=INIT, const AT &a=AT()) { return resize(m,n,Noinit()).init0(a); }

      Matrix<General,Ref,Ref,AT>& resize(int m, int n, Eye ini, const AT &a=1) { return resize(m,n,Noinit()).init(ini,a); } 

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * */
      Matrix<General,Ref,Ref,AT>& operator=(const Matrix<General,Ref,Ref,AT> &A) {
        FMATVEC_ASSERT(m == A.rows(), AT);
        FMATVEC_ASSERT(n == A.cols(), AT);
        return copy(A);
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
      Matrix<General,Ref,Ref,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
        FMATVEC_ASSERT(m == A.rows(), AT);
        FMATVEC_ASSERT(n == A.cols(), AT);
        return copy(A);
      }

      /*! \brief Reference operator
       *
       * References the matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Ref,Ref,AT>& operator&=(Matrix<General,Ref,Ref,AT> &A) {
        m=A.m; 
        n=A.n;
        memory = A.memory;
        ele = A.ele;
        lda = A.lda;
        return *this;
      }

      /*! \brief Matrix assignment
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
        inline Matrix<General,Ref,Ref,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
          if(m!=A.rows() || n!=A.cols()) resize(A.rows(),A.cols(),NONINIT);
          return copy(A);
        }

      template <class AT2>
      operator Matrix<General,Ref,Ref,AT2>() const {
        Matrix<General,Ref,Ref,AT2> ret(rows(), cols());
        for(size_t r=0; r<rows(); ++r)
          for(size_t c=0; c<cols(); ++c)
            ret(r,c) = (*this)(r,c);
        return ret;
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
	FMATVEC_ASSERT(i<m, AT);
	FMATVEC_ASSERT(j<n, AT);

	return e(i,j);
      }

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(j>=0, AT);
	FMATVEC_ASSERT(i<m, AT);
	FMATVEC_ASSERT(j<n, AT);

	return e(i,j);
      }

      AT& e(int i, int j) {
        return ele[i+j*lda];
      }

      const AT& e(int i, int j) const {
	return ele[i+j*lda];
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

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix.
       * */
      int rows() const {return m;}

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix.
       * */
      int cols() const {return n;}

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix.
       * */
      int ldim() const {return lda;}

      /*! \brief Transposed status.
       *
       * Returns the blas-conform transposed status.
       * \return CblasTrans if the matrix is in transposed form, CblasNoTrans
       * otherwise. 
       * */
      const CBLAS_TRANSPOSE blasTrans() const {
        return CblasNoTrans;
      }

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

      /*! \brief Submatrix operator.
       *
       * Returns a submatrix of the calling matrix. 
       * For example
       * \code 
       * B = A(Range<Var,Var>(1,2),Range<Var,Var>(2,4));
       * \endcode
       * yields
       * \f[ 
       * A=\begin{pmatrix}
       * 	a_{00} & a_{01} & a_{02} & a_{03} & a_{04}\\
       * 	a_{10} & a_{11} & a_{12} & a_{13} & a_{14}\\
       * 	a_{20} & a_{21} & a_{22} & a_{23} & a_{24}\\
       * 	a_{30} & a_{31} & a_{32} & a_{33} & a_{34}
       * \end{pmatrix}\quad \Rightarrow \quad
       * B=\begin{pmatrix}
       * 	 a_{12} & a_{13} & a_{14}\\
       * 	 a_{22} & a_{23} & a_{24}
       * \end{pmatrix}
       * \f]
       * \attention The submatrix and the
       * calling matrix will share the same physical memory.
       * \param I Range containing the starting and the ending row.
       * \param J Range containing the starting and the ending column.
       * \return A submatrix of the calling matrix.
       * */

      /*! \brief Submatrix operator.
       *
       * See operator()(const Range<Var,Var>&, const Range<Var,Var>&)
       * */
      inline const Matrix<General,Ref,Ref,AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const;

      /*! \brief Submatrix operator.
       * */
      inline const SquareMatrix<Ref,AT> operator()(const Range<Var,Var> &I) const;

      /*! \brief Row operator.
       * */
      inline const RowVector<Ref,AT> row(int i) const;

      /*! \brief Column operator.
       * */
      inline const Vector<Ref,AT> col(int j) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Ref,Ref,AT>& init(const AT &val=AT());
      inline Matrix<General,Ref,Ref,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<General,Ref,Ref,AT>& init(Eye, const AT &val=1);
      inline Matrix<General,Ref,Ref,AT>& init(Noinit, const AT &a=AT()) { return *this; }
      inline Matrix<General,Ref,Ref,AT>& init0(const AT &val=AT());
      inline Matrix<General,Ref,Ref,AT>& init0(Init, const AT &a=AT()) { return init0(a); }

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit inline operator std::vector<std::vector<AT>>() const;

      /*! \brief std::vector<std::vector<AT>> Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT>> object.
       * An FMATVEC_ASSERT checks for constant length of each row.
       * \param m The std::vector<std::vector<AT>> the matrix will be initialized with.
       * */
      explicit inline Matrix(const std::vector<std::vector<AT>> &m);

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the matrix
//       * */
//      explicit operator AT() const {
//        FMATVEC_ASSERT(m==1, AT);
//        FMATVEC_ASSERT(n==1, AT);
//        return e(0);
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a matrix with a AT object.
//       * \param x The AT the matrix will be initialized with.
//       * */
//      explicit Matrix(const AT &x) : memory(1), ele((AT*)memory.get()), m(1), n(1), lda(1) { ele[0] = x; }

      const Matrix<General,Ref,Ref,AT> T() const;

      template<class Row> inline void set(int j, const Vector<Row,AT> &x);

      template<class Col> inline void set(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void set(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      template<class Row> inline void add(int j, const Vector<Row,AT> &x);

      template<class Col> inline void add(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void add(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      template<class Type, class Row, class Col> inline void sub(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      inline const Matrix<General,Ref,Ref,AT> operator()(const Indices &I, const Indices &J) const;

      template<class Type, class Row, class Col> inline void set(const Indices &I, const Indices &J, const Matrix<Type,Row,Col,AT> &A);

      template<class Row> inline void set(const Indices &I, int j, const Vector<Row,AT> &x);

      inline void ref(Matrix<General,Ref,Ref,AT> &A, const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J);

      inline void ref(Matrix<Symmetric,Ref,Ref,AT> &A, const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J);

      /*! \brief Count nonzero elements.
       *
       * This function counts the nonzero elements of a matrix. ALL diagonal
       * elements are treated as NONZERO!!! (See the storage format for sparse matrix)
       * \return The number of nonzero or diagonal elemets.
       * */
      int countElements() const;
  };

  template <class AT> 
    Matrix<General,Ref,Ref,AT>::Matrix(const std::string &strs) : memory(), ele(nullptr) {
      std::istringstream iss(strs);
      iss>>*this;

      // check end of stream
      iss>>std::ws;
      if(!iss.eof())
        throw std::runtime_error("Input not fully read.");
    }
  template <class AT> Matrix<General,Ref,Ref,AT>::Matrix(const char *strs) : Matrix<General,Ref,Ref,AT>::Matrix(std::string(strs)) {}

  template <class AT>
    inline Matrix<General,Ref,Ref,AT>&  Matrix<General,Ref,Ref,AT>::init0(const AT &val) {
      for(int i=0; i<m*n; i++) ele[i]=val;
      return *this;
    }

  template <class AT>
    inline Matrix<General,Ref,Ref,AT>&  Matrix<General,Ref,Ref,AT>::init(const AT &val) {
      for(int i=0; i<rows(); i++)
        for(int j=0; j<cols(); j++)
          e(i,j) = val;
      return *this;
    }

  template <class AT>
    inline Matrix<General,Ref,Ref,AT>&  Matrix<General,Ref,Ref,AT>::init(Eye, const AT &val) {
      for(int i=0; i<m; i++)
        for(int j=0; j<n; j++)
          e(i,j) = (i==j) ? val : 0;
      return *this;
    }

  template <class AT> 
    inline const Matrix<General,Ref,Ref,AT> Matrix<General,Ref,Ref,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
      FMATVEC_ASSERT(I.end()<m, AT);
      FMATVEC_ASSERT(J.end()<n, AT);
      Matrix<General,Ref,Ref,AT> A(I.end()-I.start()+1,J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

  template <class AT> 
    inline const SquareMatrix<Ref,AT> Matrix<General,Ref,Ref,AT>::operator()(const Range<Var,Var> &I) const {
      FMATVEC_ASSERT(I.end()<m, AT);
      SquareMatrix<Ref,AT> A(I.end()-I.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,I.start()+j);

      return A;
    }

  template <class AT>
    inline const RowVector<Ref,AT>  Matrix<General,Ref,Ref,AT>::row(int i) const {

      FMATVEC_ASSERT(i>=0, AT);
      FMATVEC_ASSERT(i<m, AT);

      RowVector<Ref,AT> x(n,NONINIT);

      for(int j=0; j<n; j++)
        x.e(j) = e(i,j);

      return x;
    }

  template <class AT>
    inline const Vector<Ref,AT>  Matrix<General,Ref,Ref,AT>::col(int j) const {

      FMATVEC_ASSERT(j>=0, AT);
      FMATVEC_ASSERT(j<n, AT);

      Vector<Ref,AT> x(m,NONINIT);

      for(int i=0; i<m; i++)
        x.e(i) = e(i,j);

      return x;
    }

  template <class AT>
    inline Matrix<General,Ref,Ref,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]= e(r,c);
      }
      return ret;
    }

  template <class AT>
    inline Matrix<General,Ref,Ref,AT>::Matrix(const std::vector<std::vector<AT>> &m) : memory(m.size()*m[0].size()), ele((AT*)memory.get()), m(static_cast<int>(m.size())), n(static_cast<int>(m[0].size())), lda(static_cast<int>(m.size())) {
      for(int r=0; r<rows(); r++) {
        if(static_cast<int>(m[r].size())!=cols())
          throw std::runtime_error("The rows of the input have different length.");
        for(int c=0; c<cols(); c++)
          e(r,c)=m[r][c];
      }
    }

  /// @cond NO_SHOW

  template <class AT> template <class Type, class Row, class Col>
    inline Matrix<General,Ref,Ref,AT>& Matrix<General,Ref,Ref,AT>::copy(const Matrix<Type,Row,Col,AT> &A) {
      for(int i=0; i<m; i++)
        for(int j=0; j<n; j++)
          e(i,j) = A.e(i,j);
      return *this;
    }

  template <class AT>
    inline const Matrix<General,Ref,Ref,AT> Matrix<General,Ref,Ref,AT>::T() const {
      Matrix<General,Ref,Ref,AT> A(n,m,NONINIT);
      for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <class AT> template <class Row>
    inline void Matrix<General,Ref,Ref,AT>::set(int j, const Vector<Row,AT> &x) {
      FMATVEC_ASSERT(j<cols(), AT);
      FMATVEC_ASSERT(rows()==x.size(), AT);
      for(int i=0; i<rows(); i++)
        e(i,j) = x.e(i);
    }

  template <class AT> template <class Col>
    inline void Matrix<General,Ref,Ref,AT>::set(int i, const RowVector<Col,AT> &x) {
      FMATVEC_ASSERT(i<rows(), AT);
      FMATVEC_ASSERT(cols()==x.size(), AT);
      for(int j=0; j<cols(); j++)
        e(i,j) = x.e(j);
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Ref,Ref,AT>::set(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      FMATVEC_ASSERT(I.end()<rows(), AT);
      FMATVEC_ASSERT(J.end()<cols(), AT);
      FMATVEC_ASSERT(I.size()==A.rows(), AT);
      FMATVEC_ASSERT(J.size()==A.cols(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <class AT> template <class Row>
    inline void Matrix<General,Ref,Ref,AT>::add(int j, const Vector<Row,AT> &x) {
      FMATVEC_ASSERT(j<cols(), AT);
      FMATVEC_ASSERT(rows()==x.size(), AT);
      for(int i=0; i<rows(); i++)
        e(i,j) += x.e(i);
    }

  template <class AT> template <class Col>
    inline void Matrix<General,Ref,Ref,AT>::add(int i, const RowVector<Col,AT> &x) {
      FMATVEC_ASSERT(i<rows(), AT);
      FMATVEC_ASSERT(cols()==x.size(), AT);
      for(int j=0; j<cols(); j++)
        e(i,j) += x.e(j);
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Ref,Ref,AT>::add(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      FMATVEC_ASSERT(I.end()<rows(), AT);
      FMATVEC_ASSERT(J.end()<cols(), AT);
      FMATVEC_ASSERT(I.size()==A.rows(), AT);
      FMATVEC_ASSERT(J.size()==A.cols(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Ref,Ref,AT>::sub(const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      FMATVEC_ASSERT(I.end()<rows(), AT);
      FMATVEC_ASSERT(J.end()<cols(), AT);
      FMATVEC_ASSERT(I.size()==A.rows(), AT);
      FMATVEC_ASSERT(J.size()==A.cols(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) -= A.e(ii,jj);
    }

  template <class AT>
    inline const Matrix<General,Ref,Ref,AT> Matrix<General,Ref,Ref,AT>::operator()(const Indices &I, const Indices &J) const {
      FMATVEC_ASSERT(I.max()<rows(), AT);
      FMATVEC_ASSERT(J.max()<cols(), AT);

      Matrix<General,Ref,Ref,AT> A(I.size(),J.size(),NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
	A.e(i,j) = e(I[i],J[j]);

      return A;
    }

  template <class AT> template <class Type, class Row, class Col>
    inline void Matrix<General,Ref,Ref,AT>::set(const Indices &I, const Indices &J, const Matrix<Type,Row,Col,AT> &A) {
      FMATVEC_ASSERT(I.max()<rows(), AT);
      FMATVEC_ASSERT(J.max()<cols(), AT);
      FMATVEC_ASSERT(I.size()==A.rows(), AT);
      FMATVEC_ASSERT(J.size()==A.cols(), AT);
      for(int i=0; i<I.size(); i++)
	for(int j=0; j<J.size(); j++)
	  e(I[i],J[j]) = A.e(i,j);
    }

  template <class AT> template<class Row>
    inline void Matrix<General,Ref,Ref,AT>::set(const Indices &I, int j, const Vector<Row,AT> &x) {
      FMATVEC_ASSERT(I.max()<rows(), AT);
      FMATVEC_ASSERT(j<cols(), AT);
      FMATVEC_ASSERT(I.size()==x.size(), AT);
      for(int i=0; i<I.size(); i++)
	e(I[i],j) = x.e(i);
    }

  template <class AT>
    inline void Matrix<General,Ref,Ref,AT>::ref(Matrix<General,Ref,Ref,AT> &A, const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J) {
      FMATVEC_ASSERT(I.end()<A.rows(), AT);
      FMATVEC_ASSERT(J.end()<A.cols(), AT);
      m=I.size();
      n=J.size();
      memory = A.memory;
      ele = A.elePtr(I.start(),J.start());
      lda = A.lda;
    }

  template <class AT>
    inline void Matrix<General,Ref,Ref,AT>::ref(Matrix<Symmetric,Ref,Ref,AT> &A, const fmatvec::Range<Var,Var> &I, const fmatvec::Range<Var,Var> &J) {
      FMATVEC_ASSERT(I.end()<A.size(), AT);
      FMATVEC_ASSERT(J.end()<A.size(), AT);
      FMATVEC_ASSERT(I.start()>=J.start(), AT);
      FMATVEC_ASSERT(J.end()<=I.start(), AT);
      m=I.size();
      n=J.size();
      memory = A.memory;
      ele = A.elePtr(I.start(),J.start());
      lda = A.lda;
    }

  template <class AT>
    int Matrix<General,Ref,Ref,AT>::countElements() const {
    int k = 0;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if (fabs(e(i, j)) > 1e-16 || i==j) {
          k++;
        }
      }
    }
    return k;
  }

  /// @endcond

}

#endif
