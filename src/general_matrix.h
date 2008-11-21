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

#ifndef general_matrix_h
#define general_matrix_h

#include "memory.h"
#include "types.h"
#include "index.h"

namespace fmatvec {

  template <class AT> class Vector;
  template <class AT> class RowVector;
  template <class AT> class SquareMatrix;
  //template <class AT> class Matrix<Symmetric, AT>;

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type General and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <class AT> class Matrix<General, AT> {

    public:

 /// @cond NO_SHOW

      friend class Matrix<Symmetric, AT>;
      
      template <class T> friend const Matrix<General, T>  trans(const Matrix<General, T> &A);

      template <class T> friend Matrix<General, T>  trans(Matrix<General, T> &A);

    protected:

      Memory<AT> memory;
      AT *ele;
      int m;
      int n;
      int lda;
      bool tp;

      template <class Type> void deepCopy(const Matrix<Type, AT> &A); 

      const AT* elePtr(int i, int j) const {
	return tp ? ele+i*lda+j : ele+i+j*lda; 
      };

      AT* elePtr(int i, int j) {
	return tp ? ele+i*lda+j : ele+i+j*lda; 
      };

      Matrix(int m_, int n_, int lda_, bool tp_, Memory<AT> memory_, const AT* ele_) : memory(memory_), ele((AT*)ele_), m(m_), n(n_), lda(lda_), tp(tp_) {
      }

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      Matrix() : memory(), ele(0), m(0), n(0), lda(0) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x n.
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \remark The matrix will be initialised to
       * zero by default. This default behavior can be changed by defining 
       * FMATVEC_NO_INITIALIZATION.
       * */
      Matrix(int m_, int n_) : memory(m_*n_), ele((AT*)memory.get()), m(m_), n(n_), lda(m_), tp(false) {  

#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x n. The matrix will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * */
      Matrix(int m_, int n_, Initialization ini, const AT &a=0) : memory(m_*n_), ele((AT*)memory.get()), m(m_), n(n_), lda(m_), tp(false) {  

	if(ini == INIT) {
	  AT *el=ele;
	  for(int i=0; i<m*n; i++)
	    *el++=a;
	} else if(ini == EYE ) {
	  for(int i=0; i<m; i++) {
	    for(int j=0; j<n; j++) {
	      if (i==j) operator()(i,i) = 1;
	      else operator()(i,j) = 0;
	    }
	  }
	}  
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      Matrix(const Matrix<General, AT> &A) : memory(A.memory), ele(A.ele) ,m(A.m), n(A.n), lda(A.lda), tp(A.tp) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a matrix of size m x n with the pyhsical memory given by \em ele_.
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \param ele_ The physical memory the matrix will point to.
       * */
      Matrix(int m_, int n_, AT* ele_) : memory(), ele(ele_), m(m_), n(n_), lda(m_), tp(false) { 
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<General,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const char *str);

      /*! \brief Destructor. 
       * */
      ~Matrix() {
      }

      /*! \brief Matrix resizing. 
       *
       * Resizes the matrix to size m x n.    
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \return A reference to the calling matrix.
       * \remark The matrix will be initialised to
       * zero by default. To change this behavior, define
       * FMATVEC_NO_INITIALIZATION.
       * */
      Matrix<General, AT>& resize(int m_, int n_) {
	m=m_;n=n_;
	lda=m;
	tp = false;
	memory.resize(m*n);
	ele = (AT*)memory.get();

#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
	return *this;
      }

      /*! \brief Matrix resizing. 
       *
       * Resizes the matrix to size m x n. The matrix will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param m_ The number of rows.
       * \param n_ The number of columns.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * \return A reference to the calling matrix.
       * */
      Matrix<General, AT>& resize(int m_, int n_, Initialization ini, const AT &a=0) {
	m=m_;n=n_;
	lda=m;
	tp = false;
	memory.resize(m*n);
	ele = (AT*)memory.get();

	if(ini == INIT)
	  init(a);
	else if(ini == EYE ) {
	  for(int i=0; i<m; i++) {
	    for(int j=0; j<n; j++) {
	      if (i==j) operator()(i,i) = 1;
	      else operator()(i,j) = 0;
	    }
	  }
	}

	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      Matrix<General, AT>& operator=(const Matrix<General, AT> &A) {
#ifndef FMATVEC_NO_DEEP_ASSIGNMENT 
	return operator<<(A);
#else
	return operator>>(A);
#endif
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> Matrix<General, AT>& operator<<(const Matrix<T, AT> &A);

      /*! \brief Reference operator
       *
       * References the matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      Matrix<General, AT>& operator>>(const Matrix<General, AT> &A);

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
	assert(i<m);
	assert(j<n);
#endif

	return tp ? ele[i*lda+j] : ele[i+j*lda];
      };

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<m);
	assert(j<n);
#endif

	return tp ? ele[i*lda+j] : ele[i+j*lda];//  return ele[i*lda+j*ldb];
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

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix.
       * */
      int ldim() const {return lda;};

      /*! \brief Transposed status.
       *
       * Returns boolean transposed status.
       * \return True, if the matrix is in transposed form, false otherwise.
       * */
      bool transposed() const {return tp;};

      /*! \brief Transposed status.
       *
       * Returns the blas-conform transposed status.
       * \return CblasTrans if the matrix is in transposed form, CblasNoTrans
       * otherwise. 
       * */
      const enum CBLAS_TRANSPOSE blasTrans() const {
	return (tp)? CblasTrans : CblasNoTrans;
      };

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasColMajor.
       * */
      const enum CBLAS_ORDER blasOrder() const {
	return  CblasColMajor;
      };

      /*! \brief Submatrix operator.
       *
       * Returns a submatrix of the calling matrix. 
       * \attention The submatrix and the
       * calling matrix will share the same physical memory.
       * \param i1 The starting row. 
       * \param j1 The starting column.
       * \param i2 The ending row.
       * \param j2 The ending column.
       * \return A submatrix of the calling matrix.
       * */
      Matrix<General, AT> operator()(int i1, int j1, int i2, int j2);

      /*! \brief Submatrix operator.
       *
       * See operator()(int,int,int,int);
       * */
      const Matrix<General, AT> operator()(int i1, int j1, int i2, int j2) const;

      /*! \brief Submatrix operator.
       *
       * Returns a submatrix of the calling matrix. 
       * For example
       * \code 
       * B = A(Index(1,2),Index(2,4));
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
       * \param I Index containing the starting and the ending row. 
       * \param J Index containing the starting and the ending column. 
       * \return A submatrix of the calling matrix.
       * */
      Matrix<General, AT> operator()(const Index &I, const Index &J);

      /*! \brief Submatrix operator.
       *
       * See operator()(const Index&, const Index&)
       * */
      const Matrix<General, AT> operator()(const Index &I, const Index &J) const;

      /*! \brief Submatrix operator.
       *
       * Returns a square submatrix of the calling matrix. 
       * \attention The submatrix and the
       * calling matrix will share the same physical memory.
       * \param I Index containing the starting and the ending row. 
       * \return A submatrix of the calling matrix.
       * */
      SquareMatrix<AT> operator() (const Index &I);

      /*! \brief Submatrix operator.
       *
       * See operator()(const Index&)
       * */
      const SquareMatrix<AT> operator()(const Index &I) const;

      /*! \brief Column operator.
       *
       * Returns a vector containing the i-th column of the calling matrix. 
       * \attention The vector and the calling matrix will share the same physical memory.
       * \param i The column, that will be returned.  
       * \return A vector containing the i-th column of the calling matrix.
       * */
      Vector<AT> col(int i);

      /*! \brief Column operator.
       *
       * see col(int)
       * */
      const Vector<AT> col(int i) const;

      /*! \brief Row operator.
       *
       * Returns a RowVector containing the i-th row of the calling matrix. 
       * \attention The rowvector and the calling matrix will share the same physical memory.
       * \param i The row, that will be returned. 
       * \return A rowvector containing the i-th row of the calling matrix.
       * */
      RowVector<AT> row(int i); 

      /*! \brief Row operator.
       *
       * see row(int)
       * */
      const RowVector<AT> row(int i) const;

      /*! \brief Matrix unsizing.
       *
       * Resizes the matrix to size zero.  
       * \return A reference to the calling matrix.
       * */
      Matrix<General, AT>& resize() {m=0;n=0;return *this;};

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      Matrix<General, AT> copy() const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      Matrix<General, AT>& init(const AT &a);

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      operator vector<vector<AT> >();

      /*! \brief std::vector<std::vector<AT> > Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT> > object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT> > the matrix will be initialized with. 
       * */
      Matrix(vector<vector<AT> > m);

  };
  // ------------------------- Constructors -------------------------------------
  // ----------------------------------------------------------------------------

  template <class AT>
    Matrix<General, AT>& Matrix<General, AT>::operator>>(const Matrix<General, AT> &A) { 

      if(m==0) {
	m=A.m; 
	n=A.n;
      } else {
#ifdef FMATVEC_SIZE_CHECK
	assert(m == A.m);
	assert(n == A.n);
#endif
      }

      memory = A.memory;
      ele = A.ele;
      lda = A.lda;
      tp = A.tp; 

      return *this;
    }

   template <class AT> template< class Type>
    Matrix<General, AT>& Matrix<General, AT>::operator<<(const Matrix<Type, AT> &A) { 

      if(A.rows() == 0 || A.cols() == 0)
	return *this;

      if(m==0) {
	m = A.rows(); 
	n = A.cols();
	lda = m;
	tp = false;
	memory.resize(m*n);
	ele = (AT*)memory.get();
      } else {
#ifdef FMATVEC_SIZE_CHECK
	assert(m == A.rows());
	assert(n == A.cols());
#endif
      }

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<General, AT>&  Matrix<General, AT>::init(const AT& val) {

      for(int i=0; i<rows(); i++) 
	for(int j=0; j<cols(); j++) 
	  operator()(i,j) = val;

      return *this;
    }

  template <class AT>
    Matrix<General, AT>  Matrix<General, AT>::operator()(int i1, int j1, int i2, int j2) {
      return operator()(Index(i1,i2),Index(j1,j2));
    }

  template <class AT>
    const Matrix<General, AT>  Matrix<General, AT>::operator()(int i1, int j1, int i2, int j2) const {
      return operator()(Index(i1,i2),Index(j1,j2));
    }

  template <class AT> Matrix<General, AT> Matrix<General, AT>::operator()(const Index &I, const Index &J) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
    assert(I.end()<m);
    assert(J.end()<n);
#endif
    return Matrix<General, AT>(I.end()-I.start()+1,J.end()-J.start()+1,lda,tp,memory,elePtr(I.start(),J.start()));
  }

  template <class AT> const Matrix<General, AT> Matrix<General, AT>::operator()(const Index &I, const Index &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
    assert(I.end()<m);
    assert(J.end()<n);
#endif
    return Matrix<General, AT>(I.end()-I.start()+1,J.end()-J.start()+1,lda,tp,memory,elePtr(I.start(),J.start()));
  }
  
 template <class AT> const SquareMatrix<AT> Matrix<General, AT>::operator()(const Index &I) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
    assert(I.end()<m);
#endif
    return SquareMatrix<AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start(),I.start()));
  }
 template <class AT> SquareMatrix<AT> Matrix<General, AT>::operator()(const Index &I) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
    assert(I.end()<m);
#endif
    return SquareMatrix<AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start(),I.start()));
  }

  template <class AT>
    Vector<AT>  Matrix<General, AT>::col(int i) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<n);
#endif

      return Vector<AT>(m,lda,tp,memory,elePtr(0,i));
    }

  template <class AT>
    const Vector<AT>  Matrix<General, AT>::col(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<n);
#endif

      return Vector<AT>(m,lda,tp,memory,elePtr(0,i));
    }

  template <class AT>
    RowVector<AT>  Matrix<General, AT>::row(int i) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<m);
#endif

      return RowVector<AT>(n,lda,tp,memory,elePtr(i,0));
    }

  template <class AT>
    const RowVector<AT>  Matrix<General, AT>::row(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<m);
#endif

      return RowVector<AT>(n,lda,tp,memory,elePtr(i,0));
    }

  template <class AT>
    Matrix<General, AT> Matrix<General, AT>::copy() const {

      Matrix<General, AT> A(m,n);
      A.deepCopy(*this);

      return A;
    }

  template <class AT> template <class Type>
    void Matrix<General, AT>::deepCopy(const Matrix<Type, AT> &A) { 
      for(int i=0; i<m; i++) 
	for(int j=0; j<n; j++)
          operator()(i,j) = A.operator()(i,j);
    }

  template <class AT>
    Matrix<General, AT>::operator vector<vector<AT> >() {
      vector<vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=operator()(r,c);
      }
      return ret;
    }

  template <class AT>
    Matrix<General, AT>::Matrix(vector<vector<AT> > m) : memory(m.size()*m[0].size()), ele((AT*)memory.get()), m(m.size()), n(m[0].size()), lda(m.size()), tp(false) {
#ifndef FMATVEC_NO_INITIALIZATION 
      init(0);
#endif
      for(int r=0; r<rows(); r++) {
        assert(m[r].size()==cols());
        for(int c=0; c<cols(); c++)
          operator()(r,c)=m[r][c];
      }
    }

   /// @cond NO_SHOW
   
  template<> template<>
  void Matrix<General, double>::deepCopy(const Matrix<General, double> &A);

  template<> template<>
  void Matrix<General, double>::deepCopy(const Matrix<Symmetric, double> &A);

   /// @endcond

}

#endif
