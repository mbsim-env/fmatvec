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

#ifndef fixed_general_matrix_h
#define fixed_general_matrix_h

#include "types.h"
#include <cstdlib>
#include <stdexcept>

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type General, Fixed and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <int M, int N, class AT> class Matrix<General,Fixed<M>,Fixed<N>,AT> {

    public:

      typedef AT value_type;

 /// @cond NO_SHOW

    protected:

      AT ele[M][N];

      template <class Type, class Row, class Col> inline void deepCopy(const Matrix<Type,Row,Col,AT> &A); 
      template<class Row> inline void deepCopy(const Matrix<Symmetric,Row,Row,AT> &A);

 /// @endcond
 
    public:

// Works with -std=gnu++0x only
//      template<class Ini=All<AT> >
//      Matrix(Ini ini=All<AT>()) {
//        init(ini);
//      }
//      template<class Ini=All<AT> >
//      Matrix(int m_, int n_, Ini ini=All<AT>()) {
//        init(ini);
//      }

      Matrix(Noinit) { }
      Matrix(Init ini=INIT, const AT &a=AT()) { init(a); }
      Matrix(Eye ini, const AT &a=1) { init(ini,a); }
      Matrix(int m, int n, Noinit) { }
      Matrix(int m, int n, Init ini=INIT, const AT &a=AT()) { init(a); }
      Matrix(int m, int n, Eye ini, const AT &a=1) { init(ini,a); }

      template<class Row, class Col>
      Matrix(const Matrix<General,Row,Col,AT> &A) {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == M); 
	assert(A.cols() == N);
#endif

	deepCopy(A);
      }

      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) {

#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == M); 
	assert(A.cols() == N);
#endif

	deepCopy(A);
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<General, Fixed,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const std::string &strs);
      Matrix(const char *strs);

      template <class Type, class Row, class Col>
      inline Matrix<General,Fixed<M>,Fixed<N>,AT>& operator=(const Matrix<Type,Row,Col,AT> &A);
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
	assert(j<N);
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
	assert(j<N);
#endif

	return e(i,j);
      };

      AT& e(int i, int j) {
	return ele[i][j];
      };

      /*! \brief Element operator
       *
       * See e(int,int) 
       * */
      const AT& e(int i, int j) const {
	return ele[i][j];
      };

      /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return &(ele[0][0]);};

      /*! \brief Pointer operator
       *
       * See operator()() 
       * */
      const AT* operator()() const {return &(ele[0][0]);};

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix.
       * */
      int rows() const {return M;};

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix.
       * */
      int cols() const {return N;};

      /*! \brief Leading dimension.
       *
       * \return The leading dimension of the matrix
       * */
      int ldim() const {return N;};

      //! A fixed matrix is stored in c-storage order -> transposed is always true
      bool transposed() const {
	return true;
      };

      //! Resize a fixed matrix.
      //! Do nothing for the fixed dimension and throw on any other dimension.
      void resize(int m, int n) {
        if(m!=M || n!=N)
          throw std::runtime_error("A fixed matrix cannot be resized.");
      }

      /*! \brief Transposed status.
       *
       * Returns the blas-conform transposed status.
       * \return CblasTrans if the matrix is in transposed form, CblasNoTrans
       * otherwise. 
       * */
      const CBLAS_TRANSPOSE blasTrans() const {
	return CblasNoTrans;
      };

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasColMajor.
       * */
      const CBLAS_ORDER blasOrder() const {
	return  CblasRowMajor;
      };

//      /*! \brief Submatrix operator.
//       *
//       * Returns a submatrix of the calling matrix. 
//       * \attention The submatrix and the
//       * calling matrix will share the same physical memory.
//       * \param i1 The starting row. 
//       * \param j1 The starting column.
//       * \param i2 The ending row.
//       * \param j2 The ending column.
//       * \return A submatrix of the calling matrix.
//       * */
//      inline Matrix<General, AT> operator()(int i1, int j1, int i2, int j2);
//
//      /*! \brief Submatrix operator.
//       *
//       * See operator()(int,int,int,int);
//       * */
//      inline const Matrix<General, AT> operator()(int i1, int j1, int i2, int j2) const;

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
       * \param I Range<Var,Var> containing the starting and the ending row. 
       * \param J Range<Var,Var> containing the starting and the ending column. 
       * \return A submatrix of the calling matrix.
       * */
//      inline Matrix<General, Var, AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J);

      /*! \brief Submatrix operator.
       *
       * See operator()(const Range<Var,Var>&, const Range<Var,Var>&)
       * */
      inline const Matrix<General,Var,Var,AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const;
      
      template <int M1, int M2, int N1, int N2>
      inline const Matrix<General,Fixed<M2-M1+1>,Fixed<N2-N1+1>,AT> operator()(const Range<Fixed<M1>,Fixed<M2> > &I, const Range<Fixed<N1>,Fixed<N2> > &J) const;
      template <int M1, int M2>
      inline const Matrix<General,Fixed<M2-M1+1>,Var,AT> operator()(const Range<Fixed<M1>,Fixed<M2> > &I, const Range<Var,Var > &J) const;

      template <int N1, int N2>
      inline const Matrix<General,Var,Fixed<N2-N1+1>,AT> operator()(const Range<Var,Var> &I, const Range<Fixed<N1>,Fixed<N2> > &J) const;

      inline const RowVector<Fixed<N>,AT> row(int i) const;
      inline const Vector<Fixed<M>,AT> col(int j) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Fixed<M>,Fixed<N>,AT>& init(const AT &val=AT()); 
      inline Matrix<General,Fixed<M>,Fixed<N>,AT>& init(Init, const AT &a=AT()) { return init(a); };
      inline Matrix<General,Fixed<M>,Fixed<N>,AT>& init(Eye, const AT &val=1);
      inline Matrix<General,Fixed<M>,Fixed<N>,AT>& init(Noinit, const AT &a=AT()) { return *this; }

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >() const;

      /*! \brief std::vector<std::vector<AT> > Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT> > object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT> > the matrix will be initialized with. 
       * */
      Matrix(std::vector<std::vector<AT> > m);

      inline const Matrix<General,Fixed<N>,Fixed<M>,AT> T() const;

      /*!
       * \brief set Column j of matrix to given Vector
       */
      template<class Row> inline void set(int j, const Vector<Row,AT> &x);

      /*!
       * \brief set Row i of matrix to given RowVector
       */
      template<class Col> inline void set(int i, const RowVector<Col,AT> &x);

      /*!
       * \brief set the submatrix - specified by the Range operators - to the values of A
       * \param I Range of starting and ending row
       * \param J Range of starting and ending column
       * \param A Matrix with the values the submatrix should take
       */
      template<class Type, class Row, class Col> inline void set(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      /*!
       * \brief add to Column j of the matrix the given Vector
       */
      template<class Row> inline void add(int j, const Vector<Row,AT> &x);

      /*!
       * \brief add to Row i of the matrix the given RowVector
       */
      template<class Col> inline void add(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void add(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

  };

  template <int M, int N, class AT> 
    Matrix<General,Fixed<M>,Fixed<N>,AT>::Matrix(const std::string &strs) {
      std::istringstream iss(strs);
      iss>>*this;

      // check end of stream
      iss>>std::ws;
      if(!iss.eof())
        throw std::runtime_error("Input not fully read.");
    }
  template <int M, int N, class AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::Matrix(const char *strs) :
    Matrix<General,Fixed<M>,Fixed<N>,AT>::Matrix(std::string(strs)) {}

  template <int M, int N, class AT> template<class Type, class Row, class Col>
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>& Matrix<General,Fixed<M>,Fixed<N>,AT>::operator=(const Matrix<Type,Row,Col,AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.rows() == M); 
      assert(A.cols() == N);
#endif

      deepCopy(A);

      return *this;
    }

  template <int M, int N, class AT>
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>& Matrix<General,Fixed<M>,Fixed<N>,AT>::init(const AT &val) {
       for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++) 
          e(i,j) = val;
      return *this;
    }

  template <int M, int N, class AT>
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>& Matrix<General,Fixed<M>,Fixed<N>,AT>::init(Eye, const AT &val) {
      for(int i=0; i<M; i++)
        for(int j=0; j<N; j++)
          e(i,j) = (i==j) ? val : 0;
      return *this;
    }

  template <int M, int N, class AT> 
    inline const Matrix<General,Var,Var,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(J.end()<N);
#endif
      Matrix<General,Var,Var,AT> A(I.end()-I.start()+1,J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

 template <int M, int N, class AT> template <int M1, int M2, int N1, int N2>
    inline const Matrix<General,Fixed<M2-M1+1>,Fixed<N2-N1+1>,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::operator()(const Range<Fixed<M1>,Fixed<M2> > &I, const Range<Fixed<N1>,Fixed<N2> > &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(M2<M);
      assert(N2<N);
#endif
      Matrix<General,Fixed<M2-M1+1>,Fixed<N2-N1+1>,AT> A(NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(M1+i,J.start()+j);

      return A;
    }

 template <int M, int N, class AT> template <int M1, int M2>
    inline const Matrix<General,Fixed<M2-M1+1>,Var,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::operator()(const Range<Fixed<M1>,Fixed<M2> > &I, const Range<Var,Var> &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(M2<M);
      assert(J.end()<N);
#endif
      Matrix<General,Fixed<M2-M1+1>,Var,AT> A(J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

  template <int M, int N, class AT> template <int N1, int N2>
    inline const Matrix<General,Var,Fixed<N2-N1+1>,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::operator()(const Range<Var,Var> &I, const Range<Fixed<N1>,Fixed<N2> > &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(N2<N);
#endif
      Matrix<General,Var,Fixed<N2-N1+1>,AT> A(I.end()-I.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

  template <int M, int N, class AT>
    inline const RowVector<Fixed<N>,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::row(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<M);
#endif

      RowVector<Fixed<N>,AT> x(NONINIT);

      for(int j=0; j<N; j++)
        x.e(j) = e(i,j);

      return x;

    }

  template <int M, int N, class AT>
    inline const Vector<Fixed<M>,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::col(int j) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

      Vector<Fixed<M>,AT> x(NONINIT);

      for(int i=0; i<M; i++)
        x.e(i) = e(i,j);

      return x;

    }

  template <int M, int N, class AT>
    inline const Matrix<General,Fixed<N>,Fixed<M>,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::T() const {
      Matrix<General,Fixed<N>,Fixed<M>,AT> A(NONINIT);
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int M, int N, class AT> template <class Row>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::set(int j, const Vector<Row,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j<cols());
      assert(rows()==x.size());
#endif
      for(int i=0; i<rows(); i++)
        e(i,j) = x.e(i);
    }

  template <int M, int N, class AT> template <class Col>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::set(int i, const RowVector<Col,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i<rows());
      assert(cols()==x.size());
#endif
      for(int j=0; j<cols(); j++)
        e(i,j) = x.e(j);
    }

  template <int M, int N, class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::set(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<rows());
      assert(J.end()<cols());
      assert(I.size()==A.rows());
      assert(J.size()==A.cols());
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <int M, int N, class AT> template <class Row>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::add(int j, const Vector<Row,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j<cols());
      assert(rows()==x.size());
#endif
      for(int i=0; i<rows(); i++)
        e(i,j) += x.e(i);
    }

  template <int M, int N, class AT> template <class Col>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::add(int i, const RowVector<Col,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i<rows());
      assert(cols()==x.size());
#endif
      for(int j=0; j<cols(); j++)
        e(i,j) += x.e(j);
    }

  template <int M, int N, class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::add(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<rows());
      assert(J.end()<cols());
      assert(I.size()==A.rows());
      assert(J.size()==A.cols());
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <int M, int N, class AT>
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>::operator std::vector<std::vector<AT> >() const {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <int M, int N, class AT>
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>::Matrix(std::vector<std::vector<AT> > m) {
      if(m.size() != M)
        throw std::runtime_error("The input has "+std::to_string(m.size())+" rows but "+std::to_string(M)+" rows are required.");
      if(m[0].size() != N)
        throw std::runtime_error("The input has "+std::to_string(m[0].size())+" columne but "+std::to_string(N)+" columns are required.");
      for(int r=0; r<rows(); r++) {
        if(static_cast<int>(m[r].size())!=cols())
          throw std::runtime_error("The rows of the input have different length.");
        for(int c=0; c<cols(); c++)
          e(r,c)=m[r][c];
      }
    }

  /// @cond NO_SHOW

  template <int M, int N, class AT> template <class Type, class Row, class Col>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::deepCopy(const Matrix<Type,Row,Col,AT> &A) { 
      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++)
          e(i,j) = A.e(i,j);
    }

  template<int M, int N, class AT> template<class Row>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::deepCopy(const Matrix<Symmetric,Row,Row,AT> &A) {
      for(int i=0; i<A.size(); i++) {
        e(i,i) = A.ej(i,i);
        for(int j=i+1; j<A.size(); j++)
          e(i,j) = e(j,i) = A.ej(i,j);
      }
    }

  /// @endcond

  /*! A matrix representing a rotation matrix.
   * This class is simply derived from general square matrix from which all member are inherited.
   * We must use this to be able to distinguish it from a general square matrix at compile time.
   * The constructors of general matrix are redefined here since these are not inherited.
   * Moreover we can overload some member functions here which capitalize the special
   * properties of rotation matrices, link "RotMat inv() { return trans(*this); }."
   */
  template<int M, class AT>
  class Matrix<Rotation,Fixed<M>,Fixed<M>,AT> : public Matrix<General,Fixed<M>,Fixed<M>,AT> {
    public:
      typedef AT value_type;
      // Constructors are not inherited. Hence we must redefine all ctors here.

      Matrix(Noinit ini) { }
      Matrix(Init ini=INIT, const AT &a=0) { this->init(a); }
      Matrix(Eye ini, const AT &a=1) { this->init(ini,a); }
      Matrix(int m, int n, Noinit ini) { }
      Matrix(int m, int n, Init ini, const AT &a=0) { this->init(a); }
      Matrix(int m, int n, Eye ini, const AT &a=1) { this->init(ini,a); }
      explicit Matrix(const Matrix<General,Fixed<M>,Fixed<M>,AT>&  A) { this->deepCopy(A); }

      template<class Row>
      Matrix(const Matrix<Symmetric,Row,Row,AT> &A) {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.size() == M); 
#endif
	this->deepCopy(A);
      }

      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == M); 
	assert(A.cols() == M);
#endif
	this->deepCopy(A);
      }

      // TODO: we should add some overloaded member functions here which capitalize the special
      // properties of rotation matrices, link "RotMat inv() { return trans(*this); }"
  };


}

#endif
