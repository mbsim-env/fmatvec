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

#ifndef fixed_var_general_matrix_h
#define fixed_var_general_matrix_h

#include "types.h"
#include "range.h"
#include <cstdlib>
#include <stdexcept>

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type GeneralFixed and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <int M, class AT> class Matrix<General,Fixed<M>,Var,AT> {

    public:

      typedef AT value_type;

 /// @cond NO_SHOW

    protected:

      int N{0};

      AT *ele;

      template <class Type, class Row, class Col> inline void deepCopy(const Matrix<Type,Row,Col,AT> &A); 

 /// @endcond
 
    public:

     /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      explicit Matrix() :  ele(nullptr) { }

//      template<class Ini=All<AT> >
//      Matrix(int n, Ini ini=All<AT>()) :  N(n), ele(new AT[M*N]) {
//        init(ini);
//      }
//      template<class Ini=All<AT> >
//      Matrix(int m, int n, Ini ini=All<AT>()) :  N(n), ele(new AT[M*N]) {
//        init(ini);
//      }

      explicit Matrix(int n, Noinit) : N(n), ele(new AT[M*N]) { }
      explicit Matrix(int n, Init ini=INIT, const AT &a=AT()) : N(n), ele(new AT[M*N]) { init(a); }
      explicit Matrix(int n, Eye ini, const AT &a=1) : N(n), ele(new AT[M*N]) { init(ini,a); }
      explicit Matrix(int m, int n, Noinit) : N(n), ele(new AT[M*N]) { assert(m==M); }
      explicit Matrix(int m, int n, Init ini=INIT, const AT &a=AT()) : N(n), ele(new AT[M*N]) { assert(m==M); init(a); }
      explicit Matrix(int m, int n, Eye ini, const AT &a=1) : N(n), ele(new AT[M*N]) { assert(m==M); init(ini,a); }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      Matrix(const Matrix<General,Fixed<M>,Var,AT> &A) : N(A.N), ele(new AT[M*N]) {
	deepCopy(A);
      }

      template<class Row, class Col>
      Matrix(const Matrix<General,Row,Col,AT> &A) : N(A.cols()), ele(new AT[M*N]) {

	deepCopy(A);
      }

      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : N(A.cols()), ele(new AT[M*N]) {

#ifndef FMATVEC_NO_SIZE_CHECK
	assert(A.rows() == M); 
#endif

	deepCopy(A);
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<GeneralFixed,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const std::string &strs);
      Matrix(const char *strs);

      /*! \brief Destructor. 
       * */
      ~Matrix() {
	delete[] ele;
      }

      Matrix<General,Fixed<M>,Var,AT>& resize() { 
	delete[] ele;
	N = 0;
	ele = nullptr;
        return *this;
      }

      Matrix<General,Fixed<M>,Var,AT>& resize(int n, Noinit) { 
	delete[] ele;
	N = n;
	ele = new AT[M*N];
        return *this;
      }

      Matrix<General,Fixed<M>,Var,AT>& resize(int n, Init ini=INIT, const AT &a=AT()) { return resize(n,Noinit()).init(a); }

      Matrix<General,Fixed<M>,Var,AT>& resize(int n, Eye ini, const AT &a=1) { return resize(n,Noinit()).init(ini,a); }

      //! Resize a fixed-var matrix.
      //! Throw if the fixed dimension is different and resize the var dimension
      void resize(int m, int n) {
        if(m!=M)
          throw std::runtime_error("A fixed-var matrix can only be resized in the second dimension.");
        resize(n);
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Fixed<M>,Var,AT>& operator=(const Matrix<General,Fixed<M>,Var,AT> &A);

      template <class Type, class Row, class Col>
      inline Matrix<General,Fixed<M>,Var,AT>& operator=(const Matrix<Type,Row,Col,AT> &A);

      template<class Type, class Row, class Col>
      inline Matrix<General,Fixed<M>,Var,AT>& operator<<(const Matrix<Type,Row,Col,AT> &A);

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
	return ele[i*N+j];
      };

      /*! \brief Element operator
       *
       * See e(int,int) 
       * */
      const AT& e(int i, int j) const {
	return ele[i*N+j];
      };

      AT& e(int i) {
	return ele[i];
      };

      /*! \brief Element operator
       *
       * See e(int,int) 
       * */
      const AT& e(int i) const {
	return ele[i];
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
      int ldim() const {return M;};

      //! The storage format for a fixed-var matrix is fortran-storage order -> transposed is always false
      bool transposed() const {
	return false;
      };

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
	return  CblasColMajor;
      };

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
      //inline Matrix<GeneralVar, AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J);

      /*! \brief Submatrix operator.
       *
       * See operator()(const Range<Var,Var>&, const Range<Var,Var>&)
       * */
      inline const Matrix<General,Var,Var,AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const;

      template <int M1, int M2>
      inline const Matrix<General,Fixed<M2-M1+1>,Var,AT> operator()(const Range<Fixed<M1>,Fixed<M2> > &I, const Range<Var,Var> &J) const;

      inline const RowVector<Var,AT> row(int i) const;
      inline const Vector<Fixed<M>,AT> col(int j) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Fixed<M>,Var,AT>& init(const AT &val=AT()); 
      inline Matrix<General,Fixed<M>,Var,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<General,Fixed<M>,Var,AT>& init(Eye, const AT &val=1);
      inline Matrix<General,Fixed<M>,Var,AT>& init(Noinit, const AT &a=AT()) { return *this; }

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
      inline Matrix(std::vector<std::vector<AT> > m);

      inline const Matrix<General,Var,Fixed<M>,AT> T() const;

      template<class Row> inline void set(int j, const Vector<Row,AT> &x);

      template<class Col> inline void set(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void set(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      template<class Row> inline void add(int j, const Vector<Row,AT> &x);

      template<class Col> inline void add(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void add(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

  };

  template <int M, class AT> 
    Matrix<General,Fixed<M>,Var,AT>::Matrix(const std::string &strs) :  ele(nullptr) {
      std::istringstream iss(strs);
      iss>>*this;

      // check end of stream
      iss>>std::ws;
      if(!iss.eof())
        throw std::runtime_error("Input not fully read.");
    }
  template <int M, class AT> Matrix<General,Fixed<M>,Var,AT>::Matrix(const char *strs) :
    Matrix<General,Fixed<M>,Var,AT>::Matrix(std::string(strs)) {}

  template <int M, class AT> template< class Type, class Row, class Col>
    inline Matrix<General,Fixed<M>,Var,AT>& Matrix<General,Fixed<M>,Var,AT>::operator=(const Matrix<Type,Row,Col,AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(M == A.rows()); 
#endif
      if(!ele) {
        delete[] ele;
        N = A.cols(); 
        ele = new AT[M*N];
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(N == A.cols());
#endif
      }

      deepCopy(A);

      return *this;
    }

  template <int M, class AT>
    inline Matrix<General,Fixed<M>,Var,AT>& Matrix<General,Fixed<M>,Var,AT>::operator=(const Matrix<General,Fixed<M>,Var,AT> &A) { 

      if(!ele) {
        delete[] ele;
        N = A.cols(); 
        ele = new AT[M*N];
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(N == A.cols());
#endif
      }

      deepCopy(A);

      return *this;
    }

  template <int M, class AT> template< class Type, class Row, class Col>
    inline Matrix<General,Fixed<M>,Var,AT>& Matrix<General,Fixed<M>,Var,AT>::operator<<(const Matrix<Type,Row,Col,AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(M == A.rows());
#endif
      if(N!=A.cols()) {
        delete[] ele;
        N = A.cols();
        ele = new AT[M*N];
      }

      deepCopy(A);

      return *this;
    }

  template <int M, class AT>
    inline Matrix<General,Fixed<M>,Var,AT>&  Matrix<General,Fixed<M>,Var,AT>::init(const AT &val) {
      for(int i=0; i<M*N; i++) 
        e(i) = val;
      return *this;
    }

  template <int M, class AT>
    inline Matrix<General,Fixed<M>,Var,AT>&  Matrix<General,Fixed<M>,Var,AT>::init(Eye eye, const AT &val) {
      for(int i=0; i<M; i++)
        for(int j=0; j<N; j++)
          e(i,j) = (i==j) ? val : 0;
      return *this;
    }

  template <int M, class AT> template <int M1, int M2>
    inline const Matrix<General,Fixed<M2-M1+1>,Var,AT> Matrix<General,Fixed<M>,Var,AT>::operator()(const Range<Fixed<M1>,Fixed<M2> > &I, const Range<Var,Var> &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(M2<M);
      assert(J.end()<N);
#endif
      Matrix<General,Fixed<M2-M1+1>,Var,AT> A(J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(M1+i,J.start()+j);

      return A;
    }

  template <int M, class AT> 
    inline const Matrix<General,Var,Var,AT> Matrix<General,Fixed<M>,Var,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
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

  template <int M, class AT>
    inline const RowVector<Var,AT> Matrix<General,Fixed<M>,Var,AT>::row(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<M);
#endif

      RowVector<Var,AT> x(N,NONINIT);

      for(int j=0; j<N; j++)
        x.e(j) = e(i,j);

      return x;

    }

  template <int M, class AT>
    inline const Vector<Fixed<M>,AT> Matrix<General,Fixed<M>,Var,AT>::col(int j) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

      Vector<Fixed<M>,AT> x(NONINIT);

      for(int i=0; i<M; i++)
        x.e(i) = e(i,j);

      return x;

    }

  template <int M, class AT>
    inline const Matrix<General,Var,Fixed<M>,AT> Matrix<General,Fixed<M>,Var,AT>::T() const {
      Matrix<General,Var,Fixed<M>,AT> A(cols(),NONINIT);
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int M, class AT> template <class Row>
    inline void Matrix<General,Fixed<M>,Var,AT>::set(int j, const Vector<Row,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j<cols());
      assert(rows()==x.size());
#endif
      for(int i=0; i<rows(); i++)
        e(i,j) = x.e(i);
    }

  template <int M, class AT> template <class Col>
    inline void Matrix<General,Fixed<M>,Var,AT>::set(int i, const RowVector<Col,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i<rows());
      assert(cols()==x.size());
#endif
      for(int j=0; j<cols(); j++)
        e(i,j) = x.e(j);
    }

  template <int M, class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Fixed<M>,Var,AT>::set(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

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

  template <int M, class AT> template <class Row>
    inline void Matrix<General,Fixed<M>,Var,AT>::add(int j, const Vector<Row,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j<cols());
      assert(rows()==x.size());
#endif
      for(int i=0; i<rows(); i++)
        e(i,j) += x.e(i);
    }

  template <int M, class AT> template <class Col>
    inline void Matrix<General,Fixed<M>,Var,AT>::add(int i, const RowVector<Col,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i<rows());
      assert(cols()==x.size());
#endif
      for(int j=0; j<cols(); j++)
        e(i,j) += x.e(j);
    }

  template <int M, class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Fixed<M>,Var,AT>::add(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

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

  template <int M, class AT>
    inline Matrix<General,Fixed<M>,Var,AT>::operator std::vector<std::vector<AT> >()  const{
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <int M, class AT>
    inline Matrix<General,Fixed<M>,Var,AT>::Matrix(std::vector<std::vector<AT> > m) :
      Matrix<General,Fixed<M>,Var,AT>(!m.empty() ? m[0].size() : 0) {
      if(m.size() != M)
        throw std::runtime_error("The input has "+std::to_string(m.size())+" rows but "+std::to_string(M)+" rows are required.");
      if(static_cast<int>(m[0].size()) != N)
        throw std::runtime_error("The input has "+std::to_string(m[0].size())+" columns but "+std::to_string(N)+" columns are required.");
      for(int r=0; r<rows(); r++) {
        if(static_cast<int>(m[r].size())!=cols())
          throw std::runtime_error("The rows of the input have different length.");
        for(int c=0; c<cols(); c++)
          e(r,c)=m[r][c];
      }
    }

  /// @cond NO_SHOW

  template <int M, class AT> template <class Type, class Row, class Col>
    inline void Matrix<General,Fixed<M>,Var,AT>::deepCopy(const Matrix<Type,Row,Col,AT> &A) { 
      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++)
          e(i,j) = A.e(i,j);
    }

  /// @endcond

}

#endif
