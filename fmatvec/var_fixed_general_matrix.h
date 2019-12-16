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

#ifndef var_fixed_general_matrix_h
#define var_fixed_general_matrix_h

#include "types.h"
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
  template <int N, class AT> class Matrix<fmatvec::General,fmatvec::Var,fmatvec::Fixed<N>,AT> {

    public:
      typedef AT value_type;

 /// @cond NO_SHOW

    protected:

      int M{0};

      AT *ele;

      template <class Type, class Row, class Col> inline void deepCopy(const Matrix<Type,Row,Col,AT> &A); 
      //inline void deepCopy(const Matrix<General,Var,Fixed<N>,AT> &A); 

 /// @endcond
 
    public:

      explicit Matrix() :  ele(nullptr) { }

//      template<class Ini=All<AT> >
//      Matrix(int m, Ini ini=All<AT>()) :  M(m), ele(new AT[M*N]) {
//        init(ini);
//      }
//      template<class Ini=All<AT> >
//      Matrix(int m, int n, Ini ini=All<AT>()) :  M(m), ele(new AT[M*N]) {
//        init(ini);
//      }

      explicit Matrix(int m, Noinit) : M(m), ele(new AT[M*N]) { }
      explicit Matrix(int m, Init ini=INIT, const AT &a=AT()) : M(m), ele(new AT[M*N]) { init(a); }
      explicit Matrix(int m, Eye ini, const AT &a=1) : M(m), ele(new AT[M*N]) { init(ini,a); }
      explicit Matrix(int m, int n, Noinit) : M(m), ele(new AT[M*N]) { assert(n==N); }
      explicit Matrix(int m, int n, Init ini=INIT, const AT &a=AT()) : M(m), ele(new AT[M*N]) {  assert(n==N); init(a); }
      explicit Matrix(int m, int n, Eye ini, const AT &a=1) : M(m), ele(new AT[M*N]) {  assert(n==N); init(ini,a); }

      /*! \brief Copy Constructor
       *
       * Constructs a copy to the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      Matrix(const Matrix<General,Var,Fixed<N>,AT> &A) : M(A.M), ele(new AT[M*N]) {
	deepCopy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy to the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Row, class Col>
      Matrix(const Matrix<General,Row,Col,AT> &A) : M(A.rows()), ele(new AT[M*N]) {
	assert(A.cols() == N); 
	deepCopy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy to the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : M(A.rows()), ele(new AT[M*N]) {
	assert(A.cols() == N); 
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

      Matrix<General,Var,Fixed<N>,AT>& resize(int m, Noinit) {
	delete[] ele;
	M = m;
	ele = new AT[M*N];
        return *this;
      }

      Matrix<General,Var,Fixed<N>,AT>& resize(int m, Init ini=INIT, const AT &a=AT()) { return resize(m,Noinit()).init(a); }

      Matrix<General,Var,Fixed<N>,AT>& resize(int m, Eye ini, const AT &a=1) { return resize(m,Noinit()).init(ini,a); }

      //! Resize a var-fixed matrix.
      //! Throw if the fixed dimension is different and resize the var dimension
      void resize(int m, int n) {
        if(n!=N)
          throw std::runtime_error("A var-fixed matrix can only be resized in the first dimension.");
        resize(m);
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Var,Fixed<N>,AT>& operator=(const Matrix<General,Var,Fixed<N>,AT> &A) {
        assert(M == A.rows());
        deepCopy(A);
        return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      template <class Type, class Row, class Col>
      inline Matrix<General,Var,Fixed<N>,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
        assert(N == A.cols()); 
        assert(M == A.rows());
        deepCopy(A);
        return *this;
      }

      /*! \brief Matrix replacement
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template <class Type, class Row, class Col>
      inline Matrix<General,Var,Fixed<N>,AT>& replace(const Matrix<Type,Row,Col,AT> &A); 

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
	assert(i<M);
	assert(j<N);

	return e(i,j);
      };

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
	assert(i>=0);
	assert(j>=0);
	assert(i<M);
	assert(j<N);

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

      //! A var-fixed matrix is stored in fortran-storage order -> return always false
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

      inline const Matrix<General,Var,Var,AT> operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const;

      inline const RowVector<Fixed<N>,AT> row(int i) const;
      inline const Vector<Var,AT> col(int j) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Var,Fixed<N>,AT>& init(const AT &val=AT()); 
      inline Matrix<General,Var,Fixed<N>,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<General,Var,Fixed<N>,AT>& init(Eye, const AT &val=1);
      inline Matrix<General,Var,Fixed<N>,AT>& init(Noinit, const AT &a=AT()) { return *this; }

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

      inline const Matrix<General,Fixed<N>,Var,AT> T() const;

      template<class Row> inline void set(int j, const Vector<Row,AT> &x);

      template<class Col> inline void set(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void set(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

      template<class Row> inline void add(int j, const Vector<Row,AT> &x);

      template<class Col> inline void add(int i, const RowVector<Col,AT> &x);

      template<class Type, class Row, class Col> inline void add(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A);

  };

  template <int N, class AT> 
    Matrix<General,Var,Fixed<N>,AT>::Matrix(const std::string &strs) :  ele(0) {
      std::istringstream iss(strs);
      iss>>*this;

      // check end of stream
      iss>>std::ws;
      if(!iss.eof())
        throw std::runtime_error("Input not fully read.");
    }
  template <int N, class AT> Matrix<General,Var,Fixed<N>,AT>::Matrix(const char * strs) :
    Matrix<General,Var,Fixed<N>,AT>::Matrix(std::string(strs)) {}

  template <int N, class AT> template< class Type, class Row, class Col>
    inline Matrix<General,Var,Fixed<N>,AT>& Matrix<General,Var,Fixed<N>,AT>::replace(const Matrix<Type,Row,Col,AT> &A) { 

      assert(N == A.cols());

      if(M!=A.rows()) {
        delete[] ele;
        M = A.rows();
        ele = new AT[M*N];
      }

      deepCopy(A);

      return *this;
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>&  Matrix<General,Var,Fixed<N>,AT>::init(const AT &val) {
      for(int i=0; i<M*N; i++) 
        e(i) = val;
      return *this;
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>&  Matrix<General,Var,Fixed<N>,AT>::init(Eye eye, const AT &val) {
      for(int i=0; i<M; i++)
        for(int j=0; j<N; j++)
          e(i,j) = (i==j) ? val : 0;
      return *this;
    }

  template <int N, class AT>
    inline const Matrix<General,Var,Var,AT> Matrix<General,Var,Fixed<N>,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
      assert(I.end()<M);
      assert(J.end()<N);
      Matrix<General,Var,Var,AT> A(I.end()-I.start()+1,J.end()-J.start()+1,NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          A.e(i,j) = e(I.start()+i,J.start()+j);

      return A;
    }

  template <int N, class AT>
    inline const RowVector<Fixed<N>,AT> Matrix<General,Var,Fixed<N>,AT>::row(int i) const {

      assert(i>=0);
      assert(i<M);

      RowVector<Fixed<N>,AT> x(NONINIT);

      for(int j=0; j<N; j++)
        x.e(j) = e(i,j);

      return x;

    }

  template <int N, class AT>
    inline const Vector<Var,AT> Matrix<General,Var,Fixed<N>,AT>::col(int j) const {

      assert(j>=0);
      assert(j<N);

      Vector<Var,AT> x(M,NONINIT);

      for(int i=0; i<M; i++)
        x.e(i) = e(i,j);

      return x;

    }

  template <int N, class AT>
    inline const Matrix<General,Fixed<N>,Var,AT> Matrix<General,Var,Fixed<N>,AT>::T() const {
      Matrix<General,Fixed<N>,Var,AT> A(rows(),NONINIT);
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int N, class AT> template <class Row>
    inline void Matrix<General,Var,Fixed<N>,AT>::set(int j, const Vector<Row,AT> &x) {
      assert(j<cols());
      assert(rows()==x.size());
      for(int i=0; i<rows(); i++)
        e(i,j) = x.e(i);
    }

  template <int N, class AT> template <class Col>
    inline void Matrix<General,Var,Fixed<N>,AT>::set(int i, const RowVector<Col,AT> &x) {
      assert(i<rows());
      assert(cols()==x.size());
      for(int j=0; j<cols(); j++)
        e(i,j) = x.e(j);
    }

  template <int N, class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Var,Fixed<N>,AT>::set(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      assert(I.end()<rows());
      assert(J.end()<cols());
      assert(I.size()==A.rows());
      assert(J.size()==A.cols());

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <int N, class AT> template <class Row>
    inline void Matrix<General,Var,Fixed<N>,AT>::add(int j, const Vector<Row,AT> &x) {
      assert(j<cols());
      assert(rows()==x.size());
      for(int i=0; i<rows(); i++)
        e(i,j) += x.e(i);
    }

  template <int N, class AT> template <class Col>
    inline void Matrix<General,Var,Fixed<N>,AT>::add(int i, const RowVector<Col,AT> &x) {
      assert(i<rows());
      assert(cols()==x.size());
      for(int j=0; j<cols(); j++)
        e(i,j) += x.e(j);
    }

  template <int N, class AT> template<class Type, class Row, class Col>
    inline void Matrix<General,Var,Fixed<N>,AT>::add(const Range<Var,Var> &I, const Range<Var,Var> &J, const Matrix<Type,Row,Col,AT> &A) {

      assert(I.end()<rows());
      assert(J.end()<cols());
      assert(I.size()==A.rows());
      assert(J.size()==A.cols());

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>::operator std::vector<std::vector<AT> >() const {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>::Matrix(std::vector<std::vector<AT> > m) : M(m.size()), ele(new AT[M*N]) {
      if(m.size()<1)
        throw std::runtime_error("The input must have at least 1 row.");
      if(m[0].size() != N)
        throw std::runtime_error("The input has "+std::to_string(m[0].size())+" columns but "+std::to_string(N)+" columns are required.");
      for(int r=0; r<rows(); r++) {
        if(static_cast<int>(m[r].size())!=cols())
          throw std::runtime_error("The rows of the input have different length.");
        for(int c=0; c<cols(); c++)
          e(r,c)=m[r][c];
      }
    }

  /// @cond NO_SHOW

  template <int N, class AT> template <class Type, class Row, class Col>
    inline void Matrix<General,Var,Fixed<N>,AT>::deepCopy(const Matrix<Type,Row,Col,AT> &A) { 
      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++)
          e(i,j) = A.e(i,j);
    }

//  template<int N, class AT>
//    inline void Matrix<General,Var,Fixed<N>,AT>::deepCopy(const Matrix<General,Var,Fixed<N>,AT> &A) {
//      for(int i=0; i<M*N; i++) 
//        e(i) = A.e(i);
//    }

  /// @endcond

}

#endif
