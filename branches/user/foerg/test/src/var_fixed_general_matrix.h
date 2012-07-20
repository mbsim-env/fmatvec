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
#include <stdlib.h>

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type GeneralFixed and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <int N, class AT> class Matrix<General,Var,Fixed<N>,AT> {

    public:

 /// @cond NO_SHOW

    protected:

      int M;

      AT *ele;

      template <class Type, class Row, class Col> inline void deepCopy(const Matrix<Type,Row,Col,AT> &A); 
      //inline void deepCopy(const Matrix<General,Var,Fixed<N>,AT> &A); 

 /// @endcond
 
    public:

      Matrix() : M(0), ele(0) { }

//      template<class Ini=All<AT> >
//      Matrix(int m, Ini ini=All<AT>()) :  M(m), ele(new AT[M*N]) {
//        init(ini);
//      }
//      template<class Ini=All<AT> >
//      Matrix(int m, int n, Ini ini=All<AT>()) :  M(m), ele(new AT[M*N]) {
//        init(ini);
//      }

      Matrix(int m) : M(m), ele(new AT[M*N]) { init(0); }
      Matrix(int m, const Noinit &ini) : M(m), ele(new AT[M*N]) { }
      Matrix(int m, const All<AT> &ini) : M(m), ele(new AT[M*N]) { init(ini); }
      Matrix(int m, const Eye<AT> &ini) : M(m), ele(new AT[M*N]) { init(ini); }
      Matrix(int m, int n, const Noinit &ini) : M(m), ele(new AT[M*N]) { }

      // For compatibility
      Matrix(int m, const All<AT> &ini, const AT &a) :  M(m), ele(new AT[M*N]) { init(a); }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      Matrix(const Matrix<General,Var,Fixed<N>,AT> &A) : M(A.M), ele(new AT[M*N]) {
	deepCopy(A);
      }

      template<class Row, class Col>
      Matrix(const Matrix<General,Row,Col,AT> &A) : M(A.rows()), ele(new AT[M*N]) {

	deepCopy(A);
      }

      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : M(A.rows()), ele(new AT[M*N]) {

#ifndef FMATVEC_NO_SIZE_CHECK
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
       * Matrix<GeneralFixed,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const char *str);

      /*! \brief Destructor. 
       * */
      ~Matrix() {
	delete[] ele;
      }

//      template<class Ini=All<AT> >
//      Matrix<General,Var,Fixed<N>,AT>& resize(int m=0, Ini ini=All<AT>()) {
//	delete[] ele;
//	M=m;
//	ele = new AT[M*N];
//        init(ini);
//        return *this;
//      }

      Matrix<General,Var,Fixed<N>,AT>& resize(int m=0) { return resize(m,All<AT>()); }

      template<class Ini>
      Matrix<General,Var,Fixed<N>,AT>& resize(int m, const Ini &ini) {
	delete[] ele;
	M=m;
	ele = new AT[M*N];
        init(ini);
        return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Var,Fixed<N>,AT>& operator=(const Matrix<General,Var,Fixed<N>,AT> &A);

      template <class Type, class Row, class Col>
      inline Matrix<General,Var,Fixed<N>,AT>& operator=(const Matrix<Type,Row,Col,AT> &A);
      
      template <class Type, class Row, class Col>
      inline Matrix<General,Var,Fixed<N>,AT>& operator<<(const Matrix<Type,Row,Col,AT> &A);

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

      inline const RowVector<Fixed<N>,AT> row(int j) const;
      inline const Vector<Var,AT> col(int j) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Var,Fixed<N>,AT>& init(const AT &a=0); 
      inline Matrix<General,Var,Fixed<N>,AT>& init(const All<AT> &all) { return init(all.a); }
      inline Matrix<General,Var,Fixed<N>,AT>& init(Eye<AT> eye);
      inline Matrix<General,Var,Fixed<N>,AT>& init(Noinit) { return *this; }

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      inline operator std::vector<std::vector<AT> >();

      /*! \brief std::vector<std::vector<AT> > Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT> > object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT> > the matrix will be initialized with. 
       * */
      inline Matrix(std::vector<std::vector<AT> > m);

      inline const Matrix<General,Fixed<N>,Var,AT> T() const;

      inline void set(int j, const RowVector<Fixed<N>,AT> &x);

      template<int K> 
	inline void set(const Index &I, const Index &J, const Matrix<General,Var,Fixed<K>,AT> &A);

      template<int K> 
	inline void add(const Index &I, const Index &J, const Matrix<General,Var,Fixed<K>,AT> &A);


  };

  template <int N, class AT> 
    Matrix<General,Var,Fixed<N>,AT>::Matrix(const char *strs) {
      // if 'strs' is a single scalar, surround it first with '[' and ']'.
      // This is more Matlab-like, because e.g. '5' and '[5]' is just the same.
      // (This functionallitiy is needed e.g. by MBXMLUtils (OpenMBV,MBSim))
      std::istringstream iss(strs);
      char c;
      iss>>c;
      if(c=='[') iss.str(strs);
      else iss.str(std::string("[")+strs+"]");

      M=0;
      int n=0;
      int buf=0;
      iss >> c;
      AT x;
      do {
        iss >> x;
        iss >> c;
        if(c==';') {
          if(buf)
            assert(buf == n);

          buf=n;
          n=0;
          M++;
        }
        else if(c==',')
          n++;
        c='0';
      } while(iss);

      n++; M++;
      assert(n==N);
      ele = new AT[M*N];
      iss.clear();
      iss.seekg(0);
      iss >> c;
      for(int i=0; i<M; i++)
        for(int j=0; j<N; j++) {
          iss >> e(i,j);
          iss >> c;
        }
    }

  template <int N, class AT> template< class Type, class Row, class Col>
    inline Matrix<General,Var,Fixed<N>,AT>& Matrix<General,Var,Fixed<N>,AT>::operator=(const Matrix<Type,Row,Col,AT> &A) { 

#ifndef FMATVEC_RESIZE_VOID
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(M == A.rows()); 
      assert(N == A.cols());
#endif
#else
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(N == A.cols()); 
#endif
      if(M==0) {
        delete[] ele;
        M = A.rows(); 
        ele = new AT[M*N];
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(M == A.rows());
#endif
      }
#endif

      deepCopy(A);

      return *this;
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>& Matrix<General,Var,Fixed<N>,AT>::operator=(const Matrix<General,Var,Fixed<N>,AT> &A) { 

#ifndef FMATVEC_RESIZE_VOID
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(M == A.rows()); 
#endif
#else
      if(M==0) {
        delete[] ele;
        M = A.rows(); 
        ele = new AT[M*N];
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(M == A.rows());
#endif
      }
#endif

      deepCopy(A);

      return *this;
    }

  template <int N, class AT> template< class Type, class Row, class Col>
    inline Matrix<General,Var,Fixed<N>,AT>& Matrix<General,Var,Fixed<N>,AT>::operator<<(const Matrix<Type,Row,Col,AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(N == A.cols());
#endif
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
    inline Matrix<General,Var,Fixed<N>,AT>&  Matrix<General,Var,Fixed<N>,AT>::init(Eye<AT> eye) {
      for(int i=0; i<M; i++) {
        e(i,i) = eye.a;
        for(int j=0; j<i; j++) {
          e(i,j) = e(j,i) = 0; 
        }
      }
      return *this;
    }

  template <int N, class AT>
    inline const Matrix<General,Var,Var,AT> Matrix<General,Var,Fixed<N>,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
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

  template <int N, class AT>
    inline const RowVector<Fixed<N>,AT> Matrix<General,Var,Fixed<N>,AT>::row(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i>=0);
      assert(i<M);
#endif

      RowVector<Fixed<N>,AT> x(NONINIT);

      for(int j=0; j<N; j++)
        x.e(j) = e(i,j);

      return x;

    }

  template <int N, class AT>
    inline const Vector<Var,AT> Matrix<General,Var,Fixed<N>,AT>::col(int j) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(j>=0);
      assert(j<N);
#endif

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

  template <int N, class AT>
    inline void Matrix<General,Var,Fixed<N>,AT>::set(int j, const RowVector<Fixed<N>,AT> &x) {
      for(int i=0; i<N; i++)
        e(i,j) = x.e(i);
    }

  template <int N, class AT> template<int K>
    inline void Matrix<General,Var,Fixed<N>,AT>::set(const Index &I, const Index &J, const Matrix<General,Var,Fixed<K>,AT> &A) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(J.end()<N);
      assert(I.size()==A.rows());
      assert(J.size()==K);
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) = A.e(ii,jj);
    }

  template <int N, class AT> template<int K>
    inline void Matrix<General,Var,Fixed<N>,AT>::add(const Index &I, const Index &J, const Matrix<General,Var,Fixed<K>,AT> &A) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(J.end()<N);
      assert(I.size()==A.rows());
      assert(J.size()==K);
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        for(int j=J.start(), jj=0; j<=J.end(); j++, jj++)
          e(i,j) += A.e(ii,jj);
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

  template <int N, class AT>
    inline Matrix<General,Var,Fixed<N>,AT>::Matrix(std::vector<std::vector<AT> > m) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(m.size() == M);
      assert(m[0].size() == N);
#endif
      for(int r=0; r<rows(); r++) {
        assert(m[r].size()==cols());
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
