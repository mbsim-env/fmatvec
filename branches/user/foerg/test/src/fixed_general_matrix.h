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
#include <stdlib.h>

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

 /// @cond NO_SHOW

    protected:

      AT ele[M][N];

      template <class Type, class Row, class Col> inline void deepCopy(const Matrix<Type,Row,Col,AT> &A); 
      template<class Row> inline void deepCopy(const Matrix<Symmetric,Row,Row,AT> &A);

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * */
      Matrix() {
#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      Matrix(Initialization ini, const AT &a=0) {  

	if(ini == INIT) {
          init(a);
	} else if(ini == EYE ) {
	  for(int i=0; i<M; i++) {
	    for(int j=0; j<N; j++) {
	      if (i==j) e(i,j) = 1;
	      else e(i,j) = 0;
	    }
	  }
	}
      }

      Matrix(NOINIT) {}
      Matrix(int m, int n, NOINIT) {}
      Matrix(SCALAR, const AT &a=0) { init(a); }
      Matrix(int m, int n, SCALAR, const AT &a=0) { init(a); }

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
      Matrix(const char *str);

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

      inline const RowVector<Fixed<N>,AT> row(int j) const;
      inline const Vector<Fixed<M>,AT> col(int i) const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<General,Fixed<M>,Fixed<N>,AT>& init(const AT &a);

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
      Matrix(std::vector<std::vector<AT> > m);

      inline const Matrix<General,Fixed<M>,Fixed<N>,AT> T() const;

      inline void set(int j, const Vector<Fixed<M>,AT> &x);

      inline void set(int i, const RowVector<Fixed<N>,AT> &x);

  };

  template <int M, int N, class AT> 
    Matrix<General,Fixed<M>,Fixed<N>,AT>::Matrix(const char *strs) {
      // if 'strs' is a single scalar, surround it first with '[' and ']'.
      // This is more Matlab-like, because e.g. '5' and '[5]' is just the same.
      // (This functionallitiy is needed e.g. by MBXMLUtils (OpenMBV,MBSim))
      std::istringstream iss(strs);
      char c;
      iss>>c;
      if(c=='[') iss.str(strs);
      else iss.str(std::string("[")+strs+"]");

      int m=0, n=0;
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
          m++;
        }
        else if(c==',')
          n++;
        c='0';
      } while(iss);

      n++; m++;
      assert(m==M);
      assert(n==N);
      iss.clear();
      iss.seekg(0);
      iss >> c;
      for(int i=0; i<M; i++)
        for(int j=0; j<N; j++) {
          iss >> e(i,j);
          iss >> c;
        }
    }

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
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>& Matrix<General,Fixed<M>,Fixed<N>,AT>::init(const AT& val) {

      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++) 
          e(i,j) = val;

      return *this;
    }

  template <int M, int N, class AT> 
    inline const Matrix<General,Var,Var,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::operator()(const Range<Var,Var> &I, const Range<Var,Var> &J) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<M);
      assert(J.end()<N);
#endif
      Matrix<General,Var,Var,AT> A(I.end()-I.start()+1,J.end()-J.start()+1,NOINIT());

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

      RowVector<Fixed<N>,AT> x(0,NOINIT());

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

      Vector<Fixed<M>,AT> x(0,NOINIT());

      for(int i=0; i<M; i++)
        x.e(i) = e(i,j);

      return x;

    }

  template <int M, int N, class AT>
    inline const Matrix<General,Fixed<M>,Fixed<N>,AT> Matrix<General,Fixed<M>,Fixed<N>,AT>::T() const {
      Matrix<General,Fixed<N>,Fixed<M>,AT> A(0,0,NOINIT());
      for(int i=0; i<N; i++)
        for(int j=0; j<M; j++)
          A.e(i,j) = e(j,i);
      return A;
    }

  template <int M, int N, class AT>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::set(int j, const Vector<Fixed<M>,AT> &x) {
      for(int i=0; i<M; i++)
        e(i,j) = x.e(i);
    }

  template <int M, int N, class AT>
    inline void Matrix<General,Fixed<M>,Fixed<N>,AT>::set(int i, const RowVector<Fixed<N>,AT> &x) {
      for(int j=0; j<N; j++)
        e(i,j) = x.e(j);
    }

  template <int M, int N, class AT>
    inline Matrix<General,Fixed<M>,Fixed<N>,AT>::operator std::vector<std::vector<AT> >() {
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

}

#endif
