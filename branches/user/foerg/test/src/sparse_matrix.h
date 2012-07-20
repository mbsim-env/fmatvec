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

#ifndef sparse_matrix_h
#define sparse_matrix_h

#include "matrix.h"
#include "square_matrix.h"
#include "types.h"
#include "_memory.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for sparse quadratic matrices.
   *
   * Template class Matrix with shape type Sparse and atomic type
   * AT. The matrix ist stored in compressed row-wise skyline format, BUT the
   * diagonal elements are ALWAYS stored (even if they are zero) as the
   * FIRST entry of each row.
   * The template parameter AT defines the atomic type of the
   * matrix. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT>
    class Matrix<Sparse,Ref,Ref,AT> {

      protected:

    /// @cond NO_SHOW
    
	Memory<AT> memEle;
	Memory<int> memI, memJ;
	AT *ele;
	int *I, *J;
	int m, n, k;

	void deepCopy(const Matrix<Sparse,Ref,Ref,AT> &x);

	void deepCopy(const SquareMatrix<Ref,AT> &A);

    /// @endcond


      public:

	/*! \brief Standard constructor
	 *
	 * Constructs a matrix with no size. 
	 * */
	Matrix() : memEle(), memI(), memJ(), ele(0), I(0), J(0), m(0), n(0), k(0) { }

//        template<class Ini=All<AT> >
//          Matrix(int n_, Ini ini=All<AT>()) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) {
//            init(ini);
//          }
//        template<class Ini=All<AT> >
//          Matrix(int m_, int n_, Ini ini) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) {  
//            init(ini);
//          }

        Matrix(int n_) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { init(0); }
        Matrix(int n_, int k_) : memEle(k_), memI(n_+1), memJ(k_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(k_) {  init(0); }
        Matrix(int n_, const Noinit &ini) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { }
        Matrix(int n_, const All<AT> &ini) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { init(ini); }
        Matrix(int n_, const Eye<AT> &ini) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { init(ini); }

        Matrix(int m_, int n_, const Noinit &ini) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { }

	/*! \brief Copy Constructor
	 *
	 * Constructs a reference to the matrix \em A.
	 * \attention The physical memory of the matrix \em A will not be copied, only
	 * referenced.
	 * \param A The matrix that will be referenced.
	 * */
	Matrix(const Matrix<Sparse,Ref,Ref,AT> &A) : memEle(A.memEle), memI(A.memI), memJ(A.memJ), ele(A.ele), I(A.I), J(A.J), m(A.m), n(A.n), k(A.k) {
	}

	/*! \brief Destructor. 
	 * */
	~Matrix() {
	}

	/*! \brief Matrix resizing. 
	 *
	 * Resizes the matrix to size n x n. The matrix will be initialized to the value given by \em a
	 * (default 0), if ini is set to INIT. If init is set to NONINIT, the
	 * matrix will not be initialized.
	 * \param n_ The number of columns/rows.
	 * \param k_ The number of nonzero elements.
	 * \return A reference to the calling matrix.
	 * */
//        template<class Ini=All<AT> >
//          Matrix<Sparse,Ref,Ref,AT>& resize(int n_=0, int k_=0, Ini ini=All<AT>()) {
//            m=n_;n=n_;k=k_;
//            memEle.resize(k);
//            memI.resize(m+1);
//            memJ.resize(k);
//            ele = (AT*)memEle.get();
//            I = (int*)memI.get();
//            J = (int*)memJ.get();
//            init(ini);
//
//            return *this;
//          }

        Matrix<Sparse,Ref,Ref,AT>& resize(int n=0, int k=0) { return resize(n,k,All<AT>()); }

        template<class Ini>
          Matrix<Sparse,Ref,Ref,AT>& resize(int n_, int k_, const Ini &ini) {
            m=n_;n=n_;k=k_;
            memEle.resize(k);
            memI.resize(m+1);
            memJ.resize(k);
            ele = (AT*)memEle.get();
            I = (int*)memI.get();
            J = (int*)memJ.get();
            init(ini);

            return *this;
          }

	/*! \brief Copy operator
	 *
	 * Copies the sparse matrix given by \em A.
	 * \param A The matrix to be copied. 
	 * \return A reference to the calling matrix.
	 * */
	inline Matrix<Sparse,Ref,Ref,AT>& operator<<(const Matrix<Sparse,Ref,Ref,AT> &A);

	/*! \brief Reference operator
	 *
	 * References the sparse matrix given by \em A.
	 * \param A The matrix to be referenced. 
	 * \return A reference to the calling matrix.
	 * */
	inline Matrix<Sparse,Ref,Ref,AT>& operator>>(const Matrix<Sparse,Ref,Ref,AT> &A);

	/*! \brief Element operator
	 *
	 * See operator<<(const Matrix<Sparse,Ref,Ref,T>&) 
	 * */
	inline Matrix<Sparse,Ref,Ref,AT>& operator<<(const SquareMatrix<Ref,AT> &A);

	/*! \brief Assignment operator
	 *
	 * Copies the sparse matrix given by \em A by calling operator<<().
	 * \param A The matrix to be assigned. 
	 * \return A reference to the calling matrix.
	 * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
	 * \sa operator<<(), operator>>()
	 * */
	inline Matrix<Sparse,Ref,Ref,AT>& operator=(const Matrix<Sparse,Ref,Ref,AT> &A);

	/*! \brief Pointer operator
	 *
	 * See operator()() 
	 * */
	const AT* operator()() const {
	  return ele;
	};

	/*! \brief Pointer operator.
	 *
	 * Returns the pointer to the first element.
	 * \return The pointer to the first element.
	 * */
	AT* operator()() {
	  return ele;
	};

	/*! \brief Pointer operator.
	 *
	 * See Ip() 
	 * */
	const int* Ip() const {
	  return I;
	};

	/*! \brief Pointer operator.
	 *
	 * \todo Docu
	 * */
	int* Ip() {
	  return I;
	};

	/*! \brief Pointer operator.
	 *
	 * See Jp() 
	 * */
	const int* Jp() const {
	  return J;
	};

	/*! \brief Pointer operator.
	 *
	 * \todo Docu
	 * */
	int* Jp() {
	  return J;
	};

	/*! \brief Number of rows.
	 *
	 * \return The number of rows of the matrix
	 * */
	int rows() const {return m;};

	/*! \brief Number of columns.
	 *
	 * \return The number of columns of the matrix
	 * */
	int cols() const {return n;};

	/*! \brief Matrix duplicating.
	 *
	 * The calling matrix returns a \em deep copy of itself.  
	 * \return The duplicate.
	 * */
	Matrix<Sparse,Ref,Ref,AT> copy() const;

	/*! \brief Initialization.
	 *
	 * Initializes all elements of the calling matrix with 
	 * the value given by \em a.
	 * \param a Value all elements will be initialized with.
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Sparse,Ref,Ref,AT>& init(const AT &a);
        inline Matrix<Sparse,Ref,Ref,AT>& init(const All<AT> &all);
        inline Matrix<Sparse,Ref,Ref,AT>& init(Eye<AT> eye);
        inline Matrix<Sparse,Ref,Ref,AT>& init(Noinit) { return *this; }

    };
  // ------------------------- Constructors -------------------------------------
  // ----------------------------------------------------------------------------

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::operator>>(const Matrix<Sparse,Ref,Ref,AT> &A) { 
      m=A.m;
      n=A.n;
      k=A.k;
      memEle = A.memEle;
      memI = A.memI;
      memJ = A.memJ;
      ele = A.ele;
      I = A.I;
      J = A.J;
      k = A.k;

      return *this;
    }

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::operator=(const Matrix<Sparse,Ref,Ref,AT> &A) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(m == A.m);
      assert(n == A.n);
      assert(k == A.k);
#endif

      deepCopy(A);

      return *this;
    }

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::operator<<(const Matrix<Sparse,Ref,Ref,AT> &A) { 

      if(m!=A.m || n!=A.n || k!=A.k) {
	m = A.m;
	n = A.n;
	k = A.k;
	memEle.resize(k);
	memI.resize(m+1);
	memJ.resize(k);
	ele = (AT*)memEle.get();
	I = (int*)memI.get();
	J = (int*)memJ.get();
      }

      deepCopy(A);

      return *this;
    }

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::operator<<(const SquareMatrix<Ref,AT> &A) { 

      if(m!=A.rows() || n!=A.cols()) {
        m = A.rows();
        n = A.cols();
        k = countElements(A);
        memEle.resize(k);
        memI.resize(m+1);
        memJ.resize(k);
        ele = (AT*)memEle.get();
        I = (int*)memI.get();
        J = (int*)memJ.get();
      } else {
        int k_ = countElements(A);
        if(k != k_) {
          k = k_;
          memEle.resize(k);
          memJ.resize(k);
          ele = (AT*)memEle.get();
          J = (int*)memJ.get();
        }
      }

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<Sparse,Ref,Ref,AT> Matrix<Sparse,Ref,Ref,AT>::copy() const {

      Matrix<Sparse,Ref,Ref,AT> A(m,NONINIT);

      A.deepCopy(*this);

      return A;
    }

  template <class AT>
    void Matrix<Sparse,Ref,Ref,AT>::deepCopy(const Matrix<Sparse,Ref,Ref,AT> &A) { 
      for(int i=0; i<=m; i++) {
	I[i] = A.I[i];
      }
      for(int i=0; i<k; i++) {
	ele[i] = A.ele[i];
	J[i] = A.J[i];
      }
    }

  template <class AT> void Matrix<Sparse,Ref,Ref,AT>::deepCopy(const SquareMatrix<Ref,AT> &A) { 
      int k=0;
      int i;
      for(i=0; i<A.size(); i++) {
	ele[k]=A(i,i);
	J[k]=i;
	I[i]=k++;
	for(int j=0; j<i; j++) {
	  // TODO eps
	  if(fabs(A(i,j))>1e-16) {
	    ele[k]=A(i,j);
	    J[k++]=j;
	  }
	}
	for(int j=i+1; j<A.size(); j++) {
	  // TODO eps
	  if(fabs(A(i,j))>1e-16) {
	    ele[k]=A(i,j);
	    J[k++]=j;
	  }
	}
      }
      if(n) I[i]=k;
    }

  template <class AT>
    Matrix<Sparse,Ref,Ref,AT>&  Matrix<Sparse,Ref,Ref,AT>::init(const AT& val) {
      return init(All<AT>(val));
    }

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>&  Matrix<Sparse,Ref,Ref,AT>::init(const All<AT> &all) {
      for(int i=0; i<k; i++) {
	ele[i] = all.a;
      }
      return *this;
    }

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>&  Matrix<Sparse,Ref,Ref,AT>::init(Eye<AT> eye) {
      throw;
      return *this;
    }

/* template <class AT> template <class Type> void Matrix<Sparse,Ref,Ref,AT>::deepCopy(const Matrix<Type, AT> &A) { 
      k=0;
      AT *eleBuf = new AT[m*n];
      int *JBuf = new int[m*n];
      int i;
      for(i=0; i<A.rows(); i++) {
	eleBuf[k]=A(i,i);
	JBuf[k]=i;
	I[i]=k++;
	for(int j=0; j<i; j++) {
	  // TODO eps
	  if(fabs(A(i,j))>1e-16) {
	    eleBuf[k]=A(i,j);
	    JBuf[k++]=j;
	  }
	}
	for(int j=i+1; j<A.cols(); j++) {
	  // TODO eps
	  if(fabs(A(i,j))>1e-16) {
	    eleBuf[k]=A(i,j);
	    JBuf[k++]=j;
	  }
	}
      }
      I[i]=k;
      memEle.resize(k);
      memJ.resize(k);
      ele = (AT*)memEle.get();
      J = (int*)memJ.get();
      for(int i=0; i<k; i++) {
	ele[i] = eleBuf[i];
	J[i] = JBuf[i];
      }

      delete [] eleBuf;
      delete [] JBuf;
    }

    */

  //template <> extern void Matrix<Sparse,Ref,Ref,double>::deepCopy(const Matrix<Sparse,Ref,Ref,double> &A);

}

#endif
