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
	int *I{nullptr}, *J{nullptr};
	int m{0}, n{0}, k{0};

        template <class Type, class Row, class Col> inline Matrix<Sparse,Ref,Ref,AT>& copy(const Matrix<Type,Row,Col,AT> &A);
	template <class Row> inline Matrix<Sparse,Ref,Ref,AT>& copy(const Matrix<Symmetric,Row,Row,AT> &A);
	inline Matrix<Sparse,Ref,Ref,AT>& copy(const Matrix<Sparse,Ref,Ref,AT> &A);

    /// @endcond


      public:

        typedef AT value_type;

	/*! \brief Standard constructor
	 *
	 * Constructs a matrix with no size. 
	 * */
	explicit Matrix() : memEle(),  ele(0) { }

//        template<class Ini=All<AT>>
//          Matrix(int n_, Ini ini=All<AT>()) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) {
//            init(ini);
//          }
//        template<class Ini=All<AT>>
//          Matrix(int m_, int n_, Ini ini) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) {  
//            init(ini);
//          }

        explicit Matrix(int m_, int n_, int k_, Noinit) : memEle(k_), memI(n_+1), memJ(k_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(k_) { }
        explicit Matrix(int m_, int n_, int k_, Init ini=INIT, const AT &a=AT()) : memEle(k_), memI(n_+1), memJ(k_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(k_) {  init(a); }

        explicit Matrix(int n_, Noinit) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { }
        explicit Matrix(int n_, Init ini=INIT, const AT &a=AT()) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { init(a); }

        explicit Matrix(int m_, int n_, Noinit) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { }
        explicit Matrix(int m_, int n_, Init ini=INIT, const AT &a=AT()) : memEle(n_*n_), memI(n_+1), memJ(n_*n_), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(n_), n(n_), k(n_*n_) { init(a); }


	/*! \brief Copy Constructor
	 *
	 * Constructs a reference to the matrix \em A.
	 * \param A The matrix that will be referenced.
	 * */
	Matrix(const Matrix<Sparse,Ref,Ref,AT> &A) : memEle(A.k), memI(A.n+1), memJ(A.k), ele((AT*)memEle.get()), I((int*)memI.get()), J((int*)memJ.get()), m(A.m), n(A.n), k(A.k) {
          copy(A);
	}

	/*! \brief Destructor. 
	 * */
	~Matrix() = default;

        Matrix<Sparse,Ref,Ref,AT>& resize(int n_, int k_, Noinit) {
            m = n_; n = n_; k = k_;
            memEle.resize(k);
            memI.resize(m+1);
            memJ.resize(k);
            ele = (AT*)memEle.get();
            I = (int*)memI.get();
            J = (int*)memJ.get();
            return *this;
        }

        Matrix<Sparse,Ref,Ref,AT>& resize(int n, int k, Init ini=INIT, const AT &a=AT()) { return resize(m,k,Noinit()).init(a); }

	/*! \brief Assignment operator
	 *
	 * Copies the sparse matrix given by \em A.
	 * \param A The matrix to be assigned. 
	 * \return A reference to the calling matrix.
	 * */
	inline Matrix<Sparse,Ref,Ref,AT>& operator=(const Matrix<Sparse,Ref,Ref,AT> &A) {
          assert(m == A.m);
          assert(k == A.k);
          return copy(A);
        }

	/*! \brief Assignment operator
	 *
	 * Copies the sparse matrix given by \em A.
	 * \param A The matrix to be assigned. 
	 * \return A reference to the calling matrix.
	 * */
        template<class Type, class Row, class Col>
	inline Matrix<Sparse,Ref,Ref,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
          assert(A.rows() == A.cols());
          assert(m == A.cols());
          assert(k == A.countElements());
          return copy(A);
        }

	/*! \brief Reference operator
	 *
	 * References the sparse matrix given by \em A.
	 * \param A The matrix to be referenced. 
	 * \return A reference to the calling matrix.
	 * */
        inline Matrix<Sparse,Ref,Ref,AT>& operator&=(Matrix<Sparse,Ref,Ref,AT> &A) {
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

        /*! \brief Matrix assignment
	 *
	 * Copies the sparse matrix given by \em A.
	 * \param A The matrix to be copied.
	 * \return A reference to the calling matrix.
	 * */
        template<class Type, class Row, class Col>
	inline Matrix<Sparse,Ref,Ref,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
          assert(A.rows() == A.cols());
          int k_ = A.countElements();
          if(m!=A.rows() || k!=k_) resize(A.rows(),k_,NONINIT);
          return copy(A);
        }

	/*! \brief Pointer operator
	 *
	 * See operator()() 
	 * */
	const AT* operator()() const {
	  return ele;
	}

	/*! \brief Pointer operator.
	 *
	 * Returns the pointer to the first element.
	 * \return The pointer to the first element.
	 * */
	AT* operator()() {
	  return ele;
	}

        AT& operator()(int i, int j) {
          throw std::runtime_error("Matrix<Sparse, Ref, Ref, AT>::operator(int i, int j) is not implemented.");
        }

        const AT& operator()(int i, int j) const {
          throw std::runtime_error("Matrix<Sparse, Ref, Ref, AT>::operator(int i, int j) const is not implemented.");
        }

        int ldim() {
          throw std::runtime_error("Matrix<Sparse, Ref, Ref, AT>::ldim() cannot be called.");
        }

	/*! \brief Pointer operator.
	 *
	 * See Ip() 
	 * */
	const int* Ip() const {
	  return I;
	}

	/*! \brief Pointer operator.
	 *
	 * \todo Docu
	 * */
	int* Ip() {
	  return I;
	}

	/*! \brief Pointer operator.
	 *
	 * See Jp() 
	 * */
	const int* Jp() const {
	  return J;
	}

	/*! \brief Pointer operator.
	 *
	 * \todo Docu
	 * */
	int* Jp() {
	  return J;
	}

	/*! \brief Number of rows.
	 *
	 * \return The number of rows of the matrix
	 * */
	int rows() const {return m;}

	/*! \brief Number of columns.
	 *
	 * \return The number of columns of the matrix
	 * */
	int cols() const {return n;}

        /*! \brief Storage convention.
         * */
        const CBLAS_ORDER blasOrder() const {
          throw std::runtime_error("Matrix<Sparse, Ref, Ref, AT>::blasOrder() cannot be called.");
        }

	/*! \brief Initialization.
	 *
	 * Initializes all elements of the calling matrix with 
	 * the value given by \em a.
	 * \param a Value all elements will be initialized with.
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Sparse,Ref,Ref,AT>& init(const AT &val);
        inline Matrix<Sparse,Ref,Ref,AT>& init(Init, const AT &a=AT()) { return init(a); }
        inline Matrix<Sparse,Ref,Ref,AT>& init(Noinit, const AT &a=AT()) { return *this; }

        int countElements() const { return k; }
    };
  // ------------------------- Constructors -------------------------------------
  // ----------------------------------------------------------------------------

  template <class AT> template <class Type, class Row, class Col>
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::copy(const Matrix<Type,Row,Col,AT> &A) {
      int k=0;
      int i;
      for(i=0; i<A.rows(); i++) {
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
	for(int j=i+1; j<A.cols(); j++) {
	  // TODO eps
	  if(fabs(A(i,j))>1e-16) {
	    ele[k]=A(i,j);
	    J[k++]=j;
	  }
	}
      }
      if(n) I[i]=k;
      return *this;
    }

  template <class AT> template <class Row> 
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::copy(const Matrix<Symmetric,Row,Row,AT> &A) {
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
      return *this;
  }

  template <class AT>
    inline Matrix<Sparse,Ref,Ref,AT>& Matrix<Sparse,Ref,Ref,AT>::copy(const Matrix<Sparse,Ref,Ref,AT> &A) {
      for(int i=0; i<=m; i++) {
	I[i] = A.I[i];
      }
      for(int i=0; i<k; i++) {
	ele[i] = A.ele[i];
	J[i] = A.J[i];
      }
      return *this;
    }

  template <class AT>
    Matrix<Sparse,Ref,Ref,AT>&  Matrix<Sparse,Ref,Ref,AT>::init(const AT& val) {
      for(int i=0; i<k; i++) {
	ele[i] = val;
      }
      return *this;
    }

}

#endif
