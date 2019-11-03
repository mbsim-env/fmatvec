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

#ifndef diagonal_matrix_h
#define diagonal_matrix_h

#include "matrix.h"
#include "types.h"
#include "_memory.h"

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for diagonal matrices.
   *
   * Template class Matrix of shape type Diagonal.
   * The template parameter AT defines the
   * atomic type of the matrix. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class Matrix<Diagonal,Ref,Ref,AT> {

      protected:

    /// @cond NO_SHOW

	Memory<AT> memory;
	AT *ele;
	int n{0};

	void deepCopy(const Matrix<Diagonal,Ref,Ref,AT> &A);

	explicit Matrix(int n_, Memory<AT> memory_, const AT* ele_) : memory(memory_), ele((AT*)ele_), n(n_) {
	}

	const AT* elePtr(int i) const {
	  return ele+i;
	}

	AT* elePtr(int i) {
	  return ele+i;
	}

    /// @endcond

      public:

        typedef AT value_type;

	/*! \brief Standard constructor
	 *
	 * Construckts a matrix with no size.
	 * */
	explicit Matrix() : memory(), ele(0) { }

//        template<class Ini=All<AT> >
//          Matrix(int n_, Ini ini=All<AT>()) : memory(n_), ele((AT*)memory.get()), n(n_) { init(ini); }
//        template<class Ini=All<AT> >
//          Matrix(int m_, int n_, Ini ini=All<AT>()) : memory(n_), ele((AT*)memory.get()), n(n_) { init(ini); }

        explicit Matrix(int n_, Noinit) : memory(n_), ele((AT*)memory.get()), n(n_) { }
        explicit Matrix(int n_, Init ini=INIT, const AT &a=AT()) : memory(n_), ele((AT*)memory.get()), n(n_) { init(a); }
        explicit Matrix(int n_, Eye ini, const AT &a=1) : memory(n_), ele((AT*)memory.get()), n(n_) { init(ini,a); }
        explicit Matrix(int m_, int n_, Noinit) : memory(n_), ele((AT*)memory.get()), n(n_) { }

	/*! \brief Copy Constructor
	 *
	 * Constructs a reference to the diagonal matrix \em A.
	 * \attention The physical memory of the matrix \em A 
	 * will not be copied, only referenced.
	 * \param A The matrix that will be referenced.
	 * */
	Matrix(const Matrix<Diagonal,Ref,Ref,AT> &A) : memory(A.memory), ele(A.ele) ,n(A.n) {
	}

	/*! \brief Destructor. 
	 * */
	~Matrix() = default;

        Matrix<Diagonal,Ref,Ref,AT>& resize() { 
          n = 0;
          memory.resize(0);
          ele = 0;
          return *this;
        }

        Matrix<Diagonal,Ref,Ref,AT>& resize(int n_, Noinit) { 
          n = n_;
          memory.resize(n);
          ele = (AT*)memory.get();
          return *this;
        }
 
        Matrix<Diagonal,Ref,Ref,AT>& resize(int n, Init ini=INIT, const AT &a=AT()) { resize(n,Noinit()).init(a); }
        Matrix<Diagonal,Ref,Ref,AT>& resize(int n, Eye ini, const AT &a=1) { resize(n,Noinit()).init(ini,a); }
 
	/*! \brief Copy operator
	 *
	 * Copies the diagonal matrix given by \em A.
	 * \param A The matrix to be copied. 
	 * \return A reference to the calling matrix.
	 * */
	inline Matrix<Diagonal,Ref,Ref,AT>& operator<<(const Matrix<Diagonal,Ref,Ref,AT> &A);

	/*! \brief Reference operator
	 *
	 * References the diagoanl matrix given by \em A.
	 * \param A The matrix to be referenced. 
	 * \return A reference to the calling matrix.
	 * */
	inline Matrix<Diagonal,Ref,Ref,AT>& operator>>(const Matrix<Diagonal,Ref,Ref,AT> &A);

	/*! \brief Assignment operator
	 *
	 * Copies the square matrix given by \em A by calling operator<<().
	 * \param A The matrix to be assigned. 
	 * \return A reference to the calling matrix.
	 * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
	 * \sa operator<<(), operator>>()
	 * */
	inline Matrix<Diagonal,Ref,Ref,AT>& operator=(const Matrix<Diagonal,Ref,Ref,AT> &A);

	Matrix<Diagonal,Ref,Ref,AT>(const std::vector<std::vector<AT>> &A);

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
	const AT& operator()(int i, int j) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
	  assert(i>=0);
	  assert(j>=0);
	  assert(i<n);
	  assert(j<n);
#endif

	  return e(i,j);
	};

	/*! \brief Element operator
	 *
	 * See operator()(int) 
	 * */
	const AT& operator()(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
	  assert(i>=0);
	  assert(i<n);
#endif

	  return e(i);
	};

	/*! \brief Element operator
	 *
	 * Returns a reference to the element in the i-th row and the i-th column. 
	 * \param i The i-th row and column of the matrix
	 * \return A reference to the element A(i,i).
	 * \remark The bounds are checked by default. 
	 * To change this behavior, define
	 * FMATVEC_NO_BOUNDS_CHECK.
	 * \sa operator()(int) const
	 * */
	AT& operator()(int i) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
	  assert(i>=0);
	  assert(i<n);
#endif

	  return e(i);
	};

	const AT& e(int i, int j) const {
	  static AT zero=0;
	  return i==j ? ele[i] : zero;
	};

	const AT& e(int i) const {
	  return ele[i];
	};

	AT& e(int i) {
	  return ele[i];
	};

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

	/*! \brief Number of rows.
	 *
	 * \return The number of rows of the matrix
	 * */
	int rows() const {return n;};

	/*! \brief Number of columns.
	 *
	 * \return The number of columns of the matrix
	 * */
	int cols() const {return n;};

	/*! \brief Number of rows and columns.
	 *
	 * \return The number of rows and columns of the matrix
	 * */
	int size() const {return n;};

	/*! \brief Storage convention.
	 *
	 * Returns the blas-conform storage convention. 
	 * The elements are stored * in columnmajor form, 
	 * i.e. the elements are stored columnwise. 
	 * \return CblasColMajor.
	 * */
  const CBLAS_ORDER blasOrder() const {
    return  CblasColMajor;
	};

	/*! \brief Matrix duplicating.
	 *
	 * The calling matrix returns a \em deep copy of itself.  
	 * \return The duplicate.
	 * */
	Matrix<Diagonal,Ref,Ref,AT> copy() const;

	/*! \brief Initialization.
	 *
	 * Initializes all elements of the calling matrix with 
	 * the value given by \em a.
	 * \param a Value all elements will be initialized with.
	 * \return A reference to the calling matrix.
	 * */
	Matrix<Diagonal,Ref,Ref,AT>& init(const AT &val=AT());
        inline Matrix<Diagonal,Ref,Ref,AT>& init(Init, const AT &a=AT()) { return init(a); }
        inline Matrix<Diagonal,Ref,Ref,AT>& init(Eye, const AT &a=1) { return init(a); }
        inline Matrix<Diagonal,Ref,Ref,AT>& init(Noinit, const AT &a=AT()) { return *this; }

        /*! \brief Cast to std::vector<std::vector<AT> >.
         *
         * \return The std::vector<std::vector<AT> > representation of the matrix
         * */
        operator std::vector<std::vector<AT> >() const;
    };

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>& Matrix<Diagonal,Ref,Ref,AT>::operator>>(const Matrix<Diagonal,Ref,Ref,AT> &A) { 
      n=A.n;
      memory = A.memory;
      ele = A.ele;

      return *this;
    }

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>& Matrix<Diagonal,Ref,Ref,AT>::operator=(const Matrix<Diagonal,Ref,Ref,AT> &A) { 
      if(!ele) {
	n = A.n;
	memory.resize(n);
	ele = (AT*)memory.get();
      } 
      else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(n == A.n);
#endif
      }

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>::Matrix(const std::vector<std::vector<AT>> &A) : Matrix<Diagonal,Ref,Ref,AT>(A.size()) { 
      for(size_t r=0; r<A.size(); ++r) {
        if(n!=static_cast<int>(A[r].size()))
          throw std::runtime_error("The rows of the input have different lenght.");
        for(size_t c=0; c<A[r].size(); ++c) {
          if(r==c)
            e(r)=A[r][r];
          if(r!=c && (fabs(A[r][c])>1e-13 || fabs(A[r][c])>1e-13))
            throw std::runtime_error("The input is not diagonal.");
        }
      }
    }


  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>& Matrix<Diagonal,Ref,Ref,AT>::operator<<(const Matrix<Diagonal,Ref,Ref,AT> &A) { 

      if(n!=A.n) {
	n = A.n;
	memory.resize(n);
	ele = (AT*)memory.get();
      } 

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>&  Matrix<Diagonal,Ref,Ref,AT>::init(const AT &val) {
      for(int i=0; i<rows(); i++) 
        e(i) = val;
      return *this;
    }

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT> Matrix<Diagonal,Ref,Ref,AT>::copy() const {

      Matrix<Diagonal,Ref,Ref,AT> A(n);

      A.deepCopy(*this);

      return A;
    }

  template <class AT>
    void Matrix<Diagonal,Ref,Ref,AT>::deepCopy(const Matrix<Diagonal,Ref,Ref,AT> &A) { 
      for(int i=0; i<n; i++)
        e(i) = A.e(i);
    }

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>::operator std::vector<std::vector<AT> >() const {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

}

#endif
