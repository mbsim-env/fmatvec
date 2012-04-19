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

#ifndef general_band_matrix_h
#define general_band_matrix_h

#include "vector.h"
#include "types.h"
#include "_memory.h"
#include <cstdlib>

namespace fmatvec {

  /*! 
   *  \brief This is a matrix class for general band matrices.
   *  
   * Template class Matrix with shape type GeneralBand and atomic type AT. 
   * The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <class AT> class Matrix<GeneralBand, AT> {

   /// @cond NO_SHOW

    protected:

      Memory<AT> memory;
      AT *ele;
      int n;
      int kl;
      int ku;

      void deepCopy(const Matrix<GeneralBand, AT> &A); 

      Matrix(int n_, int kl_, int ku_, Memory<AT> memory_, const AT* ele_) : memory(memory_), ele((AT*)ele_), n(n_), kl(kl_), ku(ku_) {
      }

      /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a band matrix with no size. 
       * */
      Matrix() : memory(), ele(0), n(0), kl(0), ku(0) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a band matrix of size n x n.
       * \param n_ The number of rows and columns.
       * \param kl_ The number of subdiagonals.
       * \param ku_ The number of superdiagonals.
       * \remark The matrix will be initialised to
       * zero by default. This default behavior can be changed by defining 
       * FMATVEC_NO_INITIALIZATION.
       * */
      Matrix(int n_, int kl_, int ku_) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) {  

#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      /*! \brief Regular Constructor
       *
       * Constructs a band matrix of size n x n. 
       * The matrix will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * matrix will not be initialized.
       * \param n_ The number of rows and columns.
       * \param kl_ The number of subdiagonals.
       * \param ku_ The number of superdiagonals.
       * \remark The matrix will be initialised to
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * */
      Matrix(int n_, int kl_, int ku_, Initialization ini, const AT &a=0) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) {  

	if(ini == INIT)
	  init(a);
	else if(ini == EYE){
	  init(0);
	  for(int i=0; i<n; i++)
	    ele[ku+i+i*(kl+ku)] = 1;
	} 

      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the band matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      Matrix(const Matrix<GeneralBand, AT> &A) : memory(A.memory), ele(A.ele) , n(A.n), kl(A.kl), ku(A.ku) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a band matrix of size n x n with the pyhsical memory 
       * given by \em ele_.
       * \param n_ The number of rows and columns.
       * \param kl_ The number of subdiagonals.
       * \param ku_ The number of superdiagonals.
       * \param ele_ The physical memory the matrix will point to.
       * */
      Matrix(int n_, int kl_, int ku_, AT* ele_) : memory(), ele(ele_), n(n_), kl(kl_), ku(ku_) { 
      }

      /*! \brief Destructor. 
       * */
      ~Matrix() {
      }

      /*! \brief Matrix resizing. 
       *
       * Resizes the matrix to size n x n.    
       * \param n_ The number of rows and columns.
       * \param kl_ The number of subdiagonals.
       * \param ku_ The number of superdiagonals.
       * \return A reference to the calling matrix.
       * \remark The matrix will be initialised to
       * zero by default. To change this behavior, define
       * FMATVEC_NO_INITIALIZATION.
       * */
      Matrix<GeneralBand, AT>& resize(int n_, int kl_, int ku_) {
	n=n_;
	kl=kl_;
	ku=ku_;
	memory.resize(n*(kl+ku+1));
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
       * \param n_ The number of rows and columns.
       * \param kl_ The number of subdiagonals.
       * \param ku_ The number of superdiagonals.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the matrix will be initialized with (default 0)
       * \return A reference to the calling matrix.
       * */
      Matrix<GeneralBand, AT>& resize(int n_, int kl_, int ku_, Initialization ini, const AT &a=0) {
	n=n_;
	kl=kl_;
	ku=ku_;
	memory.resize(n*(kl+ku+1));
	ele = (AT*)memory.get();

	if(ini == INIT)
	  init(a);
	else if(ini == EYE){
	  init(0);
	  for(int i=0; i<n; i++)
	    ele[ku+i+i*(kl+ku)] = 1;
	}

	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the band matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      Matrix<GeneralBand, AT>& operator=(const Matrix<GeneralBand, AT> &A) {
#ifndef FMATVEC_NO_DEEP_ASSIGNMENT 
	return operator<<(A);
#else
	return operator>>(A);
#endif
      }

      /*! \brief Copy operator
       *
       * Copies the band matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      Matrix<GeneralBand, AT>& operator<<(const Matrix<GeneralBand, AT> &A);

      /*! \brief Reference operator
       *
       * References the band matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      Matrix<GeneralBand, AT>& operator>>(const Matrix<GeneralBand, AT> &A);

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
	//     assert(i-j<=kl);
	//      assert(i-j>=-ku);
#endif
	static AT zero=0;

	return ((i-j>kl) || (i-j<-ku)) ? zero : ele[ku+i+j*(kl+ku)];
      };

      /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return ele;};

      /*! \brief Pointer operator.
       *
       * see operator()()
       * */
      const AT* operator()() const {return ele;};

      /*! \brief Size.
       *
       * \return The size of the matrix
       * */
      int size() const {return n;};

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

      /*! \brief Diagonal operator
       *
       * See operator()(int) 
       * */
      const Vector<AT> operator()(int i) const;

      /*! \brief Diagonal operator
       *
       * If i<0, the i-th subdiagonal is returned. Otherwise the i-th
       * superdiagonal is returned.  
       * \param i The i-th super- and subdiagonal,
       * respectively.
       * */
      Vector<AT> operator()(int i);

      /*! \brief Matrix unsizing.
       *
       * Resizes the matrix to size zero.  
       * \return A reference to the calling matrix.
       * */
      Matrix<GeneralBand, AT>& resize() {n=0;return *this;};

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      Matrix<GeneralBand, AT> copy() const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      Matrix<GeneralBand, AT>& init(const AT &a);

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      operator std::vector<std::vector<AT> >();
  };

  template <class AT>
    Matrix<GeneralBand, AT>& Matrix<GeneralBand, AT>::operator>>(const Matrix<GeneralBand, AT> &A) { 

      if(n==0) {
	n=A.n;
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(n == A.n);
#endif
      }

      memory = A.memory;
      ele = A.ele;

      return *this;
    }

  template <class AT>
    Matrix<GeneralBand, AT>& Matrix<GeneralBand, AT>::operator<<(const Matrix<GeneralBand, AT> &A) { 

      if(A.rows() == 0 || A.cols() == 0)
	return *this;

      if(n==0) {
	n = A.n;
	memory.resize(n*n);
	ele = (AT*)memory.get();
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(n == A.n);
#endif
      }

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<GeneralBand, AT>&  Matrix<GeneralBand, AT>::init(const AT& val) {

      for(int i=0; i<kl+ku+1; i++)
	for(int j=0; j<n; j++)
	  ele[i+j*(kl+ku+1)] = val;

      return *this;
    }

  template <class AT>
    Vector<AT> Matrix<GeneralBand, AT>::operator()(int i) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i<=ku);
      assert(i>=-kl);
#endif

      //return Vector<AT>(n-i,ku+kl+1,memory,ele[ku-i + (i>0?i:0)*(kl+ku+1)]);
      return Vector<AT>(n-abs(i),ku+kl+1,true,memory,ele+ku-i + (i>0?i:0)*(kl+ku+1));
      //(ku-i,i>0?i:0));
    }

  template <class AT>
    const Vector<AT> Matrix<GeneralBand, AT>::operator()(int i) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(i<=ku);
      assert(i>=-kl);
#endif

      return Vector<AT>(n-abs(i),ku+kl+1,true,memory,ele+ku-i + (i>0?i:0)*(kl+ku+1));
      //(ku-i,i>0?i:0));
    }

  template <class AT>
    Matrix<GeneralBand, AT> Matrix<GeneralBand, AT>::copy() const {

      Matrix<GeneralBand, AT> A(n,kl,ku);
      A.deepCopy(*this);

      return A;
    }

  template <class AT>
    void Matrix<GeneralBand, AT>::deepCopy(const Matrix<GeneralBand, AT> &A) { 

      for(int i=0; i<kl+ku+1; i++)
	for(int j=0; j<n; j++)
	  ele[i+j*(kl+ku+1)] = A.ele[i+j*(kl+ku+1)];
    }

  template <class AT>
    Matrix<GeneralBand, AT>::operator std::vector<std::vector<AT> >() {
      std::vector<std::vector<AT> > ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=operator()(r,c);
      }
      return ret;
    }

  //extern void Matrix<GeneralBand, double>::deepCopy(const Matrix<GeneralBand, double> &A);
}

#endif
