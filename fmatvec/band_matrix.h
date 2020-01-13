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
  template <class AT> class Matrix<GeneralBand,Ref,Ref,AT> {

   /// @cond NO_SHOW

    protected:

      Memory<AT> memory;
      AT *ele;
      int n{0};
      int kl{0};
      int ku{0};

      inline Matrix<GeneralBand,Ref,Ref,AT>& copy(const Matrix<GeneralBand,Ref,Ref,AT> &A);

      /// @endcond

    public:

      typedef AT value_type;

      /*! \brief Standard constructor
       *
       * Constructs a band matrix with no size. 
       * */
      explicit Matrix() : memory(), ele(0) {
      }

//      template<class Ini=All<AT>>
//        Matrix(int n_, int kl_, int ku_, Ini ini=All<AT>()) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) {  
//          init(ini);
//        }
//      template<class Ini=All<AT>>
//        Matrix(int m_, int n_, Ini ini=All<AT>()) : memory(n_*(n_+n_+1)), ele((AT*)memory.get()), n(n_), kl(n_), ku(n_) {  
//          init(ini);
//        }

      explicit Matrix(int n_, int kl_, int ku_, Noinit) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) { }
      explicit Matrix(int n_, int kl_, int ku_, Init ini=INIT, const AT &a=AT()) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) { init(a); }
      explicit Matrix(int n_, int kl_, int ku_, Eye ini, const AT &a=1) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) { init(ini,a); }
      explicit Matrix(int m_, int n_, int kl_, int ku_, Noinit) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) { }
      explicit Matrix(int m_, int n_, int kl_, int ku_, Init ini=INIT, const AT &a=AT()) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) { init(a); }
      explicit Matrix(int m_, int n_, int kl_, int ku_, Eye ini, const AT &a=1) : memory(n_*(kl_+ku_+1)), ele((AT*)memory.get()), n(n_), kl(kl_), ku(ku_) { init(ini,a); }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the band matrix \em A.
       * \param A The matrix that will be copied.
       * */
      Matrix(const Matrix<GeneralBand,Ref,Ref,AT> &A) : memory(A.n*(A.kl+A.ku+1)), ele((AT*)memory.get()), n(A.n), kl(A.kl), ku(A.ku) {
        copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Type, class Row, class Col>
      explicit Matrix(const Matrix<Type,Row,Col,AT> &A) : memory(A.rows()*(2*A.rows()-1)), ele((AT*)memory.get()), n(A.rows()), kl(n-1), ku(n-1) {
        assert(A.rows() == A.cols());
	copy(A);
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
      explicit Matrix(int n_, int kl_, int ku_, AT* ele_) : memory(), ele(ele_), n(n_), kl(kl_), ku(ku_) { 
      }

      /*! \brief Destructor. 
       * */
      ~Matrix() = default;

      Matrix<GeneralBand,Ref,Ref,AT>& resize(int n_, int kl_, int ku_, Noinit) { 
          n = n_; kl = kl_; ku = ku_;
          memory.resize(n*(kl+ku+1));
          ele = (AT*)memory.get();
          return *this;
      }

      Matrix<GeneralBand,Ref,Ref,AT>& resize(int n, int kl, int ku=0, Init ini=INIT, const AT &a=AT()) { return resize(n,kl,ku,Noinit()).init(a); }

      Matrix<GeneralBand,Ref,Ref,AT>& resize(int n, int kl, int ku, Eye ini, const AT &a=1) { return resize(n,kl,ku,Noinit()).init(a); }

      /*! \brief Assignment operator
       *
       * Copies the band matrix given by \em A.
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<GeneralBand,Ref,Ref,AT>& operator=(const Matrix<GeneralBand,Ref,Ref,AT> &A) {
        assert(n == A.n);
        return copy(A);
      }

      /*! \brief Assignment operator
       *
       * Copies the band matrix given by \em A.
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
      inline Matrix<GeneralBand,Ref,Ref,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
        assert(n == A.rows());
        assert(n == A.cols());
        return copy(A);
      }

      /*! \brief Reference operator
       *
       * References the band matrix given by \em A.
       * \param A The matrix to be referenced. 
       * \return A reference to the calling matrix.
       * */
      inline Matrix<GeneralBand,Ref,Ref,AT>& operator&=(Matrix<GeneralBand,Ref,Ref,AT> &A) {
        n = A.n;
        memory = A.memory;
        ele = A.ele;
        return *this;
      }

      /*! \brief Matrix assignment
       *
       * Copies the band matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
      inline Matrix<GeneralBand,Ref,Ref,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
        if(n!=A.n) resize(A.rows(),NONINIT);
        return copy(A);
      }

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \remark The bounds are checked in debug mode. 
       * \sa operator()(int,int) const
       * */
      const AT& operator()(int i, int j) const {
	assert(i>=0);
	assert(j>=0);
	assert(i<n);
	assert(j<n);
	//     assert(i-j<=kl);
	//      assert(i-j>=-ku);
	static AT zero=0;

	return ((i-j>kl) || (i-j<-ku)) ? zero : ele[ku+i+j*(kl+ku)];
      }

      /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return ele;}

      /*! \brief Pointer operator.
       *
       * see operator()()
       * */
      const AT* operator()() const {return ele;}

      /*! \brief Size.
       *
       * \return The size of the matrix
       * */
      int size() const {return n;}

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix
       * */
      int rows() const {return n;}

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix
       * */
      int cols() const {return n;}

      /*! \brief Storage convention.
       *
       * Returns the blas-conform storage convention. 
       * The elements are stored in columnmajor form,
       * i.e. the elements are stored columnwise. 
       * \return CblasColMajor.
       * */
      const CBLAS_ORDER blasOrder() const {
        return CblasColMajor;
      }

//      /*! \brief Diagonal operator
//       *
//       * See operator()(int)
//       * */
//      inline const Vector<Ref,AT> operator()(int i) const;
//
//      /*! \brief Diagonal operator
//       *
//       * If i<0, the i-th subdiagonal is returned. Otherwise the i-th
//       * superdiagonal is returned.
//       * \param i The i-th super- and subdiagonal,
//       * respectively.
//       * */
//      inline Vector<Ref,AT> operator()(int i);

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<GeneralBand,Ref,Ref,AT>& init(const AT &val=AT());
      inline Matrix<GeneralBand,Ref,Ref,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline Matrix<GeneralBand,Ref,Ref,AT>& init(Eye, const AT &val=1);
      inline Matrix<GeneralBand,Ref,Ref,AT>& init(Noinit, const AT &a=AT()) { return *this; }

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit inline operator std::vector<std::vector<AT>>() const;
  };

  template <class AT>
    inline Matrix<GeneralBand,Ref,Ref,AT>&  Matrix<GeneralBand,Ref,Ref,AT>::init(const AT& val) {
      for(int i=0; i<kl+ku+1; i++)
	for(int j=0; j<n; j++)
	  ele[i+j*(kl+ku+1)] = val;
      return *this;
    }

  template <class AT>
    inline Matrix<GeneralBand,Ref,Ref,AT>&  Matrix<GeneralBand,Ref,Ref,AT>::init(Eye, const AT &val) {
      for(int i=0; i<n; i++)
        ele[i] = val;
      return *this;
    }

//  template <class AT>
//    inline Vector<Ref,AT> Matrix<GeneralBand,Ref,Ref,AT>::operator()(int i) {
//      assert(i<=ku);
//      assert(i>=-kl);
//
//      //return Vector<Ref,AT>(n-i,ku+kl+1,memory,ele[ku-i + (i>0?i:0)*(kl+ku+1)]);
//      return Vector<Ref,AT>(n-abs(i),ku+kl+1,true,memory,ele+ku-i + (i>0?i:0)*(kl+ku+1));
//      //(ku-i,i>0?i:0));
//    }
//
//  template <class AT>
//    inline const Vector<Ref,AT> Matrix<GeneralBand,Ref,Ref,AT>::operator()(int i) const {
//      assert(i<=ku);
//      assert(i>=-kl);
//
//      return Vector<Ref,AT>(n-abs(i),ku+kl+1,true,memory,ele+ku-i + (i>0?i:0)*(kl+ku+1));
//      //(ku-i,i>0?i:0));
//    }

  template <class AT>
    inline Matrix<GeneralBand,Ref,Ref,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=operator()(r,c);
      }
      return ret;
    }

  /// @cond NO_SHOW

  template <class AT>
    inline Matrix<GeneralBand,Ref,Ref,AT>& Matrix<GeneralBand,Ref,Ref,AT>::copy(const Matrix<GeneralBand,Ref,Ref,AT> &A) {
      for(int i=0; i<kl+ku+1; i++)
	for(int j=0; j<n; j++)
	  ele[i+j*(kl+ku+1)] = A.ele[i+j*(kl+ku+1)];
      return *this;
    }

  /// @endcond

}

#endif
