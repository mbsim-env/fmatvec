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

      inline Matrix<Diagonal,Ref,Ref,AT>& copy(const Matrix<Diagonal,Ref,Ref,AT> &A);

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

      //        template<class Ini=All<AT>>
      //          Matrix(int n_, Ini ini=All<AT>()) : memory(n_), ele((AT*)memory.get()), n(n_) { init(ini); }
      //        template<class Ini=All<AT>>
      //          Matrix(int m_, int n_, Ini ini=All<AT>()) : memory(n_), ele((AT*)memory.get()), n(n_) { init(ini); }

      explicit Matrix(int n_, Noinit) : memory(n_), ele((AT*)memory.get()), n(n_) { }
      explicit Matrix(int n_, Init ini=INIT, const AT &a=AT()) : memory(n_), ele((AT*)memory.get()), n(n_) { init(a); }
      explicit Matrix(int n_, Eye ini, const AT &a=1) : memory(n_), ele((AT*)memory.get()), n(n_) { init(ini,a); }
      explicit Matrix(int m_, int n_, Noinit) : memory(n_), ele((AT*)memory.get()), n(n_) { }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the diagonal matrix \em A.
       * \param A The matrix that will be copied.
       * */
      Matrix(const Matrix<Diagonal,Ref,Ref,AT> &A) : memory(A.n), ele((AT*)memory.get()) ,n(A.n) {
        copy(A);
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the diagonal matrix \em A.
       * \param A The matrix that will be copied.
       * */
      template<class Type, class Row, class Col>
        Matrix(const Matrix<Type,Row,Col,AT> &A) : memory(A.rows()), ele((AT*)memory.get()) ,n(A.rows()) {
          assert(A.rows() == A.cols());
          copy(A);
        }

      /*! \brief Destructor.
       * */
      ~Matrix() = default;

      Matrix<Diagonal,Ref,Ref,AT>& resize(int n_, Noinit) {
        n = n_;
        memory.resize(n);
        ele = (AT*)memory.get();
        return *this;
      }

      Matrix<Diagonal,Ref,Ref,AT>& resize(int n, Init ini=INIT, const AT &a=AT()) { resize(n,Noinit()).init(a); }
      Matrix<Diagonal,Ref,Ref,AT>& resize(int n, Eye ini, const AT &a=1) { resize(n,Noinit()).init(ini,a); }

      /*! \brief Assignment operator
       *
       * Copies the diagonal matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Diagonal,Ref,Ref,AT>& operator=(const Matrix<Diagonal,Ref,Ref,AT> &A) {
        assert(n == A.n);
        return copy(A);
      }

      /*! \brief Assignment operator
       *
       * Copies the diagonal matrix given by \em A.
       * \param A The matrix to be assigned.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
        inline Matrix<Diagonal,Ref,Ref,AT>& operator=(const Matrix<Type,Row,Col,AT> &A) {
          assert(n == A.rows());
          assert(n == A.cols());
          return copy(A);
        }

      /*! \brief Reference operator
       *
       * References the diagoanl matrix given by \em A.
       * \param A The matrix to be referenced.
       * \return A reference to the calling matrix.
       * */
      inline Matrix<Diagonal,Ref,Ref,AT>& operator&=(const Matrix<Diagonal,Ref,Ref,AT> &A) {
        n=A.n;
        memory = A.memory;
        ele = A.ele;
        return *this;
      }

      /*! \brief Matrix assignment
       *
       * Copies the diagonal matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template<class Type, class Row, class Col>
        inline Matrix<Diagonal,Ref,Ref,AT>& operator<<=(const Matrix<Type,Row,Col,AT> &A) {
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

        return e(i,j);
      };

      /*! \brief Element operator
       *
       * See operator()(int)
       * */
      const AT& operator()(int i) const {

        assert(i>=0);
        assert(i<n);

        return e(i);
      };

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the i-th column.
       * \param i The i-th row and column of the matrix
       * \return A reference to the element A(i,i).
       * \remark The bounds are checked in debug mode.
       * \sa operator()(int) const
       * */
      AT& operator()(int i) {

        assert(i>=0);
        assert(i<n);

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

      /*! \brief std::vector<std::vector<AT>> Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT>> object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT>> the matrix will be initialized with.
       * */
      explicit Matrix<Diagonal,Ref,Ref,AT>(const std::vector<std::vector<AT>> &m);

      /*! \brief Cast to std::vector<std::vector<AT>>.
       *
       * \return The std::vector<std::vector<AT>> representation of the matrix
       * */
      explicit operator std::vector<std::vector<AT>>() const;

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the matrix
//       * */
//      explicit operator AT() const {
//        assert(n==1);
//        return ele[0];
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a matrix with a AT object.
//       * \param x The AT the matrix will be initialized with.
//       * */
//      explicit Matrix(const AT &x) : memory(1), ele((AT*)memory.get()), n(1) { ele[0] = x; }
  };

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>::Matrix(const std::vector<std::vector<AT>> &m) : Matrix<Diagonal,Ref,Ref,AT>(m.size()) {
      for(size_t r=0; r<m.size(); ++r) {
        if(n!=static_cast<int>(m[r].size()))
          throw std::runtime_error("The rows of the input have different lenght.");
        for(size_t c=0; c<m[r].size(); ++c) {
          if(r==c)
            e(r)=m[r][r];
          if(r!=c && (fabs(m[r][c])>1e-13 || fabs(m[r][c])>1e-13))
            throw std::runtime_error("The input is not diagonal.");
        }
      }
    }


  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>&  Matrix<Diagonal,Ref,Ref,AT>::init(const AT &val) {
      for(int i=0; i<rows(); i++) 
        e(i) = val;
      return *this;
    }

  template <class AT>
    inline Matrix<Diagonal,Ref,Ref,AT>& Matrix<Diagonal,Ref,Ref,AT>::copy(const Matrix<Diagonal,Ref,Ref,AT> &A) {
      for(int i=0; i<n; i++)
        e(i) = A.e(i);
      return *this;
    }

  template <class AT>
    Matrix<Diagonal,Ref,Ref,AT>::operator std::vector<std::vector<AT>>() const {
      std::vector<std::vector<AT>> ret(rows());
      for(int r=0; r<rows(); r++) {
        ret[r].resize(cols());
        for(int c=0; c<cols(); c++)
          ret[r][c]=e(r,c);
      }
      return ret;
    }

}

#endif
