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

#ifndef fixed_vector_h
#define fixed_vector_h

#include "fixed_general_matrix.h"
#include <vector>
#include <cstring>
#include <stdexcept>

namespace fmatvec {

  /*! 
   *  \brief This is a vector class of general shape in dense storage format.
   *
   * Template class Vector of type General, id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the vector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <int M, class AT>
  class Vector<Fixed<M>,AT> : public Matrix<General,Fixed<M>,Fixed<1>,AT> {
    using Matrix<General,Fixed<M>,Fixed<1>,AT>::ele;

    public:

    typedef AT* iterator;
    typedef const AT* const_iterator;

    typedef AT value_type;

    /// @cond NO_SHOW

    protected:

    template<class Row> inline void deepCopy(const Vector<Row,AT> &x);

    /// @endcond
    
    public:

      explicit Vector(Noinit ini) : Matrix<General,Fixed<M>,Fixed<1>,AT>(ini) { }
      explicit Vector(Init ini=INIT, const AT &a=AT()) : Matrix<General,Fixed<M>,Fixed<1>,AT>(ini,a) { }
      explicit Vector(int m, Noinit ini) : Matrix<General,Fixed<M>,Fixed<1>,AT>(ini) { assert(m==M); }
      explicit Vector(int m, Init ini=INIT, const AT &a=AT()) : Matrix<General,Fixed<M>,Fixed<1>,AT>(ini,a) { assert(m==M); }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      Vector(const Vector<Fixed<M>,AT> &x) : Matrix<General,Fixed<M>,Fixed<1>,AT>(x) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General, AT>&) 
       * */
      template<class Row>
      Vector(const Vector<Row,AT> &A) : Matrix<General,Fixed<M>,Fixed<1>,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General, AT>&) 
       * */
      template<class Type, class Row, class Col>
      explicit Vector(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Fixed<M>,Fixed<1>,AT>(A) {
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a vector with a string in a matlab-like
       * notation. The entries are seperated by semicolons.
       * For example
       * \code 
       * Vector<double> x("[3;1;2]");
       * \endcode
       * constructs the vector
       * \f[ x=\begin{pmatrix}3\\ 1 \\ 2\end{pmatrix}  \f]
       * \param str The string the vector will be initialized with. 
       * */
      Vector(const char *str) : Matrix<General,Fixed<M>,Fixed<1>,AT>(str) {
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be assigned.
       * \return A reference to the calling vector.
       * */
      inline Vector<Fixed<M>,AT>& operator=(const Vector<Fixed<M>,AT> &x) {
        deepCopy(x);
        return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be assigned.
       * \return A reference to the calling vector.
       * */
      template <class Row>
      inline Vector<Fixed<M>,AT>& operator=(const Vector<Row,AT> &x) {
        assert(x.size() == M);
        deepCopy(x);
        return *this;
      }

      /*! \brief Vector replacement
       *
       * Copies the vector given by \em x.
       * \param x The vector to be copied. 
       * \return A reference to the calling vector.
       * */
      template <class Row>
      inline Vector<Fixed<M>,AT>& replace(const Vector<Row,AT> &x) { return operator=(x); }

      template <class AT2>
      operator Vector<Fixed<M>,AT2>() const {
        Vector<Fixed<M>,AT2> ret;
        for(size_t i=0; i<M; ++i)
          ret(i) = (*this)(i);
        return ret;
      }

      /*! \brief Element operator
       *
       * Returns a reference to the i-th element. 
       * \param i The i-th element.
       * \return A reference to the element x(i).
       * \remark The bounds are checked in debug mode.
       * \sa operator()(int) const
       * */
      AT& operator()(int i) {

	assert(i>=0);
	assert(i<M);

	return e(i);
      };

      /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {

	assert(i>=0);
	assert(i<M);

	return e(i);
      };

      iterator begin() { return &ele[0][0]; }
      iterator end() { return &ele[M][0]; }
      const_iterator begin() const { return &ele[0][0]; }
      const_iterator end() const { return &ele[M][0]; }
      const_iterator cbegin() const noexcept { return &ele[0][0]; }
      const_iterator cend() const noexcept { return &ele[M][0]; }

      AT& e(int i) {
	return ele[i][0];
      };

      /*! \brief Element operator
       *
       * See e(int) 
       * */
      const AT& e(int i) const {
	return ele[i][0];
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling vector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling vector.
       * */
      inline Vector<Fixed<M>,AT>& init(const AT& val=AT());
      inline Vector<Fixed<M>,AT>& init(Init, const AT& a=AT()) { return init(a); }
      inline Vector<Fixed<M>,AT>& init(Noinit, const AT& a=AT()) { return *this; }

      /*! \brief Size.
       *
       * \return The size of the vector.
       * */
      int size() const {return M;};

      //! Resize a fixed vector
      //! Do nothing for the fixed dimension and throw for any other dimension.
      void resize(int m) {
        if(m!=M)
          throw std::runtime_error("A fixed vector cannot be resized.");
      }

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return 1;};

      template <int M1, int M2>
      inline const Vector<Fixed<M2-M1+1>,AT> operator()(const Range<Fixed<M1>,Fixed<M2> > &I) const;

      using Matrix<General,Fixed<M>,Fixed<1>,AT>::operator();

      /*! \brief Cast to std::vector<AT>.
       *
       * \return The std::vector<AT> representation of the vector
       * */
      inline operator std::vector<AT>() const;

      /*! \brief std::vector<AT> Constructor.
       * Constructs and initializes a vector with a std::vector<AT> object.
       * \param v The std::vector<AT> the vector will be initialized with. 
       * */
      inline Vector(std::vector<AT> v);

      /*!
       * \brief return the transpose of the Vector, i.e. a RowVector
       */
      inline const RowVector<Fixed<M>,AT> T() const;

      /*!
       * \brief set a subvector - specified by the range - to the given vector
       */
      template<class Type, class Row, class Col>
      inline void set(const Range<Var,Var> &I, const Matrix<Type,Row,Col,AT> &A);

  };

  template <int M, class AT>
    inline Vector<Fixed<M>,AT>& Vector<Fixed<M>,AT>::init(const AT &val) {
      for(int i=0; i<M; i++) 
	e(i) = val; 
      return *this;
    }

  template <int M, class AT> template <int M1, int M2>
    inline const Vector<Fixed<M2-M1+1>,AT> Vector<Fixed<M>,AT>::operator()(const Range<Fixed<M1>,Fixed<M2> > &I) const {
      assert(M2<M);
      Vector<Fixed<M2-M1+1>,AT> x(NONINIT);

      for(int i=0; i<x.size(); i++) 
        x.e(i) = e(M1+i);

      return x;
    }

  template <int M, class AT>
    inline const RowVector<Fixed<M>,AT> Vector<Fixed<M>,AT>::T() const {
      RowVector<Fixed<M>,AT> x(NONINIT);
      for(int i=0; i<M; i++)
	x.e(i) = e(i);
      return x;
    }

  template <int M, class AT>
   inline Vector<Fixed<M>,AT>::operator std::vector<AT>() const {
      std::vector<AT> ret(size());
      if(size()>0) memcpy(&ret[0], &operator()(0), sizeof(AT)*size());
      return ret;
    }

  template <int M, class AT>
    inline Vector<Fixed<M>,AT>::Vector(std::vector<AT> v) : Matrix<General,Fixed<M>,Fixed<1>,AT>() {
      if(size()>0) memcpy(&operator()(0), &v[0], sizeof(AT)*size());
    }

  /// @cond NO_SHOW

  template <int M, class AT> template<class Row>
    inline void Vector<Fixed<M>,AT>::deepCopy(const Vector<Row,AT> &x) {
      for(int i=0; i<M; i++)
	e(i) = x.e(i);
    }

  /// @endcond

  template<int M, class AT> template<class Type, class Row, class Col>
  inline void Vector<Fixed<M>,AT>::set(const Range<Var,Var> &I, const Matrix<Type,Row,Col,AT> &A) {
      Matrix<General,Fixed<M>,Fixed<1>,AT>::set(I, Range<Var,Var>(0,0), A);
  }

}

#endif
