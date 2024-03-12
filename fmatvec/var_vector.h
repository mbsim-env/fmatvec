/* Copyright (C) 2003-2012 Martin Förg

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

#ifndef var_vector_h
#define var_vector_h

#include "var_general_matrix.h"
#include <vector>
#include <cstring>

namespace fmatvec {

  /*! 
   *  \brief This is a vector class of general shape in dense storage format.
   *
   * Template class Vector of type General, id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the vector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class Vector<Var,AT> : public Matrix<General,Var,Fixed<1>,AT> {
    using Matrix<General,Var,Fixed<1>,AT>::M;
    using Matrix<General,Var,Fixed<1>,AT>::ele;
    friend class RowVector<Var,AT>;

    public:
    static constexpr bool isVector {true};

    using iterator = AT *;
    using const_iterator = const AT *;

    using value_type = AT;

    /// @cond NO_SHOW

    protected:

    template<class Row> inline Vector<Var,AT>& copy(const Vector<Row,AT> &x);

    /// @endcond
    
    public:

      /*! \brief Standard constructor
       *
       * Constructs a vector with no size. 
       * */
      explicit Vector() : Matrix<General,Var,Fixed<1>,AT>() { }

      // move
      Vector(Vector<Var,AT> &&src)  noexcept : Matrix<General,Var,Fixed<1>,AT>(std::move(src)) {}
      Vector<Var,AT>& operator=(Vector<Var,AT> &&src)  noexcept {
        FMATVEC_ASSERT(M == src.size(), AT);
        src.M=0;
        delete[]ele;
        ele=src.ele;
        src.ele=nullptr;
        return *this;
      }
      Vector(Transpose, RowVector<Var, AT> &&src) : Matrix<General,Var,Fixed<1>,AT>() {
        M=src.N;
        src.N=0;
        ele=src.ele;
        src.ele=nullptr;
      }

//      template<class Ini=All<AT>>
//        Vector(int m, Ini ini=All<AT>()) : Matrix<General,Var,Fixed<1>,AT>(m,ini) { } 

      explicit Vector(int m, Noinit ini) : Matrix<General,Var,Fixed<1>,AT>(m,ini) { } 
      explicit Vector(int m, Init ini=INIT, const AT &a=AT()) : Matrix<General,Var,Fixed<1>,AT>(m,ini,a) { } 

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the vector \em x.
       * \param x The vector that will be copied.
       * */
      Vector(const Vector<Var,AT> &x) : Matrix<General,Var,Fixed<1>,AT>(x) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the vector \em x.
       * \param x The vector that will be copied.
       * */
      template<class Row>
      Vector(const Vector<Row,AT> &A) : Matrix<General,Var,Fixed<1>,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the vector \em x.
       * \param x The vector that will be copied.
       * */
      template<class Type, class Row, class Col>
      explicit Vector(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Var,Fixed<1>,AT>(A) {
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
      Vector(const char *str) : Matrix<General,Var,Fixed<1>,AT>(str) {
      }

      Vector<Var,AT>& resize(int m, Noinit) {
        Matrix<General,Var,Fixed<1>,AT>::resize(m,Noinit());
        return *this;
      }

      Vector<Var,AT>& resize(int m, Init ini=INIT, const AT &a=AT()) {
        Matrix<General,Var,Fixed<1>,AT>::resize(m,ini,a);
        return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be assigned.
       * \return A reference to the calling vector.
       * */
      inline Vector<Var,AT>& operator=(const Vector<Var,AT> &x) {
        FMATVEC_ASSERT(M == x.size(), AT);
        return copy(x);
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be assigned.
       * \return A reference to the calling vector.
       * */
      template <class Row>
      inline Vector<Var,AT>& operator=(const Vector<Row,AT> &x) {
        FMATVEC_ASSERT(M == x.size(), AT);
        return copy(x);
      }

      /*! \brief Vector assignment
       *
       * Copies the vector given by \em x.
       * \param x The vector to be copied.
       * \return A reference to the calling vector.
       * */
      template <class Row>
      inline Vector<Var,AT>& operator<<=(const Vector<Row,AT> &x) {
        if(M!=x.size()) resize(x.size(),NONINIT);
        return copy(x);
      }
      // move
      inline Vector<Var,AT>& operator<<=(Vector<Var,AT> &&src) {
        M=src.M;
        src.M=0;
        delete[]ele;
        ele=src.ele;
        src.ele=nullptr;
        return *this;
      }

      template <class AT2>
      operator Vector<Var,AT2>() const {
        Vector<Var,AT2> ret(size());
        for(size_t i=0; i<size(); ++i)
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

	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(i<M, AT);

	return e(i);
      }

      /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {

	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(i<M, AT);

	return e(i);
      }

      iterator begin() { return &ele[0]; }
      iterator end() { return &ele[M]; }
      const_iterator begin() const { return &ele[0]; }
      const_iterator end() const { return &ele[M]; }
      const_iterator cbegin() const noexcept { return &ele[0]; }
      const_iterator cend() const noexcept { return &ele[M]; }

      AT& e(int i) {
	return ele[i];
      }

      /*! \brief Element operator
       *
       * See e(int) 
       * */
      const AT& e(int i) const {
	return ele[i];
      }

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling vector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling vector.
       * */
      inline Vector<Var,AT>& init(const AT& val=AT());
      inline Vector<Var,AT>& init(Init, const AT& a=AT()) { return init(a); }
      inline Vector<Var,AT>& init(Noinit, const AT& a=AT()) { return *this; }

      /*! \brief Size.
       *
       * \return The size of the vector.
       * */
      constexpr int size() const {return M;}

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return 1;}

      using Matrix<General,Var,Fixed<1>,AT>::operator();

      /*! \brief Cast to std::vector<AT>.
       *
       * \return The std::vector<AT> representation of the vector
       * */
      explicit inline operator std::vector<AT>() const;

      /*! \brief std::vector<AT> Constructor.
       * Constructs and initializes a vector with a std::vector<AT> object.
       * \param v The std::vector<AT> the vector will be initialized with. 
       * */
      explicit inline Vector(const std::vector<AT> &v);

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the vector
//       * */
//      explicit operator AT() const {
//        FMATVEC_ASSERT(M==1, AT);
//        return ele[0];
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a vector with a AT object.
//       * \param x The AT the vector will be initialized with.
//       * */
//      explicit Vector(const AT &x) : Matrix<General,Var,Fixed<1>,AT>(x) { }

#if defined(SWIG) && SWIG_VERSION < 0x040000
      inline const RowVector<Var,AT> T() const;
      inline RowVector<Var,AT> T();
#else
      inline const RowVector<Var,AT> T() const &;
      inline RowVector<Var,AT> T() &;
#endif
      // move
#ifndef SWIG
      inline RowVector<Var,AT> T() &&;
#endif

      inline const Vector<Var,AT> operator()(const Range<Var,Var> &I) const;

      template <class Row> inline void set(const Range<Var,Var> &I, const Vector<Row,AT> &x);

      template <class Row> inline void add(const Range<Var,Var> &I, const Vector<Row,AT> &x);

      inline const Vector<Var,AT> operator()(const Indices &I) const;

      template<class Row> inline void set(const Indices &I, const Vector<Row,AT> &x);
  };

  template <class AT>
    inline Vector<Var,AT>& Vector<Var,AT>::init(const AT &val) {
        for(int i=0; i<M; i++) 
          e(i) = val;
      return *this;
    }

  template <class AT>
#if defined(SWIG) && SWIG_VERSION < 0x040000
    inline const RowVector<Var,AT> Vector<Var,AT>::T() const {
#else
    inline const RowVector<Var,AT> Vector<Var,AT>::T() const & {
#endif
      RowVector<Var,AT> x(M,NONINIT);
      for(int i=0; i<M; i++)
        x.e(i) = e(i);
      return x;
    }
  template <class AT>
#if defined(SWIG) && SWIG_VERSION < 0x040000
    inline RowVector<Var,AT> Vector<Var,AT>::T() {
#else
    inline RowVector<Var,AT> Vector<Var,AT>::T() & {
#endif
      RowVector<Var,AT> x(M,NONINIT);
      for(int i=0; i<M; i++)
        x.e(i) = e(i);
      return x;
    }
  // move
#ifndef SWIG
  template <class AT>
    inline RowVector<Var,AT> Vector<Var,AT>::T() && {
      return RowVector<Var,AT>(Transpose(), std::move(*this));
    }
#endif

  template <class AT>
    inline Vector<Var,AT>::operator std::vector<AT>() const {
      std::vector<AT> ret(size());
      for(int i=0; i<size(); ++i)
        ret[i] = e(i);
      return ret;
    }

  template <class AT>
    inline Vector<Var,AT>::Vector(const std::vector<AT> &v) : Matrix<General,Var,Fixed<1>,AT>(static_cast<int>(v.size()),1,NONINIT) {
      for(int i=0; i<size(); ++i)
        e(i) = v[i];
    }

  template <class AT>
    inline const Vector<Var,AT> Vector<Var,AT>::operator()(const Range<Var,Var> &I) const {
      FMATVEC_ASSERT(I.end()<M, AT);
      Vector<Var,AT> x(I.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        x.e(i) = e(I.start()+i);

      return x;
    }

  template <class AT> template <class Row>
    inline void Vector<Var,AT>::set(const Range<Var,Var> &I, const Vector<Row,AT> &x) {
      FMATVEC_ASSERT(I.end()<size(), AT);
      FMATVEC_ASSERT(I.size()==x.size(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) = x.e(ii);
    }

  template <class AT> template <class Row>
    inline void Vector<Var,AT>::add(const Range<Var,Var> &I, const Vector<Row,AT> &x) {
      FMATVEC_ASSERT(I.end()<size(), AT);
      FMATVEC_ASSERT(I.size()==x.size(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) += x.e(ii);
    }

  template <class AT>
    inline const Vector<Var,AT> Vector<Var,AT>::operator()(const Indices &I) const {
      FMATVEC_ASSERT(I.max()<size(), AT);

      Vector<Var,AT> x(I.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
	x.e(i) = e(I[i]);

      return x;
    }

  template <class AT> template <class Row>
    inline void Vector<Var,AT>::set(const Indices &I, const Vector<Row,AT> &x) {
      FMATVEC_ASSERT(I.max()<size(), AT);
      FMATVEC_ASSERT(I.size()==x.size(), AT);
      for(int i=0; i<I.size(); i++)
	e(I[i]) = x.e(i);
    }

  /// @cond NO_SHOW

  template <class AT> template <class Row>
    inline Vector<Var,AT>& Vector<Var,AT>::copy(const Vector<Row,AT> &x) {
      for(int i=0; i<M; i++)
        e(i) = x.e(i);
      return *this;
    }

  /// @endcond

}

#endif
