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
  template <class AT> class Vector<fmatvec::Var,AT> : public Matrix<fmatvec::General,fmatvec::Var,fmatvec::Fixed<1>,AT> {
    using Matrix<General,Var,Fixed<1>,AT>::M;
    using Matrix<General,Var,Fixed<1>,AT>::ele;

    public:

    typedef AT* iterator;
    typedef const AT* const_iterator;

    typedef AT value_type;

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

//      template<class Ini=All<AT> >
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
        assert(M == x.size());
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
        assert(M == x.size());
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

      iterator begin() { return &ele[0]; }
      iterator end() { return &ele[M]; }
      const_iterator begin() const { return &ele[0]; }
      const_iterator end() const { return &ele[M]; }
      const_iterator cbegin() const noexcept { return &ele[0]; }
      const_iterator cend() const noexcept { return &ele[M]; }

      AT& e(int i) {
	return ele[i];
      };

      /*! \brief Element operator
       *
       * See e(int) 
       * */
      const AT& e(int i) const {
	return ele[i];
      };

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
      int size() const {return M;};

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return 1;};

      using Matrix<General,Var,Fixed<1>,AT>::operator();

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

      inline const RowVector<Var,AT> T() const;

      inline const Vector<Var,AT> operator()(const fmatvec::Range<Var,Var> &I) const;

      template <class Row>
      inline void set(const fmatvec::Range<Var,Var> &I, const Vector<Row,AT> &x);

      template <class Row>
      inline void add(const fmatvec::Range<Var,Var> &I, const Vector<Row,AT> &x);
  };

  template <class AT>
    inline Vector<Var,AT>& Vector<Var,AT>::init(const AT &val) {
        for(int i=0; i<M; i++) 
          e(i) = val;
      return *this;
    }

  template <class AT>
    inline const RowVector<Var,AT> Vector<Var,AT>::T() const {
      RowVector<Var,AT> x(M,NONINIT);
      for(int i=0; i<M; i++)
        x.e(i) = e(i);
      return x;
    }

  template <class AT>
    inline Vector<Var,AT>::operator std::vector<AT>() const {
      std::vector<AT> ret(size());
      if(size()>0) memcpy(&ret[0], &operator()(0), sizeof(AT)*size());
      return ret;
    }

  template <class AT>
    inline Vector<Var,AT>::Vector(std::vector<AT> v) : Matrix<General,Var,Fixed<1>,AT>(v.size(),1) {
      if(size()>0) memcpy(&operator()(0), &v[0], sizeof(AT)*size());
    }

  template <class AT>
    inline const Vector<Var,AT> Vector<Var,AT>::operator()(const fmatvec::Range<Var,Var> &I) const {
      assert(I.end()<M);
      Vector<Var,AT> x(I.end()-I.start()+1,NONINIT);

      for(int i=0; i<x.size(); i++)
        x.e(i) = e(I.start()+i);

      return x;
    }

  template <class AT> template <class Row>
    inline void Vector<Var,AT>::set(const fmatvec::Range<Var,Var> &I, const Vector<Row,AT> &x) {
      assert(I.end()<size());
      assert(I.size()==x.size());

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) = x.e(ii);
    }

  template <class AT> template <class Row>
    inline void Vector<Var,AT>::add(const fmatvec::Range<Var,Var> &I, const Vector<Row,AT> &x) {
      assert(I.end()<size());
      assert(I.size()==x.size());

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) += x.e(ii);
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
