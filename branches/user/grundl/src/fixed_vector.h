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
  class Vector<Fixed<M>, AT> : public Matrix<General, Fixed<M>, Fixed<1>, AT> {
      using Matrix<General, Fixed<M>, Fixed<1>, AT>::ele;

    public:

      /// @cond NO_SHOW

    protected:

      template <class Row> inline void deepCopy(const Vector<Row, AT> &x);

      /// @endcond

    public:

//      template<class Ini=All<AT> >
//        Vector(Ini ini=All<AT>()) : Matrix<General,Fixed<M>,Fixed<1>,AT>(ini) { }
//      template<class Ini=All<AT> >
//        Vector(int m, Ini ini=All<AT>()) : Matrix<General,Fixed<M>,Fixed<1>,AT>(ini) { }

      Vector(Noinit ini) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(ini) {
      }
      Vector(Init ini = INIT, const AT &a = 0) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(ini, a) {
      }
      Vector(int m, Noinit ini) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(ini) {
      }
      Vector(int m, Init ini = INIT, const AT &a = 0) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(ini, a) {
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
      Vector(const char *str) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(str) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      template <class Row>
      Vector(const Vector<Row, AT> &x) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(x) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General, AT>&) 
       * */
      template <class Type, class Row, class Col>
      explicit Vector(const Matrix<Type, Row, Col, AT> &A) :
          Matrix<General, Fixed<M>, Fixed<1>, AT>(A) {
      }

      /*!
       * \brief standard virtual destructor
       */
      virtual ~Vector() {
      }

      template <class Row>
      inline Vector<Fixed<M>, AT>& operator=(const Vector<Row, AT> &x);

      /*! \brief Element operator
       *
       * Returns a reference to the i-th element. 
       * \param i The i-th element.
       * \return A reference to the element x(i).
       * \remark The bounds are checked by default. 
       * To change this behavior, define
       * FMATVEC_NO_BOUNDS_CHECK.
       * \sa operator()(int) const
       * */
      AT& operator()(int i) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
        assert(i >= 0);
        assert(i < M);
#endif

        return e(i);
      }
      ;

      /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
        assert(i >= 0);
        assert(i < M);
#endif

        return e(i);
      }
      ;

      AT& e(int i) {
        return ele[i][0];
      }
      ;

      /*! \brief Element operator
       *
       * See e(int) 
       * */
      const AT& e(int i) const {
        return ele[i][0];
      }
      ;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling vector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling vector.
       * */
      inline Vector<Fixed<M>, AT>& init(const AT& a = 0);
      inline Vector<Fixed<M>, AT>& init(Init, const AT& a = 0) {
        return init(a);
      }
      inline Vector<Fixed<M>, AT>& init(Noinit, const AT& a = 0) {
        return *this;
      }

      /*! \brief Size.
       *
       * \return The size of the vector.
       * */
      int size() const {
        return M;
      }
      ;

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {
        return 1;
      }
      ;

      template <int M1, int M2>
      inline const Vector<Fixed<M2 - M1 + 1>, AT> operator()(const Range<Fixed<M1>, Fixed<M2> > &I) const;

      inline const Vector<Var, AT> operator()(const Index &I) const;

      using Matrix<General, Fixed<M>, Fixed<1>, AT>::operator();

      /*! \brief Cast to std::vector<AT>.
       *
       * \return The std::vector<AT> representation of the vector
       * */
      inline operator std::vector<AT>();

      /*! \brief std::vector<AT> Constructor.
       * Constructs and initializes a vector with a std::vector<AT> object.
       * \param v The std::vector<AT> the vector will be initialized with. 
       * */
      inline Vector(std::vector<AT> v);

      /*!
       * \brief return the transpose of the Vector, i.e. a RowVector
       */
      inline const RowVector<Fixed<M>, AT> T() const;

      /*!
       * \brief set a subvector - specified by the range - to the given vector
       */
      template <class Type, class Row, class Col>
      inline void set(const Range<Var, Var> &I, const Matrix<Type, Row, Col, AT> &A);

      /*!
       * \brief set a subvector - specified by the range - to the given vector
       */
      template <int M1, int M2, class Row>
      inline void set(const Range<Fixed<M1>, Fixed<M2> > &I, const Vector<Row, AT> &A);

  };
  // END Class definition

  template <int M, class AT>
  inline Vector<Fixed<M>, AT>& Vector<Fixed<M>, AT>::init(const AT &val) {
    for (int i = 0; i < M; i++)
      e(i) = val;
    return *this;
  }

  template <int M, class AT> template <int M1, int M2>
  inline const Vector<Fixed<M2 - M1 + 1>, AT> Vector<Fixed<M>, AT>::operator()(const Range<Fixed<M1>, Fixed<M2> > &I) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
    assert(M2 < M);
#endif
    Vector<Fixed<M2 - M1 + 1>, AT> x(NONINIT);

    for (int i = 0; i < x.size(); i++)
      x.e(i) = e(M1 + i);

    return x;
  }

  template <int M, class AT>
  inline const Vector<Var, AT> Vector<Fixed<M>, AT>::operator()(const Index &I) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
    assert(I.end() < M);
#endif
    Vector<Var, AT> x(I.end() - I.start() +1, NONINIT);

    for (int i = 0; i < x.size(); i++)
      x.e(i) = e(I.start() + i);

    return x;
  }

  template <int M, class AT> template <class Row>
  inline Vector<Fixed<M>, AT>& Vector<Fixed<M>, AT>::operator=(const Vector<Row, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
    assert(x.size() == M);
#endif

    deepCopy(x);

    return *this;
  }

  template <int M, class AT>
  inline const RowVector<Fixed<M>, AT> Vector<Fixed<M>, AT>::T() const {
    RowVector<Fixed<M>, AT> x(NONINIT);
    for (int i = 0; i < M; i++)
      x.e(i) = e(i);
    return x;
  }

  template <int M, class AT>
  inline Vector<Fixed<M>, AT>::operator std::vector<AT>() {
    std::vector < AT > ret(size());
    if (size() > 0)
      memcpy(&ret[0], &operator()(0), sizeof(AT) * size());
    return ret;
  }

  template <int M, class AT>
  inline Vector<Fixed<M>, AT>::Vector(std::vector<AT> v) :
      Matrix<General, Fixed<M>, Fixed<1>, AT>() {
    if (size() > 0)
      memcpy(&operator()(0), &v[0], sizeof(AT) * size());
  }

  /// @cond NO_SHOW

  template <int M, class AT> template <class Row>
  inline void Vector<Fixed<M>, AT>::deepCopy(const Vector<Row, AT> &x) {
    for (int i = 0; i < M; i++)
      e(i) = x.e(i);
  }

  /// @endcond

  template <int M, class AT> template <class Type, class Row, class Col>
  inline void Vector<Fixed<M>, AT>::set(const Range<Var, Var> &I, const Matrix<Type, Row, Col, AT> &A) {
    Matrix<General, Fixed<M>, Fixed<1>, AT>::set(I, Index(0, 0), A);
  }

  template <int M, class AT>
  template <int M1, int M2, class Row>
  inline void Vector<Fixed<M>, AT>::set(const Range<Fixed<M1>, Fixed<M2> > &I, const Vector<Row, AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(M2-M1+1 == x.size());
#endif
      for(int i = M1; i < M2 + 1; i++)
        e(i) = x.e(i-M1);

  }

}

#endif
