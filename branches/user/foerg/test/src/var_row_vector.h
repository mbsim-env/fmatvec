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

#ifndef var_row_vector_h
#define var_row_vector_h

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
  template <class AT> class RowVector<Var,AT> : public Matrix<General,Fixed<1>,Var,AT> {
    using Matrix<General,Fixed<1>,Var,AT>::N;
    using Matrix<General,Fixed<1>,Var,AT>::ele;

    public:

    /// @cond NO_SHOW

    protected:

    template<class Row> inline void deepCopy(const RowVector<Row,AT> &x);

    /// @endcond
    
    public:

      /*! \brief Standard constructor
       *
       * Constructs a vector with no size. 
       * */
      RowVector() : Matrix<General,Fixed<1>,Var,AT>() {
      }

//      template<class Ini=All<AT> >
//        RowVector(int n, Ini ini=All<AT>()) : Matrix<General,Fixed<1>,Var,AT>(n,ini) { } 

      RowVector(int n) : Matrix<General,Fixed<1>,Var,AT>(n) { } 
      RowVector(int n, Noinit &ini) : Matrix<General,Fixed<1>,Var,AT>(n,ini) { } 
      RowVector(int n, const All<AT> &ini) : Matrix<General,Fixed<1>,Var,AT>(n,ini) { } 
      RowVector(int n, const Eye<AT> &ini) : Matrix<General,Fixed<1>,Var,AT>(n,ini) { } 

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a vector with a string in a matlab-like
       * notation. The entries are seperated by semicolons.
       * For example
       * \code 
       * RowVector<double> x("[3;1;2]");
       * \endcode
       * constructs the vector
       * \f[ x=\begin{pmatrix}3\\ 1 \\ 2\end{pmatrix}  \f]
       * \param str The string the vector will be initialized with. 
       * */
      RowVector(const char *str) : Matrix<General,Fixed<1>,Var,AT>(str) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      RowVector(const RowVector<Var,AT> &x) : Matrix<General,Fixed<1>,Var,AT>(x) {
      }

      template<class Row>
      RowVector(const RowVector<Row,AT> &x) : Matrix<General,Fixed<1>,Var,AT>(x) {
      }

      template<class Type, class Row, class Col>
      explicit RowVector(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Fixed<1>,Var,AT>(A) {
      }

//      template<class Ini=All<AT> >
//        RowVector<Var,AT>& resize(int n=0, Ini ini=All<AT>()) {
//          Matrix<General,Fixed<1>,Var,AT>::resize(1,n,ini);
//          return *this;
//      }

      RowVector<Var,AT>& resize(int n=0) {
        Matrix<General,Fixed<1>,Var,AT>::resize(1,n);
        return *this;
      }

      template<class Ini>
        RowVector<Var,AT>& resize(int n, const Ini &ini) {
          Matrix<General,Fixed<1>,Var,AT>::resize(1,n,ini);
          return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x .
       * \param x The vector to be assigned. 
       * \return A reference to the calling vector.
       * */
      inline RowVector<Var,AT>& operator=(const RowVector<Var,AT> &x);

      template <class Row>
      inline RowVector<Var,AT>& operator=(const RowVector<Row,AT> &x);

      /*! \brief Copy operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be copied. 
       * \return A reference to the calling vector.
       * */
      template <class Row>
        inline RowVector<Var,AT>& operator<<(const RowVector<Row,AT> &x);

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
	assert(i>=0);
	assert(i<N);
#endif

	return e(i);
      };

      /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(i<N);
#endif

	return e(i);
      };

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
      inline RowVector<Var,AT>& init(const AT& a);
      inline RowVector<Var,AT>& init(const All<AT> &all) { return init(all.a); }
      inline RowVector<Var,AT>& init(Eye<AT> eye);
      inline RowVector<Var,AT>& init(Noinit) { return *this; }

      /*! \brief Size.
       *
       * \return The size of the vector.
       * */
      int size() const {return N;};

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return 1;};

      using Matrix<General,Fixed<1>,Var,AT>::operator();

      /*! \brief Cast to std::vector<AT>.
       *
       * \return The std::vector<AT> representation of the vector
       * */
      inline operator std::vector<AT>();

      /*! \brief std::vector<AT> Constructor.
       * Constructs and initializes a vector with a std::vector<AT> object.
       * \param v The std::vector<AT> the vector will be initialized with. 
       * */
      inline RowVector(std::vector<AT> v);

      inline const Vector<Var,AT> T() const;

      template <class Row>
      inline void set(const Range<Var,Var> &I, const RowVector<Row,AT> &x);
  };

  template <class AT>
    inline RowVector<Var,AT>& RowVector<Var,AT>::init(const AT &val) {
      for(int i=0; i<N; i++) 
        e[i] = val; 
      return *this;
    }

  template <class AT>
    inline RowVector<Var,AT>& RowVector<Var,AT>::operator=(const RowVector<Var,AT> &x) { 

#ifndef FMATVEC_RESIZE_VOID
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(N == x.size());
#endif
#else
      if(N==0) {
        delete[] ele;
        N = x.size(); 
        ele = new AT[N];
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(N == x.size());
#endif
      }
#endif

      deepCopy(x);

      return *this;
    }

  template <class AT> template<class Row>
    inline RowVector<Var,AT>& RowVector<Var,AT>::operator=(const RowVector<Row,AT> &x) { 

#ifndef FMATVEC_RESIZE_VOID
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(N == x.size());
#endif
#else
      if(N==0) {
        delete[] ele;
        N = x.size(); 
        ele = new AT[N];
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(N == x.size());
#endif
      }
#endif

      deepCopy(x);

      return *this;
    }

  template <class AT> template<class Row>
    inline RowVector<Var,AT>& RowVector<Var,AT>::operator<<(const RowVector<Row,AT> &x) { 

      if(N!=x.size()) {
        delete[] ele;
        N = x.size(); 
        ele = new AT[N];
      }

      deepCopy(x);

      return *this;
    }

  template <class AT>
    inline const Vector<Var,AT> RowVector<Var,AT>::T() const {
      Vector<Var,AT> x(NONINIT);
      for(int i=0; i<N; i++)
        x.e(i) = e(i);
      return x;
    }

  template <class AT>
    inline RowVector<Var,AT>::operator std::vector<AT>() {
      std::vector<AT> ret(size());
      if(size()>0) memcpy(&ret[0], &operator()(0), sizeof(AT)*size());
      return ret;
    }

  template <class AT>
    inline RowVector<Var,AT>::RowVector(std::vector<AT> v) : Matrix<General,Fixed<1>,Var,AT>(v.size(),1) {
      if(size()>0) memcpy(&operator()(0), &v[0], sizeof(AT)*size());
    }

  template <class AT> template <class Row>
    inline void RowVector<Var,AT>::set(const Range<Var,Var> &I, const RowVector<Row,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<size());
      assert(I.size()==x.size());
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) = x.e(ii);
    }

  /// @cond NO_SHOW

  template <class AT> template <class Row>
    inline void RowVector<Var,AT>::deepCopy(const RowVector<Row,AT> &x) {
      for(int i=0; i<N; i++)
        e(i) = x.e(i);
    }

  /// @endcond

}

#endif
