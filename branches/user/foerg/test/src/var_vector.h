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
  template <class AT> class Vector<General,Var,Fixed<1>,AT> : public Matrix<General,Var,Fixed<1>,AT> {
    using Matrix<General,Var,Fixed<1>,AT>::M;
    using Matrix<General,Var,Fixed<1>,AT>::ele;

    public:

    /// @cond NO_SHOW

    protected:

    template<class Type, class Row, class Col> inline void deepCopy(const Vector<Type,Row,Col,AT> &x);

    /// @endcond
    
    public:

      /*! \brief Standard constructor
       *
       * Constructs a vector with no size. 
       * */
      Vector(int m=0) : Matrix<General,Var,Fixed<1>,AT>(m) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a vector of size m. The vector will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * vector will not be initialized.
       * \param m The size.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the vector will be initialized with (default 0)
       * */
      Vector(int M, Initialization ini, const AT &a=0) : Matrix<General,Var,Fixed<1>,AT>(M,ini,a) {
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

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      Vector(const Vector<General,Var,Fixed<1>,AT> &x) : Matrix<General,Var,Fixed<1>,AT>(x) {
      }

      template<class Row>
      Vector(const Vector<General,Row,Fixed<1>,AT> &x) : Matrix<General,Var,Fixed<1>,AT>(x) {
      }

      template<class Type, class Row, class Col>
      explicit Vector(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Var,Fixed<1>,AT>(A) {
      }

      Vector<General,Var,Fixed<1>,AT>& resize(int n) {
	Matrix<General,Var,Fixed<1>,AT>::resize(n,1);
	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x .
       * \param x The vector to be assigned. 
       * \return A reference to the calling vector.
       * */
      inline Vector<General,Var,Fixed<1>,AT>& operator=(const Vector<General,Var,Fixed<1>,AT> &x);

      template <class Type, class Row, class Col>
      inline Vector<General,Var,Fixed<1>,AT>& operator=(const Vector<Type,Row,Col,AT> &x);

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
	assert(i<M);
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
	assert(i<M);
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
      Vector<General,Var,Fixed<1>,AT>& init(const AT& a); 

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
      inline operator std::vector<AT>();

      /*! \brief std::vector<AT> Constructor.
       * Constructs and initializes a vector with a std::vector<AT> object.
       * \param v The std::vector<AT> the vector will be initialized with. 
       * */
      inline Vector(std::vector<AT> v);

      inline const RowVector<General,Var,Fixed<1>,AT> T() const;

      template <class Type, class Row, class Col>
      inline void set(const Range<Var,Var> &I, const Vector<Type,Row,Col,AT> &x);
  };

  template <class AT>
    inline Vector<General,Var,Fixed<1>,AT>& Vector<General,Var,Fixed<1>,AT>::init(const AT& val) {

      for(int i=0; i<M; i++) 
        e(i) = val; 

      return *this;
    }

  template <class AT>
    inline Vector<General,Var,Fixed<1>,AT>& Vector<General,Var,Fixed<1>,AT>::operator=(const Vector<General,Var,Fixed<1>,AT> &x) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(M == x.M);
#endif

      deepCopy(x);

      return *this;
    }

  template <class AT> template<class Type, class Row, class Col>
    inline Vector<General,Var,Fixed<1>,AT>& Vector<General,Var,Fixed<1>,AT>::operator=(const Vector<Type,Row,Col,AT> &x) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(M == x.size());
#endif

      deepCopy(x);

      return *this;
    }

  template <class AT>
    inline const RowVector<General,Var,Fixed<1>,AT> Vector<General,Var,Fixed<1>,AT>::T() const {
      RowVector<General,Var,Fixed<1>,AT> x(NONINIT);
      for(int i=0; i<M; i++)
        x.e(i) = e(i);
      return x;
    }

  template <class AT>
    inline Vector<General,Var,Fixed<1>,AT>::operator std::vector<AT>() {
      std::vector<AT> ret(size());
      if(size()>0) memcpy(&ret[0], &operator()(0), sizeof(AT)*size());
      return ret;
    }

  template <class AT>
    inline Vector<General,Var,Fixed<1>,AT>::Vector(std::vector<AT> v) : Matrix<General,Var,Fixed<1>,AT>(v.size(),1) {
      if(size()>0) memcpy(&operator()(0), &v[0], sizeof(AT)*size());
    }

  template <class AT> template <class Type, class Row, class Col>
    inline void Vector<General,Var,Fixed<1>,AT>::set(const Range<Var,Var> &I, const Vector<Type,Row,Col,AT> &x) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<size());
      assert(I.size()==x.size());
#endif

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) = x.e(ii);
    }

  /// @cond NO_SHOW

  template <class AT> template <class Type, class Row, class Col>
    inline void Vector<General,Var,Fixed<1>,AT>::deepCopy(const Vector<Type,Row,Col,AT> &x) {
      for(int i=0; i<M; i++)
        e(i) = x.e(i);
    }

  /// @endcond

}

#endif
