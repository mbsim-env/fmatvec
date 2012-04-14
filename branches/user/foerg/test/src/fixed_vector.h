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

#include "general_fixed_matrix.h"
#include <vector>
#include <cstring>

namespace fmatvec {

  template <int N, class AT> class FixedRowVector;

  /*! 
   *  \brief This is a vector class of general shape in dense storage format.
   *
   * Template class Vector of type General, id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the vector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <int M, class AT> class FixedVector : public Matrix<GeneralFixed<M,1>, AT> {
    using Matrix<GeneralFixed<M,1>, AT>::ele;

    public:

    /// @cond NO_SHOW

    //template <class T> friend RowVector<T> trans(const Vector<T> &x); 
    //template <class T> friend Vector<T> trans(const RowVector<T> &x);     
    
    //friend class RowVector<AT>;

    //friend Vector<AT> Matrix<General, AT>::col(int i);
    //friend const Vector<AT> Matrix<General, AT>::col(int i) const;

    protected:

    void deepCopy(const FixedVector<M,AT> &x);
    void deepCopy(const Vector<AT> &x);

    /// @endcond
    
    public:

      /*! \brief Standard constructor
       *
       * Constructs a vector with no size. 
       * */
      FixedVector() : Matrix<GeneralFixed<M,1>, AT>() {
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
      FixedVector(Initialization ini, const AT &a=0) : Matrix<GeneralFixed<M,1>, AT>(ini,a) {
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
      FixedVector(const char *str) : Matrix<GeneralFixed<M,1>, AT>(str) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      FixedVector(const FixedVector<M,AT> &x) : Matrix<GeneralFixed<M,1>, AT>(x) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      explicit FixedVector(const Vector<AT> &x) : Matrix<GeneralFixed<M,1>, AT>(x) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<AT>&) 
       * */
      explicit FixedVector(const Matrix<GeneralFixed<M,1>, AT> &A) : Matrix<GeneralFixed<M,1>, AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<AT>&) 
       * */
      explicit FixedVector(const Matrix<General, AT> &A) : Matrix<GeneralFixed<M,1>, AT>(A) {
      }

      /*! \brief Copy operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be copied. 
       * \return A reference to the calling vector.
       * */
      FixedVector<M,AT>& operator<<(const FixedVector<M,AT> &x);

      FixedVector<M,AT>& operator<<(const Vector<AT> &x);

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x by calling operator<<().
       * \param x The vector to be assigned. 
       * \return A reference to the calling vector.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      FixedVector<M,AT>& operator=(const FixedVector<M,AT> &x) {
	return operator<<(x);
      }

      FixedVector<M,AT>& operator=(const Vector<AT> &x) {
	return operator<<(x);
      }

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
      FixedVector<M,AT>& init(const AT& a); 

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

      /*! \brief Vector duplicating.
       *
       * The calling vector returns a \em deep copy of itself.  
       * \return The duplicate of the calling vector.
       * */
      FixedVector<M,AT> copy() const;

      using Matrix<GeneralFixed<M,1>, AT>::operator();

      /*! \brief Cast to std::vector<AT>.
       *
       * \return The std::vector<AT> representation of the vector
       * */
      operator std::vector<AT>();

      /*! \brief std::vector<AT> Constructor.
       * Constructs and initializes a vector with a std::vector<AT> object.
       * \param v The std::vector<AT> the vector will be initialized with. 
       * */
      FixedVector(std::vector<AT> v);

      FixedRowVector<M,AT> T() {
	FixedRowVector<M,AT> x(NONINIT);
	for(int i=0; i<M; i++)
	  x.e(i) = e(i);
	return x;
      };

      const FixedRowVector<M,AT> T() const {
	FixedRowVector<M,AT> x(NONINIT);
	for(int i=0; i<M; i++)
	  x.e(i) = e(i);
	return x;
      }

  };

  template <int M, class AT>
    FixedVector<M,AT>& FixedVector<M,AT>::init(const AT& val) {

      for(int i=0; i<M; i++) 
	e(i) = val; 

      return *this;
    }

  template <int M, class AT>
    FixedVector<M,AT>& FixedVector<M,AT>::operator<<(const FixedVector<M,AT> &x) { 

      deepCopy(x);

      return *this;
    }

  template <int M, class AT>
    FixedVector<M,AT>& FixedVector<M,AT>::operator<<(const Vector<AT> &x) { 

#ifdef FMATVEC_SIZE_CHECK
       assert(x.size() == M);
#endif

      deepCopy(x);

      return *this;
    }

  template <int M, class AT>
    FixedVector<M,AT> FixedVector<M,AT>::copy() const {

      FixedVector<M,AT> x(NONINIT);
      x.deepCopy(*this);

      return x;
    }

  template <int M, class AT>
    void FixedVector<M,AT>::deepCopy(const Vector<AT> &x) {
      for(int i=0; i<M; i++)
	e(i) = x.operator()(i);
    }

  template <int M, class AT>
    void FixedVector<M,AT>::deepCopy(const FixedVector<M,AT> &x) {
      for(int i=0; i<M; i++)
	e(i) = x.e(i);
    }

//  template <class AT>
//    Vector<AT>::operator std::vector<AT>() {
//      std::vector<AT> ret(size());
//      if(size()>0) memcpy(&ret[0], &operator()(0), sizeof(AT)*size());
//      return ret;
//    }
//
//  template <class AT>
//    Vector<AT>::Vector(std::vector<AT> v) : Matrix<General, AT>(v.size(),1) {
//      if(size()>0) memcpy(&operator()(0), &v[0], sizeof(AT)*size());
//    }
}

#endif
