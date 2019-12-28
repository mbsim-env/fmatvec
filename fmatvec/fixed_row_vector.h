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

#ifndef fixed_row_vector_h
#define fixed_row_vector_h

#include "fixed_general_matrix.h"
#include <stdexcept>

namespace fmatvec {

  /*! 
   *  \brief This is a rowvector class of general shape in dense storage format.
   *
   * Template class RowVector of type General,Fixed<1>,Fixed<N>,id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the rowvector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <int N, class AT> class RowVector<Fixed<N>,AT> : public Matrix<General,Fixed<1>,Fixed<N>,AT> {

    using Matrix<General,Fixed<1>,Fixed<N>,AT>::ele;

    public:

      typedef AT* iterator;
      typedef const AT* const_iterator;

      typedef AT value_type;

    /// @cond NO_SHOW
    
    protected:

    template<class Col> inline RowVector<Fixed<N>,AT>& copy(const RowVector<Col,AT> &x);

    /// @endcond

    public:

//      template<class Ini=All<AT>>
//        RowVector(Ini ini=All<AT>()) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini) { }
//      template<class Ini=All<AT>>
//        RowVector(int n, Ini ini=All<AT>()) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini) { }

      explicit RowVector(Noinit ini) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini) { }
      explicit RowVector(Init ini=INIT, const AT &a=AT()) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini,a) { }
      explicit RowVector(int n, Noinit ini) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini) { assert(n==N); }
      explicit RowVector(int n, Init ini=INIT, const AT &a=AT()) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini,a) { assert(n==N); }

      /*! \brief Copy Constructor
       *
       * See RowVector(const RowVector<Fixed<N>,AT> 
       * */
      RowVector(const RowVector<Fixed<N>,AT> &A) = default;

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General, AT>&) 
       * */
      template<class Row>
      RowVector(const RowVector<Row,AT> &A) : Matrix<General,Fixed<1>,Fixed<N>,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General, AT>&) 
       * */
      template<class Type, class Row, class Col>
      explicit RowVector(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Fixed<1>,Fixed<N>,AT>(A) {
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a vector with a string in a matlab-like
       * notation. The entries are seperated by semicolons.
       * For example
       * \code 
       * RowVector<Fixed<N>,AT>e> x("[3,1,2]");
       * \endcode
       * constructs the vector
       * \f[ x=\begin{pmatrix}3 1 2\end{pmatrix}  \f]
       * \param str The string the vector will be initialized with. 
       * */
      RowVector(const char *str) : Matrix<General,Fixed<1>,Fixed<N>,AT>(str) {
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be assigned.
       * \return A reference to the calling vector.
       * */
      inline RowVector<Fixed<N>,AT>& operator=(const RowVector<Fixed<N>,AT> &x) = default;

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be assigned.
       * \return A reference to the calling vector.
       * */
      template <class Row>
      inline RowVector<Fixed<N>,AT>& operator=(const RowVector<Row,AT> &x) {
        assert(x.size() == N);
        return copy(x);
      }

      /*! \brief Matrix assignment
       *
       * Copies the band matrix given by \em A.
       * \param A The matrix to be copied.
       * \return A reference to the calling matrix.
       * */
      template <class Row>
      inline RowVector<Fixed<N>,AT>& operator<<=(const RowVector<Row,AT> &x) { return operator=(x); }

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
	assert(i<N);
	return e(i);
      };

     /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {
	assert(i>=0);
	assert(i<N);
	return e(i);
      };

      iterator begin() { return &ele[0][0]; }
      iterator end() { return &ele[0][N]; }
      const_iterator begin() const { return &ele[0][0]; }
      const_iterator end() const { return &ele[0][N]; }
      const_iterator cbegin() const noexcept { return &ele[0][0]; }
      const_iterator cend() const noexcept { return &ele[0][N]; }

      AT& e(int i) {
	return ele[0][i];
      };

      /*! \brief Element operator
       *
       * See e(int) 
       * */
      const AT& e(int i) const {
	return ele[0][i];
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling rowvector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Fixed<N>,AT>& init(const AT &val=AT()); 
      inline RowVector<Fixed<N>,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline RowVector<Fixed<N>,AT>& init(Noinit, const AT &a=AT()) { return *this; }

      /*! \brief Size.
       *
       * \return The size of the rowvector.
       * */
      int size() const {return N;};

      //! Resize a fixed vector
      //! Do nothing for the fixed dimension and throw on any other dimension.
      void resize(int n) {
        if(n!=N)
          throw std::runtime_error("A fixed row vector cannot be resized.");
      }

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return 1;};

      using Matrix<General,Fixed<1>,Fixed<N>,AT>::operator();

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the vector
//       * */
//      explicit operator AT() const {
//        assert(N==1);
//        return ele[0];
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a vector with a AT object.
//       * \param x The AT the vector will be initialized with.
//       * */
//      explicit RowVector(const AT &x) : Matrix<General,Fixed<1>,Fixed<N>,AT>(x) { }

      inline const Vector<Fixed<N>,AT> T() const;
  };

  template <int N, class AT>
    inline RowVector<Fixed<N>,AT>& RowVector<Fixed<N>,AT>::init(const AT &val) {
      for(int i=0; i<N; i++) 
	e(i) = val; 
      return *this;
    }

  template <int N, class AT>
    inline const Vector<Fixed<N>,AT> RowVector<Fixed<N>,AT>::T() const {
      Vector<Fixed<N>,AT> x(NONINIT);
      for(int i=0; i<N; i++)
        x.e(i) = e(i);
      return x;
    }

   /// @cond NO_SHOW
 
  template <int N, class AT> template <class Col>
    inline RowVector<Fixed<N>,AT>& RowVector<Fixed<N>,AT>::copy(const RowVector<Col,AT> &x) {
      for(int i=0; i<N; i++)
        e(i) = x.e(i);
      return *this;
    }

  /// @endcond

}
#endif
