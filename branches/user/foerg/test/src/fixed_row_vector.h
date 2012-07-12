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

    /// @cond NO_SHOW
    
    protected:

    template<class Col> inline void deepCopy(const RowVector<Col,AT> &x);

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a rowvector with no size. 
       * */
      RowVector() : Matrix<General,Fixed<1>,Fixed<N>,AT>() {
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

      /*! \brief Regular Constructor
       *
       * Constructs a vector of size n. The rowvector will be initialized 
       * to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * rowvector will not be initialized.
       * \param n The size.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the rowvector will be initialized with (default 0)
       * */
      RowVector(Initialization ini, const AT &a=0) : Matrix<General,Fixed<1>,Fixed<N>,AT>(ini,a) {
      }

      RowVector(int n, Initialization ini, const AT &a=0) : Matrix<General,Fixed<1>,Fixed<N>,AT>(1,n,ini,a) {
      }

      /*! \brief Copy Constructor
       *
       * See RowVector(const RowVector<Fixed<N>,AT> 
       * */
      template<class Col>
      RowVector(const RowVector<Col,AT> &A) : Matrix<General,Fixed<1>,Fixed<N>,AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General, AT>&) 
       * */
      template<class Type, class Row, class Col>
      explicit RowVector(const Matrix<Type,Row,Col,AT> &A) : Matrix<General,Fixed<1>,Fixed<N>,AT>(A) {
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
      inline RowVector<Fixed<N>,AT>& init(const AT& a);

      /*! \brief Size.
       *
       * \return The size of the rowvector.
       * */
      int size() const {return N;};

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return 1;};

      using Matrix<General,Fixed<1>,Fixed<N>,AT>::operator();

      inline const Vector<Fixed<N>,AT> T() const;
  };

  template <int N, class AT>
    inline RowVector<Fixed<N>,AT>& RowVector<Fixed<N>,AT>::init(const AT& val) {

      for(int i=0; i<N; i++) 
	ele[i] = val;

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
    inline void RowVector<Fixed<N>,AT>::deepCopy(const RowVector<Col,AT> &x) {
      for(int i=0; i<N; i++)
        e(i) = x.e(i);
    }

  /// @endcond

}
#endif
