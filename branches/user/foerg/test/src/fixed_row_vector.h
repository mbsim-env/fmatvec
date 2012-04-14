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

#include "general_fixed_matrix.h"

namespace fmatvec {

   template <int M, class AT> class FixedVector;

  /*! 
   *  \brief This is a rowvector class of general shape in dense storage format.
   *
   * Template class FixedRowVector of type GeneralFixed<1,N>, id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the rowvector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <int N, class AT> class FixedRowVector : public Matrix<GeneralFixed<1,N>, AT> {

    using Matrix<GeneralFixed<1,N>, AT>::ele;

    public:

    /// @cond NO_SHOW
    
    //friend class Vector<AT>;

    //friend RowVector<AT> Matrix<General, AT>::row(int i);
    //friend const RowVector<AT> Matrix<General, AT>::row(int i) const;

    protected:

    void deepCopy(const FixedRowVector<N,AT> &x);
    void deepCopy(const RowVector<AT> &x);

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a rowvector with no size. 
       * */
      FixedRowVector() : Matrix<GeneralFixed<1,N>, AT>() {
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a vector with a string in a matlab-like
       * notation. The entries are seperated by semicolons.
       * For example
       * \code 
       * FixedRowVector<double> x("[3,1,2]");
       * \endcode
       * constructs the vector
       * \f[ x=\begin{pmatrix}3 1 2\end{pmatrix}  \f]
       * \param str The string the vector will be initialized with. 
       * */
      FixedRowVector(const char *str) : Matrix<GeneralFixed<1,N>, AT>(str) {
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
      FixedRowVector(Initialization ini, const AT &a=0) : Matrix<GeneralFixed<1,N>, AT>(ini,a) {
      }

      /*! \brief Copy Constructor
       *
       * See FixedRowVector(const FixedRowVector<AT>&) 
       * */
      explicit FixedRowVector(const Matrix<GeneralFixed<1,N>, AT> &A) : Matrix<GeneralFixed<1,N>, AT>(A) {
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the rowvector \em x.
       * \attention The physical memory of the rowvector
       * \em x will not be copied, only referenced.
       * \param x The rowvector that will be referenced.
       * */
      FixedRowVector(const FixedRowVector<N,AT> &x) : Matrix<GeneralFixed<1,N>, AT>(x) {
      }

      /*! \brief Copy operator
       *
       * Copies the rowvector given by \em x.
       * \param x The rowvector to be copied. 
       * \return A reference to the calling rowvector.
       * */
      FixedRowVector<N,AT>& operator<<(const FixedRowVector<N,AT> &x);

      /*! \brief Assignment operator
       *
       * Copies the rowvector given by \em x by calling operator<<().
       * \param x The rowvector to be assigned. 
       * \return A reference to the calling rowvector.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      FixedRowVector<N,AT>& operator=(const FixedRowVector<N,AT> &x) {
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
	assert(i<N);
#endif
	return ele[i];
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
	return ele[i];
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling rowvector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling rowvector.
       * */
      FixedRowVector<N,AT>& init(const AT& a);

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

      /*! \brief Rowvector duplicating.
       *
       * The calling rowvector returns a \em deep copy of itself.  
       * \return The duplicate of the calling rowvector.
       * */
      FixedRowVector<N,AT> copy() const;

      using Matrix<GeneralFixed<1,N>, AT>::operator();

      FixedVector<N,AT> T() {
	FixedVector<N,AT> x(NONINIT);
	for(int i=0; i<N; i++)
	  x()[i] = ele[i];
      };

      const FixedVector<N,AT> T() const {
	FixedVector<N,AT> x(NONINIT);
	for(int i=0; i<N; i++)
	  x()[i] = ele[i];
      }
  };

  /////////////////////////////////////////////////////////////////////

  template <int N, class AT>
    FixedRowVector<N,AT>& FixedRowVector<N,AT>::init(const AT& val) {

      for(int i=0; i<N; i++) 
	ele[i] = val;

      return *this;
    }

  template <int N, class AT>
    FixedRowVector<N,AT>& FixedRowVector<N,AT>::operator<<(const FixedRowVector<N,AT> &x) { 

      deepCopy(x);

      return *this;
    }

  template <int N, class AT>
    FixedRowVector<N,AT> FixedRowVector<N,AT>::copy() const {

      FixedRowVector<N,AT> x(NONINIT);
      x.deepCopy(*this);

      return x;
    }

  template <int N, class AT>
    void FixedRowVector<N,AT>::deepCopy(const FixedRowVector<N,AT> &x) {
      for(int i=0; i<N; i++)
	ele[i] = x.ele[i];
    }

}
#endif