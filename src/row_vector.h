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
 *   mfoerg@users.berlios.de
 *
 */

#ifndef row_vector_h
#define row_vector_h

#include "general_matrix.h"

namespace fmatvec {
  /*! 
   *  \brief This is a rowvector class of general shape in dense storage format.
   *
   * Template class RowVector of type General, id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the rowvector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class RowVector : public Matrix<General, AT> {
    using Matrix<General, AT>::m;
    using Matrix<General, AT>::n;
    using Matrix<General, AT>::lda;
    using Matrix<General, AT>::ele;
    using Matrix<General, AT>::tp;
    using Matrix<General, AT>::memory;
    using Matrix<General, AT>::elePtr;
    public:

    /*! \brief Rowvector transposing. 
     *
     * Returns the transpose of the rowvector \em x. 
     \param x The rowvector, that will be transposed.
     \return The transpose of rowvector \em x.
     * */
    template <class T> friend const Vector<T> trans(const RowVector<T> &x); 

    /*! \brief Rowvector transposing. 
     *
     * See trans(const RowVector<T> &)
     * */
    template <class T> friend Vector<T> trans(RowVector<T> &x); 

    /*! \brief Vector transposing. 
     *
     * Returns the transpose of the vector \em x. 
     \param x The vector, that will be transposed.
     \return The transpose of vector \em x.
     * */
    template <class T> friend RowVector<T> trans(Vector<T> &x);

    /*! \brief Vector transposing. 
     *
     * See trans(const Vector<T> &)
     * */
    template <class T> friend const RowVector<T> trans(const Vector<T> &x);

    /// @cond NO_SHOW
    friend RowVector<AT> Matrix<General, AT>::row(int i);
    friend const RowVector<AT> Matrix<General, AT>::row(int i) const;
    /// @endcond

    protected:

    /*! \internal 
     * */
    void deepCopy(const RowVector<AT> &x);

    /*! \internal 
     * */
    const AT* elePtr(int i) const {
      return tp ? ele+i : ele+lda*i;
    };

    /*! \internal 
     * */
    AT* elePtr(int i) {
      return tp ? ele+i : ele+lda*i;
    };

    /*! \internal 
     * */
    RowVector(int n_, int lda_, bool tp, Memory<AT> memory, const AT* ele_) : Matrix<General, AT>(1,n_, lda_, tp, memory, ele_) {
    }

    public:

      /*! \brief Standard constructor
       *
       * Constructs a rowvector with no size. 
       * */
      RowVector() : Matrix<General, AT>() {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a rowvector of size n. 
       * \param n The size.
       * \remark The rowvector will be initialised to
       * zero by default. This default behavior can be
       * changed by defining FMATVEC_NO_INITIALIZATION.
       * */
      RowVector(int n) : Matrix<General, AT>(1,n) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a rowvector of size m with the pyhsical memory given by \em ele_;
       * \param n The size.
       * \param ele The physical memory the rowvector will point to.
       * */
      RowVector(int n, AT* ele) : Matrix<General, AT>(1,n,ele) {
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
      RowVector(int n, Initialization ini, const AT &a=0) : Matrix<General, AT>(1,n,ini,a) {
      }

      /*! \brief Copy Constructor
       *
       * See RowVector(const RowVector<AT>&) 
       * */
      explicit RowVector(const Matrix<General, AT> &x) : Matrix<General, AT>(x) {
	LogicAssert(x.rows()==1,"Number of cols() must be 1 in RowVector<AT>::RowVector(const Matrix<General, AT> &x)");
      }

      /*! \brief Rowvector resizing. 
       *
       * Resizes the rowvector to size n.    
       * \param n The size.
       * \return A reference to the calling rowvector.
       * \remark The rowvector will be initialised to
       * zero by default. To change this behavior, define
       * FMATVEC_NO_INITIALIZATION.
       * */
      RowVector& resize(int n_) {
	Matrix<General, AT>::resize(1,n);
	return *this;
      }

      /*! \brief Rowvector resizing. 
       *
       * Resizes the rowvector to size n. 
       * The vector will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * vector will not be initialized.
       * \param n The size.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the rowvector will be initialized with (default 0)
       * \return A reference to the calling rowvector.
       * */
      RowVector& resize(int n, Initialization ini, const AT &a=0) {
	Matrix<General, AT>::resize(1,n,ini,a);
	return *this;
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the rowvector \em x.
       * \attention The physical memory of the rowvector
       * \em x will not be copied, only referenced.
       * \param x The rowvector that will be referenced.
       * */
      RowVector(const RowVector<AT> &x) : Matrix<General, AT>(x) {
      }

      /*! \brief Copy operator
       *
       * Copies the rowvector given by \em x.
       * \param A The rowvector to be copied. 
       * \return A reference to the calling rowvector.
       * */
      RowVector<AT>& operator<<(const RowVector<AT> &x);

      /*! \brief Reference operator
       *
       * References the rowvector given by \em x.
       * \param A The rowvector to be referenced. 
       * \return A reference to the calling rowvector.
       * */
      RowVector<AT>& operator>>(const RowVector<AT> &x);

      /*! \brief Assignment operator
       *
       * Copies the rowvector given by \em x by calling operator<<().
       * \param x The rowvector to be assigned. 
       * \return A reference to the calling rowvector.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      RowVector<AT>& operator=(const RowVector<AT> &x) {
#ifndef FMATVEC_NO_DEEP_ASSIGNMENT 
	return operator<<(x);
#else
	return operator>>(x);
#endif
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
	assert(i<n);
#endif
	return tp ? ele[i] : ele[i*lda];
      };

     /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(i<n);
#endif
	return tp ? ele[i] : ele[i*lda];
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling rowvector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling rowvector.
       * */
      RowVector<AT>& init(const AT& val);

      /*! \brief Size.
       *
       * \return The size of the rowvector.
       * */
      int size() const {return n;};

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return tp?1:lda;};

      /*! \brief Rowvector duplicating.
       *
       * The calling rowvector returns a \em deep copy of itself.  
       * \return The duplicate of the calling rowvector.
       * */
      RowVector<AT> copy() const;

      /*! \brief Subrowvector operator.
       *
       * Returns a subrowvector of the calling rowvector. 
       * \attention The subrowvector and the
       * calling rowvector will share the same physical memory.
       * \param i1 The starting element. 
       * \param i2 The ending element.
       * \return A subrowvector of the calling rowvector.
       * */
      RowVector<AT> operator()(int i1, int i2);

      /*! \brief Subvector operator.
       *
       * See operator()(int,int);
       * */
      const RowVector<AT> operator()(int i1, int i2) const;

      /*! \brief Subrowvector operator.
       *
       * Returns a subrowvector of the calling rowvector. 
       * \attention The subrowvector and the
       * calling rowvector will share the same physical memory.
       * \param I Index containing the starting and the ending element. 
       * \return A subrowvector of the calling rowvector.
       * */
      RowVector<AT> operator()(const Index &I);

      /*! \brief Subrowvector operator.
       *
       * See operator()(const Index&)
       * */
      const RowVector<AT> operator()(const Index &I) const;

      using Matrix<General, AT>::operator();
      using Matrix<General, AT>::resize;
  };

  /////////////////////////////////////////////////////////////////////

  template <class AT>
    RowVector<AT>& RowVector<AT>::init(const AT& val) {

      for(int i=0; i<n; i++) 
	operator()(i) = val;

      return *this;
    }

  template <class AT>
    RowVector<AT>& RowVector<AT>::operator<<(const RowVector<AT> &x) { 

      if(x.size() == 0)
	return *this;

      if(n==0) {
	m = x.m; 
	n = x.n;
	lda = m;
	tp = false;
	memory.resize(n);
	ele = (AT*)memory.get();
      } else {
#ifdef FMATVEC_SIZE_CHECK
	assert(n == x.n);
#endif
      }

      deepCopy(x);

      return *this;
    }

  template <class AT>
    RowVector<AT>& RowVector<AT>::operator>>(const RowVector<AT> &x) { 

      if(n==0) {
	m = x.m; 
	n = x.n;
      } else {
#ifdef FMATVEC_SIZE_CHECK
	assert(n == x.n);
#endif
      }

      memory = x.memory;
      ele = x.ele;
      lda = x.lda;
      tp = x.tp; 

      return *this;
    }

  template <class AT>
    RowVector<AT> RowVector<AT>::copy() const {

      RowVector<AT> x(m);
      x.deepCopy(*this);

      return x;
    }

  template <class AT>
    void RowVector<AT>::deepCopy(const RowVector<AT> &x) {
      for(int i=0; i<size(); i++)
	operator()(i) = x(i);
    }

  template <> extern void RowVector<double>::deepCopy(const RowVector<double> &A);

  template <class AT> RowVector<AT> RowVector<AT>::operator()(int i1, int i2) {
    return operator()(Index(i1,i2));
  }

  template <class AT> const RowVector<AT> RowVector<AT>::operator()(int i1, int i2) const {
    return operator()(Index(i1,i2));
  }

  template <class AT>
    const RowVector<AT> RowVector<AT>::operator()(const Index &I) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<n);
#endif

      return RowVector<AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  template <class AT>
    RowVector<AT> RowVector<AT>::operator()(const Index &I) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<n);
#endif

      return RowVector<AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }



}
#endif
