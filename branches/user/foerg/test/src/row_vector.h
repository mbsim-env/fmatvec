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

#ifndef row_vector_h
#define row_vector_h

#include "general_matrix.h"

namespace fmatvec {
  /*! 
   *  \brief This is a rowvector class of general shape in dense storage format.
   *
   * Template class RowVector of type General,Ref,Ref,id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the rowvector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class RowVector<Ref,AT> : public Matrix<General,Ref,Ref,AT> {
    using Matrix<General,Ref,Ref,AT>::m;
    using Matrix<General,Ref,Ref,AT>::n;
    using Matrix<General,Ref,Ref,AT>::lda;
    using Matrix<General,Ref,Ref,AT>::ele;
    using Matrix<General,Ref,Ref,AT>::tp;
    using Matrix<General,Ref,Ref,AT>::memory;
    using Matrix<General,Ref,Ref,AT>::elePtr;

    public:

    /// @cond NO_SHOW
    
    template <class T> friend Vector<Ref,T> trans(const RowVector<Ref,T> &x); 
    template <class T> friend RowVector<Ref,T> trans(const Vector<Ref,T> &x);
    //friend RowVector<Ref,double> trans(const Vector<Ref,double> &x);

    friend class Vector<Ref,AT>;

    friend RowVector<Ref,AT> Matrix<General,Ref,Ref,AT>::row(int i);
    friend const RowVector<Ref,AT> Matrix<General,Ref,Ref,AT>::row(int i) const;

    protected:

    inline void deepCopy(const RowVector<Ref,AT> &x);

    const AT* elePtr(int i) const {
      return tp ? ele+i : ele+lda*i;
    };

    AT* elePtr(int i) {
      return tp ? ele+i : ele+lda*i;
    };

    RowVector(int n_, int lda_, bool tp, Memory<AT> memory, const AT* ele_) : Matrix<General,Ref,Ref,AT>(1,n_, lda_, tp, memory, ele_) {
    }

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a rowvector with no size. 
       * */
      RowVector() : Matrix<General,Ref,Ref,AT>(1,0) {
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a vector with a string in a matlab-like
       * notation. The entries are seperated by semicolons.
       * For example
       * \code 
       * RowVector<double> x("[3,1,2]");
       * \endcode
       * constructs the vector
       * \f[ x=\begin{pmatrix}3 1 2\end{pmatrix}  \f]
       * \param str The string the vector will be initialized with. 
       * */
      RowVector(const char *str) : Matrix<General,Ref,Ref,AT>(str) {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(m==1);
#endif
      }

      /*! \brief Regular Constructor
       *
       * Constructs a rowvector of size n. 
       * \param n The size.
       * \remark The rowvector will be initialised to
       * zero by default. This default behavior can be
       * changed by defining FMATVEC_NO_INITIALIZATION.
       * */
      RowVector(int n) : Matrix<General,Ref,Ref,AT>(1,n) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a rowvector of size m with the pyhsical memory given by \em ele_;
       * \param n The size.
       * \param ele The physical memory the rowvector will point to.
       * */
      RowVector(int n, AT* ele) : Matrix<General,Ref,Ref,AT>(1,n,ele) {
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
      RowVector(int n, Initialization ini, const AT &a=0) : Matrix<General,Ref,Ref,AT>(1,n,ini,a) {
      }

      /*! \brief Copy Constructor
       *
       * See RowVector(const RowVector<Ref,AT>&) 
       * */
      explicit RowVector(const Matrix<General,Ref,Ref,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {
	LogicAssert(x.rows()==1,"Number of cols() must be 1 in RowVector<Ref,AT>::RowVector(const Matrix<General,Ref,Ref,AT> &x)");
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
      RowVector& resize(int n, Initialization ini=INIT, const AT &a=0) {
	Matrix<General,Ref,Ref,AT>::resize(1,n,ini,a);
	return *this;
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the rowvector \em x.
       * \attention The physical memory of the rowvector
       * \em x will not be copied, only referenced.
       * \param x The rowvector that will be referenced.
       * */
      RowVector(const RowVector<Ref,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {
      }

      /*! \brief Copy operator
       *
       * Copies the rowvector given by \em x.
       * \param x The rowvector to be copied. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator<<(const RowVector<Ref,AT> &x);

      /*! \brief Reference operator
       *
       * References the rowvector given by \em x.
       * \param x The rowvector to be referenced. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator>>(const RowVector<Ref,AT> &x);

      /*! \brief Assignment operator
       *
       * Copies the rowvector given by \em x by calling operator<<().
       * \param x The rowvector to be assigned. 
       * \return A reference to the calling rowvector.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      inline RowVector<Ref,AT>& operator=(const RowVector<Ref,AT> &x);

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
	return e(i);
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
	return e(i);
      };

      AT& er(int i) {
	return ele[i*lda];
      };

      const AT& er(int i) const {
	return ele[i*lda];
      };

      AT& et(int i) {
	return ele[i];
      };

      const AT& et(int i) const {
	return ele[i];
      };

       AT& e(int i) {
	return tp ? et(i) : er(i);
      };

      const AT& e(int i) const {
	return tp ? et(i) : er(i);
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling rowvector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& init(const AT& a);

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
      inline RowVector<Ref,AT> copy() const;

      /*! \brief Subrowvector operator.
       *
       * Returns a subrowvector of the calling rowvector. 
       * \attention The subrowvector and the
       * calling rowvector will share the same physical memory.
       * \param i1 The starting element. 
       * \param i2 The ending element.
       * \return A subrowvector of the calling rowvector.
       * */
      inline RowVector<Ref,AT> operator()(int i1, int i2);

      /*! \brief Subvector operator.
       *
       * See operator()(int,int);
       * */
      inline const RowVector<Ref,AT> operator()(int i1, int i2) const;

      /*! \brief Subrowvector operator.
       *
       * Returns a subrowvector of the calling rowvector. 
       * \attention The subrowvector and the
       * calling rowvector will share the same physical memory.
       * \param I Index containing the starting and the ending element. 
       * \return A subrowvector of the calling rowvector.
       * */
      inline RowVector<Ref,AT> operator()(const Index &I);

      /*! \brief Subrowvector operator.
       *
       * See operator()(const Index&)
       * */
      inline const RowVector<Ref,AT> operator()(const Index &I) const;

      using Matrix<General,Ref,Ref,AT>::operator();

      Vector<Ref,AT> T() {
	return Vector<Ref,AT>(n,lda,tp?false:true,memory,ele);
      }

      const Vector<Ref,AT> T() const {
	return Vector<Ref,AT>(n,lda,tp?false:true,memory,ele);
      }
  };

  /////////////////////////////////////////////////////////////////////

  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::init(const AT& val) {

      if(tp) {
	for(int i=0; i<n; i++) 
	  ele[i] = val; // operator()(i) = val;
      }
      else {
	for(int i=0; i<n; i++) 
	  ele[i*lda] = val;
      }

      return *this;
    }

  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::operator=(const RowVector<Ref,AT> &x) { 

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(n == x.size());
#endif

      deepCopy(x);

      return *this;
    }

  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::operator<<(const RowVector<Ref,AT> &x) { 

      if(n!=x.size()) {
	m = 1; 
	n = x.size();
	lda = m;
	tp = false;
	memory.resize(n);
	ele = (AT*)memory.get();
      } 

      deepCopy(x);

      return *this;
    }

  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::operator>>(const RowVector<Ref,AT> &x) { 

      m = 1;
      n = x.size();
      memory = x.memory;
      ele = x.ele;
      lda = x.lda;
      tp = x.tp; 

      return *this;
    }

  template <class AT>
    inline RowVector<Ref,AT> RowVector<Ref,AT>::copy() const {

      RowVector<Ref,AT> x(n,NONINIT);
      x.deepCopy(*this);

      return x;
    }

  template <class AT> 
    inline RowVector<Ref,AT> RowVector<Ref,AT>::operator()(int i1, int i2) {
      return operator()(Index(i1,i2));
    }

  template <class AT> 
    inline const RowVector<Ref,AT> RowVector<Ref,AT>::operator()(int i1, int i2) const {
      return operator()(Index(i1,i2));
    }

  template <class AT>
    inline const RowVector<Ref,AT> RowVector<Ref,AT>::operator()(const Index &I) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<n);
#endif

      return RowVector<Ref,AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  template <class AT>
    inline RowVector<Ref,AT> RowVector<Ref,AT>::operator()(const Index &I) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<n);
#endif

      return RowVector<Ref,AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  /// @cond NO_SHOW

  template <class AT>
    inline void RowVector<Ref,AT>::deepCopy(const RowVector<Ref,AT> &x) {
      if(tp) {
        if(x.tp) 
          for(int i=0; i<size(); i++)
            et(i) = x.et(i); 
        else 
          for(int i=0; i<size(); i++)
            et(i) = x.er(i);
      }
      else {
        if(x.tp)
          for(int i=0; i<size(); i++)
            er(i) = x.et(i); 
        else 
          for(int i=0; i<size(); i++)
            er(i) = x.er(i);
      }
    }

  /// @endcond


}
#endif
