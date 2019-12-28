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

    typedef AT* iterator;
    typedef const AT* const_iterator;

    typedef AT value_type;

    /// @cond NO_SHOW
    
    template <class T> friend Vector<Ref,T> trans(const RowVector<Ref,T> &x); 
    template <class T> friend RowVector<Ref,T> trans(const Vector<Ref,T> &x);
    //friend RowVector<Ref,double> trans(const Vector<Ref,double> &x);

    friend class Vector<Ref,AT>;

    friend RowVector<Ref,AT> Matrix<General,Ref,Ref,AT>::row(int i);
    friend const RowVector<Ref,AT> Matrix<General,Ref,Ref,AT>::row(int i) const;

    protected:

    inline RowVector<Ref,AT>& copy(const RowVector<Ref,AT> &x);
    template<class Col> inline RowVector<Ref,AT>& copy(const RowVector<Col,AT> &x);

    const AT* elePtr(int i) const {
      return tp ? ele+i : ele+lda*i;
    };

    AT* elePtr(int i) {
      return tp ? ele+i : ele+lda*i;
    };

    explicit RowVector(int n_, int lda_, bool tp, Memory<AT> memory, const AT* ele_) : Matrix<General,Ref,Ref,AT>(1,n_, lda_, tp, memory, ele_) {
    }

    /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs a rowvector with no size. 
       * */
      explicit RowVector() : Matrix<General,Ref,Ref,AT>() { m=1; }

      explicit RowVector(int n, Noinit ini) : Matrix<General,Ref,Ref,AT>(1,n,ini) { } 
      explicit RowVector(int n, Init ini=INIT, const AT &a=AT()) : Matrix<General,Ref,Ref,AT>(1,n,ini,a) { } 

      /*! \brief Copy Constructor
       *
       * Constructs a copy of the vector \em x.
       * \param x The vector that will be copied.
       * */
      RowVector(const RowVector<Ref,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {
      }

      template<class Row>
      RowVector(const RowVector<Row,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {
      }

      template<class Type, class Row, class Col>
      explicit RowVector(const Matrix<Type,Row,Col,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {
	assert(x.rows()==1);
      }

      /*! \brief Regular Constructor
       *
       * Constructs a rowvector of size m with the pyhsical memory given by \em ele_;
       * \param n The size.
       * \param ele The physical memory the rowvector will point to.
       * */
      explicit RowVector(int n, AT* ele) : Matrix<General,Ref,Ref,AT>(1,n,ele) {
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
	assert(m==1);
      }

      RowVector<Ref,AT>& resize(int n, Noinit) {
        Matrix<General,Ref,Ref,AT>::resize(1,n,Noinit());
        return *this;
      }

      RowVector<Ref,AT>& resize(int n, Init ini=INIT, const AT &a=AT()) {
        Matrix<General,Ref,Ref,AT>::resize(1,n,ini,a);
        return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the rowvector given by \em x.
       * \param x The rowvector to be assigned. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator=(const RowVector<Ref,AT> &x) {
        assert(n == x.size());
        return copy(x);
      }

      /*! \brief Assignment operator
       *
       * Copies the rowvector given by \em x.
       * \param x The rowvector to be assigned. 
       * \return A reference to the calling rowvector.
       * */
      template <class Row>
      inline RowVector<Ref,AT>& operator=(const RowVector<Row,AT> &x) {
        assert(n == x.size());
        return copy(x);
      }

      /*! \brief Reference operator
       *
       * References the rowvector given by \em x.
       * \param x The rowvector to be referenced. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator&=(const RowVector<Ref,AT> &x) {
        m = 1;
        n = x.n;
        memory = x.memory;
        ele = x.ele;
        lda = x.lda;
        tp = x.tp; 
        return *this;
      }

      /*! \brief Reference operator
       *
       * References the rowvector given by \em x.
       * \param x The rowvector to be referenced. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator&=(const Matrix<General,Ref,Ref,AT> &A) {
        assert(A.rows() == 1);
        m = 1;
        n = A.cols();
        memory = A.memory;
        ele = A.ele;
        lda = A.lda;
        tp = A.tp; 
        return *this;
      }

      /*! \brief Rowvector assignment
       *
       * Copies the rowvector given by \em x.
       * \param x The rowvector to be copied.
       * \return A reference to the calling rowvector.
       * */
      template <class Row>
        inline RowVector<Ref,AT>& operator<<=(const RowVector<Row,AT> &x) {
          if(n!=x.size()) resize(x.size(),NONINIT);
          return copy(x);
        }

      template <class AT2>
      operator Vector<Ref,AT2>() const {
        Vector<Ref,AT2> ret(size());
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
	assert(i<n);
	return e(i);
      };

     /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {
	assert(i>=0);
	assert(i<n);
	return e(i);
      };

      AT& er(int i) {
	return ele[i*lda];
      };

      iterator begin() { assert(lda==1); return &ele[0]; }
      iterator end() { assert(lda==1); return &ele[n*lda]; }
      const_iterator begin() const { assert(lda==1); return &ele[0]; }
      const_iterator end() const { assert(lda==1); return &ele[n*lda]; }
      const_iterator cbegin() const noexcept { assert(lda==1); return &ele[0]; }
      const_iterator cend() const noexcept { assert(lda==1); return &ele[n*lda]; }

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
      inline RowVector<Ref,AT>& init(const AT& val); 
      inline RowVector<Ref,AT>& init(Init, const AT &a=AT()) { return init(a); }
      inline RowVector<Ref,AT>& init(Noinit, const AT &a=AT()) { return *this; }

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
       * \param I Range containing the starting and the ending element.
       * \return A subrowvector of the calling rowvector.
       * */
      inline RowVector<Ref,AT> operator()(const Range<Var,Var> &I);

      /*! \brief Subrowvector operator.
       *
       * See operator()(const Range<Var,Var>&)
       * */
      inline const RowVector<Ref,AT> operator()(const Range<Var,Var> &I) const;

      using Matrix<General,Ref,Ref,AT>::operator();

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the vector
//       * */
//      explicit operator AT() const {
//        assert(n==1);
//        return e(0);
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a vector with a AT object.
//       * \param x The AT the vector will be initialized with.
//       * */
//      explicit RowVector(const AT &x) : Matrix<General,Ref,Ref,AT>(x) { }

      Vector<Ref,AT> T() {
	return Vector<Ref,AT>(n,lda,tp?false:true,memory,ele);
      }

      const Vector<Ref,AT> T() const {
	return Vector<Ref,AT>(n,lda,tp?false:true,memory,ele);
      }
  };

  /////////////////////////////////////////////////////////////////////


  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::init(const AT &val) {

      if(tp) {
	for(int i=0; i<n; i++) 
	  ele[i] = val; 
      }
      else {
	for(int i=0; i<n; i++) 
	  ele[i*lda] = val;
      }

      return *this;
    }

  template <class AT> 
    inline RowVector<Ref,AT> RowVector<Ref,AT>::operator()(int i1, int i2) {
      return operator()(Range<Var,Var>(i1,i2));
    }

  template <class AT> 
    inline const RowVector<Ref,AT> RowVector<Ref,AT>::operator()(int i1, int i2) const {
      return operator()(Range<Var,Var>(i1,i2));
    }

  template <class AT>
    inline const RowVector<Ref,AT> RowVector<Ref,AT>::operator()(const Range<Var,Var> &I) const {

      assert(I.end()<n);

      return RowVector<Ref,AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  template <class AT>
    inline RowVector<Ref,AT> RowVector<Ref,AT>::operator()(const Range<Var,Var> &I) {

      assert(I.end()<n);

      return RowVector<Ref,AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  /// @cond NO_SHOW

  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::copy(const RowVector<Ref,AT> &x) {
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
      return *this;
    }

  template <class AT> template <class Col>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::copy(const RowVector<Col,AT> &x) {
      if(tp) 
	for(int i=0; i<size(); i++)
	  et(i) = x.e(i);
      else 
	for(int i=0; i<size(); i++)
	  er(i) = x.e(i);
      return *this;
    }

  /// @endcond

}

#endif
