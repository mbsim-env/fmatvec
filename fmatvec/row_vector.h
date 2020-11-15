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
    using Matrix<General,Ref,Ref,AT>::memory;
    using Matrix<General,Ref,Ref,AT>::elePtr;

    public:
    static constexpr bool isVector {true};

    typedef AT* iterator;
    typedef const AT* const_iterator;

    typedef AT value_type;

    /// @cond NO_SHOW
    
    friend class Vector<Ref,AT>;

    protected:

    template<class Col> inline RowVector<Ref,AT>& copy(const RowVector<Col,AT> &x);

    const AT* elePtr(int i) const {
      return ele+lda*i;
    }

    AT* elePtr(int i) {
      return ele+lda*i;
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
	FMATVEC_ASSERT(x.rows()==1, AT);
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
	FMATVEC_ASSERT(m==1, AT);
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
        FMATVEC_ASSERT(n == x.size(), AT);
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
        FMATVEC_ASSERT(n == x.size(), AT);
        return copy(x);
      }

      /*! \brief Reference operator
       *
       * References the rowvector given by \em x.
       * \param x The rowvector to be referenced. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator&=(RowVector<Ref,AT> &x) {
        m = 1;
        n = x.n;
        memory = x.memory;
        ele = x.ele;
        lda = x.lda;
        return *this;
      }

      /*! \brief Reference operator
       *
       * References the rowvector given by \em x.
       * \param x The rowvector to be referenced. 
       * \return A reference to the calling rowvector.
       * */
      inline RowVector<Ref,AT>& operator&=(Matrix<General,Ref,Ref,AT> &A) {
        FMATVEC_ASSERT(A.rows() == 1, AT);
        m = 1;
        n = A.cols();
        memory = A.memory;
        ele = A.ele;
        lda = A.lda;
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
	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(i<n, AT);
	return e(i);
      }

     /*! \brief Element operator
       *
       * See operator()(int) 
       * */
      const AT& operator()(int i) const {
	FMATVEC_ASSERT(i>=0, AT);
	FMATVEC_ASSERT(i<n, AT);
	return e(i);
      }

      iterator begin() { FMATVEC_ASSERT(lda==1, AT); return &ele[0]; }
      iterator end() { FMATVEC_ASSERT(lda==1, AT); return &ele[n*lda]; }
      const_iterator begin() const { FMATVEC_ASSERT(lda==1, AT); return &ele[0]; }
      const_iterator end() const { FMATVEC_ASSERT(lda==1, AT); return &ele[n*lda]; }
      const_iterator cbegin() const noexcept { FMATVEC_ASSERT(lda==1, AT); return &ele[0]; }
      const_iterator cend() const noexcept { FMATVEC_ASSERT(lda==1, AT); return &ele[n*lda]; }

      AT& e(int i) {
	return ele[i*lda];
      }

      const AT& e(int i) const {
	return ele[i*lda];
      }

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
      int size() const {return n;}

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return lda;}

      /*! \brief Subrowvector operator.
       * */
      inline const RowVector<Ref,AT> operator()(const Range<Var,Var> &I) const;

      using Matrix<General,Ref,Ref,AT>::operator();

//      /*! \brief Cast to AT.
//       *
//       * \return The AT representation of the vector
//       * */
//      explicit operator AT() const {
//        FMATVEC_ASSERT(n==1, AT);
//        return e(0);
//      }
//
//      /*! \brief AT Constructor.
//       * Constructs and initializes a vector with a AT object.
//       * \param x The AT the vector will be initialized with.
//       * */
//      explicit RowVector(const AT &x) : Matrix<General,Ref,Ref,AT>(x) { }

      template<class Row> inline void set(const fmatvec::Range<Var,Var> &I, const RowVector<Row,AT> &x);

      template<class Row> inline void add(const fmatvec::Range<Var,Var> &I, const RowVector<Row,AT> &x);

      inline void ref(RowVector<Ref,AT> &x, const fmatvec::Range<Var,Var> &I);

      inline void ref(Matrix<General,Ref,Ref,AT> &A, int i);

      inline void ref(Matrix<General,Ref,Ref,AT> &A, int i, const fmatvec::Range<Var,Var> &J);

      const Vector<Ref,AT> T() const;
  };

  /////////////////////////////////////////////////////////////////////

  template <class AT>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::init(const AT &val) {
      for(int i=0; i<n; i++)
        ele[i*lda] = val;
      return *this;
    }

  template <class AT>
    inline const RowVector<Ref,AT> RowVector<Ref,AT>::operator()(const Range<Var,Var> &I) const {

      FMATVEC_ASSERT(I.end()<n, AT);
      RowVector<Ref,AT> x(I.end()-I.start()+1,NONINIT);

      for(int i=0; i<x.size(); i++)
        x.e(i) = e(I.start()+i);

      return x;
    }

  /// @cond NO_SHOW

  template <class AT> template <class Col>
    inline RowVector<Ref,AT>& RowVector<Ref,AT>::copy(const RowVector<Col,AT> &x) {
      for(int i=0; i<size(); i++)
        e(i) = x.e(i);
      return *this;
    }

  template <class AT> template <class Row>
    inline void RowVector<Ref,AT>::set(const Range<Var,Var> &I, const RowVector<Row,AT> &x) {

      FMATVEC_ASSERT(I.end()<size(), AT);
      FMATVEC_ASSERT(I.size()==x.size(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) = x.e(ii);
    }

  template <class AT> template <class Row>

    inline void RowVector<Ref,AT>::add(const Range<Var,Var> &I, const RowVector<Row,AT> &x) {
      FMATVEC_ASSERT(I.end()<size(), AT);
      FMATVEC_ASSERT(I.size()==x.size(), AT);

      for(int i=I.start(), ii=0; i<=I.end(); i++, ii++)
        e(i) += x.e(ii);
    }

  template <class AT>
    inline void RowVector<Ref,AT>::ref(RowVector<Ref,AT> &x, const fmatvec::Range<Var,Var> &I) {

      FMATVEC_ASSERT(I.end()<x.size(), AT);

      m=1;
      n=I.size();
      memory = x.memory;
      ele = x.elePtr(I.start());
      lda = x.lda;
    }

  template <class AT>
    inline void RowVector<Ref,AT>::ref(Matrix<General,Ref,Ref,AT> &A, int i) {

      FMATVEC_ASSERT(i<A.rows(), AT);

      m=1;
      n=A.cols();
      memory = A.memory;
      ele = A.elePtr(i,0);
      lda = A.lda;
    }

  template <class AT>
    inline void RowVector<Ref,AT>::ref(Matrix<General,Ref,Ref,AT> &A, int i, const fmatvec::Range<Var,Var> &J) {

      FMATVEC_ASSERT(i<A.rows(), AT);
      FMATVEC_ASSERT(J.end()<A.cols(), AT);

      m=1;
      n=J.size();
      memory = A.memory;
      ele = A.elePtr(i,J.start());
      lda = A.lda;
    }

  template <class AT>
    inline const Vector<Ref,AT> RowVector<Ref,AT>::T() const {
      Vector<Ref,AT> x(n,NONINIT);
      for(int i=0; i<n; i++)
        x.e(i) = e(i);
      return x;
    }

  /// @endcond

}

#endif
