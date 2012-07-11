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

#ifndef vector_h
#define vector_h

#include "general_matrix.h"
#include <vector>
#include <cstring>

namespace fmatvec {

  /*! 
   *  \brief This is a vector class of general shape in dense storage format.
   *
   * Template class Vector of type General<Ref,Ref>, id. shape is general
   * and storage form is dense. The template parameter AT defines the
   * atomic type of the vector. Valid types are int, float,
   * double, complex<float> and complex<double> 
   * */
  template <class AT> class Vector<General,Ref,Fixed<1>,AT> : public Matrix<General,Ref,Ref,AT> {
    using Matrix<General,Ref,Ref,AT>::m;
    using Matrix<General,Ref,Ref,AT>::n;
    using Matrix<General,Ref,Ref,AT>::lda;
    using Matrix<General,Ref,Ref,AT>::ele;
    using Matrix<General,Ref,Ref,AT>::tp;
    using Matrix<General,Ref,Ref,AT>::memory;
    using Matrix<General,Ref,Ref,AT>::elePtr;

    public:

    /// @cond NO_SHOW

    template <class T> friend RowVector<General,Fixed<1>,Ref,T> trans(const Vector<General,Ref,Fixed<1>,T> &x); 
    template <class T> friend Vector<General,Ref,Fixed<1>,T> trans(const RowVector<General,Fixed<1>,Ref,T> &x);     
    
    friend class RowVector<General,Fixed<1>,Ref,AT>;

    friend Vector<General,Ref,Fixed<1>,AT> Matrix<General,Ref,Ref,AT>::col(int i);
    friend const Vector<General,Ref,Fixed<1>,AT> Matrix<General,Ref,Ref,AT>::col(int i) const;

    protected:

    inline void deepCopy(const Vector<General,Ref,Fixed<1>,AT> &x);
    template<class Type, class Row, class Col> inline void deepCopy(const Vector<Type,Row,Col,AT> &x);

    AT* elePtr(int i) {
      return tp ? ele+lda*i : ele+i;
    };

    const AT* elePtr(int i) const {
      return tp ? ele+lda*i : ele+i;
    };

    Vector(int n_, int lda_, bool tp, Memory<AT> memory, const AT* ele_) : Matrix<General,Ref,Ref,AT>(n_, 1, lda_, tp, memory, ele_) {
    }

    /// @endcond
    
    public:

      /*! \brief Standard constructor
       *
       * Constructs a vector with no size. 
       * */
      Vector() : Matrix<General,Ref,Ref,AT>() {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a vector of size m. 
       * \param m The size.
       * \remark The vector will be initialised to
       * zero by default. This default behavior can be
       * changed by defining FMATVEC_NO_INITIALIZATION.
       * */
      Vector(int m) : Matrix<General,Ref,Ref,AT>(m,1) {
      }

      /*! \brief Regular Constructor
       *
       * Constructs a vector of size m with the pyhsical memory given by \em ele_;
       * \param m The size.
       * \param ele The physical memory the vector will point to.
       * */
      Vector(int m, AT* ele) : Matrix<General,Ref,Ref,AT>(m,1,ele) { 
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
      Vector(int m, Initialization ini, const AT &a=0) : Matrix<General,Ref,Ref,AT>(m,1,ini,a) {
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
      Vector(const char *str) : Matrix<General,Ref,Ref,AT>(str) {
#ifndef FMATVEC_NO_SIZE_CHECK
	assert(n==1);
#endif
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the vector \em x.
       * \attention The physical memory of the vector
       * \em x will not be copied, only referenced.
       * \param x The vector that will be referenced.
       * */
      Vector(const Vector<General,Ref,Fixed<1>,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {
      }

      /*! \brief Copy Constructor
       *
       * See Vector(const Vector<General,Ref,Fixed<1>,AT>&) 
       * */
      explicit Vector(const Matrix<General,Ref,Ref,AT> &x) : Matrix<General,Ref,Ref,AT>(x) {

#ifndef FMATVEC_NO_SIZE_CHECK
	assert(x.cols()==1);
#endif
      }

     template<class Type, class Row, class Col>
      explicit Vector(const Matrix<Type,Row,Col,AT> &x) : Matrix<General,Ref,Ref,AT>(x)  {
      }

      /*! \brief Copy operator
       *
       * Copies the vector given by \em x.
       * \param x The vector to be copied. 
       * \return A reference to the calling vector.
       * */
      template <class Type, class Row, class Col>
        inline Vector<General,Ref,Fixed<1>,AT>& operator<<(const Vector<Type,Row,Col,AT> &x);

      /*! \brief Reference operator
       *
       * References the vector given by \em x.
       * \param x The vector to be referenced. 
       * \return A reference to the calling vector.
       * */
      inline Vector<General,Ref,Fixed<1>,AT>& operator>>(const Vector<General,Ref,Fixed<1>,AT> &x);

      /*! \brief Vector resizing. 
       *
       * Resizes the vector to size n.    
       * \param n The size.
       * \return A reference to the calling vector.
       * \remark The vector will be initialised to
       * zero by default. To change this behavior, define
       * FMATVEC_NO_INITIALIZATION.
       * */
      Vector<General,Ref,Fixed<1>,AT>& resize(int n) {
	Matrix<General,Ref,Ref,AT>::resize(n,1);
	return *this;
      }

      /*! \brief Vector resizing. 
       *
       * Resizes the vector to size n. The vector will be initialized to the value given by \em a
       * (default 0), if ini is set to INIT. If init is set to NONINIT, the
       * vector will not be initialized.
       * \param n The size.
       * \param ini INIT means initialization, NONINIT means no initialization.
       * \param a The value, the vector will be initialized with (default 0)
       * \return A reference to the calling vector.
       * */
      Vector<General,Ref,Fixed<1>,AT>& resize(int n, Initialization ini, const AT &a=0) {
	Matrix<General,Ref,Ref,AT>::resize(n,1,ini,a);
	return *this;
      }

      /*! \brief Assignment operator
       *
       * Copies the vector given by \em x by calling operator<<().
       * \param x The vector to be assigned. 
       * \return A reference to the calling vector.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      Vector<General,Ref,Fixed<1>,AT>& operator=(const Vector<General,Ref,Fixed<1>,AT> &x) {
#ifndef FMATVEC_NO_DEEP_ASSIGNMENT 
	return operator<<(x);
#else
	return operator>>(x);
#endif
      }

      template <class Type, class Row, class Col>
      Vector<General,Ref,Fixed<1>,AT>& operator=(const Vector<Type,Row,Col,AT> &x) {
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
	assert(i<m);
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
	assert(i<m);
#endif

	return e(i);
      };

      AT& er(int i) {
	return ele[i];
      };

      const AT& er(int i) const {
	return ele[i];
      };

      AT& et(int i) {
	return ele[i*lda];
      };

      const AT& et(int i) const {
	return ele[i*lda];
      };

       AT& e(int i) {
	return tp ? et(i) : er(i);
      };

      const AT& e(int i) const {
	return tp ? et(i) : er(i);
      };

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling vector
       * with the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling vector.
       * */
      inline Vector<General,Ref,Fixed<1>,AT>& init(const AT& a); 

      /*! \brief Size.
       *
       * \return The size of the vector.
       * */
      int size() const {return m;};

      /*! \brief Increment.
       *
       * \todo Docu
       *
       * \return The increment.
       * */
      int inc() const {return tp?lda:1;};

      /*! \brief Vector duplicating.
       *
       * The calling vector returns a \em deep copy of itself.  
       * \return The duplicate of the calling vector.
       * */
      inline Vector<General,Ref,Fixed<1>,AT> copy() const;

      /*! \brief Subvector operator.
       *
       * Returns a subvector of the calling vector. 
       * \attention The subvector and the
       * calling vector will share the same physical memory.
       * \param i1 The starting element. 
       * \param i2 The ending element.
       * \return A subvector of the calling vector.
       * */
      inline Vector<General,Ref,Fixed<1>,AT> operator()(int i1, int i2);

      /*! \brief Subvector operator.
       *
       * See operator()(int,int);
       * */
      inline const Vector<General,Ref,Fixed<1>,AT> operator()(int i1, int i2) const;

      /*! \brief Subvector operator.
       *
       * Returns a subvector of the calling vector. 
       * \attention The subvector and the
       * calling vector will share the same physical memory.
       * \param I Index containing the starting and the ending element. 
       * \return A subvector of the calling vector.
       * */
      inline Vector<General,Ref,Fixed<1>,AT> operator()(const Index &I);

      /*! \brief Subvector operator.
       *
       * See operator()(const Index&)
       * */
      inline const Vector<General,Ref,Fixed<1>,AT> operator()(const Index &I) const;

      using Matrix<General,Ref,Ref,AT>::operator();
      using Matrix<General,Ref,Ref,AT>::resize;

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

      RowVector<General,Fixed<1>,Ref,AT> T() {
        return RowVector<General,Fixed<1>,Ref,AT>(m,lda,tp?false:true,memory,ele);
      };

      const RowVector<General,Fixed<1>,Ref,AT> T() const {
        return RowVector<General,Fixed<1>,Ref,AT>(m,lda,tp?false:true,memory,ele);
      }

  };

  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT>& Vector<General,Ref,Fixed<1>,AT>::init(const AT& val) {

      if(tp) 
        for(int i=0; i<m; i++) 
          et(i) = val;
      else 
        for(int i=0; i<m; i++) 
          er(i) = val;

      return *this;
    }

  template <class AT> template<class Type, class Row, class Col>
    inline Vector<General,Ref,Fixed<1>,AT>& Vector<General,Ref,Fixed<1>,AT>::operator<<(const Vector<Type,Row,Col,AT> &x) { 

      if(x.size() == 0)
        return *this;

      if(m==0) {
        m = x.size(); 
        n = 1;
        lda = m;
        tp = false;
        memory.resize(m);
        ele = (AT*)memory.get();
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(m == x.size());
#endif
      }

      deepCopy(x);

      return *this;
    }

  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT>& Vector<General,Ref,Fixed<1>,AT>::operator>>(const Vector<General,Ref,Fixed<1>,AT> &x) { 

      if(m==0) {
        m = x.m; 
        n = x.n;
      } else {
#ifndef FMATVEC_NO_SIZE_CHECK
        assert(m == x.m);
#endif
      }

      memory = x.memory;
      ele = x.ele;
      lda = x.lda;
      tp = x.tp; 

      return *this;
    }

  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT> Vector<General,Ref,Fixed<1>,AT>::copy() const {

      Vector<General,Ref,Fixed<1>,AT> x(m,NONINIT);
      x.deepCopy(*this);

      return x;
    }

  template <class AT> 
    inline Vector<General,Ref,Fixed<1>,AT> Vector<General,Ref,Fixed<1>,AT>::operator()(int i1, int i2) {
    return operator()(Index(i1,i2));
  }

  template <class AT> 
    inline const Vector<General,Ref,Fixed<1>,AT> Vector<General,Ref,Fixed<1>,AT>::operator()(int i1, int i2) const {
    return operator()(Index(i1,i2));
  }

  template <class AT>
    inline const Vector<General,Ref,Fixed<1>,AT> Vector<General,Ref,Fixed<1>,AT>::operator()(const Index &I) const {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<m);
#endif

      return Vector<General,Ref,Fixed<1>,AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT> Vector<General,Ref,Fixed<1>,AT>::operator()(const Index &I) {

#ifndef FMATVEC_NO_BOUNDS_CHECK
      assert(I.end()<m);
#endif

      return Vector<General,Ref,Fixed<1>,AT>(I.end()-I.start()+1,lda,tp,memory,elePtr(I.start()));
    }

  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT>::operator std::vector<AT>() {
      std::vector<AT> ret(size());
      if(size()>0) memcpy(&ret[0], &operator()(0), sizeof(AT)*size());
      return ret;
    }

  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT>::Vector(std::vector<AT> v) : Matrix<General,Ref,Ref,AT>(v.size(),1) {
      if(size()>0) memcpy(&operator()(0), &v[0], sizeof(AT)*size());
    }

  /// @cond NO_SHOW

  template <class AT>
    inline void Vector<General,Ref,Fixed<1>,AT>::deepCopy(const Vector<General,Ref,Fixed<1>,AT> &x) {
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

  template <class AT> template <class Type, class Row, class Col>
    inline void Vector<General,Ref,Fixed<1>,AT>::deepCopy(const Vector<Type,Row,Col,AT> &x) {
      if(tp) 
	for(int i=0; i<size(); i++)
	  et(i) = x.e(i);
      else 
	for(int i=0; i<size(); i++)
	  er(i) = x.e(i);
    }

  /// @endcond

}

#endif
