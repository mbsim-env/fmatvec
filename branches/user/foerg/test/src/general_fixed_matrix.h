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

#ifndef general_fixed_matrix_h
#define general_fixed_matrix_h

#include "types.h"
#include <stdlib.h>

namespace fmatvec {

//  template <class AT> class Vector;
//  template <class AT> class RowVector;
//  template <class AT> class SquareMatrix;
  //template <class AT> class Matrix<Symmetric, AT>;

  /*! 
   *  \brief This is a matrix class for general matrices.
   *  
   * Template class Matrix with shape type ThreeByThree and atomic type AT. The
   * storage form is dense. The template parameter AT defines the atomic type
   * of the matrix. Valid types are int, float, double, complex<float> and
   * complex<double> 
   * */
  template <class AT> class Matrix<ThreeByThree, AT> {

    public:

 /// @cond NO_SHOW

      friend class Matrix<Symmetric, AT>;
      
      template <class T> friend Matrix<ThreeByThree, T>  trans(const Matrix<ThreeByThree, T> &A);

    protected:

      AT ele[9];

      template <class Type> void deepCopy(const Matrix<Type, AT> &A); 

 /// @endcond
 
    public:

      /*! \brief Standard constructor
       *
       * Constructs a matrix with no size. 
       * */
      Matrix() {
#ifndef FMATVEC_NO_INITIALIZATION 
	init(0);
#endif
      }

      Matrix(Initialization ini, const AT &a=0) {  

	if(ini == INIT) {
	  ele[0] = a; ele[1] = a; ele[2] = a; ele[3] = a; ele[4] = a; ele[5] = a; ele[6] = a; ele[7] = a; ele[8] = a;
	} else if(ini == EYE ) {
	  ele[0] = 1; ele[1] = 0; ele[2] = 0; ele[3] = 0; ele[4] = 1; ele[5] = 0; ele[6] = 0; ele[7] = 0; ele[8] = 1;
	}
      }

      /*! \brief Copy Constructor
       *
       * Constructs a reference to the matrix \em A.
       * \attention The physical memory of the matrix \em A will not be copied, only
       * referenced.
       * \param A The matrix that will be referenced.
       * */
      Matrix(const Matrix<ThreeByThree, AT> &A) {
	ele[0] =  A.ele[0] ; ele[1] =  A.ele[1] ; ele[2] =  A.ele[2] ; ele[3] =  A.ele[3] ; ele[4] =  A.ele[4] ; ele[5] =  A.ele[5] ; ele[6] =  A.ele[6] ; ele[7] =  A.ele[7] ; ele[8] =  A.ele[8] ;
      }

      /*! \brief String Constructor. 
       *
       * Constructs and initializes a matrix with a string in a matlab-like
       * notation. The rows are seperated by semicolons, the columns by commas.
       * For example
       * \code 
       * Matrix<ThreeByThree,double> A("[3,2;1,2]");
       * \endcode 
       * constructs the matrix
       * \f[ A=\begin{pmatrix}3 & 2\\ 1 & 2\end{pmatrix}  \f]
       * \param str The string the matrix will be initialized with. 
       * */
      Matrix(const char *str);

      /*! \brief Destructor. 
       * */
      ~Matrix() {
      }

      /*! \brief Assignment operator
       *
       * Copies the matrix given by \em A by calling operator<<().
       * \param A The matrix to be assigned. 
       * \return A reference to the calling matrix.
       * \remark To call operator>>() by default, define FMATVEC_NO_DEEP_ASSIGNMENT
       * \sa operator<<(), operator>>()
       * */
      Matrix<ThreeByThree, AT>& operator=(const Matrix<ThreeByThree, AT> &A) {
	return operator<<(A);
      }

      template <class Type>
      Matrix<ThreeByThree, AT>& operator=(const Matrix<Type, AT> &A) {
	return operator<<(A);
      }

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      template<class T> Matrix<ThreeByThree, AT>& operator<<(const Matrix<T, AT> &A);

      /*! \brief Copy operator
       *
       * Copies the matrix given by \em A.
       * \param A The matrix to be copied. 
       * \return A reference to the calling matrix.
       * */
      Matrix<ThreeByThree, AT>& operator<<(const Matrix<ThreeByThree, AT> &A);

      /*! \brief Element operator
       *
       * Returns a reference to the element in the i-th row and the j-th column. 
       * \param i The i-th row of the matrix
       * \param j The j-th column of the matrix
       * \return A reference to the element A(i,j).
       * \remark The bounds are checked by default. 
       * To change this behavior, define
       * FMATVEC_NO_BOUNDS_CHECK.
       * \sa operator()(int,int) const
       * */
      AT& operator()(int i, int j) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<3);
	assert(j<3);
#endif

	return ele[i+j*3];
      };

      /*! \brief Element operator
       *
       * See operator()(int,int) 
       * */
      const AT& operator()(int i, int j) const {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i>=0);
	assert(j>=0);
	assert(i<3);
	assert(j<3);
#endif

	return ele[i+j*3];//  return ele[i*lda+j*ldb];
      };

      /*! \brief Pointer operator.
       *
       * Returns the pointer to the first element.
       * \return The pointer to the first element.
       * */
      AT* operator()() {return ele;};

      /*! \brief Pointer operator
       *
       * See operator()() 
       * */
      const AT* operator()() const {return ele;};

      /*! \brief Number of rows.
       *
       * \return The number of rows of the matrix.
       * */
      int rows() const {return 3;};

      /*! \brief Number of columns.
       *
       * \return The number of columns of the matrix.
       * */
      int cols() const {return 3;};

      /*! \brief Matrix duplicating.
       *
       * The calling matrix returns a \em deep copy of itself.  
       * \return The duplicate.
       * */
      Matrix<ThreeByThree, AT> copy() const;

      /*! \brief Initialization.
       *
       * Initializes all elements of the calling matrix with 
       * the value given by \em a.
       * \param a Value all elements will be initialized with.
       * \return A reference to the calling matrix.
       * */
      Matrix<ThreeByThree, AT>& init(const AT &a);

      /*! \brief Cast to std::vector<std::vector<AT> >.
       *
       * \return The std::vector<std::vector<AT> > representation of the matrix
       * */
      // operator std::vector<std::vector<AT> >();

      /*! \brief std::vector<std::vector<AT> > Constructor.
       * Constructs and initializes a matrix with a std::vector<std::vector<AT> > object.
       * An assert checks for constant length of each row.
       * \param m The std::vector<std::vector<AT> > the matrix will be initialized with. 
       * */
      Matrix(std::vector<std::vector<AT> > m);

      Matrix<ThreeByThree, AT> T() {
	Matrix<ThreeByThree, AT> A(NONINIT);
	A.ele[0] = ele[0]; A.ele[1] = ele[3]; A.ele[2] = ele[6]; A.ele[3] = ele[1]; A.ele[4] = ele[4]; A.ele[5] = ele[7]; A.ele[6] = ele[2]; A.ele[7] = ele[5]; A.ele[8] = ele[8];
	return A;
      }

      const Matrix<ThreeByThree, AT> T() const {
	Matrix<ThreeByThree, AT> A(NONINIT);
	A.ele[0] = ele[0]; A.ele[1] = ele[3]; A.ele[2] = ele[6]; A.ele[3] = ele[1]; A.ele[4] = ele[4]; A.ele[5] = ele[7]; A.ele[6] = ele[2]; A.ele[7] = ele[5]; A.ele[8] = ele[8];
	return A;
      }
  };
  // ------------------------- Constructors -------------------------------------
  template <class AT> 
    Matrix<ThreeByThree, AT>::Matrix(const char *strs) {
    // if 'strs' is a single scalar, surround it first with '[' and ']'.
    // This is more Matlab-like, because e.g. '5' and '[5]' is just the same.
    // (This functionallitiy is needed e.g. by MBXMLUtils (OpenMBV,MBSim))
    std::istringstream iss(strs);
    char c;
    iss>>c;
    if(c=='[') iss.str(strs);
    else iss.str(std::string("[")+strs+"]");

    int m = 0,n=0;
    int buf=0;
    iss >> c;
    AT x;
    do {
      iss >> x;
      iss >> c;
      if(c==';') {
	if(buf)
	  assert(buf == n);

	buf=n;
	n=0;
	m++;
      }
      else if(c==',')
	n++;
      c='0';
    } while(iss);

    n++; m++;
    assert(m==3);
    assert(n==3);
    iss.clear();
    iss.seekg(0);
    iss >> c;
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++) {
	iss >> ele[i+j*3];
	iss >> c;
      }
  }
  // ----------------------------------------------------------------------------

   template <class AT> template< class Type>
    Matrix<ThreeByThree, AT>& Matrix<ThreeByThree, AT>::operator<<(const Matrix<Type, AT> &A) { 

      if(A.rows() == 0 || A.cols() == 0)
	return *this;

#ifdef FMATVEC_SIZE_CHECK
      assert(m == A.rows());
      assert(n == A.cols());
#endif

      deepCopy(A);

      return *this;
    }

   template <class AT>
    Matrix<ThreeByThree, AT>& Matrix<ThreeByThree, AT>::operator<<(const Matrix<ThreeByThree, AT> &A) { 

      deepCopy(A);

      return *this;
    }

  template <class AT>
    Matrix<ThreeByThree, AT>&  Matrix<ThreeByThree, AT>::init(const AT& val) {

      ele[0] = val; ele[1] = val; ele[2] = val; ele[3] = val; ele[4] = val; ele[5] = val; ele[6] = val; ele[7] = val; ele[8] = val;

      return *this;
    }

  template <class AT>
    Matrix<ThreeByThree, AT> Matrix<ThreeByThree, AT>::copy() const {

      Matrix<ThreeByThree, AT> A(NONINIT);
      A.deepCopy(*this);

      return A;
    }

  template <class AT> template <class Type>
    void Matrix<ThreeByThree, AT>::deepCopy(const Matrix<Type, AT> &A) { 
      for(int i=0; i<3; i++) 
	for(int j=0; j<3; j++)
          ele[i+j*3] = A.operator()(i,j);
    }

//  template <class AT>
//    Matrix<ThreeByThree, AT>::operator std::vector<std::vector<AT> >() {
//      std::vector<std::vector<AT> > ret(rows());
//      for(int r=0; r<rows(); r++) {
//        ret[r].resize(cols());
//        for(int c=0; c<cols(); c++)
//          ret[r][c]=operator()(r,c);
//      }
//      return ret;
//    }
//
//  template <class AT>
//    Matrix<ThreeByThree, AT>::Matrix(std::vector<std::vector<AT> > m) : memory(m.size()*m[0].size()), ele((AT*)memory.get()), m(m.size()), n(m[0].size()), lda(m.size()), tp(false) {
//#ifndef FMATVEC_NO_INITIALIZATION 
//      init(0);
//#endif
//      for(int r=0; r<rows(); r++) {
//        assert(m[r].size()==cols());
//        for(int c=0; c<cols(); c++)
//          operator()(r,c)=m[r][c];
//      }
//    }

   /// @cond NO_SHOW
   
  template<> template<>
  void Matrix<ThreeByThree, double>::deepCopy(const Matrix<ThreeByThree, double> &A) {
    ele[0] = A.ele[0]; ele[1] = A.ele[1]; ele[2] = A.ele[2]; ele[3] = A.ele[3]; ele[4] = A.ele[4]; ele[5] = A.ele[5]; ele[6] = A.ele[6]; ele[7] = A.ele[7]; ele[8] = A.ele[8];
  }

//  template<> template<>
//  void Matrix<ThreeByThree, double>::deepCopy(const Matrix<Symmetric, double> &A);
//
   /// @endcond

  template <class AT>
  inline Matrix<ThreeByThree, AT> operator+(const Matrix<ThreeByThree, AT> &A1, const Matrix<ThreeByThree, AT> &A2) {
    Matrix<ThreeByThree, AT> A3;
    A3()[0] = A1()[0] + A2()[0];
    A3()[1] = A1()[1] + A2()[1];
    A3()[2] = A1()[2] + A2()[2];
    A3()[3] = A1()[3] + A2()[3];
    A3()[4] = A1()[4] + A2()[4];
    A3()[5] = A1()[5] + A2()[5];
    A3()[6] = A1()[6] + A2()[6];
    A3()[7] = A1()[7] + A2()[7];
    A3()[8] = A1()[8] + A2()[8];
    return A3;
  }

  template <class AT>
  inline Matrix<ThreeByThree, AT> trans(const Matrix<ThreeByThree, AT> &A) {
    Matrix<ThreeByThree, AT> B(NONINIT);
    B()[0] = A()[0]; B()[1] = A()[3]; B()[2] = A()[6]; B()[3] = A()[1]; B()[4] = A()[4]; B()[5] = A()[7]; B()[6] = A()[2]; B()[7] = A()[5]; B()[8] = A()[8];
    return B;
  }

}

#endif
