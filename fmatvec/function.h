#ifndef _FMATVEC_FUNCTION_H_
#define _FMATVEC_FUNCTION_H_

#include <fmatvec/fmatvec.h>
#include <fmatvec/atom.h>
#include <stdexcept>

namespace fmatvec {

/*! Just a dummy class representing a compile error if used e.g. by
 * a assign operator (=) with a double or Vec.
 */
class ErrorType {
  ErrorType() = delete;
};

/*! Defines the static size (dimension) of the template argument T.
 * Maximal two dimensions are supported.
 * See the specializations for types which have a static size.
 */
template<typename T>
struct StaticSize;

//! Defines the static size (dimension) of a double.
template<>
struct StaticSize<int> {
  enum { size1=1, size2=1 };
};

//! Defines the static size (dimension) of a double.
template<>
struct StaticSize<double> {
  enum { size1=1, size2=1 };
};

//! Defines the static size (dimension) of a vector.
template<typename Shape, typename AT>
struct StaticSize<Vector<Shape, AT> > {
  enum { size1=0, size2=1 };
};

//! Defines the static size (dimension) of a row vector.
template<typename Shape, typename AT>
struct StaticSize<RowVector<Shape, AT> > {
  enum { size1=1, size2=0 };
};

//! Defines the static size (dimension) of a fixed vector.
template<int N, typename AT>
struct StaticSize<Vector<Fixed<N>, AT> > {
  enum { size1=N, size2=1 };
};

//! Defines the static size (dimension) of a fixed rowvector.
template<int N, typename AT>
struct StaticSize<RowVector<Fixed<N>, AT> > {
  enum { size1=1, size2=N };
};

//! Defines the static size (dimension) of a variable matrix.
template<typename Storage, typename RowShape, typename ColShape, typename AT>
struct StaticSize<Matrix<Storage, RowShape, ColShape, AT> > {
  enum { size1=0, size2=0 };
};

//! Defines the static size (dimension) of a fixed matrix.
template<typename Storage, int N, int M, typename AT>
struct StaticSize<Matrix<Storage, Fixed<N>, Fixed<M>, AT> > {
  enum { size1=N, size2=M };
};

//! Defines the static size (dimension) of a partial fixed matrix.
template<typename Storage, int N, typename ColShape, typename AT>
struct StaticSize<Matrix<Storage, Fixed<N>, ColShape, AT> > {
  enum { size1=N, size2=0 };
};

//! Defines the static size (dimension) of a partial fixed matrix.
template<typename Storage, typename RowShape, int M, typename AT>
struct StaticSize<Matrix<Storage, RowShape, Fixed<M>, AT> > {
  enum { size1=0, size2=M };
};

//! Defines the static size (dimension) of a fixed square matrix.
template<int N, typename AT>
struct StaticSize<SquareMatrix<Fixed<N>, AT> > {
  enum { size1=N, size2=N };
};

//! Defines the static size (dimension) of a variable square matrix.
template<typename Shape, typename AT>
struct StaticSize<SquareMatrix<Shape, AT> > {
  enum { size1=0, size2=0 };
};

/*! Defines the resulting type of the derivative of a value of type Dep with respect to a value of type Indep.
 * See the specializations.
 */
template<typename Dep, typename Indep>
struct Der {
  typedef ErrorType type;
};

/*! Defines the type of the derivative of a value of type Dep with respect to a scalar as type Dep.
 * The partial derivative operator with respect to a scalar is define like in common mathematics.
 */
template<typename Dep>
struct Der<Dep, double> {
  typedef Dep type;
};

/*! Defines the type of the derivative of a scalar with respect to a vector as row vector.
 * The partial derivative operator with respect to a vector is define like in common mathematics:
 * a column vector where each entry corresponds to the partial derivative with respect the
 * corresponding entry of the independent vector.
 */
template<typename IndepVecShape>
struct Der<double, Vector<IndepVecShape, double> > {
  typedef RowVector<IndepVecShape, double> type;
};

/*! Defines the type of the derivative of a vector with respect to a vector as matrix.
 * The partial derivative operator with respect to a vector is define like in common mathematics:
 * a column vector where each entry corresponds to the partial derivative with respect the
 * corresponding entry of the independent vector.
 */
template<typename DepVecShape, typename IndepVecShape>
struct Der<Vector<DepVecShape, double>, Vector<IndepVecShape, double> > {
  typedef Matrix<General, DepVecShape, IndepVecShape, double> type;
};

/*! Defines the type of the derivative of a rotation matrix with respect to a scalar as vector.
 * The partial derivative operator \f$ \texttt{parDer}_x(\boldsymbol{R}_{12}) \f$ of a
 * rotation matrix \f$ \boldsymbol{R}_{12} \f$ with respect to a scalar \f$ x \f$ in defined
 * by the Function class (member functions parDerX) as follows:
 * \f[
 *   \texttt{parDer}_x(\boldsymbol{R}_{12}) = \widetilde{\left(\frac{\partial\boldsymbol{R}_{12}}{\partial x}\boldsymbol{R}^T_{12}\right)} \in \mathbb{R}^{3\times 1}
 * \f]
 * where the "inverse tilde" operator (\f$ \widetilde{o} \f$) transform a skew symmetric matrix to a vector.
 * This class spezialization define the type of \f$ \texttt{parDer}_x \f$.
 */
template<>
struct Der<Matrix<Rotation, Fixed<3>, Fixed<3>, double>, double> {
  typedef Vector<Fixed<3>, double> type;
};

/*! Defines the type of the derivative of a rotation matrix with respect to a vector as matrix.
 * The partial derivative operator \f$ \texttt{parDer}_{\boldsymbol{x}}(\boldsymbol{R}_{12}) \f$ of a
 * rotation matrix \f$ \boldsymbol{R}_{12} \f$ with respect to a vector \f$ \boldsymbol{x} \f$ in defined
 * by the Function class (member functions parDerX) as follows:
 * \f[
 *   \texttt{parDer}_{\boldsymbol{x}}(\boldsymbol{R}_{12}) = \left[
 *     \widetilde{\left(\frac{\partial\boldsymbol{R}_{12}}{\partial x_1}\boldsymbol{R}^T_{12}\right)},
 *     \widetilde{\left(\frac{\partial\boldsymbol{R}_{12}}{\partial x_1}\boldsymbol{R}^T_{12}\right)},
 *     \dots
 *   \right] \in \mathbb{R}^{3\times N} \quad \text{when}\quad \boldsymbol{x} \in \mathbb{R}^{N\times 1}
 * \f]
 * where the "inverse tilde" operator (\f$ \widetilde{o} \f$) transform a skew symmetric matrix to a vector.
 * This class spezialization define the type of \f$ \texttt{parDer}_{\boldsymbol{x}} \f$.
 */
template<typename IndepVecShape>
struct Der<Matrix<Rotation, Fixed<3>, Fixed<3>, double>, Vector<IndepVecShape, double> > {
  typedef Matrix<General, Fixed<3>, IndepVecShape, double> type;
};

/*! A function object of arbitary type (defined like in std::function).
 * The number of arguments is variable and always one value is returned.
 * The type of the arguments and the return value is also variable using templates.
 * The function value is always provided and the partial and directional derivatives with
 * respect to the arguments can be proveided if implemented by a derived class. This
 * is normally done up to derivative order 2. The exact defintion of the partial
 * derivatives depend on the type of the function value and the type of the parameters
 * used for the derivative. See the Der class template specializations for a detailed
 * description of the defintion of the partial derivative dependent on the types.
 */
template<typename Sig>
class Function;

//! A function object with 1 argument
template<typename Ret, typename Arg>
class Function<Ret(Arg)> : virtual public Atom {

  public:

    using DRetDArg = typename Der<Ret, Arg>::type;
    using DDRetDDArg = typename Der<DRetDArg, Arg>::type;

    //! Compile time size of the return value: =0 == unknown compile time size
    enum { retSize1 = StaticSize<Ret>::size1, retSize2 = StaticSize<Ret>::size2 };

    //! Compile time size of the argument: =1 == scalar; >1 == vector; =0 == unknown compile time vector size
    static constexpr int argSize = StaticSize<Arg>::size1;

    //! Return the size of the return value: =0 == unknown size
    virtual std::pair<int, int> getRetSize() const { return std::make_pair(retSize1, retSize2); }

    //! Return the size of the argument: =1 == scalar; >1 == vector; =0 == unknown vector size
    virtual int getArgSize() const { return argSize; }

    //! Function value: pure virtual (MUST be implemented by derived class)
    virtual Ret operator()(const Arg &arg)=0;

    //! First derivative: partial derivative of the function value with respect to the argument.
    virtual DRetDArg parDer(const Arg &arg) {
      throw std::runtime_error("parDer must be overloaded by derived class.");
    }

    //! First derivative: directional derivative of the function value with respect to the argument.
    virtual Ret dirDer(const Arg &argDir, const Arg &arg) {
      throw std::runtime_error("dirDer must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer with respect to the argument.
    virtual DDRetDDArg parDerParDer(const Arg &arg) {
      throw std::runtime_error("parDerParDer must be overloaded by derived class.");
    }

    //! Second derivative: directional derivative of parDer with respect to the argument.
    virtual DRetDArg parDerDirDer(const Arg &argDir, const Arg &arg) {
      throw std::runtime_error("parDerDirDer must be overloaded by derived class.");
    }

    //! Second derivative: directional derivative of dirDer with respect to the argument.
    virtual Ret dirDerDirDer(const Arg &argDir_1, const Arg &argDir_2, const Arg &arg) {
      throw std::runtime_error("dirDerDirDer must be overloaded by derived class.");
    }

    //! Returns true, if the partial derivative of the function value with respect to the argument 
    //  is constant.
    virtual bool constParDer() const { return false; }
};

//! A function object with 2 arguments
template<typename Ret, typename Arg1, typename Arg2>
class Function<Ret(Arg1, Arg2)> : virtual public Atom {

  public:

    using DRetDArg1 = typename Der<Ret, Arg1>::type;
    using DRetDArg2 = typename Der<Ret, Arg2>::type;
    using DDRetDDArg1 = typename Der<DRetDArg1, Arg1>::type;
    using DDRetDDArg2 = typename Der<DRetDArg2, Arg2>::type;
    using DDRetDArg1DArg2 = typename Der<DRetDArg1, Arg2>::type;

    //! Compile time size of the return value: =0 == unknown compile time size
    enum { retSize1 = StaticSize<Ret>::size1, retSize2 = StaticSize<Ret>::size2 };

    //! Compile time size of the first argument: =1 == scalar; >1 == vector; =0 == unknown compile time vector size
    static constexpr int arg1Size = StaticSize<Arg1>::size1;

    //! Compile time size of the second argument: =1 == scalar; >1 == vector; =0 == unknown compile time vector size
    static constexpr int arg2Size = StaticSize<Arg2>::size1;

    //! Return the size of the return value: =0 == unknown size
    virtual std::pair<int, int> getRetSize() const { return std::make_pair(retSize1, retSize2); }

    //! Return the size of the first argument: =1 == scalar; >1 == vector; =0 == unknown vector size
    virtual int getArg1Size() const { return arg1Size; }
    //! Return the size of the second argument: =1 == scalar; >1 == vector; =0 == unknown vector size
    virtual int getArg2Size() const { return arg2Size; }

    //! Function value: pure virtual (MUST be implemented by derived class)
    virtual Ret operator()(const Arg1 &arg1, const Arg2 &arg2)=0;

    //! First derivative: partial derivative of the function value with respect to the first argument.
    virtual DRetDArg1 parDer1(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1 must be overloaded by derived class.");
    }

    //! First derivative: directional derivative of the function value with respect to the first argument.
    virtual Ret dirDer1(const Arg1 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer1 must be overloaded by derived class.");
    }

    //! First derivative: partial derivative of the function value with respect to the second argument.
    virtual DRetDArg2 parDer2(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2 must be overloaded by derived class.");
    }

    //! First derivative: directional derivative of the function value with respect to the second argument.
    virtual Ret dirDer2(const Arg2 &arg2Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer2 must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer1 with respect to the first argument.
    virtual DDRetDDArg1 parDer1ParDer1(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1ParDer1 must be overloaded by derived class.");
    }

    //! Second derivative: directional derivative of parDer1 with respect to the first argument.
    virtual DRetDArg1 parDer1DirDer1(const Arg1 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1DirDer1 must be overloaded by derived class.");
    }

    //! Second derivative: directional derivative of dirDer1 with respect to the first argument.
    virtual Ret dirDer1DirDer1(const Arg1 &arg1Dir_1, const Arg1 &arg1Dir_2, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer1DirDer1 must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer2 with respect to the first argument.
    virtual DDRetDDArg2 parDer2ParDer2(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2ParDer2 must be overloaded by derived class.");
    }

    //! Second derivative: directional derivative of parDer2 with respect to the first argument.
    virtual DRetDArg2 parDer2DirDer2(const Arg2 &arg2Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2DirDer2 must be overloaded by derived class.");
    }

    //! Second derivative: directional derivative of dirDer2 with respect to the first argument.
    virtual Ret dirDer2DirDer2(const Arg2 &arg2Dir_1, const Arg2 &arg2Dir_2, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer2DirDer2 must be overloaded by derived class.");
    }

    //! Second mixed derivative: partial derivative of parDer1 with respect to the second argument.
    virtual DDRetDArg1DArg2 parDer1ParDer2(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1ParDer2 must be overloaded by derived class.");
    }

    //! Second mixed derivative: directional derivative of parDer1 with respect to the second argument.
    virtual DRetDArg1 parDer1DirDer2(const Arg2 &arg2Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1DirDer2 must be overloaded by derived class.");
    }

    //! Second mixed derivative: directional derivative of dirDer2 with respect to the first argument.
    virtual Ret dirDer2DirDer1(const Arg2 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer2DirDer1 must be overloaded by derived class.");
    }

    //! Second mixed derivative: partial derivative of dirDer1 with respect to the second argument.
    virtual DRetDArg2 parDer2DirDer1(const Arg1 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2DirDer1 must be overloaded by derived class.");
    }

    //! Returns true, if the partial derivative of the function value with respect to the first argument 
    //  is constant.
    virtual bool constParDer1() const { return false; }

    //! Returns true, if the partial derivative of the function value with respect to the second argument 
    //  is constant.
    virtual bool constParDer2() const { return false; }
};

// A function object with 3, 4, 5, ... arguments
// Not required till now!

}

#endif
