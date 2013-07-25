#ifndef _FMATVEC_FUNCTION_H_
#define _FMATVEC_FUNCTION_H_

#include <fmatvec.h>
#include <stdexcept>

namespace MBXMLUtils {
  class TiXmlElement;
  class TiXmlNode;
}

namespace fmatvec {



/*! Just a dummy class representing a compile error if used e.g. by
 * a assign operator (=) with a double or Vec.
 */
class ErrorType {};



/*! Defines the size type (dimension) of a argument of a Function object.
 * This general template defines this type to ErrorType (not defined size type).
 * See the specializations for working types.
 */
template<typename T>
struct Size {
  typedef ErrorType type;
};
//! Defines the size type of a double as int.
template<>
struct Size<double> {
  typedef int type;
};
//! Defines the size type of a (column) vector as int.
template<typename Shape>
struct Size<Vector<Shape, double> > {
  typedef int type;
};
//! Defines the size type of a row vector as int.
template<typename Shape>
struct Size<RowVector<Shape, double> > {
  typedef int type;
};
//! Defines the size type of a matrix as a fixed size vector of type int (two integers).
template<typename Type, typename RowShape, typename ColShape>
struct Size<Matrix<Type, RowShape, ColShape, double> > {
  typedef Vector<Fixed<2>, int> type;
};



/*! Defines the resulting type of the derivative of a value of type Dep with respect to a value of type Indep.
 * This general template defines this type to ErrorType (not defined derivative).
 * See the specializations for working types.
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
 * The partial derivative operator \f$ \texttt{parDer}_x(\boldsymbol{R}) \f$ of a
 * rotation matrix \f$ \boldsymbol{R} \f$ with respect to a scalar \f$ x \f$ in defined
 * by the Function class (member functions parDerX) as follows:
 * \f[
 *   \texttt{parDer}_x(\boldsymbol{R}) = \widetilde{\boldsymbol{R}^T \frac{\partial\boldsymbol{R}}{\partial x}}
 * \f]
 * where the "inverse tilde" operator (\f$ \widetilde{o} \f$) transform a skew symmetric matrix to a vector.
 * This class spezialization define the type of \f$ \texttt{parDer}_x \f$.
 */
template<typename DepMatShape>
struct Der<Matrix<Rotation, DepMatShape, DepMatShape, double>, double> {
  typedef Vector<DepMatShape, double> type;
};
/*! Defines the type of the derivative of a rotation matrix with respect to a vector as matrix.
 * The partial derivative operator \f$ \texttt{parDer}_{\boldsymbol{x}}(\boldsymbol{R}) \f$ of a
 * rotation matrix \f$ \boldsymbol{R} \f$ with respect to a vector \f$ \boldsymbol{x} \f$ in defined
 * by the Function class (member functions parDerX) as follows:
 * \f[
 *   \texttt{parDer}_{\boldsymbol{x}}(\boldsymbol{R}) = \left[
 *     \widetilde{\boldsymbol{R}^T \frac{\partial\boldsymbol{R}}{\partial x_1}},
 *     \widetilde{\boldsymbol{R}^T \frac{\partial\boldsymbol{R}}{\partial x_1}},
 *     \dots
 *   \right]
 * \f]
 * where the "inverse tilde" operator (\f$ \widetilde{o} \f$) transform a skew symmetric matrix to a vector.
 * This class spezialization define the type of \f$ \texttt{parDer}_{\boldsymbol{x}} \f$.
 */
template<typename DepMatShape, typename IndepVecShape>
struct Der<Matrix<Rotation, DepMatShape, DepMatShape, double>, Vector<IndepVecShape, double> > {
  typedef Matrix<General, DepMatShape, IndepVecShape, double> type;
};




/*! Just a base class for all template based Function classes.
 * (required to have a common base class e.g. for object factories)
 */
class FunctionBase {
  public:
    virtual void initializeUsingXML(MBXMLUtils::TiXmlElement *element) { }
    virtual MBXMLUtils::TiXmlElement *writeXMLFile(MBXMLUtils::TiXmlNode *parent) { return NULL; }
};



/*! A function object of arbitary type (defined like in boost::function).
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
class Function<Ret(Arg)> : public FunctionBase {

  public:

    //! Return the size of first argument
    virtual typename Size<Arg>::type getArgSize() const {
      throw std::runtime_error("getArgSize must be overloaded by derived class.");
    }

    //! Function value: pure virtual (MUST be implemented by derived class)
    virtual Ret operator()(const Arg &arg)=0;

    //! First derivative: partial derivative of the function value with respect to the first argument.
    virtual typename Der<Ret, Arg>::type parDer(const Arg &arg) {
      throw std::runtime_error("parDer must be overloaded by derived class.");
    }
    //! First derivative: directional derivative of the function value with respect to the first argument.
    virtual Ret dirDer(const Arg &argDir, const Arg &arg) {
      throw std::runtime_error("dirDer must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg>::type, Arg>::type parDerParDer(const Arg &arg) {
      throw std::runtime_error("parDerParDer must be overloaded by derived class.");
    }
    //! Second derivative: directional derivative of parDer with respect to the first argument.
    virtual typename Der<Ret, Arg>::type parDerDirDer(const Arg &argDir, const Arg &arg) {
      throw std::runtime_error("parDerDirDer must be overloaded by derived class.");
    }

};



//! A function object with 2 arguments
template<typename Ret, typename Arg1, typename Arg2>
class Function<Ret(Arg1, Arg2)> : public FunctionBase {

  public:

    //! Return the size of first argument
    virtual typename Size<Arg1>::type getArg1Size() const {
      throw std::runtime_error("getArg1Size must be overloaded by derived class.");
    }
    //! Return the size of first argument
    virtual typename Size<Arg2>::type getArg2Size() const {
      throw std::runtime_error("getArg2Size must be overloaded by derived class.");
    }

    //! Function value: pure virtual (MUST be implemented by derived class)
    virtual Ret operator()(const Arg1 &arg1, const Arg2 &arg2)=0;

    //! First derivative: partial derivative of the function value with respect to the first argument.
    virtual typename Der<Ret, Arg1>::type parDer1(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1 must be overloaded by derived class.");
    }
    //! First derivative: directional derivative of the function value with respect to the first argument.
    virtual Ret dirDer1(const Arg1 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer1 must be overloaded by derived class.");
    }
    //! First derivative: partial derivative of the function value with respect to the second argument.
    virtual typename Der<Ret, Arg2>::type parDer2(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2 must be overloaded by derived class.");
    }
    //! First derivative: directional derivative of the function value with respect to the second argument.
    virtual Ret dirDer2(const Arg2 &arg2Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("dirDer2 must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer1 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg1>::type, Arg1>::type parDer1ParDer1(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1ParDer1 must be overloaded by derived class.");
    }
    //! Second derivative: directional derivative of parDer1 with respect to the first argument.
    virtual typename Der<Ret, Arg1>::type parDer1DirDer1(const Arg1 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1DirDer1 must be overloaded by derived class.");
    }
    //! Second derivative: partial derivative of parDer2 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg2>::type, Arg2>::type parDer2ParDer2(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2ParDer2 must be overloaded by derived class.");
    }
    //! Second derivative: directional derivative of parDer2 with respect to the first argument.
    virtual typename Der<Ret, Arg2>::type parDer2DirDer2(const Arg2 &arg2Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2DirDer2 must be overloaded by derived class.");
    }

    //! Second mixed derivative: partial derivative of parDer1 with respect to the second argument.
    virtual typename Der<typename Der<Ret, Arg1>::type, Arg2>::type parDer1ParDer2(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1ParDer2 must be overloaded by derived class.");
    }
    //! Second mixed derivative: directional derivative of parDer1 with respect to the second argument.
    virtual typename Der<Ret, Arg1>::type parDer1DirDer2(const Arg2 &arg2Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer1DirDer2 must be overloaded by derived class.");
    }
    //! Second mixed derivative: partial derivative of parDer2 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg2>::type, Arg1>::type parDer2ParDer1(const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2ParDer1 must be overloaded by derived class.");
    }
    //! Second mixed derivative: directional derivative of parDer2 with respect to the first argument.
    virtual typename Der<Ret, Arg2>::type parDer2DirDer1(const Arg1 &arg1Dir, const Arg1 &arg1, const Arg2 &arg2) {
      throw std::runtime_error("parDer2DirDer1 must be overloaded by derived class.");
    }

};



// A function object with 3, 4, 5, ... arguments
// Not required till now!

}

#endif
