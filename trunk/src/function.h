#ifndef _FMATVEC_FUNCTION_H_
#define _FMATVEC_FUNCTION_H_

#include <fmatvec.h>
#include <stdexcept>

namespace fmatvec {



// Just a dummy class representing a compile error if used e.g. by
// a assign operator (=) with a double or Vec.
class ErrorType {};



/** Defines the size type (dimension) of a argument.
 * Scalar does not define this.
 * Vector has one int.
 * Matrix has two int (fmatvec int vector)
 * This general template defines a ErrorType (not defined size type, uses for all other than Vec and Mat).
 * See the specializations for working types.
 */
template<typename T>
struct Size {
  typedef ErrorType type;
};
// specialization for vector: return one int
template<typename Shape>
struct Size<Vector<Shape, double> > {
  typedef int type;
};
// specialization for row vector: return one int
template<typename Shape>
struct Size<RowVector<Shape, double> > {
  typedef int type;
};
// specialization for matrix: return two int (fmatvec int vector)
template<typename Type, typename RowShape, typename ColShape>
struct Size<Matrix<Type, RowShape, ColShape, double> > {
  typedef Vector<Fixed<2>, int> type;
};



/** Create a type, equal the type of the derivative of a type "Dep" with respect to type "Indep".
 * e.g.: Der<Vec, double>::type expands to Vec
 * e.g.: Der<Mat, double>::type expands to Mat
 * e.g.: Der<Vec3, Vec>::type expands to Mat3xV
 * e.g.: Der<VecV, VecV>::type expands to MatVxV
 * This general template defines a ErrorType (not working derivative).
 * See the specializations for working types.
 */
template<typename Dep, typename Indep>
struct Der {
  typedef ErrorType type;
};
// specialization, if the independent type is a scalar: derivative type = Dep
template<typename Dep>
struct Der<Dep, double> {
  typedef Dep type;
};
// specialization, if the depdentent type is a scalar and the independent type is a vector: derivative type = RowVector
template<typename IndepVecShape>
struct Der<double, Vector<IndepVecShape, double> > {
  typedef RowVector<IndepVecShape, double> type;
};
// specialization, if the dependent and independent type is a vector: derivative type = Matrix
template<typename DepVecShape, typename IndepVecShape>
struct Der<Vector<DepVecShape, double>, Vector<IndepVecShape, double> > {
  typedef Matrix<General, DepVecShape, IndepVecShape, double> type;
};
// specialization, if the dependent type is a rotation matrix and the independent type is a vector
template<typename DepMatShape, typename IndepVecShape>
struct Der<Matrix<Rotation, DepMatShape, DepMatShape, double>, Vector<IndepVecShape, double> > {
  typedef Matrix<General, DepMatShape, IndepVecShape, double> type;
};
// specialization, if the dependent type is a rotation matrix and the independent type is a double
template<typename DepMatShape>
struct Der<Matrix<Rotation, DepMatShape, DepMatShape, double>, double> {
  typedef Vector<DepMatShape, double> type;
};




// Just a base class for all template based Function classes
// (required to have a common base class e.g. for object factories)
class FunctionBase {};



/*! A function object of arbitary type (defined like in boost::function).
 * The number of arguments is variable and always one value is returned.
 * The type of the arguments and the return value is also variable using templates.
 */
template<typename Sig>
class Function;



//! A function object with 1 argument
template<typename Ret, typename Arg1>
class Function<Ret(Arg1)> : public FunctionBase {

  public:

    //! Return the size of first argument
    virtual typename Size<Arg1>::type getArg1Size() {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

    //! Function value: pure virtual (MUST be implemented by derived class)
    virtual Ret operator()(Arg1 arg1)=0;

    //! First derivative: partial derivative of the function value with respect to the first argument.
    virtual typename Der<Ret, Arg1>::type parDer1(Arg1 arg1) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! First derivative: directional derivative of the function value with respect to the first argument.
    virtual Ret dirDer1(Arg1 arg1Dir, Arg1 arg1) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer1 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg1>::type, Arg1>::type parDer1ParDer1(Arg1 arg1) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second derivative: directional derivative of parDer1 with respect to the first argument.
    virtual typename Der<Ret, Arg1>::type parDer1DirDer1(Arg1 arg1Dir, Arg1 arg1) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

};



//! A function object with 2 arguments
template<typename Ret, typename Arg1, typename Arg2>
class Function<Ret(Arg1, Arg2)> : public FunctionBase {

  public:

    //! Return the size of first argument
    virtual typename Size<Arg1>::type getArg1Size() {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Return the size of first argument
    virtual typename Size<Arg2>::type getArg2Size() {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

    //! Function value: pure virtual (MUST be implemented by derived class)
    virtual Ret operator()(Arg1 arg1, Arg2 arg2)=0;

    //! First derivative: partial derivative of the function value with respect to the first argument.
    virtual typename Der<Ret, Arg1>::type parDer1(Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! First derivative: directional derivative of the function value with respect to the first argument.
    virtual Ret dirDer1(Arg1 arg1Dir, Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! First derivative: partial derivative of the function value with respect to the second argument.
    virtual typename Der<Ret, Arg2>::type parDer2(Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! First derivative: directional derivative of the function value with respect to the second argument.
    virtual Ret dirDer2(Arg2 arg2Dir, Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

    //! Second derivative: partial derivative of parDer1 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg1>::type, Arg1>::type parDer1ParDer1(Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second derivative: directional derivative of parDer1 with respect to the first argument.
    virtual typename Der<Ret, Arg1>::type parDer1DirDer1(Arg1 arg1Dir, Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second derivative: partial derivative of parDer2 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg2>::type, Arg2>::type parDer2ParDer2(Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second derivative: directional derivative of parDer2 with respect to the first argument.
    virtual typename Der<Ret, Arg2>::type parDer2DirDer2(Arg2 arg2Dir, Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

    //! Second mixed derivative: partial derivative of parDer1 with respect to the second argument.
    virtual typename Der<typename Der<Ret, Arg1>::type, Arg2>::type parDer1ParDer2(Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second mixed derivative: directional derivative of parDer1 with respect to the second argument.
    virtual typename Der<Ret, Arg1>::type parDer1DirDer2(Arg2 arg2Dir, Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second mixed derivative: partial derivative of parDer2 with respect to the first argument.
    virtual typename Der<typename Der<Ret, Arg2>::type, Arg1>::type parDer2ParDer1(Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }
    //! Second mixed derivative: directional derivative of parDer2 with respect to the first argument.
    virtual typename Der<Ret, Arg2>::type parDer2DirDer1(Arg1 arg1Dir, Arg1 arg1, Arg2 arg2) {
      throw std::runtime_error("Must be overloaded by derived class.");
    }

};



// A function object with 3, 4, 5, ... arguments
// Not required till now!

}

#endif
