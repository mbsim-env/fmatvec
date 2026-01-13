#ifndef _FMATVEC_AST_H_
#define _FMATVEC_AST_H_

#include "fmatvec/function.h"
#include <map>
#include <array>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/container/small_vector.hpp>
#include <utility>
#include <fmatvec/stream.h>

namespace fmatvec {

// forward declarations

class IndependentVariable;

namespace AST {
  class Vertex;
  class Symbol;
  class Operation;
  class NativeFunction;
  template<class T> class Constant;
  FMATVEC_EXPORT SymbolicExpression substScalar(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);
}

template<class... Arg> class Eval;
class EvalHelper;

//! A symbolic expression.
//! This class represent an arbitary symbolic expression.
//! This can either be the special case of just a symbol
//! or an arbitary expression consisting of a hierarchiy
//! of operations consisting itself of expressions, symbols or constants.
class SymbolicExpression : public std::shared_ptr<const AST::Vertex> {
  template<class... Arg> friend class Eval;
  friend class EvalHelper;
  friend class AST::Operation;
  friend class AST::Constant<long>;
  friend class AST::Constant<double>;
  friend class AST::NativeFunction;
  friend FMATVEC_EXPORT SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep);
  friend SymbolicExpression AST::substScalar(const SymbolicExpression &se,
                                             const IndependentVariable& a, const SymbolicExpression &b);
  protected:
    template<class T> SymbolicExpression(const shared_ptr<T> &x);
#ifndef SWIG
    static const struct ConstructSymbol{} constructSymbol; // just used for tag dispatching
#endif
    SymbolicExpression(ConstructSymbol);

    // do not allow to call these functions inherited from std::shared_ptr
    using std::shared_ptr<const AST::Vertex>::get;
    using std::shared_ptr<const AST::Vertex>::operator*;
    using std::shared_ptr<const AST::Vertex>::operator->;
  public:
    //! Creates the value 0.
    FMATVEC_EXPORT SymbolicExpression();
    //! Creates a expression form the specified string (the string must be a serialized SymbolicExpression).
    FMATVEC_EXPORT SymbolicExpression(const std::string &str);
    //! Creates a integer constant.
    FMATVEC_EXPORT SymbolicExpression(int x);
    FMATVEC_EXPORT SymbolicExpression(long x);
    //! Creates a double constant.
    FMATVEC_EXPORT SymbolicExpression(double x);

    //! Garbage collect everything.
    //! The only "garbage" which can be left by this class are empty weak_ptr's stored in global maps.
    //! This routine removes such garbage.
    FMATVEC_EXPORT static void garbageCollect();

    // default ctors and assignment operators
    FMATVEC_EXPORT SymbolicExpression(const SymbolicExpression& x) = default;
    FMATVEC_EXPORT SymbolicExpression(SymbolicExpression&& x) = default;
    FMATVEC_EXPORT SymbolicExpression& operator=(const SymbolicExpression& x) = default;
    FMATVEC_EXPORT SymbolicExpression& operator=(SymbolicExpression&& x) = default;
    FMATVEC_EXPORT ~SymbolicExpression() = default;

    // operators
#ifndef SWIG
    FMATVEC_EXPORT SymbolicExpression operator+(const SymbolicExpression &b) const;
    FMATVEC_EXPORT SymbolicExpression operator-(const SymbolicExpression &b) const;
    FMATVEC_EXPORT SymbolicExpression operator*(const SymbolicExpression &b) const;
    FMATVEC_EXPORT SymbolicExpression operator/(const SymbolicExpression &b) const;
#endif
    FMATVEC_EXPORT SymbolicExpression operator-() const;
    FMATVEC_EXPORT SymbolicExpression& operator+=(const SymbolicExpression &b);
    FMATVEC_EXPORT SymbolicExpression& operator-=(const SymbolicExpression &b);
    FMATVEC_EXPORT SymbolicExpression& operator*=(const SymbolicExpression &b);
    FMATVEC_EXPORT SymbolicExpression& operator/=(const SymbolicExpression &b);

    // We implement operator<<= to be able to use SymbolicExpression fully equivalent to the vector/matrix classes
    // which have a assign method which redimension before assigning. Redimensioning is not required for a scalar
    // type like assign
    FMATVEC_EXPORT SymbolicExpression& operator<<=(const SymbolicExpression &src);

#if defined(FMATVEC_DEBUG) && !defined(SWIG)
    FMATVEC_EXPORT static signed long evalOperationsCount;
#endif
};

//! A independent variable.
//! Any SymbolicExpression can be partialliy differentiated with respect to a independent variable.
//! An independent varible can also be assigned a value which is used if eval is called.
class IndependentVariable : public SymbolicExpression {
  friend class AST::Symbol;
  friend FMATVEC_EXPORT std::istream& operator>>(std::istream& s, IndependentVariable &v);
  public:
    //! Creates a IndependentVariable variable (each call to this ctor creates a new independent variable)
    FMATVEC_EXPORT IndependentVariable();
    //! Creates a IndependentVariable variable from the specified string (the string is a serialized IndependentVariable).
    FMATVEC_EXPORT IndependentVariable(const std::string &str);

#ifndef SWIG
    //! Set the double value of the independent value.
    inline IndependentVariable& operator^=(double x);
#endif

  private:
    IndependentVariable(const shared_ptr<const AST::Symbol> &x);
};

// define the operator results of SymbolicExpression and IndependentVariable with other AT types.
FMATVEC_OPERATORRESULT1(SymbolicExpression, SymbolicExpression)

FMATVEC_OPERATORRESULT2(SymbolicExpression, int, SymbolicExpression)
FMATVEC_OPERATORRESULT2(SymbolicExpression, long, SymbolicExpression)
FMATVEC_OPERATORRESULT2(SymbolicExpression, double, SymbolicExpression)

FMATVEC_OPERATORRESULT1(IndependentVariable, SymbolicExpression)

FMATVEC_OPERATORRESULT2(IndependentVariable, int, SymbolicExpression)
FMATVEC_OPERATORRESULT2(IndependentVariable, long, SymbolicExpression)
FMATVEC_OPERATORRESULT2(IndependentVariable, double, SymbolicExpression)
FMATVEC_OPERATORRESULT2(IndependentVariable, SymbolicExpression, SymbolicExpression)

// Some member function definition of SymbolicExpression/IndependentVariable are moved to the end of this file
// since they need the defintion of the other class defined in this file.

FMATVEC_EXPORT SymbolicExpression operator+(double a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator-(double a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator*(double a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator/(double a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator+(int a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator-(int a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator*(int a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator/(int a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator+(long a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator-(long a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator*(long a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression operator/(long a, const SymbolicExpression &b);

//! Generate a new SymbolicExpression being the partial derivate of dep
//! with respect to indep (indep must be a symbol).
FMATVEC_EXPORT SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep);

//! Write a SymbolicExpression to a stream using serialization.
FMATVEC_EXPORT std::ostream& operator<<(std::ostream& s, const SymbolicExpression& se);
//! Create/initialize a SymbolicExpression from a stream using deserialization.
FMATVEC_EXPORT std::istream& operator>>(std::istream& s, SymbolicExpression &se);
//! Create/initialize a IndependentVariable from a stream using deserialization.
FMATVEC_EXPORT std::istream& operator>>(std::istream& s, IndependentVariable &v);

// function operations overloaded for SymbolicExpression
FMATVEC_EXPORT SymbolicExpression pow(const SymbolicExpression &a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression log(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression sqrt(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression sin(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression cos(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression tan(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression sinh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression cosh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression tanh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression asin(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression acos(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression atan(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression atan2(const SymbolicExpression &y, const SymbolicExpression &x);
FMATVEC_EXPORT SymbolicExpression asinh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression acosh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression atanh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression exp(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression sign(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression heaviside(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression abs(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression min(const SymbolicExpression &a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression max(const SymbolicExpression &a, const SymbolicExpression &b);
FMATVEC_EXPORT SymbolicExpression condition(const SymbolicExpression &c, const SymbolicExpression &gt, const SymbolicExpression &lt);

#ifndef SWIG
namespace AST { // internal namespace

/* ***** Struct for a "bytecode" for fast runtime evaluation *****
 * This class uses very low-level technices like raw-pointers to gain maximal performance 
 * during runtime evaluation.
 *
 * This class store a kind of "bytecode" which can be used to evaluate an mathematical expression very fast.
 * This class is not copyable and not moveable since it uses intern pointers between instances of this class.
 * For optimal performance this class should be used in continous memory e.g. as std::vector<ByteCode>.
 * For the same reason, to keep all data near together for best processor memory cache usage, all members
 * of these class also to not use allocated memory: std::array instead of std::vector is used.
 * The usage of std::vector<ByteCode> requires that the move-ctor must be provided at compile time
 * (it cannot be delete'd) but the move-ctor throws at runtime when called to ensure its not used.
 * Hence std::vector<ByteCode> can only be used without reallocation: std::vector::reserve must be called before use.
 * To add an element only std::vector::emplace_back() can be used.
*/
struct ByteCode {
  static constexpr size_t N { 12 };
  FMATVEC_EXPORT ByteCode(size_t n);
  ByteCode(const ByteCode &) = delete;
  FMATVEC_EXPORT ByteCode(ByteCode &&) noexcept;
  ByteCode &operator=(const ByteCode &) = delete;
  ByteCode &operator=(ByteCode &&) = delete;
  FMATVEC_EXPORT ~ByteCode() = default;
  using Arg = boost::container::small_vector<double*,N>;
  std::function<void(double*, const Arg&)> func; // the operation this byteCode entry executes
  double  retValue; // storage of the return value of the operation: retPtr may point to this value
  double* retPtr; // a pointer to which the operation can write its result to
  Arg argsPtr; // pointers from which the operation reads its arguments
};

template<class Func, class ArgS>
struct SymbolicFuncWrapArg1 {
  static SymbolicExpression call(
    const std::shared_ptr<Function<Func>> &func,
    const ArgS &arg);
};

template<>
struct SymbolicFuncWrapArg1<double(double), SymbolicExpression> {
  static SymbolicExpression call(
    const std::shared_ptr<Function<double(double)>> &func,
    const SymbolicExpression &arg);
};

template<class Type>
struct SymbolicFuncWrapArg1<double(Vector<Type, double>), Vector<Type, SymbolicExpression>> {
  static SymbolicExpression call(
    const std::shared_ptr<Function<double(Vector<Type, double>)>> &func,
    const Vector<Type, SymbolicExpression> &arg);
};
  
template<class Func, class Arg1S, class Arg2S>
struct SymbolicFuncWrapArg2 {
  static SymbolicExpression call(
    const std::shared_ptr<Function<Func>> &func,
    const Arg1S &arg1, const Arg2S &arg2);
};

template<>
struct SymbolicFuncWrapArg2<double(double,double), SymbolicExpression, SymbolicExpression> {
  static SymbolicExpression call(
    const std::shared_ptr<Function<double(double,double)>> &func,
    const SymbolicExpression &arg1, const SymbolicExpression &arg2);
};

template<class Type>
struct SymbolicFuncWrapArg2<double(Vector<Type,double>,double), Vector<Type, SymbolicExpression>, SymbolicExpression> {
  static SymbolicExpression call(
    const std::shared_ptr<Function<double(Vector<Type,double>,double)>> &func,
    const Vector<Type, SymbolicExpression> &arg1, const SymbolicExpression &arg2);
};

template<class Type>
struct SymbolicFuncWrapArg2<double(double,Vector<Type,double>), SymbolicExpression, Vector<Type, SymbolicExpression>> {
  static SymbolicExpression call(
    const std::shared_ptr<Function<double(double,Vector<Type,double>)>> &func,
    const SymbolicExpression &arg1, const Vector<Type, SymbolicExpression> &arg2);
};

template<class Type1, class Type2>
struct SymbolicFuncWrapArg2<double(Vector<Type1,double>,Vector<Type2,double>), Vector<Type1, SymbolicExpression>, Vector<Type2, SymbolicExpression>> {
  static SymbolicExpression call(
    const std::shared_ptr<Function<double(Vector<Type1,double>,Vector<Type2,double>)>> &func,
    const Vector<Type1, SymbolicExpression> &arg1, const Vector<Type2, SymbolicExpression> &arg2);
};

// ***** Vertex *****

//! A abstract class for a Vertex of the AST (abstract syntax tree).
class FMATVEC_EXPORT Vertex {
  friend Operation;
  friend NativeFunction;
  public:
    virtual ~Vertex() = default;

    //! Generate a new AST being the partial derivate of this AST with respect to the variable x.
    virtual SymbolicExpression parDer(const IndependentVariable &x) const=0;
    //! Rreturn true if this Vertex is a constant integer.
    inline virtual bool isConstantInt() const;
    //! Returns true if this Vertex is a constant with value 0.
    //! Note that only integer constants can be 0. Double constants are never treated as 0 by this function.
    bool isZero() const;
    //! Returns true if this Vertex is a constant with value 1.
    //! Note that only integer constants can be 1. Double constants are never treated as 1 by this function.
    bool isOne() const;

    virtual std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                          std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const=0;

    virtual void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const=0;

  protected:

    // helper function to make it easy to implement new expression optimizations. See ast.cc Operation::create for details.
    // Returns true if this Vertex the Vertex (SymbolicExpression) b. Every Symbol variables in this Vertex are free. This
    // means that true is also returned if the Symbols in this Vertex can be replaced by anything such that the expressions
    // are equal. All these required replacements are stored in m.
    // NOTE:
    // We define our own less operator for std::map<IndependentVariable, ...> since we abort in the overloaded function
    // inline bool operator<(const fmatvec::SymbolicExpression& a, const fmatvec::SymbolicExpression& b) noexcept,
    // see below. This abort is needed to avoid that e.g. std::min/std::max can be called with SymbolicExpression arguments.
    // If your code calls this abort function you may call somewhere std::min(SymbolicExpression, SymbolicExpression) instead of
    // fmatvec::min(SymbolicExpression, SymbolicExpression). The same applies for std::max/fmatvec::max.
    // I was not able to solve this in a better way!
    struct LessIV {
      bool operator()(const IndependentVariable& a, const IndependentVariable& b) const {
        return a.std::shared_ptr<const AST::Vertex>::get() < b.std::shared_ptr<const AST::Vertex>::get();
      }
    };
    using MapIVSE = std::map<IndependentVariable, SymbolicExpression, LessIV>;
    virtual bool equal(const SymbolicExpression &b, MapIVSE &m) const=0;
};

inline bool Vertex::isConstantInt() const {
  return false;
}

// ***** Constant *****

//! A vertex of the AST representing a constant (long or double)
template<class T>
class Constant : public Vertex, public std::enable_shared_from_this<Constant<T>> {
  friend SymbolicExpression;
  public:

    FMATVEC_EXPORT static SymbolicExpression create(const T& c_);
    FMATVEC_EXPORT SymbolicExpression parDer(const IndependentVariable &x) const override;
    inline bool isConstantInt() const override;
    //! Get the constant value of the vertex.
    inline const T& getValue() const;
    FMATVEC_EXPORT std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                  std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    FMATVEC_EXPORT void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

  private:

    Constant(const T& c_);
    const T c;
    bool equal(const SymbolicExpression &b, MapIVSE &m) const override;
    using CacheKey = T;
    static std::map<CacheKey, std::weak_ptr<const Constant>> cache;
};

template<>
inline bool Constant<double>::isConstantInt() const {
  return false;
}

template<>
inline bool Constant<long>::isConstantInt() const {
  return true;
}

template<class T>
const T& Constant<T>::getValue() const {
  return c;
}

// ***** Symbol *****

//! A vertex of the AST representing a independent variable.
class Symbol : public Vertex, public std::enable_shared_from_this<Symbol> {
  friend SymbolicExpression;
  friend IndependentVariable;
  public:

    FMATVEC_EXPORT static IndependentVariable create(const boost::uuids::uuid& uuid_=boost::uuids::random_generator()());
    FMATVEC_EXPORT SymbolicExpression parDer(const IndependentVariable &x) const override;
    //! Set the value of this independent variable.
    //! This has an influence on the evaluation of all ASTs which depend on this independent variable.
    inline void setValue(double x_) const;

    FMATVEC_EXPORT std::string getUUIDStr() const;

    FMATVEC_EXPORT std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                  std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    FMATVEC_EXPORT void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

  private:

    Symbol(const boost::uuids::uuid& uuid_);
    bool equal(const SymbolicExpression &b, MapIVSE &m) const override;
    mutable double x = 0.0;
    boost::uuids::uuid uuid; // each variable has a uuid (this is only used when the AST is serialized and for caching)
    using CacheKey = boost::uuids::uuid;
    static std::map<CacheKey, std::weak_ptr<const Symbol>> cache;
};

void Symbol::setValue(double x_) const {
  x=x_;
}

// ***** ScalarFunctionWrapArg *****

class ScalarFunctionWrapArg {
  public:
    virtual ~ScalarFunctionWrapArg() = default;
    virtual double operator()(const ByteCode::Arg &arg) = 0;
    virtual double dirDer(const ByteCode::Arg &arg) = 0;
    virtual double dirDerDirDer(const ByteCode::Arg &arg) = 0;
};

class ScalarFunctionWrapArgS : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgS(std::shared_ptr<fmatvec::Function<double(double)>> func_) : func(std::move(func_)) {}
    double operator()(const ByteCode::Arg &arg) override {
      return (*func)(*arg[0]);
    }
    double dirDer(const ByteCode::Arg &arg) override {
      return func->dirDer(*arg[1], *arg[0]);
    }
    double dirDerDirDer(const ByteCode::Arg &arg) override {
      return func->dirDerDirDer(*arg[1], *arg[2], *arg[0]);
    }
  private:
    std::shared_ptr<fmatvec::Function<double(double)>> func;
};

template<class Type>
class ScalarFunctionWrapArgV : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgV(const std::shared_ptr<fmatvec::Function<double(fmatvec::Vector<Type,double>)>> &func_, int argSize) : func(func_) {
      argV.resize(argSize);
      dir1V.resize(argSize);
      dir2V.resize(argSize);
    }
    double operator()(const ByteCode::Arg &arg) override {
      for(size_t i=0; i<argV.size(); ++i)
        argV(i) = *arg[i];
      return (*func)(argV);
    }
    double dirDer(const ByteCode::Arg &arg) override {
      size_t argVsize=argV.size();
      for(size_t i=0; i<argVsize; ++i)
        argV(i) = *arg[i];
      for(size_t i=0; i<argVsize; ++i)
        dir1V(i) = *arg[i+argVsize];
      return func->dirDer(dir1V, argV);
    }
    double dirDerDirDer(const ByteCode::Arg &arg) override {
      size_t argVsize=argV.size();
      for(size_t i=0; i<argVsize; ++i)
        argV(i) = *arg[i];
      for(size_t i=0; i<argVsize; ++i)
        dir1V(i) = *arg[i+argVsize];
      for(size_t i=0; i<argVsize; ++i)
        dir2V(i) = *arg[i+2*argVsize];
      return func->dirDerDirDer(dir1V, dir2V, argV);
    }
  private:
    std::shared_ptr<fmatvec::Function<double(fmatvec::Vector<Type,double>)>> func;
    fmatvec::Vector<Type,double> argV;
    fmatvec::Vector<Type,double> dir1V;
    fmatvec::Vector<Type,double> dir2V;
};

class ScalarFunctionWrapArgSS : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgSS(std::shared_ptr<fmatvec::Function<double(double,double)>> func_) : func(std::move(func_)) {}
    double operator()(const ByteCode::Arg &arg) override {
      return (*func)(*arg[0], *arg[1]);
    }
    double dirDer(const ByteCode::Arg &arg) override {
      return func->dirDer1(*arg[2], *arg[0], *arg[1]) +
             func->dirDer2(*arg[3], *arg[0], *arg[1]);
    }
    double dirDerDirDer(const ByteCode::Arg &arg) override {
      return func->dirDer1DirDer1(*arg[2], *arg[4], *arg[0], *arg[1]) + 
             func->dirDer2DirDer1(*arg[5], *arg[2], *arg[0], *arg[1]) +
             func->dirDer2DirDer1(*arg[3], *arg[4], *arg[0], *arg[1]) +
             func->dirDer2DirDer2(*arg[3], *arg[5], *arg[0], *arg[1]);
    }
  private:
    std::shared_ptr<fmatvec::Function<double(double,double)>> func;
};

template<class Type>
class ScalarFunctionWrapArgVS : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgVS(const std::shared_ptr<fmatvec::Function<double(Vector<Type,double>,double)>> &func_, int argSize) : func(func_) {
      argV.resize(argSize);
      dir1V.resize(argSize);
      dir2V.resize(argSize);
    }
    double operator()(const ByteCode::Arg &arg) override {
      for(size_t i=0; i<argV.size(); ++i)
        argV(i) = *arg[i];
      return (*func)(argV, *arg[argV.size()]);
    }
    double dirDer(const ByteCode::Arg &arg) override {
      auto argVsize=argV.size();
      for(size_t i=0; i<argVsize; ++i)
        argV(i) = *arg[i];
      for(size_t i=0; i<argVsize; ++i)
        dir1V(i) = *arg[i+1+argVsize];
      return func->dirDer1(dir1V             , argV, *arg[argVsize]) +
             func->dirDer2(*arg[2*argVsize+1], argV, *arg[argVsize]);
    }
    double dirDerDirDer(const ByteCode::Arg &arg) override {
      auto argVsize=argV.size();
      for(size_t i=0; i<argVsize; ++i)
        argV(i) = *arg[i];
      for(size_t i=0; i<argVsize; ++i)
        dir1V(i) = *arg[i+1+argVsize];
      for(size_t i=0; i<argVsize; ++i)
        dir2V(i) = *arg[i+2+2*argVsize];
      return func->dirDer1DirDer1(dir1V, dir2V, argV, *arg[argVsize]) + 
             func->dirDer2DirDer1(*arg[1+2*argVsize], dir2V             , argV, *arg[argVsize]) +
             func->dirDer2DirDer1(*arg[2+3*argVsize], dir1V             , argV, *arg[argVsize]) +
             func->dirDer2DirDer2(*arg[1+2*argVsize], *arg[2+3*argVsize], argV, *arg[argVsize]);
    }
  private:
    std::shared_ptr<fmatvec::Function<double(Vector<Type,double>,double)>> func;
    fmatvec::Vector<Type,double> argV;
    fmatvec::Vector<Type,double> dir1V;
    fmatvec::Vector<Type,double> dir2V;
};

template<class Type>
class ScalarFunctionWrapArgSV : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgSV(const std::shared_ptr<fmatvec::Function<double(double,Vector<Type,double>)>> &func_, int argSize) : func(func_) {
      argV.resize(argSize);
      dir1V.resize(argSize);
      dir2V.resize(argSize);
    }
    double operator()(const ByteCode::Arg &arg) override {
      for(size_t i=0; i<argV.size(); ++i)
        argV(i) = *arg[i+1];
      return (*func)(*arg[0], argV);
    }
    double dirDer(const ByteCode::Arg &arg) override {
      auto argVsize=argV.size();
      for(size_t i=0; i<argVsize; ++i)
        argV(i) = *arg[i+1];
      for(size_t i=0; i<argVsize; ++i)
        dir1V(i) = *arg[i+2+argVsize];
      return func->dirDer1(*arg[argVsize+1], *arg[0], argV) +
             func->dirDer2(dir1V           , *arg[0], argV);
    }
    double dirDerDirDer(const ByteCode::Arg &arg) override {
      auto argVsize=argV.size();
      for(size_t i=0; i<argVsize; ++i)
        argV(i) = *arg[i+1];
      for(size_t i=0; i<argVsize; ++i)
        dir1V(i) = *arg[i+2+argVsize];
      for(size_t i=0; i<argVsize; ++i)
        dir2V(i) = *arg[i+3+2*argVsize];
      return func->dirDer1DirDer1(*arg[argVsize+1], *arg[2*argVsize+2], *arg[0], argV) + 
             func->dirDer2DirDer1(dir2V           , *arg[argVsize+1]  , *arg[0], argV) +
             func->dirDer2DirDer1(dir1V           , *arg[2*argVsize+2], *arg[0], argV) +
             func->dirDer2DirDer2(dir1V           , dir2V             , *arg[0], argV);
    }
  private:
    std::shared_ptr<fmatvec::Function<double(double,Vector<Type,double>)>> func;
    fmatvec::Vector<Type,double> argV;
    fmatvec::Vector<Type,double> dir1V;
    fmatvec::Vector<Type,double> dir2V;
};

template<class Type1, class Type2>
class ScalarFunctionWrapArgVV : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgVV(const std::shared_ptr<fmatvec::Function<double(Vector<Type1,double>,Vector<Type2,double>)>> &func_,
                      int argaSize, int argbSize) : func(func_) {
      argV1.resize(argaSize);
      argV2.resize(argbSize);
      dir1V1.resize(argaSize);
      dir1V2.resize(argbSize);
      dir2V1.resize(argaSize);
      dir2V2.resize(argbSize);
    }
    double operator()(const ByteCode::Arg &arg) override {
      auto argV1size = argV1.size();
      auto argV2size = argV2.size();
      for(int i=0; i<argV1size; ++i)
        argV1(i)=*arg[i];
      for(int i=0; i<argV2size; ++i)
        argV2(i)=*arg[argV1size+i];
      return (*func)(argV1, argV2);
    }
    double dirDer(const ByteCode::Arg &arg) override {
      auto argV1size = argV1.size();
      auto argV2size = argV2.size();
      for(int i=0; i<argV1size; ++i)
        argV1(i)=*arg[i];
      for(int i=0; i<argV2size; ++i)
        argV2(i)=*arg[argV1size+i];
      for(size_t i=0; i<argV1size; ++i)
        dir1V1(i) = *arg[argV1size+argV2size+i];
      for(size_t i=0; i<argV2size; ++i)
        dir1V2(i) = *arg[2*argV1size+argV2size+i];
      return func->dirDer1(dir1V1, argV1, argV2) +
             func->dirDer2(dir1V2, argV1, argV2);
    }
    double dirDerDirDer(const ByteCode::Arg &arg) override {
      auto argV1size = argV1.size();
      auto argV2size = argV2.size();
      for(int i=0; i<argV1size; ++i)
        argV1(i)=*arg[i];
      for(int i=0; i<argV2size; ++i)
        argV2(i)=*arg[argV1size+i];
      for(size_t i=0; i<argV1size; ++i)
        dir1V1(i) = *arg[argV1size+argV2size+i];
      for(size_t i=0; i<argV2size; ++i)
        dir1V2(i) = *arg[2*argV1size+argV2size+i];
      for(size_t i=0; i<argV1size; ++i)
        dir2V1(i) = *arg[2*argV1size+2*argV2size+i];
      for(size_t i=0; i<argV2size; ++i)
        dir2V2(i) = *arg[3*argV1size+2*argV2size+i];
      return func->dirDer1DirDer1(dir1V1, dir2V1, argV1, argV2) + 
             func->dirDer2DirDer1(dir2V2, dir1V1, argV1, argV2) +
             func->dirDer2DirDer1(dir1V2, dir2V1, argV1, argV2) +
             func->dirDer2DirDer2(dir1V2, dir2V2, argV1, argV2);
    }
  private:
    std::shared_ptr<fmatvec::Function<double(Vector<Type1,double>,Vector<Type2,double>)>> func;
    fmatvec::Vector<Type1,double> argV1;
    fmatvec::Vector<Type2,double> argV2;
    fmatvec::Vector<Type1,double> dir1V1;
    fmatvec::Vector<Type2,double> dir1V2;
    fmatvec::Vector<Type1,double> dir2V1;
    fmatvec::Vector<Type2,double> dir2V2;
};

// ***** Function *****

//! A vertex of the AST representing an arbitary function.
class NativeFunction : public Vertex, public std::enable_shared_from_this<NativeFunction> {
  public:
    FMATVEC_EXPORT static SymbolicExpression create(const std::shared_ptr<ScalarFunctionWrapArg> &funcWrapper,
                                     const std::vector<SymbolicExpression> &argS,
                                     const std::vector<SymbolicExpression> &dir1S={},
                                     const std::vector<SymbolicExpression> &dir2S={});
    FMATVEC_EXPORT SymbolicExpression parDer(const IndependentVariable &x) const override;

    FMATVEC_EXPORT std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                                 std::map<const Vertex*,
                                                 std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    FMATVEC_EXPORT void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

  private:
    NativeFunction(const std::shared_ptr<ScalarFunctionWrapArg> &func_,
                   const std::vector<SymbolicExpression> &argS,
                   const std::vector<SymbolicExpression> &dir1S,
                   const std::vector<SymbolicExpression> &dir2S);
    bool equal(const SymbolicExpression &b, MapIVSE &m) const override;

    std::shared_ptr<ScalarFunctionWrapArg> funcWrapper;
    const std::vector<SymbolicExpression> argS;
    const std::vector<SymbolicExpression> dir1S;
    const std::vector<SymbolicExpression> dir2S;
    
    using CacheKey = std::tuple<std::weak_ptr<ScalarFunctionWrapArg>,
                     std::vector<std::weak_ptr<const Vertex>>,
                     std::vector<std::weak_ptr<const Vertex>>,
                     std::vector<std::weak_ptr<const Vertex>> >;
    struct CacheKeyComp {
      bool operator()(const CacheKey& l, const CacheKey& r) const;
      private:
        std::owner_less<std::weak_ptr<const Vertex>> ol;
        std::owner_less<std::weak_ptr<const ScalarFunctionWrapArg>> olf;
    };
    static std::map<CacheKey, std::weak_ptr<const NativeFunction>, CacheKeyComp> cache;
};

// ***** Operation *****

//! A vertex of the AST representing an operation.
class Operation : public Vertex, public std::enable_shared_from_this<Operation> {
  friend SymbolicExpression;
  friend SymbolicExpression fmatvec::AST::substScalar(const SymbolicExpression &se,
                                                      const IndependentVariable& a, const SymbolicExpression &b);
  friend boost::spirit::qi::rule<boost::spirit::istream_iterator, SymbolicExpression()>&
    fmatvec::getBoostSpiritQiRule<SymbolicExpression>();
  friend boost::spirit::karma::rule<std::ostream_iterator<char>, SymbolicExpression()>&
    fmatvec::getBoostSpiritKarmaRule<SymbolicExpression>();
  public:

    //! Defined operations.
    enum Operator { Plus, Minus, Mult, Div, Pow, Log, Sqrt, Neg, Sin, Cos, Tan, Sinh, Cosh, Tanh, ASin, ACos, ATan, ATan2, ASinh, ACosh, ATanh, Exp, Sign, Heaviside, Abs, Min, Max, Condition };
    FMATVEC_EXPORT static SymbolicExpression create(Operator op_, const std::vector<SymbolicExpression> &child_);
    FMATVEC_EXPORT SymbolicExpression parDer(const IndependentVariable &x) const override;

    FMATVEC_EXPORT Operator getOp() const { return op; }
    FMATVEC_EXPORT const std::vector<SymbolicExpression>& getChilds() const { return child; }
    FMATVEC_EXPORT std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                  std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    FMATVEC_EXPORT void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

  private:

    Operation(Operator op_, const std::vector<SymbolicExpression> &child_);
    bool equal(const SymbolicExpression &b, MapIVSE &m) const override;
    Operator op;
    std::vector<SymbolicExpression> child;
    using CacheKey = std::pair<Operator, std::vector<std::weak_ptr<const Vertex>>>;
    struct CacheKeyComp {
      bool operator()(const CacheKey& l, const CacheKey& r) const;
      private:
        std::owner_less<std::weak_ptr<const Vertex>> ol;
    };
    static std::map<CacheKey, std::weak_ptr<const Operation>, CacheKeyComp> cache;

    struct OpMap {
      std::string funcName; // used to dump/read the expression using boost spirit/karma: e.g.
                            // "plus"
      std::function<void(double*, const ByteCode::Arg&)> func; // used for runtime evaluation: e.g.
                            // [](double* r, const ByteCode::Arg& a){ *r = *a[0] + *a[1]; }
    };
    static const std::map<Operator, OpMap> opMap;
};

inline SymbolicExpression SymbolicFuncWrapArg1<double(double), SymbolicExpression>::call(
  const std::shared_ptr<fmatvec::Function<double(double)>> &func,
  const SymbolicExpression &arg) {
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgS>(func), std::vector<SymbolicExpression>{arg});
}

template<class Type>
SymbolicExpression SymbolicFuncWrapArg1<double(Vector<Type, double>), Vector<Type, SymbolicExpression>>::call(
  const std::shared_ptr<Function<double(Vector<Type, double>)>> &func,
  const Vector<Type, SymbolicExpression> &arg) {
  std::vector<SymbolicExpression> argV(arg.size());
  for(int i=0; i<arg.size(); ++i)
    argV[i]=arg(i);
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgV<Type>>(func, arg.size()), argV);
}

inline SymbolicExpression SymbolicFuncWrapArg2<double(double,double), SymbolicExpression, SymbolicExpression>::call(
  const std::shared_ptr<fmatvec::Function<double(double,double)>> &func,
  const SymbolicExpression &arg1, const SymbolicExpression &arg2) {
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgSS>(func), std::vector<SymbolicExpression>{arg1,arg2});
}

template<class Type>
SymbolicExpression SymbolicFuncWrapArg2<double(Vector<Type,double>,double), Vector<Type, SymbolicExpression>, SymbolicExpression>::call(
  const std::shared_ptr<Function<double(Vector<Type,double>,double)>> &func,
  const Vector<Type, SymbolicExpression> &arg1, const SymbolicExpression &arg2) {
  std::vector<SymbolicExpression> arg(arg1.size()+1);
  for(int i=0; i<arg1.size(); ++i)
    arg[i]=arg1(i);
  arg[arg1.size()]=arg2;
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgVS<Type>>(func, arg1.size()), arg);
}

template<class Type>
SymbolicExpression SymbolicFuncWrapArg2<double(double,Vector<Type,double>), SymbolicExpression, Vector<Type, SymbolicExpression>>::call(
  const std::shared_ptr<Function<double(double,Vector<Type,double>)>> &func,
  const SymbolicExpression &arg1, const Vector<Type, SymbolicExpression> &arg2) {
  std::vector<SymbolicExpression> arg(1+arg2.size());
  arg[0]=arg1;
  for(int i=0; i<arg2.size(); ++i)
    arg[i+1]=arg2(i);
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgSV<Type>>(func, arg2.size()), arg);
}

template<class Type1, class Type2>
SymbolicExpression SymbolicFuncWrapArg2<double(Vector<Type1,double>,Vector<Type2,double>), Vector<Type1, SymbolicExpression>, Vector<Type2, SymbolicExpression>>::call(
  const std::shared_ptr<Function<double(Vector<Type1,double>,Vector<Type2,double>)>> &func,
  const Vector<Type1, SymbolicExpression> &arg1, const Vector<Type2, SymbolicExpression> &arg2) {
  std::vector<SymbolicExpression> arg(arg1.size()+arg2.size());
  for(int i=0; i<arg1.size(); ++i)
    arg[i]=arg1(i);
  for(int i=0; i<arg2.size(); ++i)
    arg[arg1.size()+i]=arg2(i);
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgVV<Type1,Type2>>(func, arg1.size(), arg2.size()), arg);
}

} // end namespace AST
#endif

inline bool operator<(const fmatvec::SymbolicExpression& a, const fmatvec::SymbolicExpression& b) noexcept {
  // See LessIV why we need to overload this operator< and call abort!
  std::cerr<<"ABORT: See LessIV: have you called std::min/max instead of fmatvev::min/max with type SymbolicExpression?"<<std::endl;
  assert(((void)"See LessIV: have you called std::min/max instead of fmatvev::min/max with type SymbolicExpression?", false));
  abort();
}

IndependentVariable& IndependentVariable::operator^=(double x) {
  static_cast<const AST::Symbol*>(get())->setValue(x);
  return *this;
}

template<> boost::spirit::qi::rule<boost::spirit::istream_iterator, IndependentVariable()>& getBoostSpiritQiRule<IndependentVariable>();
template<> boost::spirit::qi::rule<boost::spirit::istream_iterator, SymbolicExpression()>& getBoostSpiritQiRule<SymbolicExpression>();

template<> boost::spirit::karma::rule<std::ostream_iterator<char>, IndependentVariable()>& getBoostSpiritKarmaRule<IndependentVariable>();
template<> boost::spirit::karma::rule<std::ostream_iterator<char>, SymbolicExpression()>& getBoostSpiritKarmaRule<SymbolicExpression>();

// for vectors/matrices of type IndependentVariable/SymbolicExpression throw exceptions instead of assert.
// these exceptions are thrown regardless whether NDEBUG is defined or not.
// this is done since vectors/matrix operations of this type are usually user input, should not abort/assert on errors,
// and it not time critical since it just build the expression once which is than evaluated at runtime.
template<> struct AssertUseException<IndependentVariable> { constexpr static bool value = true; };
template<> struct AssertUseException<SymbolicExpression>  { constexpr static bool value = true; };

} // end namespace fmatvec

#endif
