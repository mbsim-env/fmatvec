#ifndef _FMATVEC_AST_H_
#define _FMATVEC_AST_H_

#include <vector>
#include <map>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include "function.h"

// the following two lines are a workaround for a bug in boost 1.69
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>

#include <boost/uuid/uuid_generators.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <fmatvec/types.h>
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
}

//! A symbolic expression.
//! This class represent an arbitary symbolic expression.
//! This can either be the special case of just a symbol
//! or an arbitary expression consisting of a hierarchiy
//! of operations consisting itself of expressions, symbols or constants.
class FMATVEC_EXPORT SymbolicExpression : public std::shared_ptr<const AST::Vertex> {
  friend class AST::Operation;
  friend class AST::Constant<int>;
  friend class AST::Constant<double>;
  friend class AST::NativeFunction;
  friend FMATVEC_EXPORT SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep);
  friend FMATVEC_EXPORT double eval(const SymbolicExpression &x);
  friend FMATVEC_EXPORT SymbolicExpression subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);
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
    SymbolicExpression();
    //! Creates a expression form the specified string (the string must be a serialized SymbolicExpression).
    SymbolicExpression(const std::string &str);
    //! Creates a integer constant.
    SymbolicExpression(int x);
    //! Creates a double constant.
    SymbolicExpression(double x);

    //! Garbage collect everything.
    //! The only "garbage" which can be left by this class are empty weak_ptr's stored in global maps.
    //! This routine removes such garbage.
    static void garbageCollect();

    // default ctors and assignment operators
    SymbolicExpression(const SymbolicExpression& x) = default;
    SymbolicExpression(SymbolicExpression&& x) = default;
    SymbolicExpression& operator=(const SymbolicExpression& x) = default;
    SymbolicExpression& operator=(SymbolicExpression&& x) = default;
    ~SymbolicExpression() = default;

    // operators
#ifndef SWIG
    SymbolicExpression operator+(const SymbolicExpression &b) const;
    SymbolicExpression operator-(const SymbolicExpression &b) const;
    SymbolicExpression operator*(const SymbolicExpression &b) const;
    SymbolicExpression operator/(const SymbolicExpression &b) const;
#endif
    SymbolicExpression operator-() const;
    SymbolicExpression& operator+=(const SymbolicExpression &b);
    SymbolicExpression& operator-=(const SymbolicExpression &b);
    SymbolicExpression& operator*=(const SymbolicExpression &b);
    SymbolicExpression& operator/=(const SymbolicExpression &b);

    // We implement operator<<= to be able to use SymbolicExpression fully equivalent to the vector/matrix classes
    // which have a assign method which redimension before assigning. Redimensioning is not required for a scalar
    // type like assign
    SymbolicExpression& operator<<=(const SymbolicExpression &src);

#if !defined(NDEBUG) && !defined(SWIG)
    static unsigned long evalOperationsCount;
#endif
};

//! A independent variable.
//! Any SymbolicExpression can be partialliy differentiated with respect to a independent variable.
//! An independent varible can also be assigned a value which is used if eval is called.
class FMATVEC_EXPORT IndependentVariable : public SymbolicExpression {
  friend class AST::Symbol;
  friend FMATVEC_EXPORT std::istream& operator>>(std::istream& s, IndependentVariable &v);
  public:
    //! Creates a IndependentVariable variable (each call to this ctor creates a new independent variable)
    IndependentVariable();
    //! Creates a IndependentVariable variable from the specified string (the string is a serialized IndependentVariable).
    IndependentVariable(const std::string &str);

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
FMATVEC_OPERATORRESULT2(SymbolicExpression, double, SymbolicExpression)

FMATVEC_OPERATORRESULT1(IndependentVariable, SymbolicExpression)

FMATVEC_OPERATORRESULT2(IndependentVariable, int, SymbolicExpression)
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

//! Generate a new SymbolicExpression being the partial derivate of dep
//! with respect to indep (indep must be a symbol).
FMATVEC_EXPORT SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep);

//! Evaluate the SymbolicExpression.
//! The returned value depends on the symbolic expression and on the current values of all independent
//! variables this symbolic expression depends on.
//! Also see Symbol ant Vertex::getDependsOn().
inline double eval(const SymbolicExpression &x);

//! Write a SymbolicExpression to a stream using serialization.
FMATVEC_EXPORT std::ostream& operator<<(std::ostream& s, const SymbolicExpression& se);
//! Create/initialize a SymbolicExpression from a stream using deserialization.
FMATVEC_EXPORT std::istream& operator>>(std::istream& s, SymbolicExpression &se);
//! Create/initialize a IndependentVariable from a stream using deserialization.
FMATVEC_EXPORT std::istream& operator>>(std::istream& s, IndependentVariable &v);
//! Substitutes in se the independent variable a with the expression b and returns the new expression.
FMATVEC_EXPORT SymbolicExpression subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);

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
FMATVEC_EXPORT SymbolicExpression asinh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression acosh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression atanh(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression exp(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression sign(const SymbolicExpression &a);
FMATVEC_EXPORT SymbolicExpression abs(const SymbolicExpression &a);

#ifndef SWIG
namespace AST { // internal namespace

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

    //! Evaluate the AST.
    //! The returned value depends on the AST and on the current values of all independent variables this AST depends on.
    //! Also see Symbol ant Vertex::getDependsOn().
    virtual double eval() const=0;
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
    //! Return a list of all variables this AST depends on.
    const std::map<std::weak_ptr<const Symbol>, unsigned long, std::owner_less<std::weak_ptr<const Symbol>>>& getDependsOn() const;

  protected:

    // helper function to make it easy to implement new expression optimizations. See ast.cc Operation::create for details.
    // Returns true if this Vertex the Vertex (SymbolicExpression) b. Every Symbol variables in this Vertex are free. This
    // means that true is also returned if the Symbols in this Vertex can be replaced by anything such that the expressions
    // are equal. All these required replacements are stored in m.
    virtual bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const=0;

    // the value of this map (the unsigned long = version) must be mutable
    mutable std::map<std::weak_ptr<const Symbol>, unsigned long, std::owner_less<std::weak_ptr<const Symbol>>> dependsOn;
};

inline bool Vertex::isConstantInt() const {
  return false;
}

// ***** Constant *****

//! A vertex of the AST representing a constant (int or double)
template<class T>
class FMATVEC_EXPORT Constant : public Vertex, public std::enable_shared_from_this<Constant<T>> {
  friend SymbolicExpression;
  public:

    static SymbolicExpression create(const T& c_);
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    inline bool isConstantInt() const override;
    //! Get the constant value of the vertex.
    inline const T& getValue() const;

  private:

    Constant(const T& c_);
    const T c;
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;
    typedef T CacheKey;
    static std::map<CacheKey, std::weak_ptr<const Constant>> cache;
};

template<>
inline bool Constant<double>::isConstantInt() const {
  return false;
}

template<>
inline bool Constant<int>::isConstantInt() const {
  return true;
}

template<class T>
double Constant<T>::eval() const {
  return c;
}

template<class T>
const T& Constant<T>::getValue() const {
  return c;
}

// ***** Symbol *****

//! A vertex of the AST representing a independent variable.
class FMATVEC_EXPORT Symbol : public Vertex, public std::enable_shared_from_this<Symbol> {
  friend SymbolicExpression;
  friend IndependentVariable;
  public:

    static IndependentVariable create(const boost::uuids::uuid& uuid_=boost::uuids::random_generator()());
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    //! Set the value of this independent variable.
    //! This has an influence on the evaluation of all ASTs which depend on this independent variable.
    inline void setValue(double x_) const;
    //! Eeturn the "version" of this variable.
    //! Each call to Symbol::set increased this "version count". This is used for caching evaluations.
    inline unsigned long getVersion() const;

    std::string getUUIDStr() const;

  private:

    Symbol(const boost::uuids::uuid& uuid_);
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;
    mutable double x = 0.0;
    mutable unsigned long version;
    boost::uuids::uuid uuid; // each variable has a uuid (this is only used when the AST is serialized and for caching)
    typedef boost::uuids::uuid CacheKey;
    static std::map<CacheKey, std::weak_ptr<const Symbol>> cache;
};

void Symbol::setValue(double x_) const {
  version++;
  x=x_;
}

double Symbol::eval() const {
  return x;
}

unsigned long Symbol::getVersion() const {
  return version;
}

// ***** ScalarFunctionWrapArg *****

class ScalarFunctionWrapArg {
  public:
    virtual double operator()(const std::vector<double> &arg) = 0;
    virtual double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) = 0;
    virtual double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) = 0;
};

class ScalarFunctionWrapArgS : public ScalarFunctionWrapArg {
  public:
    ScalarFunctionWrapArgS(const std::shared_ptr<fmatvec::Function<double(double)>> &func_) : func(func_) {}
    double operator()(const std::vector<double> &arg) override {
      return (*func)(arg[0]);
    }
    double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) override {
      return func->dirDer(dir1[0], arg[0]);
    }
    double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) override {
      return func->dirDerDirDer(dir1[0], dir2[0], arg[0]);
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
    double operator()(const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size(); ++i)
        argV(i) = arg[i];
      return (*func)(argV);
    }
    double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size(); ++i)
        argV(i) = arg[i];
      for(size_t i=0; i<arg.size(); ++i)
        dir1V(i) = dir1[i];
      return func->dirDer(dir1V, argV);
    }
    double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size(); ++i)
        argV(i) = arg[i];
      for(size_t i=0; i<arg.size(); ++i)
        dir1V(i) = dir1[i];
      for(size_t i=0; i<arg.size(); ++i)
        dir2V(i) = dir2[i];
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
    ScalarFunctionWrapArgSS(const std::shared_ptr<fmatvec::Function<double(double,double)>> &func_) : func(func_) {}
    double operator()(const std::vector<double> &arg) override {
      return (*func)(arg[0], arg[1]);
    }
    double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) override {
      return func->dirDer1(dir1[0], arg[0], arg[1]) +
             func->dirDer2(dir1[1], arg[0], arg[1]);
    }
    double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) override {
      return func->dirDer1DirDer1(dir1[0], dir2[0], arg[0], arg[1]) + 
             func->dirDer2DirDer1(dir2[1], dir1[0], arg[0], arg[1]) +
             func->dirDer2DirDer1(dir1[1], dir2[0], arg[0], arg[1]) +
             func->dirDer2DirDer2(dir1[1], dir2[1], arg[0], arg[1]);
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
    double operator()(const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size()-1; ++i)
        argV(i) = arg[i];
      return (*func)(argV, arg[arg.size()-1]);
    }
    double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size()-1; ++i)
        argV(i) = arg[i];
      for(size_t i=0; i<arg.size()-1; ++i)
        dir1V(i) = dir1[i];
      return func->dirDer1(dir1V             , argV, arg[arg.size()-1]) +
             func->dirDer2(dir1[arg.size()-1], argV, arg[arg.size()-1]);
    }
    double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size()-1; ++i)
        argV(i) = arg[i];
      for(size_t i=0; i<arg.size()-1; ++i)
        dir1V(i) = dir1[i];
      for(size_t i=0; i<arg.size()-1; ++i)
        dir2V(i) = dir2[i];
      return func->dirDer1DirDer1(dir1V             , dir2V             , argV, arg[arg.size()-1]) + 
             func->dirDer2DirDer1(dir2[arg.size()-1], dir1V             , argV, arg[arg.size()-1]) +
             func->dirDer2DirDer1(dir1[arg.size()-1], dir2V             , argV, arg[arg.size()-1]) +
             func->dirDer2DirDer2(dir1[arg.size()-1], dir2[arg.size()-1], argV, arg[arg.size()-1]);
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
    double operator()(const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size()-1; ++i)
        argV(i) = arg[i+1];
      return (*func)(arg[0], argV);
    }
    double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size()-1; ++i)
        argV(i) = arg[i+1];
      for(size_t i=0; i<arg.size()-1; ++i)
        dir1V(i) = dir1[i+1];
      return func->dirDer1(dir1[0], arg[0], argV) +
             func->dirDer2(dir1V  , arg[0], argV);
    }
    double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) override {
      for(size_t i=0; i<arg.size()-1; ++i)
        argV(i) = arg[i+1];
      for(size_t i=0; i<arg.size()-1; ++i)
        dir1V(i) = dir1[i+1];
      for(size_t i=0; i<arg.size()-1; ++i)
        dir2V(i) = dir2[i+1];
      return func->dirDer1DirDer1(dir1[0], dir2[0], arg[0], argV) + 
             func->dirDer2DirDer1(dir2V  , dir1[0], arg[0], argV) +
             func->dirDer2DirDer1(dir1V  , dir2[0], arg[0], argV) +
             func->dirDer2DirDer2(dir1V  , dir2V  , arg[0], argV);
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
    double operator()(const std::vector<double> &arg) override {
      for(int i=0; i<argV1.size(); ++i)
        argV1(i)=arg[i];
      for(int i=0; i<argV2.size(); ++i)
        argV2(i)=arg[argV1.size()+i];
      return (*func)(argV1, argV2);
    }
    double dirDer(const std::vector<double> &dir1, const std::vector<double> &arg) override {
      for(int i=0; i<argV1.size(); ++i)
        argV1(i)=arg[i];
      for(int i=0; i<argV2.size(); ++i)
        argV2(i)=arg[argV1.size()+i];
      for(int i=0; i<argV1.size(); ++i)
        dir1V1(i)=dir1[i];
      for(int i=0; i<argV2.size(); ++i)
        dir1V2(i)=dir1[dir1V1.size()+i];
      return func->dirDer1(dir1V1, argV1, argV2) +
             func->dirDer2(dir1V2, argV1, argV2);
    }
    double dirDerDirDer(const std::vector<double> &dir1, const std::vector<double> &dir2, const std::vector<double> &arg) override {
      for(int i=0; i<argV1.size(); ++i)
        argV1(i)=arg[i];
      for(int i=0; i<argV2.size(); ++i)
        argV2(i)=arg[argV1.size()+i];
      for(int i=0; i<argV1.size(); ++i)
        dir1V1(i)=dir1[i];
      for(int i=0; i<argV2.size(); ++i)
        dir1V2(i)=dir1[dir1V1.size()+i];
      for(int i=0; i<argV1.size(); ++i)
        dir2V1(i)=dir2[i];
      for(int i=0; i<argV2.size(); ++i)
        dir2V2(i)=dir2[dir2V1.size()+i];
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
class FMATVEC_EXPORT NativeFunction : public Vertex, public std::enable_shared_from_this<NativeFunction> {
  public:
    static SymbolicExpression create(const std::shared_ptr<ScalarFunctionWrapArg> &funcWrapper, const std::vector<SymbolicExpression> &argS,
                                     const std::vector<SymbolicExpression> &dir1S={}, const std::vector<SymbolicExpression> &dir2S={});
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;

  private:
    NativeFunction(const std::shared_ptr<ScalarFunctionWrapArg> &func_, const std::vector<SymbolicExpression> &argS,
                   const std::vector<SymbolicExpression> &dir1S, const std::vector<SymbolicExpression> &dir2S);
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;

    std::shared_ptr<ScalarFunctionWrapArg> funcWrapper;
    const std::vector<SymbolicExpression> argS;
    const std::vector<SymbolicExpression> dir1S;
    const std::vector<SymbolicExpression> dir2S;
    mutable std::vector<double> argN;
    mutable std::vector<double> dir1N;
    mutable std::vector<double> dir2N;
    
    typedef std::tuple<std::weak_ptr<ScalarFunctionWrapArg>,
                       std::vector<std::weak_ptr<const Vertex>>,
                       std::vector<std::weak_ptr<const Vertex>>,
                       std::vector<std::weak_ptr<const Vertex>> > CacheKey;
    struct CacheKeyComp {
      bool operator()(const CacheKey& l, const CacheKey& r) const;
      private:
        std::owner_less<std::weak_ptr<const Vertex>> ol;
        std::owner_less<std::weak_ptr<const ScalarFunctionWrapArg>> olf;
    };
    static std::map<CacheKey, std::weak_ptr<const NativeFunction>, CacheKeyComp> cache;
    mutable double cacheValue;
};

double NativeFunction::eval() const {
  bool useCacheValue=true;
  for(auto &d : dependsOn) {
    auto dFirst=d.first.lock();
    if(dFirst->getVersion() > d.second)
      useCacheValue=false;
    d.second=dFirst->getVersion();
  }
  if(useCacheValue)
    return cacheValue;

  for(size_t i=0; i<argS.size(); ++i)
    argN[i]=fmatvec::eval(argS[i]);

  if(dir1S.size()==0) {
    // calculate function value (0th derivative)
    return cacheValue = (*funcWrapper)(argN);
  }
  if(dir2S.size()==0) {
    // calculate first derivative
    for(size_t i=0; i<dir1S.size(); ++i)
      dir1N[i]=fmatvec::eval(dir1S[i]);
    return cacheValue = funcWrapper->dirDer(dir1N, argN);
  }
  // calculate second derivative
  for(size_t i=0; i<dir1S.size(); ++i)
    dir1N[i]=fmatvec::eval(dir1S[i]);
  for(size_t i=0; i<dir2S.size(); ++i)
    dir2N[i]=fmatvec::eval(dir2S[i]);
  return cacheValue = funcWrapper->dirDerDirDer(dir1N, dir2N, argN);
}

// ***** Operation *****

//! A vertex of the AST representing an operation.
class FMATVEC_EXPORT Operation : public Vertex, public std::enable_shared_from_this<Operation> {
  friend SymbolicExpression;
  friend SymbolicExpression fmatvec::subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);
  friend boost::spirit::qi::rule<boost::spirit::istream_iterator, SymbolicExpression()>&
    fmatvec::getBoostSpiritQiRule<SymbolicExpression>();
  friend boost::spirit::karma::rule<std::ostream_iterator<char>, SymbolicExpression()>&
    fmatvec::getBoostSpiritKarmaRule<SymbolicExpression>();
  public:

    //! Defined operations.
    enum Operator { Plus, Minus, Mult, Div, Pow, Log, Sqrt, Neg, Sin, Cos, Tan, Sinh, Cosh, Tanh, ASin, ACos, ATan, ASinh, ACosh, ATanh, Exp, Sign, Abs };
    static SymbolicExpression create(Operator op_, const std::vector<SymbolicExpression> &child_);
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;

    Operator getOp() const { return op; }
    const std::vector<SymbolicExpression>& getChilds() const { return child; }

  private:

    Operation(Operator op_, const std::vector<SymbolicExpression> &child_);
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;
    Operator op;
    std::vector<SymbolicExpression> child;
    typedef std::pair<Operator, std::vector<std::weak_ptr<const Vertex>>> CacheKey;
    struct CacheKeyComp {
      bool operator()(const CacheKey& l, const CacheKey& r) const;
      private:
        std::owner_less<std::weak_ptr<const Vertex>> ol;
    };
    static std::map<CacheKey, std::weak_ptr<const Operation>, CacheKeyComp> cache;
    mutable double cacheValue;
    static const std::map<Operator, std::string> opMap;
};

double Operation::eval() const {
  if(!dependsOn.empty()) {
    // we eval expressions having no dependency always (such expressions are optimized out during Operation creation)
    bool useCacheValue=true;
    for(auto &d : dependsOn) {
      auto dFirst=d.first.lock();
      if(dFirst->getVersion() > d.second)
        useCacheValue=false;
      d.second=dFirst->getVersion();
    }
    if(useCacheValue)
      return cacheValue;
  }
#if !defined(NDEBUG) && !defined(SWIG)
  SymbolicExpression::evalOperationsCount++;
#endif

  // we do not use "double a=child.size()>=1 ? child[0]->eval() : 0;" here
  // for optimal performance, to avoid the child.size() call.
  #define a child[0]->eval() // value of first argument
  #define b child[1]->eval() // value of second argument
  switch(op) {
    case Plus:
      return cacheValue = a + b;
    case Minus:
      return cacheValue = a - b;
    case Mult:
      return cacheValue = a * b;
    case Div:
      return cacheValue = a / b;
    case Pow:
      if(child[1]->isConstantInt())
        return cacheValue = std::pow(a, static_cast<const Constant<int>*>(child[1].get())->getValue());
      else
        return cacheValue = std::pow(a, b);
    case Log:
      return cacheValue = std::log(a);
    case Sqrt:
      return cacheValue = std::sqrt(a);
    case Neg:
      return cacheValue = - a;
    case Sin:
      return cacheValue = std::sin(a);
    case Cos:
      return cacheValue = std::cos(a);
    case Tan:
      return cacheValue = std::tan(a);
    case Sinh:
      return cacheValue = std::sinh(a);
    case Cosh:
      return cacheValue = std::cosh(a);
    case Tanh:
      return cacheValue = std::tanh(a);
    case ASin:
      return cacheValue = std::asin(a);
    case ACos:
      return cacheValue = std::acos(a);
    case ATan:
      return cacheValue = std::atan(a);
    case ASinh:
      return cacheValue = std::asinh(a);
    case ACosh:
      return cacheValue = std::acosh(a);
    case ATanh:
      return cacheValue = std::atanh(a);
    case Exp:
      return cacheValue = std::exp(a);
    case Sign:
      return cacheValue = boost::math::sign(a);
    case Abs:
      return cacheValue = std::abs(a);
  }
  #undef a
  #undef b
  throw std::runtime_error("Unknown operation.");
}

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
  std::vector<SymbolicExpression> argV(arg1.size()+1);
  for(int i=0; i<arg1.size(); ++i)
    argV[i]=arg1(i);
  argV[arg1.size()]=arg2;
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgVS<Type>>(func, arg1.size()), argV);
}

template<class Type>
SymbolicExpression SymbolicFuncWrapArg2<double(double,Vector<Type,double>), SymbolicExpression, Vector<Type, SymbolicExpression>>::call(
  const std::shared_ptr<Function<double(double,Vector<Type,double>)>> &func,
  const SymbolicExpression &arg1, const Vector<Type, SymbolicExpression> &arg2) {
  std::vector<SymbolicExpression> argV(1+arg2.size());
  argV[0]=arg1;
  for(int i=0; i<arg2.size(); ++i)
    argV[i+1]=arg2(i);
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgSV<Type>>(func, arg2.size()), argV);
}

template<class Type1, class Type2>
SymbolicExpression SymbolicFuncWrapArg2<double(Vector<Type1,double>,Vector<Type2,double>), Vector<Type1, SymbolicExpression>, Vector<Type2, SymbolicExpression>>::call(
  const std::shared_ptr<Function<double(Vector<Type1,double>,Vector<Type2,double>)>> &func,
  const Vector<Type1, SymbolicExpression> &arg1, const Vector<Type2, SymbolicExpression> &arg2) {
  std::vector<SymbolicExpression> argV(arg1.size()+arg2.size());
  for(int i=0; i<arg1.size(); ++i)
    argV[i]=arg1(i);
  for(int i=0; i<arg2.size(); ++i)
    argV[arg1.size()+i]=arg2(i);
  return AST::NativeFunction::create(std::make_shared<AST::ScalarFunctionWrapArgVV<Type1,Type2>>(func, arg1.size(), arg2.size()), argV);
}

} // end namespace AST
#endif

IndependentVariable& IndependentVariable::operator^=(double x) {
  static_cast<const AST::Symbol*>(get())->setValue(x);
  return *this;
}

double eval(const SymbolicExpression &x) {
  return x->eval();
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
