#ifndef _FMATVEC_AST_H_
#define _FMATVEC_AST_H_

#include <vector>
#include <map>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/bimap.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <fmatvec/types.h>

namespace fmatvec {

// forward declarations

class IndependentVariable;

namespace AST {
  class Vertex;
  class Symbol;
  class Operation;
  template<class T> class Constant;
}

//! A symbolic expression.
//! This class represent an arbitary symbolic expression.
//! This can either be the special case of just a symbol
//! or an arbitary expression consisting of a hierarchiy
//! of operations consisting itself of expressions, symbols or constants.
class SymbolicExpression : public std::shared_ptr<const AST::Vertex> {
  friend class AST::Symbol;
  friend class AST::Operation;
  friend class AST::Constant<int>;
  friend class AST::Constant<double>;
  friend SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep);
  friend double eval(const SymbolicExpression &x);
  friend SymbolicExpression subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);
  friend std::ostream& operator<<(std::ostream& s, const SymbolicExpression& se);
  friend std::istream& operator>>(std::istream& s, SymbolicExpression &se);
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

    // We implement an assign method to be able to use SymbolicExpression fully equivalent to the vector/matrix classes
    // which have a assign method which redimension before assigning. Redimensioning is not required for a scalar
    // type like assign
    SymbolicExpression& assign(const SymbolicExpression &src);
    // We implement an reassign method to be able to use SymbolicExpression fully equivalent to the vector/matrix classes...
    SymbolicExpression& reassign(const SymbolicExpression &src);

#if !defined(NDEBUG) && !defined(SWIG)
    static unsigned long evalOperationsCount;
#endif
};

//! A independent variable.
//! Any SymbolicExpression can be partialliy differentiated with respect to a independent variable.
//! An independent varible can also be assigned a value which is used if eval is called.
class IndependentVariable : public SymbolicExpression {
  friend class AST::Symbol;
  friend std::istream& operator>>(std::istream& s, IndependentVariable &v);
  public:
    //! Creates a IndependentVariable variable (each call to this ctor creates a new independent variable)
    IndependentVariable();
    //! Creates a IndependentVariable variable from the specified string (the string is a serialized IndependentVariable).
    IndependentVariable(const std::string &str);

    //! Set the double value of the independent value.
    inline IndependentVariable& operator&=(double x);

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

SymbolicExpression operator+(double a, const SymbolicExpression &b);
SymbolicExpression operator-(double a, const SymbolicExpression &b);
SymbolicExpression operator*(double a, const SymbolicExpression &b);
SymbolicExpression operator/(double a, const SymbolicExpression &b);
SymbolicExpression operator+(int a, const SymbolicExpression &b);
SymbolicExpression operator-(int a, const SymbolicExpression &b);
SymbolicExpression operator*(int a, const SymbolicExpression &b);
SymbolicExpression operator/(int a, const SymbolicExpression &b);

//! Generate a new SymbolicExpression being the partial derivate of dep
//! with respect to indep (indep must be a symbol).
SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep);

//! Evaluate the SymbolicExpression.
//! The returned value depends on the symbolic expression and on the current values of all independent
//! variables this symbolic expression depends on.
//! Also see Symbol ant Vertex::getDependsOn().
inline double eval(const SymbolicExpression &x);

//! Write a SymbolicExpression to a stream using serialization.
std::ostream& operator<<(std::ostream& s, const SymbolicExpression& se);
//! Create/initialize a SymbolicExpression from a stream using deserialization.
std::istream& operator>>(std::istream& s, SymbolicExpression &se);
//! Create/initialize a IndependentVariable from a stream using deserialization.
std::istream& operator>>(std::istream& s, IndependentVariable &v);
//! Substitutes in se the independent variable a with the expression b and returns the new expression.
SymbolicExpression subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);

// function operations overloaded for SymbolicExpression
SymbolicExpression pow(const SymbolicExpression &a, const SymbolicExpression &b);
SymbolicExpression log(const SymbolicExpression &a);
SymbolicExpression sqrt(const SymbolicExpression &a);
SymbolicExpression sin(const SymbolicExpression &a);
SymbolicExpression cos(const SymbolicExpression &a);
SymbolicExpression tan(const SymbolicExpression &a);
SymbolicExpression sinh(const SymbolicExpression &a);
SymbolicExpression cosh(const SymbolicExpression &a);
SymbolicExpression tanh(const SymbolicExpression &a);
SymbolicExpression asin(const SymbolicExpression &a);
SymbolicExpression acos(const SymbolicExpression &a);
SymbolicExpression atan(const SymbolicExpression &a);
SymbolicExpression asinh(const SymbolicExpression &a);
SymbolicExpression acosh(const SymbolicExpression &a);
SymbolicExpression atanh(const SymbolicExpression &a);
SymbolicExpression exp(const SymbolicExpression &a);
SymbolicExpression sign(const SymbolicExpression &a);
SymbolicExpression abs(const SymbolicExpression &a);

#ifndef SWIG
namespace AST { // internal namespace

// ***** Vertex *****

//! A abstract class for a Vertex of the AST (abstract syntax tree).
class Vertex {
  friend Operation;
  public:

    //! Create a AST from a by deserialization a stream.
    static SymbolicExpression createFromStream(std::istream &s);
    //! Evaluate the AST.
    //! The returned value depends on the AST and on the current values of all independent variables this AST depends on.
    //! Also see Symbol ant Vertex::getDependsOn().
    virtual double eval() const=0;
    //! Generate a new AST being the partial derivate of this AST with respect to the variable x.
    virtual SymbolicExpression parDer(const IndependentVariable &x) const=0;
    //! Write/serailize a AST to a stream.
    virtual void serializeToStream(std::ostream &s) const=0;
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

bool Vertex::isConstantInt() const {
  return false;
}

// ***** Constant *****

//! A vertex of the AST representing a constant (int or double)
template<class T>
class Constant : public Vertex, public std::enable_shared_from_this<Constant<T>> {
  friend SymbolicExpression;
  public:

    static SymbolicExpression create(const T& c_);
    static SymbolicExpression createFromStream(std::istream &s);
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    void serializeToStream(std::ostream &s) const override;
    inline bool isConstantInt() const override;
    //! Get the constant value of the vertex.
    inline const T& getValue() const;

  private:

    Constant(const T& c_);
    T c;
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;
    typedef T CacheKey;
    static std::map<CacheKey, std::weak_ptr<const Constant>> cache;
};

template<>
bool Constant<double>::isConstantInt() const {
  return false;
}

template<>
bool Constant<int>::isConstantInt() const {
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
class Symbol : public Vertex, public std::enable_shared_from_this<Symbol> {
  friend SymbolicExpression;
  public:

    static IndependentVariable create(const boost::uuids::uuid& uuid_=boost::uuids::random_generator()());
    static IndependentVariable createFromStream(std::istream &s);
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    void serializeToStream(std::ostream &s) const override;
    //! Set the value of this independent variable.
    //! This has an influence on the evaluation of all ASTs which depend on this independent variable.
    inline void setValue(double x_) const;
    //! Eeturn the "version" of this variable.
    //! Each call to Symbol::set increased this "version count". This is used for caching evaluations.
    inline unsigned long getVersion() const;

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

// ***** Operation *****

//! A vertex of the AST representing an operation.
class Operation : public Vertex, public std::enable_shared_from_this<Operation> {
  friend SymbolicExpression;
  friend SymbolicExpression fmatvec::subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b);
  public:

    //! Defined operations.
    enum Operator { Plus, Minus, Mult, Div, Pow, Log, Sqrt, Neg, Sin, Cos, Tan, Sinh, Cosh, Tanh, ASin, ACos, ATan, ASinh, ACosh, ATanh, Exp, Sign, Abs };
    static SymbolicExpression create(Operator op_, const std::vector<SymbolicExpression> &child_);
    static SymbolicExpression createFromStream(std::istream &s);
    inline double eval() const override;
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    void serializeToStream(std::ostream &s) const override;

  private:

    Operation(Operator op_, const std::vector<SymbolicExpression> &child_);
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;
    Operator op;
    std::vector<SymbolicExpression> child;
    typedef std::pair<Operator, std::vector<std::weak_ptr<const Vertex>>> CacheKey;
    struct CacheKeyComp {
      bool operator()(const CacheKey& l, const CacheKey& r);
      private:
        std::owner_less<std::weak_ptr<const Vertex>> ol;
    };
    static std::map<CacheKey, std::weak_ptr<const Operation>, CacheKeyComp> cache;
    mutable double cacheValue;
    static boost::bimap<Operator, std::string> opMap;
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

} // end namespace AST
#endif

IndependentVariable& IndependentVariable::operator&=(double x) {
  static_cast<const AST::Symbol*>(get())->setValue(x);
  return *this;
}

double eval(const SymbolicExpression &x) {
  return x->eval();
}

} // end namespace fmatvec

#endif
