#ifndef _FMATVEC_AST_H_
#define _FMATVEC_AST_H_

#include <vector>
#include <map>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

namespace fmatvec {

// forward declarations

class Symbol;

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
//! If it represents a symbol, which can be checked by calling this->isSymbol(),
//! than any other expression can be differentiated with respect to this symbol.
//! See the public member function of fmatvec::AST::Vertex for all function being
//! callable as this->memberFunc(...).
class SymbolicExpression : private std::shared_ptr<const AST::Vertex> {
  friend class AST::Symbol;
  friend class AST::Operation;
  friend class AST::Constant<int>;
  friend class AST::Constant<double>;
  friend SymbolicExpression parDer(const SymbolicExpression &dep, const SymbolicExpression &indep);
  private:
    template<class T> SymbolicExpression(const shared_ptr<T> &x);
  public:
    //! Creates the value 0.
    SymbolicExpression();
    //! Creates a integer constant.
    SymbolicExpression(int x);
    //! Creates a double constant.
    SymbolicExpression(double x);
    //! Creates a symbol (each call to this ctor creates a new symbol)
    SymbolicExpression(const Symbol&);

    // default ctors and assignment operators
    SymbolicExpression(const SymbolicExpression& x) = default;
    SymbolicExpression(SymbolicExpression&& x) = default;
    SymbolicExpression& operator=(const SymbolicExpression& x) = default;
    SymbolicExpression& operator=(SymbolicExpression&& x) = default;
    ~SymbolicExpression() = default;

    // member functions from std::shared_ptr
    using std::shared_ptr<const AST::Vertex>::reset;
    using std::shared_ptr<const AST::Vertex>::swap;
    using std::shared_ptr<const AST::Vertex>::use_count;
    using std::shared_ptr<const AST::Vertex>::owner_before;

    //! Set the double value of the symbol stored in the SymbolicExpression.
    //! If this object does not store a symbol an exception is thrown.
    inline void setSymbolValue(double x);

    //! Evaluate the SymbolicExpression.
    //! The returned value depends on the symbolic expression and on the current values of all independent
    //! variables this symbolic expression depends on.
    //! Also see Symbol ant Vertex::getDependsOn().
    inline double eval() const;

    //! Write a SymbolicExpression to a XML representation (serialization).
    void writeXMLFile(std::ostream &parent) const;

    // operators
    SymbolicExpression operator+(const SymbolicExpression &b) const;
    SymbolicExpression operator-(const SymbolicExpression &b) const;
    SymbolicExpression operator*(const SymbolicExpression &b) const;
    SymbolicExpression operator/(const SymbolicExpression &b) const;
    SymbolicExpression& operator+=(const SymbolicExpression &b);
    SymbolicExpression& operator-=(const SymbolicExpression &b);
    SymbolicExpression& operator*=(const SymbolicExpression &b);
    SymbolicExpression& operator/=(const SymbolicExpression &b);

#ifndef NDEBUG
    static unsigned long evalOperationsCount;
#endif
};

// Some member function definition of vSymbolicExpression are moved ot the end of this file
// since they need the defintion of the other class defined in this file.

//! Generate a new SymbolicExpression being the partial derivate of dep
//! with respect to indep (indep must be a symbol).
SymbolicExpression parDer(const SymbolicExpression &dep, const SymbolicExpression &indep);

// function operations overloaded for SymbolicExpression
SymbolicExpression pow(const SymbolicExpression &a, const SymbolicExpression &b);
SymbolicExpression log(const SymbolicExpression &a);
SymbolicExpression sqrt(const SymbolicExpression &a);

#ifndef NDEBUG
unsigned long SymbolicExpression::evalOperationsCount = 0;
#endif

namespace AST { // internal namespace

// ***** Vertex *****

//! A abstract class for a Vertex of the AST (abstract syntax tree).
class Vertex {
  public:

    //! Create a AST from a XML representation (deserialization).
    static SymbolicExpression createUsingXML(const void *element);
    //! Evaluate the AST.
    //! The returned value depends on the AST and on the current values of all independent variables this AST depends on.
    //! Also see Symbol ant Vertex::getDependsOn().
    virtual double eval() const=0;
    //! Generate a new AST being the partial derivate of this AST with respect to the variable x.
    virtual SymbolicExpression parDer(const SymbolicExpression &x) const=0;
    //! Write a AST to a XML representation (serialization).
    virtual void writeXMLFile(std::ostream &parent) const=0;
    //! Rreturn true if this Vertex is a constant integer.
    inline virtual bool isConstantInt() const;
    //! Rreturn true if this Vertex is a symbol.
    inline virtual bool isSymbol() const;
    //! Returns true if this Vertex is a constant with value 0.
    //! Note that only integer constants can be 0. Double constants are never treated as 0 by this function.
    bool isZero() const;
    //! Returns true if this Vertex is a constant with value 1.
    //! Note that only integer constants can be 1. Double constants are never treated as 1 by this function.
    bool isOne() const;
    //! Return a list of all variables this AST depends on.
    const std::map<std::weak_ptr<const Symbol>, unsigned long, std::owner_less<std::weak_ptr<const Symbol>>>& getDependsOn() const;

  protected:

    // the value of this map (the unsigned long = version) must be mutable
    mutable std::map<std::weak_ptr<const Symbol>, unsigned long, std::owner_less<std::weak_ptr<const Symbol>>> dependsOn;
};

bool Vertex::isConstantInt() const {
  return false;
}

bool Vertex::isSymbol() const {
  return false;
}

// ***** Constant *****

//! A vertex of the AST representing a constant (int or double)
template<class T>
class Constant : public Vertex {
  public:

    static SymbolicExpression create(const T& c_);
    static SymbolicExpression createUsingXML(const void *element);
    inline double eval() const override;
    SymbolicExpression parDer(const SymbolicExpression &x) const override;
    void writeXMLFile(std::ostream &parent) const override;
    inline bool isConstantInt() const override;
    //! Get the constant value of the vertex.
    inline const T& getValue() const;

  private:

    Constant(const T& c_);
    T c;
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
class Symbol : public Vertex {
  public:

    static SymbolicExpression create(const boost::uuids::uuid& uuid_=boost::uuids::random_generator()());
    static SymbolicExpression createUsingXML(const void *element);
    inline double eval() const override;
    SymbolicExpression parDer(const SymbolicExpression &x) const override;
    void writeXMLFile(std::ostream &parent) const override;
    inline bool isSymbol() const override;
    //! Set the value of this independent variable.
    //! This has an influence on the evaluation of all ASTs which depend on this independent variable.
    inline void setValue(double x_);
    //! Eeturn the "version" of this variable.
    //! Each call to Symbol::set increased this "version count". This is used for caching evaluations.
    inline unsigned long getVersion() const;

  private:

    Symbol(const boost::uuids::uuid& uuid_);
    double x = 0.0;
    unsigned long version;
    boost::uuids::uuid uuid; // each variable has a uuid (this is only used when the AST is serialized to XML)
    typedef boost::uuids::uuid CacheKey;
    static std::map<CacheKey, std::weak_ptr<const Symbol>> cache;
};

void Symbol::setValue(double x_) {
  version++;
  x=x_;
}

double Symbol::eval() const {
  return x;
}

unsigned long Symbol::getVersion() const {
  return version;
}

bool Symbol::isSymbol() const {
  return true;
}

// ***** Operation *****

//! A vertex of the AST representing an operation.
class Operation : public Vertex, public std::enable_shared_from_this<Operation> {
  public:

    //! Defined operations.
    enum Operator { Plus, Minus, Mult, Div, Pow, Log, Sqrt };
    static SymbolicExpression create(Operator op_, const std::vector<SymbolicExpression> &child_);
    static SymbolicExpression createUsingXML(const void *element);
    inline double eval() const override;
    SymbolicExpression parDer(const SymbolicExpression &x) const override;
    void writeXMLFile(std::ostream &parent) const override;

  private:

    Operation(Operator op_, const std::vector<SymbolicExpression> &child_);
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
#ifndef NDEBUG
  SymbolicExpression::evalOperationsCount++;
#endif

  #define a child[0]->eval() // value of first argument
  #define b child[1]->eval() // value of second argument
  #define c child[2]->eval() // value of third argument
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
  }
  #undef a
  #undef b
  #undef c
  assert(0);
}

} // end namespace AST

void SymbolicExpression::setSymbolValue(double x) {
  if(!get()->isSymbol())
    throw std::runtime_error("Cannot call setSymbolValue on a object being not a symbol.");
  const_cast<AST::Symbol*>(static_cast<const AST::Symbol*>(get()))->setValue(x);
}

double SymbolicExpression::eval() const {
  return get()->eval();
}

} // end namespace fmatvec

#endif
