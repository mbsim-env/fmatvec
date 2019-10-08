#ifndef _FMATVEC_AST_H_
#define _FMATVEC_AST_H_

#include <vector>
#include <map>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

namespace fmatvec {

namespace AST {
  class Vertex;
}

class Expr : private std::shared_ptr<const AST::Vertex> {
  public:
    Expr();
    Expr(double x);
    Expr(const std::shared_ptr<const AST::Vertex> &x);
    const AST::Vertex* operator->();
    Expr operator+(const Expr &b) const;
    Expr operator-(const Expr &b) const;
    Expr operator*(const Expr &b) const;
    Expr operator/(const Expr &b) const;
    Expr& operator+=(const Expr &b);
    Expr& operator-=(const Expr &b);
    Expr& operator*=(const Expr &b);
    Expr& operator/=(const Expr &b);
};

namespace AST { // internal namespace

class Var;

// ***** Vertex *****

//! A abstract class for a Vertex of the AST (abstract syntax tree).
class Vertex {
  public:

    //! Create a AST from a XML representation (deserialization).
    static std::shared_ptr<const Vertex> createUsingXML(const void *element);
    //! Evaluate the AST.
    //! The returned value depends on the AST and on the current values of all independent variables this AST depends on.
    //! Also see Var ant Vertex::getDependsOn().
    virtual double eval() const=0;
    //! Generate a new AST being the partial derivate of this AST with respect to the variable x.
    virtual std::shared_ptr<const Vertex> parDer(const std::shared_ptr<const Var> &x) const=0;
    //! Write a AST to a XML representation (serialization).
    virtual void writeXMLFile(void *parent) const=0;
    //! Returns true if this Vertex is a constant with value 0.
    //! Note that only integer constants can be 0. Double constants are never treated as 0 by this function.
    bool isZero() const;
    //! Returns true if this Vertex is a constant with value 1.
    //! Note that only integer constants can be 1. Double constants are never treated as 1 by this function.
    bool isOne() const;
    //! Return a list of all variables this AST depends on.
    const std::map<std::weak_ptr<const Var>, unsigned long, std::owner_less<std::weak_ptr<const Var>>>& getDependsOn() const;

  protected:

    // the value of this map (the unsigned long = version) must be mutable
    mutable std::map<std::weak_ptr<const Var>, unsigned long, std::owner_less<std::weak_ptr<const Var>>> dependsOn;
};

// ***** Const *****

//! A vertex of the AST representing a constant (int or double)
template<class T>
class Const : public Vertex {
  public:

    static std::shared_ptr<const Const> create(const T& c_);
    static std::shared_ptr<const Const> createUsingXML(const void *element);
    inline double eval() const override;
    std::shared_ptr<const Vertex> parDer(const std::shared_ptr<const Var> &x) const override;
    void writeXMLFile(void *parent) const override;
    //! Get the constant value of the vertex.
    inline const T& getValue() const;

  private:

    Const(const T& c_);
    T c;
    typedef T CacheKey;
    static std::map<CacheKey, std::weak_ptr<const Const>> cache;
};

template<class T>
double Const<T>::eval() const {
  return c;
}

template<class T>
const T& Const<T>::getValue() const {
  return c;
}

// ***** Var *****

//! A vertex of the AST representing a independent variable.
class Var : public Vertex {
  public:

    static std::shared_ptr<Var> create(const boost::uuids::uuid& uuid_=boost::uuids::random_generator()());
    static std::shared_ptr<Var> createUsingXML(const void *element);
    inline double eval() const override;
    std::shared_ptr<const Vertex> parDer(const std::shared_ptr<const Var> &x) const override;
    void writeXMLFile(void *parent) const override;
    //! Set the value of this independent variable.
    //! This has an influence on the evaluation of all ASTs which depend on this independent variable.
    void set(double x_);
    //! Eeturn the "version" of this variable.
    //! Each call to Var::set increased this "version count". This is used for caching evaluations.
    inline unsigned long getVersion() const;

  private:

    Var(const boost::uuids::uuid& uuid_);
    double x;
    unsigned long version;
    boost::uuids::uuid uuid; // each variable has a uuid (this is only used when the AST is serialized to XML)
    typedef boost::uuids::uuid CacheKey;
    static std::map<CacheKey, std::weak_ptr<Var>> cache;
};

double Var::eval() const {
  return x;
}

unsigned long Var::getVersion() const {
  return version;
}

// ***** Op *****

//! A vertex of the AST representing an operation.
class Op : public Vertex, public std::enable_shared_from_this<Op> {
  public:

    //! Defined operations.
    enum Operator { Plus, Minus, Mult, Div, Pow, Log };
    static std::shared_ptr<const Vertex> create(Operator op_, const std::vector<std::shared_ptr<const Vertex>> &child_);
    static std::shared_ptr<const Vertex> createUsingXML(const void *element);
    inline double eval() const override;
    std::shared_ptr<const Vertex> parDer(const std::shared_ptr<const Var> &x) const override;
    void writeXMLFile(void *parent) const override;

  private:

    Op(Operator op_, const std::vector<std::shared_ptr<const Vertex>> &child_);
    Operator op;
    std::vector<std::shared_ptr<const Vertex>> child;
    typedef std::pair<Operator, std::vector<std::weak_ptr<const Vertex>>> CacheKey;
    struct CacheKeyComp {
      bool operator()(const CacheKey& l, const CacheKey& r);
      private:
        std::owner_less<std::weak_ptr<const Vertex>> ol;
    };
    static std::map<CacheKey, std::weak_ptr<const Op>, CacheKeyComp> cache;
    mutable double cacheValue;
};

double Op::eval() const {
  if(!dependsOn.empty()) {
    // we eval expressions having no dependency always (such expressions are optimized out during Op creation)
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

  switch(op) {
    case Plus: cacheValue = child[0]->eval() + child[1]->eval(); return cacheValue;
    case Minus: cacheValue = child[0]->eval() - child[1]->eval(); return cacheValue;
    case Mult: cacheValue = child[0]->eval() * child[1]->eval(); return cacheValue;
    case Div: cacheValue = child[0]->eval() / child[1]->eval(); return cacheValue;
    case Pow:
      if(auto cci=std::dynamic_pointer_cast<const Const<int>>(child[1]))
        cacheValue = std::pow(child[0]->eval(), cci->getValue());
      else
        cacheValue = std::pow(child[0]->eval(), child[1]->eval());
      return cacheValue;
    case Log: cacheValue = std::log(child[0]->eval()); return cacheValue;
  }
  throw std::runtime_error("Unknown operation.");
}

} // end namespace AST
} // end namespace fmatvec

#endif
