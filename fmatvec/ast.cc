#include "ast.h"
#include <boost/uuid/uuid_io.hpp>

using namespace std;

namespace fmatvec {

// ***** SymbolicExpression *****

#ifndef NDEBUG
unsigned long SymbolicExpression::evalOperationsCount = 0;
#endif

SymbolicExpression::SymbolicExpression() : shared_ptr<const AST::Vertex>(AST::Constant<int>::create(0)) {}
template<class T> SymbolicExpression::SymbolicExpression(const shared_ptr<T> &x) : shared_ptr<const AST::Vertex>(x) {}
SymbolicExpression::SymbolicExpression(int x) : shared_ptr<const AST::Vertex>(AST::Constant<int>::create(x)) {}
SymbolicExpression::SymbolicExpression(double x) : shared_ptr<const AST::Vertex>(AST::Constant<double>::create(x)) {}
SymbolicExpression::SymbolicExpression(ConstructSymbol) : shared_ptr<const AST::Vertex>(AST::Symbol::create()) {}

SymbolicExpression SymbolicExpression::operator+(const SymbolicExpression &b) const {
  return AST::Operation::create(AST::Operation::Plus, {*this, b});
}

SymbolicExpression SymbolicExpression::operator-(const SymbolicExpression &b) const {
  return AST::Operation::create(AST::Operation::Minus, {*this, b});
}

SymbolicExpression SymbolicExpression::operator*(const SymbolicExpression &b) const {
  return AST::Operation::create(AST::Operation::Mult, {*this, b});
}

SymbolicExpression SymbolicExpression::operator/(const SymbolicExpression &b) const {
  return AST::Operation::create(AST::Operation::Div, {*this, b});
}

SymbolicExpression& SymbolicExpression::operator+=(const SymbolicExpression &b) {
  *this=AST::Operation::create(AST::Operation::Plus, {*this, b});
  return *this;
}

SymbolicExpression& SymbolicExpression::operator-=(const SymbolicExpression &b) {
  *this=AST::Operation::create(AST::Operation::Minus, {*this, b});
  return *this;
}

SymbolicExpression& SymbolicExpression::operator*=(const SymbolicExpression &b) {
  *this=AST::Operation::create(AST::Operation::Mult, {*this, b});
  return *this;
}

SymbolicExpression& SymbolicExpression::operator/=(const SymbolicExpression &b) {
  *this=AST::Operation::create(AST::Operation::Div, {*this, b});
  return *this;
}

SymbolicExpression pow(const SymbolicExpression &a, const SymbolicExpression &b) {
  return AST::Operation::create(AST::Operation::Pow, {a, b});
}

SymbolicExpression log(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Log, {a});
}

SymbolicExpression sqrt(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Sqrt, {a});
}

ostream& operator<<(ostream& s, const SymbolicExpression& se) {
  s<<"{";
  se->serializeToStream(s);
  s<<" }";
  return s;
}

istream& operator>>(istream& s, SymbolicExpression &se) {
  char ch;
  s>>ch;
  if(ch!='{')
    throw runtime_error("The SymbolicExpression in the stream is not starting with {.");
  se=AST::Vertex::createFromStream(s);
  s>>ch;
  if(ch!='}')
    throw runtime_error("The SymbolicExpression in the stream is not ending with }.");
  return s;
}

AST::SubstHelper subst(SymbolicExpression &se, const map<IndependentVariable, IndependentVariable> &substitution) {
  return AST::SubstHelper(se, substitution);
}

istream& operator>>(istream& s, const AST::SubstHelper &sh) {
  char ch;
  s>>ch;
  if(ch!='{')
    throw runtime_error("The SymbolicExpression in the stream is not starting with {.");
  sh.se=AST::Vertex::createFromStream(s, sh.substitution);
  s>>ch;
  if(ch!='}')
    throw runtime_error("The SymbolicExpression in the stream is not ending with }.");
  return s;
}

SymbolicExpression operator+(double a, const SymbolicExpression &b) {
  return SymbolicExpression(a)+b;
}

SymbolicExpression operator-(double a, const SymbolicExpression &b) {
  return SymbolicExpression(a)-b;
}

SymbolicExpression operator*(double a, const SymbolicExpression &b) {
  return SymbolicExpression(a)*b;
}

SymbolicExpression operator/(double a, const SymbolicExpression &b) {
  return SymbolicExpression(a)/b;
}

SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep) {
  return dep->parDer(indep);
}

// ***** IndependentVariable *****

IndependentVariable::IndependentVariable() : SymbolicExpression(constructSymbol) {}

namespace AST { // internal namespace

// ***** Vertex *****

SymbolicExpression Vertex::createFromStream(istream &s,
  const std::map<IndependentVariable, IndependentVariable> &substitution) {
  string className;
  s>>className;
  if(className=="i")
    return Constant<int>::createFromStream(s, substitution);
  else if(className=="d")
    return Constant<double>::createFromStream(s, substitution);
  else if(className=="s")
    return Symbol::createFromStream(s, substitution);
  else if(className=="o")
    return Operation::createFromStream(s, substitution);
  else
    throw runtime_error("Unknown class "+className);
}

bool Vertex::isZero() const {
  if(isConstantInt() && static_cast<const Constant<int>*>(this)->getValue()==0)
    return true;
  return false;
}

bool Vertex::isOne() const {
  if(isConstantInt() && static_cast<const Constant<int>*>(this)->getValue()==1)
    return true;
  return false;
}

const map<weak_ptr<const Symbol>, unsigned long, owner_less<weak_ptr<const Symbol>>>& Vertex::getDependsOn() const {
  return dependsOn;
}

// ***** Constant *****

template<class T> map<typename Constant<T>::CacheKey, weak_ptr<const Constant<T>>> Constant<T>::cache;

template<class T>
SymbolicExpression Constant<T>::create(const T&c_) {
  auto r=cache.insert(make_pair(c_, weak_ptr<const Constant<T>>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<Constant<T>>(new Constant<T>(c_));
  r.first->second=newPtr;
  return newPtr;
}

template<class T>
SymbolicExpression Constant<T>::createFromStream(istream &s,
  const std::map<IndependentVariable, IndependentVariable> &substitution) {
  T c;
  s>>c;
  return Constant<T>::create(c);
}

template<class T>
SymbolicExpression Constant<T>::parDer(const IndependentVariable &x) const {
  return Constant<int>::create(0);
}

template<>
void Constant<int>::serializeToStream(ostream &s) const {
  s<<" i "<<c;
}

template<>
void Constant<double>::serializeToStream(ostream &s) const {
  s<<" d "<<c;
}

template<class T>
Constant<T>::Constant(const T& c_) : c(c_) {}

template SymbolicExpression Constant<int   >::create(const int   &c_);
template SymbolicExpression Constant<double>::create(const double&c_);
template SymbolicExpression Constant<int   >::createFromStream(istream &s,
  const std::map<IndependentVariable, IndependentVariable> &substitution);
template SymbolicExpression Constant<double>::createFromStream(istream &s,
  const std::map<IndependentVariable, IndependentVariable> &substitution);
template SymbolicExpression Constant<int   >::parDer(const IndependentVariable &x) const;
template SymbolicExpression Constant<double>::parDer(const IndependentVariable &x) const;
template void Constant<int   >::serializeToStream(ostream &s) const;
template void Constant<double>::serializeToStream(ostream &s) const;

// ***** Symbol *****

map<Symbol::CacheKey, weak_ptr<const Symbol>> Symbol::cache;

SymbolicExpression Symbol::create(const boost::uuids::uuid& uuid_) {
  auto r=cache.insert(make_pair(uuid_, weak_ptr<const Symbol>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<Symbol>(new Symbol(uuid_));
  newPtr->dependsOn.insert(make_pair(weak_ptr<const Symbol>(newPtr), 0)); // we cannot call this in Symbol::Symbol (it may be possible with weak_from_this() for C++17)
  r.first->second=newPtr;
  return newPtr;
}

SymbolicExpression Symbol::createFromStream(istream &s,
  const std::map<IndependentVariable, IndependentVariable> &substitution) {
  // get uuid from stream
  boost::uuids::uuid uuid;
  s>>uuid;
  // create a, maybe temporary, Symbol with this uuid
  auto indep=Symbol::create(uuid);
  // if this indep in substituted then return the substituted indep, if not return itself.
  auto it=substitution.find(static_cast<IndependentVariable&>(indep));
  if(it!=substitution.end())
    return it->second;
  return indep;
}

SymbolicExpression Symbol::parDer(const IndependentVariable &x) const {
  return this == x.get() ? Constant<int>::create(1) : Constant<int>::create(0);
}

void Symbol::serializeToStream(ostream &s) const {
  s<<" s "<<uuid;
}

Symbol::Symbol(const boost::uuids::uuid& uuid_) : version(0), uuid(uuid_) {}

// ***** Operation *****

map<Operation::CacheKey, weak_ptr<const Operation>, Operation::CacheKeyComp> Operation::cache;

SymbolicExpression Operation::create(Operator op_, const vector<SymbolicExpression> &child_) {
  // optimize Plus
  if(op_==Plus && child_[0]->isZero())
    return child_[1];
  if(op_==Plus && child_[1]->isZero())
    return child_[0];
  // optimize Minus
  if(op_==Minus && child_[1]->isZero())
    return child_[0];
  // optimize Mult
  if(op_==Mult && (child_[0]->isZero() || child_[1]->isZero()))
    return Constant<int>::create(0);
  if(op_==Mult && child_[0]->isOne())
    return child_[1];
  if(op_==Mult && child_[1]->isOne())
    return child_[0];
  // optimize Div
  if(op_==Div && child_[1]->isOne())
    return child_[0];
  if(op_==Div && child_[0]->isZero())
    return Constant<int>::create(0);
  // optimize Pow
  if(op_==Pow && child_[1]->isZero())
    return Constant<int>::create(1);
  if(op_==Pow && child_[1]->isOne())
    return child_[0];
  // optimize Constant arguments (if ALL are Constant)
  bool allConst=all_of(child_.begin(), child_.end(), [](const SymbolicExpression& c) {
    return c->isConstantInt() || dynamic_pointer_cast<const Constant<double>>(c)!=nullptr;
  });
  if(allConst) {
    Operation op(op_, child_);
    double doubleValue=op.eval();
    int intValue=lround(doubleValue);
    if(abs(intValue-doubleValue)<2*numeric_limits<double>::epsilon())
      return Constant<int>::create(intValue);
    return Constant<double>::create(doubleValue);
  }

  vector<weak_ptr<const Vertex>> weakChild;
  // we cannot just use std::copy here since std::shared_ptr is a private base of SymbolicExpression
  transform(child_.begin(), child_.end(), back_inserter(weakChild), [](const SymbolicExpression &x){
    return weak_ptr<const Vertex>(x);
  });
  auto r=cache.insert(make_pair(make_pair(op_, weakChild), weak_ptr<const Operation>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<Operation>(new Operation(op_, child_));
  r.first->second=newPtr;
  return newPtr;
}

SymbolicExpression Operation::createFromStream(istream &s,
  const std::map<IndependentVariable, IndependentVariable> &substitution) {
  int opInt;
  size_t size;
  vector<SymbolicExpression> child;
  s>>opInt>>size;
  for(size_t i=0; i<size; ++i)
    child.push_back(Vertex::createFromStream(s, substitution));
  return Operation::create(static_cast<Operator>(opInt), child);
}

SymbolicExpression Operation::parDer(const IndependentVariable &x) const {
  #define v SymbolicExpression(shared_from_this()) // expression to be differentiated
  #define a child[0] // expression of first argument
  #define b child[1] // expression of second argument
  #define c child[2] // expression of third argument
  #define ad child[0]->parDer(x) // pertial derivative of the expression of the first argument wrt the independent x
  #define bd child[1]->parDer(x) // pertial derivative of the expression of the second argument wrt the independent x
  #define cd child[2]->parDer(x) // pertial derivative of the expression of the third argument wrt the independent x
  switch(op) {
    case Plus:
      return ad + bd;
    case Minus:
      return ad - bd;
    case Mult:
      return ad * b + a * bd;
    case Div:
      return (ad * b - a * bd) / (b * b);
    case Pow:
      return v * (bd * log(a) + b / a * ad);
    case Log:
      return ad / a;
    case Sqrt:
      return ad * 0.5 / v;
  }
  #undef a
  #undef b
  #undef c
  #undef ad
  #undef bd
  #undef cd
  throw runtime_error("Unknown operation.");
}

void Operation::serializeToStream(ostream &s) const {
  s<<" o "<<op<<" "<<child.size();
  for(auto &c : child)
    c->serializeToStream(s);
}

Operation::Operation(Operator op_, const vector<SymbolicExpression> &child_) : op(op_), child(child_) {
  for(auto &c : child)
    dependsOn.insert(c->getDependsOn().begin(), c->getDependsOn().end());
}

bool Operation::CacheKeyComp::operator()(const CacheKey& l, const CacheKey& r) {
  if(l.first!=r.first)
    return l.first<r.first;
  if(l.second.size()!=r.second.size())
    return l.second.size()<r.second.size();
  for(size_t i=0; i<l.second.size(); ++i)
    if(ol(l.second[i],r.second[i]) || ol(r.second[i],l.second[i]))
      return ol(l.second[i], r.second[i]);
  return false;
}

SubstHelper::SubstHelper(SymbolicExpression &se_, const std::map<IndependentVariable, IndependentVariable> &substitution_)
  : se(se_), substitution(substitution_) {}

} // end namespace AST
} // end namespace fmatvec
