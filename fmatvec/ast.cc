#include "ast.h"
#include <boost/uuid/uuid_io.hpp>

using namespace std;

namespace fmatvec {

SymbolicExpression::SymbolicExpression() : shared_ptr<const AST::Vertex>(AST::Symbol::create()) {}
template<class T> SymbolicExpression::SymbolicExpression(const shared_ptr<T> &x) : shared_ptr<const AST::Vertex>(x) {}
SymbolicExpression::SymbolicExpression(int x) : shared_ptr<const AST::Vertex>(AST::Constant<int>::create(x)) {}
SymbolicExpression::SymbolicExpression(double x) : shared_ptr<const AST::Vertex>(AST::Constant<double>::create(x)) {}
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

namespace AST { // internal namespace

// ***** Vertex *****

SymbolicExpression Vertex::createUsingXML(const void *element) {
//  if(E(element)->getTagName()==AST%"Constant")
//    return Constant::createUsingXML(element);
//  else if(E(element)->getTagName()==AST%"Symbol")
//    return Symbol::createUsingXML(element);
//  else if(E(element)->getTagName()==AST%"Operation")
//    return Operation::createUsingXML(element);
//  throw runtime_error("Unknown vertex "+E(element)->getTagName().second);
  return SymbolicExpression();
}

bool Vertex::isZero() const {
  if(isConstInt() && static_cast<const Constant<int>*>(this)->getValue()==0)
    return true;
  return false;
}

bool Vertex::isOne() const {
  if(isConstInt() && static_cast<const Constant<int>*>(this)->getValue()==1)
    return true;
  return false;
}

bool Vertex::isSymbol() const {
  return dynamic_cast<const Symbol*>(this)!=nullptr;
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
SymbolicExpression Constant<T>::createUsingXML(const void *element) {
//  return create(E(element)->getText<double>());
  return SymbolicExpression();
}

template<class T>
SymbolicExpression Constant<T>::parDer(const SymbolicExpression &x) const {
  assert(dynamic_pointer_cast<const Symbol>(x));
  return Constant<int>::create(0);
}

template<class T>
void Constant<T>::writeXMLFile(ostream &parent) const {
  parent<<c;
//  E(parent)->addElementText(AST%"Constant", boost::lexical_cast<double>(c));
}

template<class T>
Constant<T>::Constant(const T& c_) : c(c_) {}

template SymbolicExpression Constant<int   >::create(const int   &c_);
template SymbolicExpression Constant<double>::create(const double&c_);
template SymbolicExpression Constant<int   >::createUsingXML(const void *element);
template SymbolicExpression Constant<double>::createUsingXML(const void *element);
template SymbolicExpression Constant<int   >::parDer(const SymbolicExpression &x) const;
template SymbolicExpression Constant<double>::parDer(const SymbolicExpression &x) const;
template void Constant<int   >::writeXMLFile(ostream &parent) const;
template void Constant<double>::writeXMLFile(ostream &parent) const;

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

SymbolicExpression Symbol::createUsingXML(const void *element) {
//  return create(boost::uuids::string_generator()(E(element)->getAttribute("uuid")));
  return SymbolicExpression();
}

SymbolicExpression Symbol::parDer(const SymbolicExpression &x) const {
  assert(dynamic_pointer_cast<const Symbol>(x));
  return this == x.get() ? Constant<int>::create(1) : Constant<int>::create(0);
}

void Symbol::writeXMLFile(ostream &parent) const {
  static map<boost::uuids::uuid, int> varName;
  auto it=varName.insert(make_pair(uuid, varName.size()+1)).first;
  parent<<"a"<<it->second;
//  DOMDocument *doc=parent->getOwnerDocument();
//  DOMElement *var=D(doc)->createElement(AST%"Symbol");
//  parent->insertBefore(var, nullptr);
//  E(var)->setAttribute("uuid", to_string(uuid));
}

void Symbol::set(double x_) {
  version++;
  x=x_;
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
    return c->isConstInt() || dynamic_pointer_cast<const Constant<double>>(c)!=nullptr;
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
  copy(child_.begin(), child_.end(), back_inserter(weakChild));
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

SymbolicExpression Operation::createUsingXML(const void *element) {
//  Operation::Operator op=static_cast<Operation::Operator>(boost::lexical_cast<int>(E(element)->getAttribute("op")));
//  vector<SymbolicExpression> arg;
//  for(auto a=element->getFirstElementChild(); a!=nullptr; a=a->getNextElementSibling())
//    arg.push_back(Vertex::createUsingXML(a));
//  return create(op, arg);
  return SymbolicExpression();
}

SymbolicExpression Operation::parDer(const SymbolicExpression &x) const {
  assert(dynamic_pointer_cast<const Symbol>(x));
  switch(op) {
    case Plus: return Operation::create(Plus, {child[0]->parDer(x), child[1]->parDer(x)});
    case Minus: return Operation::create(Minus, {child[0]->parDer(x), child[1]->parDer(x)});
    case Mult: {
      auto a=Operation::create(Mult, {child[0]->parDer(x), child[1]});
      auto b=Operation::create(Mult, {child[0], child[1]->parDer(x)});
      return Operation::create(Plus, {a, b});
    }
    case Div: {
      auto a=Operation::create(Mult, {child[0]->parDer(x), child[1]});
      auto b=Operation::create(Mult, {child[0], child[1]->parDer(x)});
      auto z=Operation::create(Minus, {a, b});
      auto n=Operation::create(Mult, {child[1], child[1]});
      return Operation::create(Div, {z, n});
    }
    case Pow: {
      auto v=shared_from_this();
      auto a=Operation::create(Mult, {child[1]->parDer(x), Operation::create(Log, {child[0]})});
      auto b=Operation::create(Mult, {Operation::create(Div, {child[1], child[0]}), child[0]->parDer(x)});
      auto c=Operation::create(Plus, {a, b});
      return Operation::create(Mult, {v, c});
    }
    case Log: return Operation::create(Div, {child[0]->parDer(x), child[0]});
    case Sqrt: {
      auto v=shared_from_this();
      auto a=Operation::create(Div, {0.5, v});
      return Operation::create(Mult, {a, child[0]->parDer(x)});
    }
  }
  throw runtime_error("Unknown operation.");
}

void Operation::writeXMLFile(ostream &parent) const {
  parent<<"("<<op;
  for(auto &c : child) {
    parent<<" ";
    c->writeXMLFile(parent);
  }
  parent<<")";
//  DOMDocument *doc=parent->getOwnerDocument();
//  DOMElement *opEle=D(doc)->createElement(AST%"Operation");
//  parent->insertBefore(opEle, nullptr);
//  E(opEle)->setAttribute("op", op);
//  for(auto &c : child)
//    c->writeXMLFile(opEle);
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

} // end namespace AST
} // end namespace fmatvec
