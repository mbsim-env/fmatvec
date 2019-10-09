#include "ast.h"
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace fmatvec {

Expr::Expr() : std::shared_ptr<const AST::Vertex>(AST::Var::create()) {}
Expr::Expr(double x) : std::shared_ptr<const AST::Vertex>(AST::Const<double>::create(x)) {}
Expr::Expr(const std::shared_ptr<const AST::Vertex> &x) : std::shared_ptr<const AST::Vertex>(x) {}
const AST::Vertex* Expr::operator->() { return get(); }
Expr Expr::operator+(const Expr &b) const { return AST::Op::create(AST::Op::Plus, {*this, b});}
Expr Expr::operator-(const Expr &b) const { return AST::Op::create(AST::Op::Minus, {*this, b}); }
Expr Expr::operator*(const Expr &b) const { return AST::Op::create(AST::Op::Mult, {*this, b}); }
Expr Expr::operator/(const Expr &b) const { return AST::Op::create(AST::Op::Div, {*this, b}); }
Expr& Expr::operator+=(const Expr &b) { *this=AST::Op::create(AST::Op::Plus, {*this, b}); return *this; }
Expr& Expr::operator-=(const Expr &b) { *this=AST::Op::create(AST::Op::Minus, {*this, b}); return *this; }
Expr& Expr::operator*=(const Expr &b) { *this=AST::Op::create(AST::Op::Mult, {*this, b}); return *this; }
Expr& Expr::operator/=(const Expr &b) { *this=AST::Op::create(AST::Op::Div, {*this, b}); return *this; }

namespace AST { // internal namespace

// ***** Vertex *****

shared_ptr<const Vertex> Vertex::createUsingXML(const void *element) {
//mfmf  if(E(element)->getTagName()==AST%"Const")
//mfmf    return Const::createUsingXML(element);
//mfmf  else if(E(element)->getTagName()==AST%"Var")
//mfmf    return Var::createUsingXML(element);
//mfmf  else if(E(element)->getTagName()==AST%"Op")
//mfmf    return Op::createUsingXML(element);
//mfmf  throw runtime_error("Unknown vertex "+E(element)->getTagName().second);
  return shared_ptr<const Vertex>();
}

bool Vertex::isZero() const {
  if(isConstInt() && static_cast<const Const<int>*>(this)->getValue()==0)
    return true;
  return false;
}

bool Vertex::isOne() const {
  if(isConstInt() && static_cast<const Const<int>*>(this)->getValue()==1)
    return true;
  return false;
}

const map<weak_ptr<const Var>, unsigned long, owner_less<weak_ptr<const Var>>>& Vertex::getDependsOn() const {
  return dependsOn;
}

// ***** Const *****

template<class T> map<typename Const<T>::CacheKey, weak_ptr<const Const<T>>> Const<T>::cache;

template<class T>
shared_ptr<const Const<T>> Const<T>::create(const T&c_) {
  auto r=cache.insert(make_pair(c_, weak_ptr<const Const<T>>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<const Const<T>>(new Const<T>(c_));
  r.first->second=newPtr;
  return newPtr;
}

template<class T>
shared_ptr<const Const<T>> Const<T>::createUsingXML(const void *element) {
//mfmf  return create(E(element)->getText<double>());
  return shared_ptr<const Const<T>>();
}

template<class T>
shared_ptr<const Vertex> Const<T>::parDer(const shared_ptr<const Var> &x) const {
  return Const<T>::create(0);
}

template<class T>
void Const<T>::writeXMLFile(void *parent) const {
//mfmf  E(parent)->addElementText(AST%"Const", boost::lexical_cast<double>(c));
}

template<class T>
Const<T>::Const(const T& c_) : c(c_) {}

template shared_ptr<const Const<int   >> Const<int   >::create(const int   &c_);
template shared_ptr<const Const<double>> Const<double>::create(const double&c_);
template shared_ptr<const Const<int   >> Const<int   >::createUsingXML(const void *element);
template shared_ptr<const Const<double>> Const<double>::createUsingXML(const void *element);
template shared_ptr<const Vertex> Const<int   >::parDer(const shared_ptr<const Var> &x) const;
template shared_ptr<const Vertex> Const<double>::parDer(const shared_ptr<const Var> &x) const;
template void Const<int   >::writeXMLFile(void *parent) const;
template void Const<double>::writeXMLFile(void *parent) const;

// ***** Var *****

map<Var::CacheKey, weak_ptr<Var>> Var::cache;

shared_ptr<Var> Var::create(const boost::uuids::uuid& uuid_) {
  auto r=cache.insert(make_pair(uuid_, weak_ptr<Var>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<Var>(new Var(uuid_));
  newPtr->dependsOn.insert(make_pair(newPtr, 0)); // we cannot call this in Var::Var (it may be possible with weak_from_this() for C++17)
  r.first->second=newPtr;
  return newPtr;
}

shared_ptr<Var> Var::createUsingXML(const void *element) {
//mfmf  return create(boost::uuids::string_generator()(E(element)->getAttribute("uuid")));
  return shared_ptr<Var>();
}

shared_ptr<const Vertex> Var::parDer(const shared_ptr<const Var> &x) const {
  return this == x.get() ? Const<int>::create(1) : Const<int>::create(0);
}

void Var::writeXMLFile(void *parent) const {
//mfmf  DOMDocument *doc=parent->getOwnerDocument();
//mfmf  DOMElement *var=D(doc)->createElement(AST%"Var");
//mfmf  parent->insertBefore(var, nullptr);
//mfmf  E(var)->setAttribute("uuid", to_string(uuid));
}

void Var::set(double x_) {
  version++;
  x=x_;
}

Var::Var(const boost::uuids::uuid& uuid_) : version(0), uuid(uuid_) {}

// ***** Op *****

map<Op::CacheKey, weak_ptr<const Op>, Op::CacheKeyComp> Op::cache;

shared_ptr<const Vertex> Op::create(Operator op_, const vector<shared_ptr<const Vertex>> &child_) {
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
    return Const<int>::create(0);
  if(op_==Mult && child_[0]->isOne())
    return child_[1];
  if(op_==Mult && child_[1]->isOne())
    return child_[0];
  // optimize Div
  if(op_==Div && child_[1]->isOne())
    return child_[0];
  // optimize Pow
  if(op_==Pow && child_[1]->isZero())
    return Const<int>::create(1);
  if(op_==Pow && child_[1]->isOne())
    return child_[0];
  // optimize Const arguments (if ALL are Const)
  bool allConst=all_of(child_.begin(), child_.end(), [](const shared_ptr<const Vertex>& c) {
    return c->isConstInt() || dynamic_pointer_cast<const Const<double>>(c)!=nullptr;
  });
  if(allConst) {
    Op op(op_, child_);
    double doubleValue=op.eval();
    int intValue=lround(doubleValue);
    if(abs(intValue-doubleValue)<2*numeric_limits<double>::epsilon())
      return Const<int>::create(intValue);
    return Const<double>::create(doubleValue);
  }

  vector<weak_ptr<const Vertex>> weakChild;
  copy(child_.begin(), child_.end(), back_inserter(weakChild));
  auto r=cache.insert(make_pair(make_pair(op_, weakChild), weak_ptr<const Op>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<const Op>(new Op(op_, child_));
  r.first->second=newPtr;
  return newPtr;
}

shared_ptr<const Vertex> Op::createUsingXML(const void *element) {
//mfmf  Op::Operator op=static_cast<Op::Operator>(boost::lexical_cast<int>(E(element)->getAttribute("op")));
//mfmf  vector<shared_ptr<const Vertex>> arg;
//mfmf  for(auto a=element->getFirstElementChild(); a!=nullptr; a=a->getNextElementSibling())
//mfmf    arg.push_back(Vertex::createUsingXML(a));
//mfmf  return create(op, arg);
  return shared_ptr<const Vertex>();
}

shared_ptr<const Vertex> Op::parDer(const shared_ptr<const Var> &x) const {
  switch(op) {
    case Plus: return Op::create(Plus, {child[0]->parDer(x), child[1]->parDer(x)});
    case Minus: return Op::create(Minus, {child[0]->parDer(x), child[1]->parDer(x)});
    case Mult: {
      auto a=Op::create(Mult, {child[0]->parDer(x), child[1]});
      auto b=Op::create(Mult, {child[0], child[1]->parDer(x)});
      return Op::create(Plus, {a, b});
    }
    case Div: {
      auto a=Op::create(Mult, {child[0]->parDer(x), child[1]});
      auto b=Op::create(Mult, {child[0], child[1]->parDer(x)});
      auto z=Op::create(Minus, {a, b});
      auto n=Op::create(Mult, {child[1], child[1]});
      return Op::create(Div, {z, n});
    }
    case Pow: {
      auto v=shared_from_this();
      auto a=Op::create(Mult, {child[1]->parDer(x), Op::create(Log, {child[0]})});
      auto b=Op::create(Mult, {Op::create(Div, {child[1], child[0]}), child[0]->parDer(x)});
      auto c=Op::create(Plus, {a, b});
      return Op::create(Mult, {v, c});
    }
    case Log: return Op::create(Div, {child[0]->parDer(x), child[0]});
  }
  throw runtime_error("Unknown operation.");
}

void Op::writeXMLFile(void *parent) const {
//mfmf  DOMDocument *doc=parent->getOwnerDocument();
//mfmf  DOMElement *opEle=D(doc)->createElement(AST%"Op");
//mfmf  parent->insertBefore(opEle, nullptr);
//mfmf  E(opEle)->setAttribute("op", op);
//mfmf  for(auto &c : child)
//mfmf    c->writeXMLFile(opEle);
}

Op::Op(Operator op_, const vector<shared_ptr<const Vertex>> &child_) : op(op_), child(child_) {
  for(auto &c : child)
    dependsOn.insert(c->getDependsOn().begin(), c->getDependsOn().end());
}

bool Op::CacheKeyComp::operator()(const CacheKey& l, const CacheKey& r) {
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
