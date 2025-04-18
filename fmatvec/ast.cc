#include "ast.h"
#include "stream_impl.h"
#include <boost/math_fwd.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/spirit/include/qi.hpp>
#include <regex>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/bind/bind_function.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/format.hpp>

using namespace std;

namespace fmatvec {

// ***** SymbolicExpression *****

#if defined(FMATVEC_DEBUG) && !defined(SWIG)
signed long SymbolicExpression::evalOperationsCount = 0;
#endif

namespace {
  SymbolicExpression stringCTorHelper(const string &str) {
    istringstream s(str);
    SymbolicExpression ret;
    s>>ret;
    return ret;
  }
}
SymbolicExpression::SymbolicExpression() : shared_ptr<const AST::Vertex>(AST::Constant<long>::create(0)) {}
SymbolicExpression::SymbolicExpression(const string &str) : SymbolicExpression(stringCTorHelper(str)) {}
template<class T> SymbolicExpression::SymbolicExpression(const shared_ptr<T> &x) : shared_ptr<const AST::Vertex>(x) {}
SymbolicExpression::SymbolicExpression(int x) : shared_ptr<const AST::Vertex>(AST::Constant<long>::create(x)) {}
SymbolicExpression::SymbolicExpression(long x) : shared_ptr<const AST::Vertex>(AST::Constant<long>::create(x)) {}
SymbolicExpression::SymbolicExpression(double x) : shared_ptr<const AST::Vertex>(AST::Constant<double>::create(x)) {}
SymbolicExpression::SymbolicExpression(ConstructSymbol) : shared_ptr<const AST::Vertex>(AST::Symbol::create()) {}

template<class T>
void erase(T &cache) {
  for(auto it=cache.begin(), endIt=cache.end(); it!=endIt; )
    if(it->second.expired()) it=cache.erase(it); else ++it;
}

void SymbolicExpression::garbageCollect() {
  erase(AST::Constant<long>::cache);
  erase(AST::Constant<double>::cache);
  erase(AST::Symbol::cache);
  erase(AST::Operation::cache);
}

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

SymbolicExpression& SymbolicExpression::operator<<=(const SymbolicExpression &src) {
  *this=src;
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

SymbolicExpression SymbolicExpression::operator-() const {
  return AST::Operation::create(AST::Operation::Neg, {*this});
}

SymbolicExpression sin(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Sin, {a});
}

SymbolicExpression cos(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Cos, {a});
}

SymbolicExpression tan(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Tan, {a});
}

SymbolicExpression sinh(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Sinh, {a});
}

SymbolicExpression cosh(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Cosh, {a});
}

SymbolicExpression tanh(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Tanh, {a});
}

SymbolicExpression asin(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::ASin, {a});
}

SymbolicExpression acos(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::ACos, {a});
}

SymbolicExpression atan(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::ATan, {a});
}

SymbolicExpression atan2(const SymbolicExpression &y, const SymbolicExpression &x) {
  return AST::Operation::create(AST::Operation::ATan2, {y, x});
}

SymbolicExpression asinh(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::ASinh, {a});
}

SymbolicExpression acosh(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::ACosh, {a});
}

SymbolicExpression atanh(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::ATanh, {a});
}

SymbolicExpression exp(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Exp, {a});
}

SymbolicExpression sign(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Sign, {a});
}

SymbolicExpression heaviside(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Heaviside, {a});
}

SymbolicExpression abs(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Abs, {a});
}

SymbolicExpression min(const SymbolicExpression &a, const SymbolicExpression &b) {
  return AST::Operation::create(AST::Operation::Min, {a, b});
}

SymbolicExpression max(const SymbolicExpression &a, const SymbolicExpression &b) {
  return AST::Operation::create(AST::Operation::Max, {a, b});
}

SymbolicExpression condition(const SymbolicExpression &c, const SymbolicExpression &gt, const SymbolicExpression &lt) {
  return AST::Operation::create(AST::Operation::Condition, {c, gt, lt});
}

ostream& operator<<(ostream& s, const SymbolicExpression& se) {
  namespace karma = boost::spirit::karma;
  using It = std::ostream_iterator<char>;

  karma::rule<It, SymbolicExpression()>& vertex = getBoostSpiritKarmaRule<SymbolicExpression>();

  if(!karma::generate(It(s), vertex, se))
    throw runtime_error("Failed to write SymbolicExpression to stream");
  return s;
}

#ifndef NDEBUG // FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID
// when in a debug build the envvar FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID is set then a Symbol is not
// serialized by it uuid but by a process global integer id. This way one can generate the same Symbol id
// in the serialized output. This is quite usefull to write tests.
// This envvar should NOT be set in normal program. It will generate wrong results if more than one
// process in involved.
static map<boost::uuids::uuid, size_t> mapUUIDInt;
#endif

namespace {
#ifndef NDEBUG // FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID
  IndependentVariable createSymbolByInt(int id) {
    auto itm=find_if(mapUUIDInt.begin(), mapUUIDInt.end(), [id](const pair<boost::uuids::uuid, int> &x){
      return id==x.second;
    });
    if(itm!=mapUUIDInt.end())
      return AST::Symbol::create(itm->first);
    else
      return AST::Symbol::create(boost::uuids::random_generator()());
  }
#endif

  IndependentVariable createSymbolByVec(const vector<char> &uuidStrVec) {
    return AST::Symbol::create(boost::lexical_cast<boost::uuids::uuid>(std::string(uuidStrVec.begin(), uuidStrVec.end())));
  }
}
istream& operator>>(istream& s, SymbolicExpression &se) {
  namespace qi = boost::spirit::qi;
  using It = boost::spirit::istream_iterator;

  qi::rule<It, SymbolicExpression()>& vertex = getBoostSpiritQiRule<SymbolicExpression>();

  auto savedFlags=s.flags();
  s.unsetf(std::ios::skipws);
  BOOST_SCOPE_EXIT(&s, &savedFlags) {
    s.flags(savedFlags);
  } BOOST_SCOPE_EXIT_END

  if(!qi::parse(It(s), It(), vertex, se))
    throw runtime_error("The stream does not contain a valid SymbolicExpression. Not parsed content of stream:\n"+
                        string(istreambuf_iterator<char>(s), istreambuf_iterator<char>()));
  return s;
}

istream& operator>>(istream& s, IndependentVariable &v) {
  SymbolicExpression se;
  s>>se;
  auto sym=dynamic_pointer_cast<const AST::Symbol>(se);
  if(!sym)
    throw runtime_error("The stream does not contain a IndependentVariable.");
  v=IndependentVariable(sym);
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

SymbolicExpression operator+(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)+b;
}
SymbolicExpression operator+(long a, const SymbolicExpression &b) {
  return SymbolicExpression(a)+b;
}

SymbolicExpression operator-(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)-b;
}
SymbolicExpression operator-(long a, const SymbolicExpression &b) {
  return SymbolicExpression(a)-b;
}

SymbolicExpression operator*(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)*b;
}
SymbolicExpression operator*(long a, const SymbolicExpression &b) {
  return SymbolicExpression(a)*b;
}

SymbolicExpression operator/(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)/b;
}
SymbolicExpression operator/(long a, const SymbolicExpression &b) {
  return SymbolicExpression(a)/b;
}

SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep) {
  return dep->parDer(indep);
}

#ifdef _MSC_VER
#ifndef SWIG
const SymbolicExpression::ConstructSymbol SymbolicExpression::constructSymbol{}; // just used for tag dispatching
#endif
#endif

// ***** IndependentVariable *****

IndependentVariable::IndependentVariable() : SymbolicExpression(constructSymbol) {}

IndependentVariable::IndependentVariable(const string &str) : SymbolicExpression(str) {}

IndependentVariable::IndependentVariable(const shared_ptr<const AST::Symbol> &x) : SymbolicExpression(x) {}

namespace AST { // internal namespace

SymbolicExpression substScalar(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b) {
  // if se is the indep a to be replaced -> return b
  if(se==a)
    return b;
  // if se is a constant -> return se
  auto o=dynamic_cast<const AST::Operation*>(se.get());
  if(!o)
    return se;
  // if se is a Operation -> copy the operation and call subst on its childs
  vector<SymbolicExpression> child;
  transform(o->child.begin(), o->child.end(), back_inserter(child), [&a, &b](const SymbolicExpression &x){
    return substScalar(x, a, b);
  });
  return AST::Operation::create(o->op, child);
}

// ***** ByteCode *****

ByteCode::ByteCode(size_t n) {
  retPtr = &retValue;
  argsPtr.resize(n);
}

ByteCode::ByteCode(ByteCode &&) noexcept {
  assert(0 && "ByteCode cannot be moved or copied. But the move-ctor is defined to allow in-place construction e.g. using emplace_back.");
}

// ***** Vertex *****

bool Vertex::isZero() const {
  if(isConstantInt() && static_cast<const Constant<long>*>(this)->getValue()==0)
    return true;
  return false;
}

bool Vertex::isOne() const {
  if(isConstantInt() && static_cast<const Constant<long>*>(this)->getValue()==1)
    return true;
  return false;
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
bool Constant<T>::equal(const SymbolicExpression &b, MapIVSE &m) const {
  // a constant is only equal to b if b is the same as this.
  return this->shared_from_this()==b;
}

template<class T>
SymbolicExpression Constant<T>::parDer(const IndependentVariable &x) const {
  return Constant<long>::create(0);
}

template<class T>
Constant<T>::Constant(const T& c_) : c(c_) {}

template<class T>
std::vector<ByteCode>::iterator Constant<T>::dumpByteCode(vector<ByteCode> &byteCode, map<const Vertex*, vector<ByteCode>::iterator> &existingVertex) const {
  auto [lastIt, inserted] = existingVertex.insert(make_pair(this, vector<ByteCode>::iterator()));
  if(!inserted) return lastIt->second;

  byteCode.emplace_back(0);
  auto it = --byteCode.end();
  lastIt->second = it;
  it->func = [](double* r, const ByteCode::Arg& a) {};
  it->retValue = c;
  return --byteCode.end();
}

template<class T>
void Constant<T>::walkVertex(const function<void(const shared_ptr<const Vertex>&)> &func) const {
  func(this->shared_from_this());
}

template SymbolicExpression Constant<long   >::create(const long   &c_);
template SymbolicExpression Constant<double>::create(const double&c_);
template bool Constant<long   >::equal(const SymbolicExpression &b, MapIVSE &m) const;
template bool Constant<double>::equal(const SymbolicExpression &b, MapIVSE &m) const;
template SymbolicExpression Constant<long   >::parDer(const IndependentVariable &x) const;
template SymbolicExpression Constant<double>::parDer(const IndependentVariable &x) const;
template std::vector<ByteCode>::iterator Constant<long   >::dumpByteCode(vector<ByteCode> &byteCode, map<const Vertex*, vector<ByteCode>::iterator> &existingVertex) const;
template std::vector<ByteCode>::iterator Constant<double>::dumpByteCode(vector<ByteCode> &byteCode, map<const Vertex*, vector<ByteCode>::iterator> &existingVertex) const;
template void Constant<long   >::walkVertex(const function<void(const shared_ptr<const Vertex>&)> &func) const;
template void Constant<double>::walkVertex(const function<void(const shared_ptr<const Vertex>&)> &func) const;

// ***** Symbol *****

map<Symbol::CacheKey, weak_ptr<const Symbol>> Symbol::cache;

IndependentVariable Symbol::create(const boost::uuids::uuid& uuid_) {
  auto r=cache.insert(make_pair(uuid_, weak_ptr<const Symbol>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<const Symbol>(new Symbol(uuid_));
  r.first->second=newPtr;
  return newPtr;
}

bool Symbol::equal(const SymbolicExpression &b, MapIVSE &m) const {
  // if this symbol is not already listed in m than this symbol is still free.
  // This means that the the symbol can be equal to b if this symbol is set to b...
  auto ret=m.insert(pair<IndependentVariable, SymbolicExpression>(shared_from_this(), SymbolicExpression()));
  if(ret.second) {
    // ... this is done here (store this symbol in m and map it to b.
    ret.first->second=b;
    return true;
  }
  // if this symbol is already listed in m then true is returned if the mapped expression equal b.
  return ret.first->second==b;
}

SymbolicExpression Symbol::parDer(const IndependentVariable &x) const {
  return this == x.get() ? Constant<long>::create(1) : Constant<long>::create(0);
}

Symbol::Symbol(const boost::uuids::uuid& uuid_) : uuid(uuid_) {}

string Symbol::getUUIDStr() const {
#ifndef NDEBUG // FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID
  if(getenv("FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID")) {
    auto res=mapUUIDInt.insert(make_pair(uuid, mapUUIDInt.size()+1));
    return "s"+to_string(res.first->second);
  }
  return to_string(uuid);
#else
  return to_string(uuid);
#endif
}

std::vector<ByteCode>::iterator Symbol::dumpByteCode(vector<ByteCode> &byteCode,
                                                  map<const Vertex*, vector<ByteCode>::iterator> &existingVertex) const {
  auto [lastIt, inserted] = existingVertex.insert(make_pair(this, vector<ByteCode>::iterator()));
  if(!inserted) return lastIt->second;

  byteCode.emplace_back(0);
  auto it = --byteCode.end();
  lastIt->second = it;
  double *xPtr = &x;
  it->func = [xPtr](double* r, const ByteCode::Arg& a) { *r = *xPtr; };
  return --byteCode.end();
}

void Symbol::walkVertex(const function<void(const shared_ptr<const Vertex>&)> &func) const {
  func(shared_from_this());
}

// ***** NativeFunction *****

map<NativeFunction::CacheKey, weak_ptr<const NativeFunction>, NativeFunction::CacheKeyComp> NativeFunction::cache;

bool NativeFunction::CacheKeyComp::operator()(const CacheKey& l, const CacheKey& r) const {
  if(olf(get<0>(l), get<0>(r)) || olf(get<0>(r), get<0>(l)))
    return olf(get<0>(l), get<0>(r));

#define FMATVEC_HELPER(ARG) \
  if(get<ARG>(l).size() != get<ARG>(r).size()) \
    return get<ARG>(l).size()<get<ARG>(r).size(); \
  for(size_t i=0; i<get<ARG>(l).size(); ++i) \
    if(ol(get<ARG>(l)[i],get<ARG>(r)[i]) || ol(get<ARG>(r)[i],get<ARG>(l)[i])) \
      return ol(get<ARG>(l)[i], get<ARG>(r)[i]);
  FMATVEC_HELPER(1)
  FMATVEC_HELPER(2)
  FMATVEC_HELPER(3)

  return false;
}

NativeFunction::NativeFunction(const shared_ptr<ScalarFunctionWrapArg> &funcWrapper_, const vector<SymbolicExpression> &argS_,
                               const vector<SymbolicExpression> &dir1S_, const vector<SymbolicExpression> &dir2S_) :
  funcWrapper(funcWrapper_), argS(argS_), dir1S(dir1S_), dir2S(dir2S_) {}

SymbolicExpression NativeFunction::create(const shared_ptr<ScalarFunctionWrapArg> &funcWrapper, const vector<SymbolicExpression> &argS,
                                          const vector<SymbolicExpression> &dir1S, const vector<SymbolicExpression> &dir2S) {
  vector<weak_ptr<const Vertex>> weakArgS;
  // we cannot just use std::copy here since std::shared_ptr is a private base of SymbolicExpression
  transform(argS.begin(), argS.end(), back_inserter(weakArgS), [](const SymbolicExpression &x){
    return weak_ptr<const Vertex>(x);
  });

  vector<weak_ptr<const Vertex>> weakDir1S;
  // we cannot just use std::copy here since std::shared_ptr is a private base of SymbolicExpression
  transform(dir1S.begin(), dir1S.end(), back_inserter(weakDir1S), [](const SymbolicExpression &x){
    return weak_ptr<const Vertex>(x);
  });

  vector<weak_ptr<const Vertex>> weakDir2S;
  // we cannot just use std::copy here since std::shared_ptr is a private base of SymbolicExpression
  transform(dir2S.begin(), dir2S.end(), back_inserter(weakDir2S), [](const SymbolicExpression &x){
    return weak_ptr<const Vertex>(x);
  });

  auto r=cache.insert(make_pair(make_tuple(weak_ptr<ScalarFunctionWrapArg>(funcWrapper),
                      weakArgS, weakDir1S, weakDir2S), weak_ptr<const NativeFunction>()));
  if(!r.second) {
    auto oldPtr=r.first->second.lock();
    if(oldPtr)
      return oldPtr;
  }
  auto newPtr=shared_ptr<NativeFunction>(new NativeFunction(funcWrapper, argS, dir1S, dir2S));
  r.first->second=newPtr;
  return newPtr;
}

SymbolicExpression NativeFunction::parDer(const IndependentVariable &x) const {
  if(dir1S.size()==0) {
    // create a new first derivative
    vector<SymbolicExpression> dir1S_;
    bool allZero=true;
    for(size_t i=0; i<argS.size(); ++i) {
      dir1S_.push_back(argS[i]->parDer(x));
      if(!dir1S_.back()->isZero())
        allZero=false;
    }
    if(allZero)
      return 0;
    return create(funcWrapper, argS, dir1S_);
  }
  if(dir2S.size()==0) {
    // create a new second derivative
    vector<SymbolicExpression> dir1S_;
    vector<SymbolicExpression> dir2S_;
    SymbolicExpression ret;
    bool all1Zero=true;
    for(size_t i=0; i<argS.size(); ++i) {
      dir1S_.push_back(dir1S[i]->parDer(x));
      if(!dir1S_.back()->isZero())
        all1Zero=false;
    }
    if(all1Zero)
      ret=0;
    else
      ret=create(funcWrapper, argS, dir1S_);
    bool all2Zero=true;
    for(size_t i=0; i<argS.size(); ++i) {
      dir2S_.push_back(argS[i]->parDer(x));
      if(!dir2S_.back()->isZero())
        all2Zero=false;
    }
    if(!all2Zero)
      ret=ret+create(funcWrapper, argS, dir1S, dir2S_);
    return ret;
  }
  throw std::runtime_error("Derivative higher than 2 of external function needed. "
                           "External functions provide only the second derivative.");
}

bool NativeFunction::equal(const SymbolicExpression &b, MapIVSE &m) const {
  // a native function equals b only if b is also a native function
  auto bnf=dynamic_pointer_cast<const NativeFunction>(b);
  if(!bnf)
    return false;
  // if the native functions are not the same -> return false
  if(funcWrapper.get()!=bnf->funcWrapper.get())
    return false;
  // if any of the argS's are not equal -> return false
  if(argS.size()!=bnf->argS.size())
    return false;
  for(size_t i=0; i<argS.size(); ++i)
    if(!argS[i]->equal(bnf->argS[i], m))
      return false;
  if(dir1S.size()!=bnf->dir1S.size())
    return false;
  for(size_t i=0; i<dir1S.size(); ++i)
    if(!dir1S[i]->equal(bnf->dir1S[i], m))
      return false;
  if(dir2S.size()!=bnf->dir2S.size())
    return false;
  for(size_t i=0; i<dir2S.size(); ++i)
    if(!dir2S[i]->equal(bnf->dir2S[i], m))
      return false;
  return true;
}

std::vector<ByteCode>::iterator NativeFunction::dumpByteCode(std::vector<ByteCode> &byteCode,
                                             std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const {

  auto [lastIt, inserted] = existingVertex.insert(make_pair(this, vector<ByteCode>::iterator()));
  if(!inserted) return lastIt->second;

  std::vector<vector<ByteCode>::iterator> allargSItVec;
  for(auto &A : {argS, dir1S, dir2S})
    for(auto &a : A) {
      auto it = a->dumpByteCode(byteCode, existingVertex);
      allargSItVec.push_back(it);
    }

  byteCode.emplace_back(argS.size()+dir1S.size()+dir2S.size());
  auto it = --byteCode.end();
  if(argS.size()+dir1S.size()+dir2S.size() > ByteCode::N)
    Atom::msgStatic(Atom::Warn)<<
      "A NativeFunction with more "<<argS.size()+dir1S.size()+dir2S.size()<<
      " arguments exists, which is more the the inplace size "<<ByteCode::N<<"."<<endl<<
      "Memory cache will not be optimal. Increase ByteCode::N at compile time to improve performance."<<endl;
  lastIt->second = it;
  auto func = funcWrapper; // we don't want to capture this (but need to capture by value the shared_ptr this->funcWapper)
  if(dir1S.size()==0 && dir2S.size()==0) // 0th derivative
    it->func = [func](double *ret, const ByteCode::Arg &arg) { *ret = (*func)(arg); };
  else if(dir2S.size()==0) // 1th derivative
    it->func = [func](double *ret, const ByteCode::Arg &arg) { *ret = func->dirDer(arg); };
  else // 2st derivative
    it->func = [func](double *ret, const ByteCode::Arg &arg) { *ret = func->dirDerDirDer(arg); };

  auto allargSIt = allargSItVec.begin();
  auto argsPtrIt = it->argsPtr.begin();
  for(; allargSIt<allargSItVec.end(); ++allargSIt, ++argsPtrIt)
    *argsPtrIt = (*allargSIt)->retPtr;

  return --byteCode.end();
}

void NativeFunction::walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const  {
  for(auto &A : {argS, dir1S, dir2S})
    for(auto &a : A)
      a->walkVertex(func);
  func(shared_from_this());
}

// ***** Operation *****

map<Operation::CacheKey, weak_ptr<const Operation>, Operation::CacheKeyComp> Operation::cache;

//MISSING: optimized calls for int arguments
#define FUNC(expr) [](double* r, const ByteCode::Arg& arg) { \
  *r = expr; \
}
#define _a *arg[0]
#define _b *arg[1]
#define _c *arg[2]
const std::map<Operation::Operator, Operation::OpMap> Operation::opMap {
//  Operator     Name          Lambda-Function
  { Plus,      { "plus"      , FUNC( _a + _b                           ) }},
  { Minus,     { "minus"     , FUNC( _a - _b                           ) }},
  { Mult,      { "mult"      , FUNC( _a * _b                           ) }},
  { Div,       { "div"       , FUNC( _a / _b                           ) }},
  { Pow,       { "pow"       , FUNC( std::pow(_a, _b)                  ) }},
  { Log,       { "log"       , FUNC( std::log(_a)                      ) }},
  { Sqrt,      { "sqrt"      , FUNC( std::sqrt(_a)                     ) }},
  { Neg,       { "neg"       , FUNC( - _a                              ) }},
  { Sin,       { "sin"       , FUNC( std::sin(_a)                      ) }},
  { Cos,       { "cos"       , FUNC( std::cos(_a)                      ) }}, 
  { Tan,       { "tan"       , FUNC( std::tan(_a)                      ) }}, 
  { Sinh,      { "sinh"      , FUNC( std::sinh(_a)                     ) }}, 
  { Cosh,      { "cosh"      , FUNC( std::cosh(_a)                     ) }}, 
  { Tanh,      { "tanh"      , FUNC( std::tanh(_a)                     ) }}, 
  { ASin,      { "asin"      , FUNC( std::asin(_a)                     ) }}, 
  { ACos,      { "acos"      , FUNC( std::acos(_a)                     ) }}, 
  { ATan,      { "atan"      , FUNC( std::atan(_a)                     ) }}, 
  { ATan2,     { "atan2"     , FUNC( std::atan2(_a, _b)                ) }}, 
  { ASinh,     { "asinh"     , FUNC( std::asinh(_a)                    ) }}, 
  { ACosh,     { "acosh"     , FUNC( std::acosh(_a)                    ) }}, 
  { ATanh,     { "atanh"     , FUNC( std::atanh(_a)                    ) }}, 
  { Exp,       { "exp"       , FUNC( std::exp(_a)                      ) }}, 
  { Sign,      { "sign"      , FUNC( boost::math::sign(_a)             ) }}, 
  { Heaviside, { "heaviside" , FUNC( 0.5 * boost::math::sign(_a) + 0.5 ) }}, 
  { Abs,       { "abs"       , FUNC( std::abs(_a)                      ) }},
  { Min,       { "min"       , FUNC( std::min(_a, _b)                  ) }}, 
  { Max,       { "max"       , FUNC( std::max(_a, _b)                  ) }},
  { Condition, { "condition" , FUNC( _a >= 0 ? _b : _c                 ) }}, 
};                                                         
#undef FUNC
#undef _a
#undef _b
#undef _c
                                                           
SymbolicExpression Operation::create(Operator op_, const vector<SymbolicExpression> &child_) {
  static bool optimizeExpressions=true; // this is "always" true (except for the initial call, see below)
  static bool firstCall=true;                              
  static vector<pair<SymbolicExpression,SymbolicExpression>> optExpr; // list of expression optimizations
  if(firstCall) { // on the first call build the list of expression optimizations
    firstCall=false;                                       
    IndependentVariable a;                                 
    IndependentVariable b;
    SymbolicExpression error(std::shared_ptr<const Vertex>{nullptr});
    // we need to disable the expression optimization during buildup of the expressions to optimize
    optimizeExpressions=false;
    optExpr={
      // list of expressions (the left ones) to optimize; the right ones are the optimized expressions
      // which must mathematically equal the left ones but are simpler.
      {    0 + a       ,  a},
      {  0.0 + a       ,  a},
      {    a + 0       ,  a},
      {    a + 0.0     ,  a},
      {    a - 0       ,  a},
      {    a - 0.0     ,  a},
      {    a - a       ,  0},
      {    0 - a       , -a},
      {      - 0       ,  0},
      {    0 * a       ,  0},
      {  0.0 * a       ,  0},
      {    a * 0       ,  0},
      {    a * 0.0     ,  0},
      {    1 * a       ,  a},
      {  1.0 * a       ,  a},
      { -1   * a       , -a},
      { -1.0 * a       , -a},
      {    a * 1       ,  a},
      {    a * 1.0     ,  a},
      {    a * (-1)    , -a},
      {    a * (-1.0)  , -a},
      {    a * a       , pow(a,2)},
      {    a / 0       , error},
      {    a / 0.0     , error},
      {    0 / a       ,  0},
      {  0.0 / a       ,  0},
      {    a / a       ,  1},
      {    a / 1       ,  a},
      {    a / 1.0     ,  a},
      {    a / -1      , -a},
      {    a / -1.0    , -a},
      { pow(a,1)       , a},
      { pow(a,1.0)     , a},
      { pow(a,-1)      , 1/a},
      { pow(a,-1.0)    , 1/a},
      { pow(a,0)       , 1},
      { pow(a,0.0)     , 1},
      { pow(a,b)*a     , pow(a,b+1)},
      { a*pow(a,b)     , pow(a,b+1)},
      { pow(a,b)/a     , pow(a,b-1)},
      { log(0)         , error},
      { log(0.0)       , error},
      { sqrt(-1)       , error},
      { sqrt(-1.0)     , error},
      { acosh(-1)      , error},
      { acosh(-1.0)    , error},
      { acosh(0)       , error},
      { acosh(0.0)     , error},
      { atanh(-1)      , error},
      { atanh(-1.0)    , error},
      { atanh(1)       , error},
      { atanh(1.0)     , error},
    };
    // now enable the optimizations again
    optimizeExpressions=true;
  }

  if(optimizeExpressions) { // optimizeExpressions is "always" true (except for the initial call, see above)
    // build the current operation
    SymbolicExpression curOp=shared_ptr<Operation>(new Operation(op_, child_));
    // loop over all optimization expressions
    for(auto &opt : optExpr) {
      SymbolicExpression in=opt.first;
      SymbolicExpression out=opt.second;
      MapIVSE m;
      if(in->equal(curOp, m)) {
        // if this optimization expression matches ...
        // ... and the result is an error -> throw
        if(!out.get())
        {
          stringstream str;
          str<<in;
          throw runtime_error("Illegal constant argument in operation: "+str.str());
        }
        // ... then return the optimized expressions with all symbols of m replaced.
        for(auto &x : m)
          out=substScalar(out, x.first, x.second);
        return out;
      }
    }

    // optimize Constant arguments (if ALL are Constant)
    bool allConst=all_of(child_.begin(), child_.end(), [](const SymbolicExpression& c) {
      return c->isConstantInt() || dynamic_pointer_cast<const Constant<double>>(c)!=nullptr;
    });
    if(allConst) {
      Operation op(op_, child_);
      vector<ByteCode> byteCode;
      byteCode.reserve(ByteCode::N+1);// preserve enought memory, we cannot move/copy a ByteCode.
      map<const Vertex*, vector<AST::ByteCode>::iterator> existingVertex;
      auto retIt = op.dumpByteCode(byteCode, existingVertex);
      std::for_each(byteCode.begin(), byteCode.end(), [](const AST::ByteCode &bc){
        bc.func(bc.retPtr, bc.argsPtr);
      });
      double doubleValue = retIt->retValue;
      if(doubleValue > static_cast<double>(numeric_limits<long>::min()) &&
         doubleValue < static_cast<double>(numeric_limits<long>::max())) {
        long intValue=lround(doubleValue);
        if(std::abs(intValue-doubleValue)<2*numeric_limits<double>::epsilon())
          return Constant<long>::create(intValue);
      }
      return Constant<double>::create(doubleValue);
    }
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

bool Operation::equal(const SymbolicExpression &b, MapIVSE &m) const {
  // a operation equals b only if b is also a operation
  auto bo=dynamic_pointer_cast<const Operation>(b);
  if(!bo)
    return false;
  // if the operation are not the same -> return false
  if(op!=bo->op)
    return false;
  // if any of the children vertexes are not equal -> return false
  if(child.size()!=bo->child.size())
    return false;
  for(size_t i=0; i<child.size(); ++i)
    if(!child[i]->equal(bo->child[i], m))
      return false;
  return true;
}

SymbolicExpression Operation::parDer(const IndependentVariable &x) const {
  auto v=SymbolicExpression(shared_from_this()); // expression to be differentiated
  auto a=child.size()>=1 ? child[0] : SymbolicExpression(); // expression of first argument
  auto b=child.size()>=2 ? child[1] : SymbolicExpression(); // expression of second argument
  auto c=child.size()>=3 ? child[2] : SymbolicExpression(); // expression of third argument
  auto ad=child.size()>=1 ? child[0]->parDer(x) : SymbolicExpression(); // pertial derivative of the expression of the first argument wrt the independent x
  auto bd=child.size()>=2 ? child[1]->parDer(x) : SymbolicExpression(); // pertial derivative of the expression of the second argument wrt the independent x
  auto cd=child.size()>=3 ? child[2]->parDer(x) : SymbolicExpression(); // pertial derivative of the expression of the third argument wrt the independent x
  switch(op) {
    case Plus:
      return ad + bd;
    case Minus:
      return ad - bd;
    case Mult:
      return ad * b + a * bd;
    case Div:
      if(dynamic_cast<const Constant<double>*>(b.get()) || dynamic_cast<const Constant<long>*>(b.get()))
        return ad / b;
      else
        return (ad * b - a * bd) / (b * b);
    case Pow:
      if(b->isConstantInt())
        return b * pow(a, static_cast<const Constant<long>*>(b.get())->getValue() - 1) * ad;
      else if(dynamic_cast<const Constant<double>*>(a.get()) || dynamic_cast<const Constant<long>*>(a.get()))
        return pow(a, b) * bd * log(a);
      else
        return pow(a, b-1) * (bd * log(a) * a + ad * b);
    case Log:
      return ad / a;
    case Sqrt:
      return ad * 0.5 / v;
    case Neg:
      return - ad;
    case Sin:
      return cos(a) * ad;
    case Cos:
      return - sin(a) * ad;
    case Tan:
      return ad / pow(cos(a), 2);
    case Sinh:
      return cosh(a) * ad;
    case Cosh:
      return sinh(a) * ad;
    case Tanh:
      return (1 - pow(v, 2)) * ad;
    case ASin:
      return ad / sqrt(1 - pow(a, 2));
    case ACos:
      return - ad / sqrt(1 - pow(a, 2));
    case ATan:
      return ad / (1 + pow(a, 2));
    case ATan2:
      return -a * bd / (pow(a, 2) + pow(b, 2)) + ad * b / (pow(a, 2) + pow(b, 2));
    case ASinh:
      return ad / sqrt(pow(a, 2) + 1);
    case ACosh:
      return ad / sqrt(pow(a, 2) - 1);
    case ATanh:
      return ad / (1 - pow(a, 2));
    case Exp:
      return exp(a) * ad;
    case Sign:
      return 0;
    case Heaviside:
      return 0;
    case Abs:
      return sign(a) * ad;
    case Min:
      return ad * heaviside(b - a) + bd * heaviside(a - b);
    case Max:
      return ad * heaviside(a - b) + bd * heaviside(b - a);
    case Condition:
      return condition(a, bd, cd);
  }
  throw runtime_error("Unknown operation.");
}

Operation::Operation(Operator op_, const vector<SymbolicExpression> &child_) : op(op_), child(child_) {}

bool Operation::CacheKeyComp::operator()(const CacheKey& l, const CacheKey& r) const {
  if(l.first!=r.first)
    return l.first<r.first;
  if(l.second.size()!=r.second.size())
    return l.second.size()<r.second.size();
  for(size_t i=0; i<l.second.size(); ++i)
    if(ol(l.second[i],r.second[i]) || ol(r.second[i],l.second[i]))
      return ol(l.second[i], r.second[i]);
  return false;
}

std::vector<ByteCode>::iterator Operation::dumpByteCode(vector<ByteCode> &byteCode,
                                                     map<const Vertex*, vector<ByteCode>::iterator> &existingVertex) const {
  auto [lastIt, inserted] = existingVertex.insert(make_pair(this, vector<ByteCode>::iterator()));
  if(!inserted) return lastIt->second;

  std::vector<vector<ByteCode>::iterator> childItVec;
  for(auto &c : child) {
    auto it = c->dumpByteCode(byteCode, existingVertex);
    childItVec.push_back(it);
  }

  byteCode.emplace_back(child.size());
  auto it = --byteCode.end();
  assert(child.size() <= ByteCode::N && "If more function arguments are needed you have to increase ByteCode::N at compile time.");
  lastIt->second = it;
  it->func = opMap.at(op).func;

  // optimization for pow with an interger exponent
  if(op == Pow && child[1]->isConstantInt())
    it->func = [](double* r, const ByteCode::Arg& a){ *r = std::pow(*a[0], static_cast<int>(*a[1])); };

  auto childIt = childItVec.begin();
  auto argsPtrIt = it->argsPtr.begin();
  for(; childIt != childItVec.end(); ++childIt, ++argsPtrIt)
    *argsPtrIt = (*childIt)->retPtr;

  return --byteCode.end();
}

void Operation::walkVertex(const function<void(const shared_ptr<const Vertex>&)> &func) const {
  for(auto &c : child)
    c->walkVertex(func);
  func(shared_from_this());
}

} // end namespace AST

template<>
boost::spirit::qi::rule<boost::spirit::istream_iterator, IndependentVariable()>& getBoostSpiritQiRule<IndependentVariable>() {
  namespace qi = boost::spirit::qi;
  namespace phx = boost::phoenix;

  static boost::spirit::qi::rule<boost::spirit::istream_iterator, IndependentVariable()> ret;
  static bool init=false;
  if(!init) {
    init=true;
#ifndef NDEBUG // FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID
    if(getenv("FMATVEC_DEBUG_SYMBOLICEXPRESSION_UUID"))
      ret = *qi::blank >> 's' >> qi::int_[qi::_val=phx::bind(&createSymbolByInt, qi::_1)];
    else
      ret = *qi::blank >> (qi::repeat(36)[qi::char_("a-z0-9-")])[qi::_val=phx::bind(&createSymbolByVec, qi::_1)];
#else
    ret = *qi::blank >> (qi::repeat(36)[qi::char_("a-f0-9-")])[qi::_val=phx::bind(&createSymbolByVec, qi::_1)];
#endif
  }
  return ret;
}

template<>
boost::spirit::qi::rule<boost::spirit::istream_iterator, SymbolicExpression()>& getBoostSpiritQiRule<SymbolicExpression>() {
  namespace qi = boost::spirit::qi;
  namespace phx = boost::phoenix;
  using It = boost::spirit::istream_iterator;

  static boost::spirit::qi::rule<boost::spirit::istream_iterator, SymbolicExpression()> vertex;
  static bool init=false;
  if(!init) {
    init=true;
    static qi::rule<It, SymbolicExpression()>  constLong;
    static qi::real_parser<double, qi::strict_real_policies<double>> strict_double;
    static qi::rule<It, SymbolicExpression()>  constDouble;
    static qi::rule<It, SymbolicExpression()>  operation;
    // NativeFunction cannot be serialized and hence not deserialized

    static qi::symbols<char, AST::Operation::Operator> opSym;
    for(auto &x : AST::Operation::opMap)
      opSym.add(x.second.funcName, x.first);

    constLong   = *qi::blank >> qi::long_[qi::_val=phx::bind(&AST::Constant<long>::create, qi::_1)];
    constDouble = *qi::blank >> strict_double[qi::_val=phx::bind(&AST::Constant<double>::create, qi::_1)];
    qi::rule<It, IndependentVariable()> &symbol = getBoostSpiritQiRule<IndependentVariable>();
    operation   = *qi::blank >> (opSym >> '(' >> (vertex % ',') >> ')')[qi::_val=phx::bind(&AST::Operation::create, qi::_1, qi::_2)];
    vertex      = symbol | operation | constDouble | constLong;
  }
  return vertex;
}

namespace {
  long getConstInt(const SymbolicExpression &se, bool &pass) {
    pass=true;
    auto c=dynamic_pointer_cast<const AST::Constant<long>>(se);
    if(c) return c->getValue();
    pass=false;
    return 0;
  }
  double getConstDouble(const SymbolicExpression &se, bool &pass) {
    pass=true;
    auto c=dynamic_pointer_cast<const AST::Constant<double>>(se);
    if(c) return c->getValue();
    pass=false;
    return 0;
  }
  AST::Operation::Operator getOperationOp(const SymbolicExpression &se, bool &pass) {
    pass=true;
    auto o=dynamic_pointer_cast<const AST::Operation>(se);
    if(o) return o->getOp();
    pass=false;
    return AST::Operation::Plus;
  }
  const vector<SymbolicExpression>& getOperationChilds(const SymbolicExpression &se, bool &pass) {
    pass=true;
    auto o=dynamic_pointer_cast<const AST::Operation>(se);
    if(o) return o->getChilds();
    pass=false;
    static vector<SymbolicExpression> empty;
    return empty;
  }
  string getSymbol(const SymbolicExpression &se, bool &pass) {
    pass=true;
    auto s=dynamic_pointer_cast<const AST::Symbol>(se);
    if(s) return s->getUUIDStr();
    pass=false;
    return {};
  }
  long getFunction(const SymbolicExpression &se, bool &pass) {
    throw runtime_error("Cannot serialize a symbolic expression containing a function.");
  }
}

template<>
boost::spirit::karma::rule<std::ostream_iterator<char>, SymbolicExpression()>& getBoostSpiritKarmaRule<SymbolicExpression>() {
  namespace karma = boost::spirit::karma;
  namespace phx = boost::phoenix;
  using It = std::ostream_iterator<char>;

  static boost::spirit::karma::rule<It, SymbolicExpression()> vertex;
  static bool init=false;
  if(!init) {
    init=true;
    static karma::rule<It, SymbolicExpression()> constLong;
    static karma::rule<It, SymbolicExpression()> constDouble;
    static karma::rule<It, SymbolicExpression()> symbol;
    static karma::rule<It, SymbolicExpression()> nativeFunction;
    static karma::rule<It, SymbolicExpression()> operation;

    auto &doubleBitIdentical=getBoostSpiritKarmaRule<double>();

    static karma::symbols<AST::Operation::Operator, string> opSym;
    for(auto &x : AST::Operation::opMap)
      opSym.add(x.first, x.second.funcName);

    constLong      = karma::long_[karma::_1=phx::bind(&getConstInt, karma::_val, karma::_pass)];
    constDouble    = doubleBitIdentical[karma::_1=phx::bind(&getConstDouble, karma::_val, karma::_pass)];
    operation      = opSym[karma::_1=phx::bind(&getOperationOp, karma::_val, karma::_pass)] << '(' <<
                     (vertex % ',')[karma::_1=phx::bind(&getOperationChilds, karma::_val, karma::_pass)] << ')';
    symbol         = karma::string[karma::_1=phx::bind(&getSymbol, karma::_val, karma::_pass)];
    nativeFunction = karma::int_[karma::_1=phx::bind(&getFunction, karma::_val, karma::_pass)]; // dummy, getFunction throw always
    vertex         = constLong | constDouble | symbol | operation | nativeFunction;
  }
  return vertex;
}

template<>
boost::spirit::karma::rule<std::ostream_iterator<char>, IndependentVariable()>& getBoostSpiritKarmaRule<IndependentVariable>() {
  namespace karma = boost::spirit::karma;
  namespace phx = boost::phoenix;
  using It = std::ostream_iterator<char>;

  static karma::rule<It, IndependentVariable()> symbol;
  static bool init=false;
  if(!init) {
    init=true;
    symbol = karma::string[karma::_1=phx::bind(&getSymbol, karma::_val, karma::_pass)];
  }
  return symbol;
}

} // end namespace fmatvec
