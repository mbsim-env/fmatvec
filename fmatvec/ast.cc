#include "ast.h"
#include "stream_impl.h"
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/scope_exit.hpp>
#include <sstream>
#include <boost/spirit/include/qi.hpp>
#include <string>
#include <iostream>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/bind/bind_function.hpp>
#include <boost/spirit/include/karma.hpp>

using namespace std;

namespace fmatvec {

// ***** SymbolicExpression *****

#if !defined(NDEBUG) && !defined(SWIG)
unsigned long SymbolicExpression::evalOperationsCount = 0;
#endif

namespace {
  SymbolicExpression stringCTorHelper(const string &str) {
    istringstream s(str);
    SymbolicExpression ret;
    s>>ret;
    return ret;
  }
}
SymbolicExpression::SymbolicExpression() : shared_ptr<const AST::Vertex>(AST::Constant<int>::create(0)) {}
SymbolicExpression::SymbolicExpression(const string &str) : SymbolicExpression(stringCTorHelper(str)) {}
template<class T> SymbolicExpression::SymbolicExpression(const shared_ptr<T> &x) : shared_ptr<const AST::Vertex>(x) {}
SymbolicExpression::SymbolicExpression(int x) : shared_ptr<const AST::Vertex>(AST::Constant<int>::create(x)) {}
SymbolicExpression::SymbolicExpression(double x) : shared_ptr<const AST::Vertex>(AST::Constant<double>::create(x)) {}
SymbolicExpression::SymbolicExpression(ConstructSymbol) : shared_ptr<const AST::Vertex>(AST::Symbol::create()) {}

template<class T>
void erase(T &cache) {
  for(auto it=cache.begin(), endIt=cache.end(); it!=endIt; )
    if(it->second.expired()) it=cache.erase(it); else ++it;
}

void SymbolicExpression::garbageCollect() {
  erase(AST::Constant<int>::cache);
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

SymbolicExpression abs(const SymbolicExpression &a) {
  return AST::Operation::create(AST::Operation::Abs, {a});
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

SymbolicExpression subst(const SymbolicExpression &se, const IndependentVariable& a, const SymbolicExpression &b) {
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
    return subst(x, a, b);
  });
  return AST::Operation::create(o->op, child);
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

SymbolicExpression operator-(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)-b;
}

SymbolicExpression operator*(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)*b;
}

SymbolicExpression operator/(int a, const SymbolicExpression &b) {
  return SymbolicExpression(a)/b;
}

SymbolicExpression parDer(const SymbolicExpression &dep, const IndependentVariable &indep) {
  return dep->parDer(indep);
}

// ***** IndependentVariable *****

IndependentVariable::IndependentVariable() : SymbolicExpression(constructSymbol) {}

IndependentVariable::IndependentVariable(const string &str) : SymbolicExpression(str) {}

IndependentVariable::IndependentVariable(const shared_ptr<const AST::Symbol> &x) : SymbolicExpression(x) {}

namespace AST { // internal namespace

// ***** Vertex *****

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
bool Constant<T>::equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const {
  // a constant is only equal to b if b is the same as this.
  return this->shared_from_this()==b;
}

template<class T>
SymbolicExpression Constant<T>::parDer(const IndependentVariable &x) const {
  return Constant<int>::create(0);
}

template<class T>
Constant<T>::Constant(const T& c_) : c(c_) {}

template SymbolicExpression Constant<int   >::create(const int   &c_);
template SymbolicExpression Constant<double>::create(const double&c_);
template bool Constant<int   >::equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const;
template bool Constant<double>::equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const;
template SymbolicExpression Constant<int   >::parDer(const IndependentVariable &x) const;
template SymbolicExpression Constant<double>::parDer(const IndependentVariable &x) const;

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
  newPtr->dependsOn.insert(make_pair(weak_ptr<const Symbol>(newPtr), 0)); // we cannot call this in Symbol::Symbol (it may be possible with weak_from_this() for C++17)
  r.first->second=newPtr;
  return newPtr;
}

bool Symbol::equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const {
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
  return this == x.get() ? Constant<int>::create(1) : Constant<int>::create(0);
}

Symbol::Symbol(const boost::uuids::uuid& uuid_) : version(0), uuid(uuid_) {}

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

// ***** Operation *****

map<Operation::CacheKey, weak_ptr<const Operation>, Operation::CacheKeyComp> Operation::cache;

std::map<Operation::Operator, string> Operation::opMap {
  { Plus,  "plus"},
  { Minus, "minus"},
  { Mult,  "mult"},
  { Div,   "div"},
  { Pow,   "pow"},
  { Log,   "log"},
  { Sqrt,  "sqrt"},
  { Neg,   "neg"},
  { Sin,   "sin"},
  { Cos,   "cos"},
  { Tan,   "tan"},
  { Sinh,  "sinh"},
  { Cosh,  "cosh"},
  { Tanh,  "tanh"},
  { ASin,  "asin"},
  { ACos,  "acos"},
  { ATan,  "atan"},
  { ASinh, "asinh"},
  { ACosh, "acosh"},
  { ATanh, "atanh"},
  { Exp,   "exp"},
  { Sign,  "sign"},
  { Abs,   "abs"}
};

SymbolicExpression Operation::create(Operator op_, const vector<SymbolicExpression> &child_) {
  static bool optimizeExpressions=true; // this is "always" true (except for the initial call, see below)
  static bool firstCall=true;
  static vector<pair<SymbolicExpression,SymbolicExpression>> optExpr; // list of expression optimizations
  if(firstCall) { // on the first call build the list of expression optimizations
    firstCall=false;
    IndependentVariable a;
    IndependentVariable b;
    // we need to disable the expression optimization during buildup of the expressions to optimize
    optimizeExpressions=false;
    optExpr={
      // list of expressions (the left ones) to optimize; the right ones are the optimized expressions
      // which must mathematically equal the left ones but are simpler.
      { 0+a       , a},
      { a+0       , a},
      { a-0       , a},
      { 0*a       , 0},
      { a*0       , 0},
      { 1*a       , a},
      { a*1       , a},
      { a*a       , pow(a,2)},
      { a/1       , a},
      { a/a       , 1}, // we do not care about the 0/0 case here :-(
      { 0/a       , 0}, // we do not care about the /0 case here :-(
      { pow(a,1)  , a},
      { pow(a,0)  , 1},
      { pow(a,b)*a, pow(a,b+1)},
      { a*pow(a,b), pow(a,b+1)},
      { pow(a,b)/a, pow(a,b-1)} // we do not care about the /0 case here :-(
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
      map<IndependentVariable, SymbolicExpression> m;
      if(in->equal(curOp, m)) {
        // if this optimization expression matches then return the optimized expressions with all symbols of m replaced.
        for(auto &x : m)
          out=subst(out, x.first, x.second);
        return out;
      }
    }

    // optimize Constant arguments (if ALL are Constant)
    bool allConst=all_of(child_.begin(), child_.end(), [](const SymbolicExpression& c) {
      return c->isConstantInt() || dynamic_pointer_cast<const Constant<double>>(c)!=nullptr;
    });
    if(allConst) {
      Operation op(op_, child_);
      double doubleValue=op.eval();
      int intValue=lround(doubleValue);
      if(std::abs(intValue-doubleValue)<2*numeric_limits<double>::epsilon())
        return Constant<int>::create(intValue);
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

bool Operation::equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const {
  // a operation equals b only if b is also a operation
  auto bo=dynamic_pointer_cast<const Operation>(b);
  if(!bo)
    return false;
  // if the operation are not the same -> return false
  if(op!=bo->op)
    return false;
  // if any of the children vertexes are not equal -> return false
  for(size_t i=0; i<child.size(); ++i)
    if(!child[i]->equal(bo->child[i], m))
      return false;
  return true;
}

SymbolicExpression Operation::parDer(const IndependentVariable &x) const {
  auto v=SymbolicExpression(shared_from_this()); // expression to be differentiated
  auto a=child.size()>=1 ? child[0] : SymbolicExpression(); // expression of first argument
  auto b=child.size()>=2 ? child[1] : SymbolicExpression(); // expression of second argument
  auto ad=child.size()>=1 ? child[0]->parDer(x) : SymbolicExpression(); // pertial derivative of the expression of the first argument wrt the independent x
  auto bd=child.size()>=2 ? child[1]->parDer(x) : SymbolicExpression(); // pertial derivative of the expression of the second argument wrt the independent x
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
      if(b->isConstantInt())
        return b * pow(a, static_cast<const Constant<int>*>(b.get())->getValue() - 1) * ad;
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
    case Abs:
      return sign(a) * ad;
  }
  throw runtime_error("Unknown operation.");
}

Operation::Operation(Operator op_, const vector<SymbolicExpression> &child_) : op(op_), child(child_) {
  for(auto &c : child)
    dependsOn.insert(c->getDependsOn().begin(), c->getDependsOn().end());
}

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
    static qi::rule<It, SymbolicExpression()>  constInt;
    static qi::real_parser<double, qi::strict_real_policies<double>> strict_double;
    static qi::rule<It, SymbolicExpression()>  constDouble;
    static qi::rule<It, SymbolicExpression()>  operation;

    static qi::symbols<char, AST::Operation::Operator> opSym;
    for(auto &x : AST::Operation::opMap)
      opSym.add(x.second, x.first);

    constInt    = *qi::blank >> qi::int_[qi::_val=phx::bind(&AST::Constant<int>::create, qi::_1)];
    constDouble = *qi::blank >> strict_double[qi::_val=phx::bind(&AST::Constant<double>::create, qi::_1)];
    qi::rule<It, IndependentVariable()> &symbol = getBoostSpiritQiRule<IndependentVariable>();
    operation   = *qi::blank >> (opSym >> '(' >> (vertex % ',') >> ')')[qi::_val=phx::bind(&AST::Operation::create, qi::_1, qi::_2)];
    vertex      = symbol | operation | constDouble | constInt;
  }
  return vertex;
}

namespace {
  int getConstInt(const SymbolicExpression &se, bool &pass) {
    pass=true;
    auto c=dynamic_pointer_cast<const AST::Constant<int>>(se);
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
    return string();
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
    static karma::rule<It, SymbolicExpression()> constInt;
    static karma::rule<It, SymbolicExpression()> constDouble;
    static karma::rule<It, SymbolicExpression()> symbol;
    static karma::rule<It, SymbolicExpression()> operation;

    auto &doubleBitIdentical=getBoostSpiritKarmaRule<double>();

    static karma::symbols<AST::Operation::Operator, string> opSym;
    for(auto &x : AST::Operation::opMap)
      opSym.add(x.first, x.second);

    constInt    = karma::int_[karma::_1=phx::bind(&getConstInt, karma::_val, karma::_pass)];
    constDouble = doubleBitIdentical[karma::_1=phx::bind(&getConstDouble, karma::_val, karma::_pass)];
    operation   = opSym[karma::_1=phx::bind(&getOperationOp, karma::_val, karma::_pass)] << '(' <<
                  (vertex % ',')[karma::_1=phx::bind(&getOperationChilds, karma::_val, karma::_pass)] << ')';
    symbol      = karma::string[karma::_1=phx::bind(&getSymbol, karma::_val, karma::_pass)];
    vertex      = constInt | constDouble | symbol | operation;
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
