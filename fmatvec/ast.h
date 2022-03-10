#ifndef _FMATVEC_AST_H_
#define _FMATVEC_AST_H_

#include <map>
#include <array>
#include <memory>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <fmatvec/stream.h>

namespace fmatvec {

// forward declarations

class IndependentVariable;

namespace AST {
  class Vertex;
  class Symbol;
  class Operation;
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
class FMATVEC_EXPORT SymbolicExpression : public std::shared_ptr<const AST::Vertex> {
  template<class... Arg> friend class Eval;
  friend class EvalHelper;
  friend class AST::Operation;
  friend class AST::Constant<int>;
  friend class AST::Constant<double>;
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
    static signed long evalOperationsCount;
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
struct FMATVEC_EXPORT ByteCode {
  static constexpr size_t N { 3 };
  ByteCode();
  ByteCode(const ByteCode &) = delete;
  ByteCode(ByteCode &&);
  ByteCode &operator=(const ByteCode &) = delete;
  ByteCode &operator=(ByteCode &&) = delete;
  ~ByteCode() = default;
  std::function<void(double*, const std::array<double*,N>&)> func; // the operation this byteCode entry executes
  double  retValue; // storage of the return value of the operation: retPtr may point to this value
  double* retPtr; // a pointer to which the operation can write its result to
  std::array<double ,N> argsValue; // storage for arguments of the operation: argsPtr may point to these values
  std::array<double*,N> argsPtr; // pointers from which the operation reads its arguments
};

// ***** Vertex *****

//! A abstract class for a Vertex of the AST (abstract syntax tree).
class FMATVEC_EXPORT Vertex {
  friend Operation;
  public:

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
    virtual bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const=0;
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
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    inline bool isConstantInt() const override;
    //! Get the constant value of the vertex.
    inline const T& getValue() const;
    std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                  std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

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
    SymbolicExpression parDer(const IndependentVariable &x) const override;
    //! Set the value of this independent variable.
    //! This has an influence on the evaluation of all ASTs which depend on this independent variable.
    inline void setValue(double x_) const;

    std::string getUUIDStr() const;

    std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                  std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

  private:

    Symbol(const boost::uuids::uuid& uuid_);
    bool equal(const SymbolicExpression &b, std::map<IndependentVariable, SymbolicExpression> &m) const override;
    mutable double x = 0.0;
    boost::uuids::uuid uuid; // each variable has a uuid (this is only used when the AST is serialized and for caching)
    typedef boost::uuids::uuid CacheKey;
    static std::map<CacheKey, std::weak_ptr<const Symbol>> cache;
};

void Symbol::setValue(double x_) const {
  x=x_;
}

// ***** Operation *****

//! A vertex of the AST representing an operation.
class FMATVEC_EXPORT Operation : public Vertex, public std::enable_shared_from_this<Operation> {
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
    static SymbolicExpression create(Operator op_, const std::vector<SymbolicExpression> &child_);
    SymbolicExpression parDer(const IndependentVariable &x) const override;

    Operator getOp() const { return op; }
    const std::vector<SymbolicExpression>& getChilds() const { return child; }
    std::vector<ByteCode>::iterator dumpByteCode(std::vector<ByteCode> &byteCode,
                                  std::map<const Vertex*, std::vector<AST::ByteCode>::iterator> &existingVertex) const override;

    void walkVertex(const std::function<void(const std::shared_ptr<const Vertex>&)> &func) const override;

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

    struct OpMap {
      std::string funcName; // used to dump/read the expression using boost spirit/karma: e.g.
                            // "plus"
      std::function<void(double*, const std::array<double*,ByteCode::N>&)> func; // used for runtime evaluation: e.g.
                            // [](double* r, const std::array<double*,N>& a){ *r = *a[0] + *a[1]; }
    };
    static const std::map<Operator, OpMap> opMap;
};

} // end namespace AST
#endif

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
