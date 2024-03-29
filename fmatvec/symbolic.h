#ifndef _fmatvec_symbolic_h_
#define _fmatvec_symbolic_h_

#include "ast.h"
#include <set>
#include "function.h"
#include <boost/hana/type.hpp>

namespace fmatvec {

// definitions of matrix/vector operations/functions only usefull for symbolic calculations
// like building partial derivatives, evaluation of symbolic expressions to double values, ...

template<class Dep, class ATIndep>
Dep parDer(const Dep &dep, const ATIndep &indep) {
  Dep ret;
  if constexpr (std::is_same_v<Dep, SymbolicExpression>) {
  }
  else if constexpr (Dep::isVector)
    ret.resize(dep.size());
  else
    ret.resize(dep.rows(), dep.cols());
  auto d=dep.begin();
  auto r=ret.begin();
  for(; d!=dep.end(); ++d, ++r)
    *r=parDer(*d, indep);
  return ret;
}

template<class ATDep, class ATIndep>
Vector<Fixed<3>, ATDep> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const ATIndep &x) {
  Matrix<General, Fixed<3>, Fixed<3>, ATDep> Rs;
  for(int r=0; r<3; ++r)
    for(int c=0; c<3; ++c)
      Rs(r,c)=parDer(R(r,c), x);
  auto retTilde=Rs*trans(R);
  Vector<Fixed<3>, ATDep> ret;
  ret(0)=retTilde(2,1);
  ret(1)=retTilde(0,2);
  ret(2)=retTilde(1,0);
  return ret;
}

template<class DepShape, class IndepShape, class ATDep, class ATIndep>
Matrix<General, DepShape, IndepShape, ATDep> parDer(const Vector<DepShape, ATDep> &dep, const Vector<IndepShape, ATIndep> &indep) {
  Matrix<General, DepShape, IndepShape, ATDep> ret(dep.size(), indep.size());
  for(int r=0; r<dep.size(); ++r)
    for(int c=0; c<indep.size(); ++c)
      ret(r,c)=parDer(dep(r), indep(c));
  return ret;
}

template<class IndepShape, class ATDep, class ATIndep>
RowVector<IndepShape, ATDep> parDer(const ATDep &dep, const Vector<IndepShape, ATIndep> &indep) {
  RowVector<IndepShape, ATDep> ret(indep.size());
  for(int r=0; r<indep.size(); ++r)
    ret(r)=parDer(dep, indep(r));
  return ret;
}

template<class ATDep, class Shape, class ATIndep>
Matrix<General, Fixed<3>, Shape, ATDep> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const Vector<Shape, ATIndep> &x) {
  Matrix<General, Fixed<3>, Shape, ATDep> ret(3,x.size());
  for(int i=0; i<x.size(); ++i)
    ret.set(i, parDer(R, x(i)));
  return ret;
}

template<class Dep, class ATDep, class ATIndep>
Dep dirDer(const Dep &dep, const ATDep &indepdir, const ATIndep &indep) {
  return parDer(dep, indep)*indepdir;
}

template<class ATDep, class ATIndep>
Vector<Fixed<3>, ATDep> dirDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const ATDep &indepdir, const ATIndep &indep) {
  return parDer(R, indep)*indepdir;
}

template<class Dep, class ShapeIndep, class ATDep, class ATIndep>
Dep dirDer(const Dep &dep, const Vector<ShapeIndep, ATDep> &indepdir, const Vector<ShapeIndep, ATIndep> &indep) {
  Dep ret;
  if constexpr (std::is_same_v<Dep, SymbolicExpression>) {
    ret=0.0;
  }
  else if constexpr (Dep::isVector)
    ret<<=Dep(dep.size(), INIT, 0.0);
  else
    ret<<=Dep(dep.rows(), dep.cols(), INIT, 0.0);
  for(int i=0; i<indep.size(); ++i)
    ret+=parDer(dep, indep(i))*indepdir(i);
  return ret;
}

template<class ATDep, class ShapeIndep, class ATIndep>
Vector<Fixed<3>, ATDep> dirDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &dep, const Vector<ShapeIndep, ATDep> &indepdir, const Vector<ShapeIndep, ATIndep> &indep) {
  return parDer(dep, indep)*indepdir;
}

// substitute scalar in scalar expression -> defined in ast.h
template<class Ret, class A, class B>
Ret subst(const Ret &src, const A &a, const B &b) {
  if constexpr (!std::is_same_v<A, IndependentVariable>) {
    if constexpr (A::isVector) {
      if(a.size()!=b.size())
        throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
    }
    else {
      if(a.rows()!=b.rows() || a.cols()!=b.cols())
        throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
    }
  }

  if constexpr (std::is_same_v<Ret, SymbolicExpression> && std::is_same_v<A, IndependentVariable>) {
    return AST::substScalar(src, a, b);
  }
  if constexpr (std::is_same_v<Ret, SymbolicExpression> && !std::is_same_v<A, IndependentVariable>) {
    auto ret=src;
    auto aIt=a.begin();
    auto bIt=b.begin();
    for(; aIt<a.end(); ++aIt, ++bIt)
      ret=subst(ret, *aIt, *bIt);
    return ret;
  }
  if constexpr (!std::is_same_v<Ret, SymbolicExpression> && std::is_same_v<A, IndependentVariable>) {
    auto ret=src;
    for(auto retIt=ret.begin(); retIt<ret.end(); ++retIt)
      *retIt=subst(*retIt, a, b);
    return ret;
  }
  if constexpr (!std::is_same_v<Ret, SymbolicExpression> && !std::is_same_v<A, IndependentVariable>) {
    auto ret=src;
    for(auto retIt=ret.begin(); retIt<ret.end(); ++retIt) {
      auto aIt=a.begin();
      auto bIt=b.begin();
      for(; aIt<a.end(); ++aIt, ++bIt)
        *retIt=subst(*retIt, *aIt, *bIt);
    }
    return ret;
  }
}



template<class Dst, class Src>
Dst& operator^=(Dst &dst, const Src &src) {
  auto dstIt=dst.begin();
  auto srcIt=src.begin();
  for(; dstIt<dst.end(); ++dstIt, ++srcIt)
    *dstIt ^= *srcIt;
  return dst;
}



/* Class for evaluating a symbolic expression
 * This class is not copy/move-able to since it uses ByteCode which itself uses internal raw pointers (for
 * performance reasons) which cannot be copied/moved.
*/
template<class... Arg>
class Eval {
  private:
    // Type for all symbolic arguments as a tuple
    using SymTuple = std::tuple<Arg...>;
    // Type for all numeric return values as a tuple
    using NumTuple = std::tuple<typename ReplaceAT<Arg, double>::Type...>;
    // The return type: equals NumTuple if more than 1 arg is given in the ctor; for only 1 arg the tuple is skipped
    using NumRetType = std::conditional_t<std::tuple_size_v<NumTuple> == 1, std::tuple_element_t<0,NumTuple>, NumTuple>;
  public:
    ~Eval();
    Eval(const Eval &src) = delete;
    Eval(Eval &&src) = delete;
    Eval<Arg...>& operator=(const Eval &src) = delete;
    Eval<Arg...>& operator=(Eval &&src) = delete;

    // construct an evaluation object for all symbolic args
    // An arg can be a symbolic scalar, vector or matrix
    Eval(const Arg&... arg);
    // evaluate all symbolic args given by the ctor and return a tuple of corrsponding evaluated numeric values.
    // Note that the return values can be get easily using "structured binding".
    inline const NumRetType& operator()() const;
  private:
    // the numeric values as a tuple
    NumTuple numTuple;
    // hold a reference to Symbols used by the evaluation to avoid deleting these symbols (since the code has pointers to these)
    std::set<std::shared_ptr<const AST::Symbol>> symbolStore;
    // helper function to walk over all SymbolicExpression (AT=atomic type) and corresponding double value in parallel.
    // For each pair the function is called.
    template<int I=0>
    void walkAT(const SymTuple &symTuple, NumTuple &numTuple,
                const std::function<void(const SymbolicExpression&, double&)> &func);

    // members for bytecode evaluation

    std::vector<AST::ByteCode> byteCode;

    // the constructor and operator() for runtime evaluation
    void ctorByteCode(const Arg&... arg);
    inline void callByteCode() const;
};

template<class... Arg>
Eval<Arg...>::~Eval() = default;

template<class... Arg>
Eval<Arg...>::Eval(const Arg&... arg) {
  ctorByteCode(arg...);
}

template<class... Arg>
auto Eval<Arg...>::operator()() const -> const NumRetType& {
  callByteCode();
  // return the numeric values: the above call has written to its addresses
  if constexpr (std::tuple_size_v<NumTuple> == 1)
    return std::get<0>(numTuple);
  else
    return numTuple;
}

template<class... Arg>
template<int I>
void Eval<Arg...>::walkAT(const SymTuple &symTuple, NumTuple &numTuple,
                          const std::function<void(const SymbolicExpression&, double&)> &func) {
  // get the I-th symbolic and numeric variable form the tuples (this function handles the I-th tuple entry)
  // (this function is called with I=0 from the external caller)
  auto &sym = std::get<I>(symTuple);
  auto &num = std::get<I>(numTuple);
  // get the type of the I-th tuple entry
  using Sym = std::tuple_element_t<I,SymTuple>;

  if constexpr (std::is_same_v<Sym, SymbolicExpression>)
    func(sym, num);
  else {
    // fill bytecode for vector or matrix value
    // resize the vector or matrix
    if constexpr (Sym::isVector)
      num.resize(sym.size());
    else
      num.resize(sym.rows(), sym.cols());
    // get begin iterator of the vector/matrix and iterator both iterator simultaniously
    auto symIt=sym.begin();
    auto numIt=num.begin();
    for(; symIt!=sym.end(); ++symIt, ++numIt)
      func(*symIt, *numIt);
  }

  // call this function recursively for the next (I+1)-th entry in the tuple (only if there is a next entry)
  if constexpr (I+1<std::tuple_size_v<SymTuple>)
    walkAT<I+1>(symTuple, numTuple, func);
}

template<class... Arg>
void Eval<Arg...>::ctorByteCode(const Arg&... arg) {
  // first walk through all vertices in all args and count the number of operations needed
  size_t byteCodeCount=0;
  std::set<const AST::Vertex*> existingVertex1;
  walkAT(SymTuple(arg...), numTuple, [&byteCodeCount, &existingVertex1](auto &sym, auto &num) {
    sym->walkVertex([&byteCodeCount, &existingVertex1](const std::shared_ptr<const AST::Vertex>& v) {
      if(!existingVertex1.insert(v.get()).second) return;

      if(auto s=std::dynamic_pointer_cast<const AST::Symbol>(v); s)
        byteCodeCount++;
      byteCodeCount++;
    });
    byteCodeCount++;
  });
  // pre-allocate byteCode (ByteCode has not move-ctor)
  byteCode.reserve(byteCodeCount);
  std::map<const AST::Vertex*, std::vector<AST::ByteCode>::iterator> existingVertex;
  // than save all symbols in the symbolStore (the code uses a pointer to Symbol::x)
  // and add code to copy from &Symbol::x to byteCode. This copy is a "slow" operation since &Symbol::x may be far away.
  walkAT(SymTuple(arg...), numTuple, [this, &existingVertex](auto &sym, auto &num) {
    sym->walkVertex([this, &existingVertex](const std::shared_ptr<const AST::Vertex>& v) {
      if(auto s=std::dynamic_pointer_cast<const AST::Symbol>(v); s) {
        symbolStore.insert(s);
        s->dumpByteCode(byteCode, existingVertex);
      }
    });
  });
  // now add the code of the symbolic expression: this is a "fast" operation since only addresses inside of byteCode are used
  // Also store a interator to each AT (atomic type) result of the symbolic expression
  std::vector< std::vector<AST::ByteCode>::iterator > exprRet;
  walkAT(SymTuple(arg...), numTuple, [this, &existingVertex, &exprRet](auto &sym, auto &num) {
    exprRet.emplace_back(sym->dumpByteCode(byteCode, existingVertex));
  });
  // lastly add code to copy the result of each AT (the iterators from above) to the address of the return value.
  // This is again a "slow" operation since the address of the return value may be far away fromo byteCode.
  auto exprRetIt = exprRet.begin();
  walkAT(SymTuple(arg...), numTuple, [this, &exprRet, &exprRetIt](auto &sym, auto &num) {
    byteCode.emplace_back(1);
    auto it = --byteCode.end();
    it->func = [](double* r, const AST::ByteCode::Arg& a) { *r = *a[0]; };
    it->argsPtr = { (*(exprRetIt++))->retPtr };
    it->retPtr = &num;
  });
}

template<class... Arg>
void Eval<Arg...>::callByteCode() const {
#if defined(FMATVEC_DEBUG) && !defined(SWIG)
  SymbolicExpression::evalOperationsCount = byteCode.size();
#endif
  std::for_each(byteCode.begin(), byteCode.end(), [](const AST::ByteCode &bc){
    bc.func(bc.retPtr, bc.argsPtr);
  });
}

template<class RetN, class ArgN>
class FunctionWrap1VecRetToScalar : public Function<double(ArgN)>  {
  public:
    FunctionWrap1VecRetToScalar(const std::shared_ptr<Function<RetN(ArgN)>> &func_, int idx_) :
      func(func_), idx(idx_) {}
    double operator()(const ArgN &arg) override {
      return (*func)(arg)(idx);
    }
    double dirDer(const ArgN &argDir, const ArgN &arg) override {
      return func->dirDer(argDir, arg)(idx);
    }
    double dirDerDirDer(const ArgN &argDir_1, const ArgN &argDir_2, const ArgN &arg) override {
      return func->dirDerDirDer(argDir_1, argDir_2, arg)(idx);
    }
  private:
    std::shared_ptr<Function<RetN(ArgN)>> func;
    int idx;
};

template<class RetN, class Arg1N, class Arg2N>
class FunctionWrap2VecRetToScalar : public Function<double(Arg1N,Arg2N)>  {
  public:
    FunctionWrap2VecRetToScalar(const std::shared_ptr<Function<RetN(Arg1N,Arg2N)>> &func_, int idx_) :
      func(func_), idx(idx_) {}
    double operator()(const Arg1N &arg1, const Arg2N &arg2) override {
      return (*func)(arg1,arg2)(idx);
    }
    double dirDer1(const Arg1N &dir1, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer1(dir1, arg1, arg2)(idx);
    }
    double dirDer2(const Arg2N &dir2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer2(dir2, arg1, arg2)(idx);
    }
    double dirDer1DirDer1(const Arg1N &dir1_1, const Arg1N &dir1_2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer1DirDer1(dir1_1, dir1_2, arg1, arg2)(idx);
    }
    double dirDer2DirDer1(const Arg2N &dir2_1, const Arg1N &dir1_2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer2DirDer1(dir2_1, dir1_2, arg1, arg2)(idx);
    }
    double dirDer2DirDer2(const Arg2N &dir2_1, const Arg2N &dir2_2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer2DirDer2(dir2_1, dir2_2, arg1, arg2)(idx);
    }
  private:
    std::shared_ptr<Function<RetN(Arg1N,Arg2N)>> func;
    int idx;
};

template<class RetN, class ArgN>
class FunctionWrap1MatRetToScalar : public Function<double(ArgN)>  {
  public:
    FunctionWrap1MatRetToScalar(const std::shared_ptr<Function<RetN(ArgN)>> &func_, int row_, int col_) :
      func(func_), row(row_), col(col_) {}
    double operator()(const ArgN &arg) override {
      return (*func)(arg)(row,col);
    }
    double dirDer(const ArgN &argDir, const ArgN &arg) override {
      return func->dirDer(argDir, arg)(row,col);
    }
    double dirDerDirDer(const ArgN &argDir_1, const ArgN &argDir_2, const ArgN &arg) override {
      return func->dirDerDirDer(argDir_1, argDir_2, arg)(row,col);
    }
  private:
    std::shared_ptr<Function<RetN(ArgN)>> func;
    int row;
    int col;
};

template<class RetN, class Arg1N, class Arg2N>
class FunctionWrap2MatRetToScalar : public Function<double(Arg1N,Arg2N)>  {
  public:
    FunctionWrap2MatRetToScalar(const std::shared_ptr<Function<RetN(Arg1N,Arg2N)>> &func_, int row_, int col_) :
      func(func_), row(row_), col(col_) {}
    double operator()(const Arg1N &arg1, const Arg2N &arg2) override {
      return (*func)(arg1,arg2)(row,col);
    }
    double dirDer1(const Arg1N &dir1, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer1(dir1, arg1, arg2)(row,col);
    }
    double dirDer2(const Arg2N &dir2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer2(dir2, arg1, arg2)(row,col);
    }
    double dirDer1DirDer1(const Arg1N &dir1_1, const Arg1N &dir1_2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer1DirDer1(dir1_1, dir1_2, arg1, arg2)(row,col);
    }
    double dirDer2DirDer1(const Arg2N &dir2_1, const Arg1N &dir1_2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer2DirDer1(dir2_1, dir1_2, arg1, arg2)(row,col);
    }
    double dirDer2DirDer2(const Arg2N &dir2_1, const Arg2N &dir2_2, const Arg1N &arg1, const Arg2N &arg2) override {
      return func->dirDer2DirDer2(dir2_1, dir2_2, arg1, arg2)(row,col);
    }
  private:
    std::shared_ptr<Function<RetN(Arg1N,Arg2N)>> func;
    int row;
    int col;
};

template<class Func, class ArgS>
typename ReplaceAT<typename std::function<Func>::result_type,SymbolicExpression>::Type symbolicFuncWrapVecAndMatRet(
  const std::shared_ptr<Function<Func>> &func, const ArgS &arg, int size1=0, int size2=0) {
  using RetN = typename std::function<Func>::result_type;
  using RetS = typename ReplaceAT<RetN,SymbolicExpression>::Type;
  using ArgN = typename ReplaceAT<ArgS,double>::Type;

  auto hasSize = boost::hana::is_valid([](auto&& vec) -> decltype(vec.size()) {});

  RetS ret;
  if constexpr (hasSize(ret)) {
    // vector
    int size=func->getRetSize().first;
    if(size==0 && size1!=0)
      size=size1;
    else if(size==0)
      size=(*func)(ArgN()).size();
    ret.resize(size);
    for(int i=0; i<size; ++i)
      ret(i) = AST::SymbolicFuncWrapArg1<double(ArgN), ArgS>::
               call(std::make_shared<FunctionWrap1VecRetToScalar<RetN,ArgN>>(func, i), arg);
    return ret;
  }
  else {
    // matrix
    auto size=func->getRetSize();
    if(size.first==0 && size1!=0)
      size.first=size1;
    else if(size.first==0)
      size.first=(*func)(ArgN()).rows();
    if(size.second==0 && size2!=0)
      size.second=size2;
    else if(size.second==0)
      size.second=(*func)(ArgN()).cols();
    ret.resize(size.first, size.second);
    for(int r=0; r<size.first; ++r)
      for(int c=0; c<size.second; ++c)
        ret(r,c) = AST::SymbolicFuncWrapArg1<double(ArgN), ArgS>::call(std::make_shared<FunctionWrap1MatRetToScalar<RetN,ArgN>>(
                   func, r, c), arg);
    return ret;
  }
}

template<class Func, class Arg1S, class Arg2S>
typename ReplaceAT<typename std::function<Func>::result_type,SymbolicExpression>::Type symbolicFuncWrapVecAndMatRet(
  const std::shared_ptr<Function<Func>> &func, const Arg1S &arg1, const Arg2S &arg2, int size1=0, int size2=0) {
  using RetN = typename std::function<Func>::result_type;
  using RetS = typename ReplaceAT<RetN,SymbolicExpression>::Type;
  using Arg1N = typename ReplaceAT<Arg1S,double>::Type;
  using Arg2N = typename ReplaceAT<Arg2S,double>::Type;

  auto hasSize = boost::hana::is_valid([](auto&& vec) -> decltype(vec.size()) {});

  RetS ret;
  if constexpr (hasSize(ret)) {
    // vector
    int size=func->getRetSize().first;
    if(size==0 && size1!=0)
      size=size1;
    else if(size==0)
      size=(*func)(Arg1N(), Arg2N()).size();
    ret.resize(size);
    for(int i=0; i<size; ++i)
      ret(i) = AST::SymbolicFuncWrapArg2<double(Arg1N,Arg2N), Arg1S, Arg2S>::
               call(std::make_shared<FunctionWrap2VecRetToScalar<RetN,Arg1N,Arg2N>>(func, i), arg1, arg2);
    return ret;
  }
  else {
    // matrix
    auto size=func->getRetSize();
    if(size.first==0 && size1!=0)
      size.first=size1;
    else if(size.first==0)
      size.first=(*func)(Arg1N(),Arg2N()).rows();
    if(size.second==0 && size2!=0)
      size.second=size2;
    else if(size.second==0)
      size.second=(*func)(Arg1N(),Arg2N()).cols();
    ret.resize(size.first, size.second);
    for(int r=0; r<size.first; ++r)
      for(int c=0; c<size.second; ++c)
        ret(r,c) = AST::SymbolicFuncWrapArg2<double(Arg1N,Arg2N), Arg1S, Arg2S>::
                   call(std::make_shared<FunctionWrap2MatRetToScalar<RetN,Arg1N,Arg2N>>(func, r, c), arg1, arg2);
    return ret;
  }
}

/*! Create a symbolic function from the native function func.
 * If func provides derivatives this symbolic function can also be differentiated using parDir(...).
 * Note that func should use a LRUCache, at least if return value of func is a vector (R^n) or a matrix (R^rxc) since
 * func will be called n or r*c times with the same arguments.
 */
template<class Func, class ArgS>
typename ReplaceAT<typename std::function<Func>::result_type,SymbolicExpression>::Type symbolicFunc(
  const std::shared_ptr<Function<Func>> &func, const ArgS &arg, int size1=0, int size2=0) {
  using RetN = typename std::function<Func>::result_type;
  using ArgN = typename ReplaceAT<ArgS,double>::Type;
  if constexpr (std::is_same_v<RetN, double>)
    return AST::SymbolicFuncWrapArg1<double(ArgN), ArgS>::call(func, arg);
  else
    return symbolicFuncWrapVecAndMatRet<Func, ArgS>(func, arg, size1, size2);
}

/*! Create a symbolic function from the native function func.
 * If func provides derivatives this symbolic function can also be differentiated using parDir(...).
 * Note that func should use a LRUCache, at least if return value of func is a vector (R^n) or a matrix (R^rxc) since
 * func will be called n or r*c times with the same arguments.
 */
template<class Func, class Arg1S, class Arg2S>
typename ReplaceAT<typename std::function<Func>::result_type,SymbolicExpression>::Type symbolicFunc(
  const std::shared_ptr<Function<Func>> &func, const Arg1S &arg1, const Arg2S &arg2, int size1=0, int size2=0) {
  using RetN = typename std::function<Func>::result_type;
  using Arg1N = typename ReplaceAT<Arg1S,double>::Type;
  using Arg2N = typename ReplaceAT<Arg2S,double>::Type;
  if constexpr (std::is_same_v<RetN, double>)
    return AST::SymbolicFuncWrapArg2<double(Arg1N, Arg2N), Arg1S, Arg2S>::call(func, arg1, arg2);
  else
    return symbolicFuncWrapVecAndMatRet<Func, Arg1S, Arg2S>(func, arg1, arg2, size1, size2);
}

}

#endif
