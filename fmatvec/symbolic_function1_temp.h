template<TEMPLATE>
class SymbolicFunction<RET(ARG)> : public virtual Function<RET(ARG)> {
  public:
    using DRetDArg = typename Function<RET(ARG)>::DRetDArg;
    using DRetDDir = typename Function<RET(ARG)>::DRetDDir;
    using DDRetDDArg = typename Function<RET(ARG)>::DDRetDDArg;
    using ArgS = typename ReplaceAT<ARG, IndependentVariable>::Type;
    using RetS = typename ReplaceAT<RET, SymbolicExpression>::Type;

    SymbolicFunction();
    SymbolicFunction(const ArgS &argS_, const RetS &retS_); // calls init() at the end
    void setIndependentVariable(const ArgS &argS_);
    void setDependentFunction(const RetS &retS_);
    ArgS& getIndependentVariable();
    RetS& getDependentFunction();
    void init(); // must be called after setIndependentVariable/setDependentFunction.

    std::pair<int, int> getRetSize() const override;
    int getArgSize() const override;

    RET operator()(const ARG &arg) override;
#ifdef PARDER
    DRetDArg parDer(const ARG &arg) override;
#endif
    DRetDDir dirDer(const ARG &argDir, const ARG &arg) override;
#ifdef PARDERPARDER
    DDRetDDArg parDerParDer(const ARG &arg) override;
#endif
#ifdef PARDER
    DRetDArg parDerDirDer(const ARG &argDir, const ARG &arg) override;
#endif
    DRetDDir dirDerDirDer(const ARG &argDir_1, const ARG &argDir_2, const ARG &arg) override;
    bool constParDer() const override;

  protected:

    ArgS argS;
    ArgS argDirS;
    ArgS argDirS_2;
    RetS retS;
    std::unique_ptr<Eval<decltype(retS)>> retSEval;
#ifdef PARDER
    std::unique_ptr<Eval<typename ReplaceAT<DRetDArg, SymbolicExpression>::Type>> pdEval;
#endif
    std::unique_ptr<Eval<typename ReplaceAT<DRetDDir, SymbolicExpression>::Type>> ddEval;
#ifdef PARDERPARDER
    std::unique_ptr<Eval<typename ReplaceAT<DDRetDDArg, SymbolicExpression>::Type>> pdpdEval;
#endif
#ifdef PARDER
    std::unique_ptr<Eval<typename ReplaceAT<DRetDArg, SymbolicExpression>::Type>> pdddEval;
#endif
    std::unique_ptr<Eval<typename ReplaceAT<DRetDDir, SymbolicExpression>::Type>> ddddEval;
    bool isParDerConst = false;
};

template<TEMPLATE>
SymbolicFunction<RET(ARG)>::SymbolicFunction() {
#ifdef PARDER
  isParDerConst=false;
#endif
}

template<TEMPLATE>
SymbolicFunction<RET(ARG)>::SymbolicFunction(const ArgS &argS_, const RetS &retS_) : argS(argS_), retS(retS_) {
#ifdef PARDER
  isParDerConst=false;
#endif
  init();
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG)>::setIndependentVariable(const ArgS &argS_) {
  argS<<=argS_;
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG)>::setDependentFunction(const RetS &retS_) {
  retS<<=retS_;
}

template<TEMPLATE>
typename SymbolicFunction<RET(ARG)>::ArgS& SymbolicFunction<RET(ARG)>::getIndependentVariable() {
  return argS;
}

template<TEMPLATE>
typename SymbolicFunction<RET(ARG)>::RetS& SymbolicFunction<RET(ARG)>::getDependentFunction() {
  return retS;
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG)>::init() {
  Helper<ArgS>::initIndep(argDirS, Helper<ArgS>::size1(argS));
  Helper<ArgS>::initIndep(argDirS, Helper<ArgS>::size1(argS));
  Helper<ArgS>::initIndep(argDirS_2, Helper<ArgS>::size1(argS));
  retSEval.reset(new Eval<decltype(retS)>{retS});
#ifdef PARDER
  pdEval.reset(new Eval<typename ReplaceAT<DRetDArg, SymbolicExpression>::Type>{fmatvec::parDer(retS, argS)});
#endif
  ddEval.reset(new Eval<typename ReplaceAT<DRetDDir, SymbolicExpression>::Type>{fmatvec::dirDer(retS, argDirS*1, argS)});
#ifdef PARDERPARDER
  pdpdEval.reset(new Eval<typename ReplaceAT<DDRetDDArg, SymbolicExpression>::Type>{fmatvec::parDer(fmatvec::parDer(retS, argS), argS)});
#endif
#ifdef PARDER
  pdddEval.reset(new Eval<typename ReplaceAT<DRetDArg, SymbolicExpression>::Type>{fmatvec::dirDer(fmatvec::parDer(retS, argS), argDirS*1, argS)});
#endif
  ddddEval.reset(new Eval<typename ReplaceAT<DRetDDir, SymbolicExpression>::Type>{fmatvec::dirDer(fmatvec::dirDer(retS, argDirS*1, argS), argDirS_2*1, argS)});
}

template<TEMPLATE>
std::pair<int, int> SymbolicFunction<RET(ARG)>::getRetSize() const {
  return std::make_pair(Helper<RetS>::size1(retS), Helper<RetS>::size2(retS));
}

template<TEMPLATE>
int SymbolicFunction<RET(ARG)>::getArgSize() const {
  return Helper<ArgS>::size1(argS);
}

template<TEMPLATE>
RET SymbolicFunction<RET(ARG)>::operator()(const ARG &arg) {
  argS^=arg;
  return (*retSEval)();
}

#ifdef PARDER
template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::parDer(const ARG &arg) -> DRetDArg {
  argS^=arg;
  return (*pdEval)();
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::dirDer(const ARG &argDir, const ARG &arg) -> DRetDDir {
  argDirS^=argDir;
  argS^=arg;
  return (*ddEval)();
}

#ifdef PARDERPARDER
template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::parDerParDer(const ARG &arg) -> DDRetDDArg {
  argS^=arg;
  return (*pdpdEval)();
}
#endif

#ifdef PARDER
template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::parDerDirDer(const ARG &argDir, const ARG &arg) -> DRetDArg {
  argDirS^=argDir;
  argS^=arg;
  return (*pdddEval)();
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::dirDerDirDer(const ARG &argDir_1, const ARG &argDir_2, const ARG &arg) -> DRetDDir {
  argDirS^=argDir_1;
  argDirS_2^=argDir_2;
  argS^=arg;
  return (*ddddEval)();
}

template<TEMPLATE>
bool SymbolicFunction<RET(ARG)>::constParDer() const {
  return isParDerConst;
}
