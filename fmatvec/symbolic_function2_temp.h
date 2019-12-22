template<TEMPLATE>
class SymbolicFunction<RET(ARG1, ARG2)> : public virtual Function<RET(ARG1, ARG2)> {
  public:
    using DRetDArg1 = typename Function<RET(ARG1, ARG2)>::DRetDArg1;
    using DRetDArg2 = typename Function<RET(ARG1, ARG2)>::DRetDArg2;
    using DRetDDir1 = typename Function<RET(ARG1, ARG2)>::DRetDDir1;
    using DRetDDir2 = typename Function<RET(ARG1, ARG2)>::DRetDDir2;
    using DDRetDDArg1 = typename Function<RET(ARG1, ARG2)>::DDRetDDArg1;
    using DDRetDDArg2 = typename Function<RET(ARG1, ARG2)>::DDRetDDArg2;
    using DDRetDArg1DArg2 = typename Function<RET(ARG1, ARG2)>::DDRetDArg1DArg2;
    using Arg1S = typename ReplaceAT<ARG1, IndependentVariable>::Type;
    using Arg2S = typename ReplaceAT<ARG2, IndependentVariable>::Type;
    using RetS = typename ReplaceAT<RET, SymbolicExpression>::Type;

    SymbolicFunction();
    SymbolicFunction(const Arg1S &arg1S_, const Arg2S &arg2S_, const RetS &retS_); // calls init() at the end
    void setIndependentVariable1(const Arg1S &arg1S_);
    void setIndependentVariable2(const Arg2S &arg2S_);
    void setDependentFunction(const RetS &retS_);
    void init(); // must be called after setIndependentVariable/setDependentFunction.

    std::pair<int, int> getRetSize() const override;
    int getArg1Size() const override;
    int getArg2Size() const override;

    RET operator()(const ARG1 &arg1, const ARG2 &arg2) override;
#ifdef PARDER1
    DRetDArg1 parDer1(const ARG1 &arg1, const ARG2 &arg2) override;
#endif
    DRetDDir1 dirDer1(const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#ifdef PARDER2
    DRetDArg2 parDer2(const ARG1 &arg1, const ARG2 &arg2) override;
#endif
    DRetDDir2 dirDer2(const ARG2 &arg2Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#ifdef PARDER1PARDER1
    DDRetDDArg1 parDer1ParDer1(const ARG1 &arg1, const ARG2 &arg2) override;
#endif
#ifdef PARDER1
    DRetDArg1 parDer1DirDer1(const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#endif
    DRetDDir1 dirDer1DirDer1(const ARG1 &arg1Dir_1, const ARG1 &arg1Dir_2, const ARG1 &arg1, const ARG2 &arg2) override;
#ifdef PARDER2PARDER2
    DDRetDDArg2 parDer2ParDer2(const ARG1 &arg1, const ARG2 &arg2) override;
#endif
#ifdef PARDER2
    DRetDArg2 parDer2DirDer2(const ARG2 &arg2Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#endif
    DRetDDir2 dirDer2DirDer2(const ARG2 &arg2Dir_1, const ARG2 &arg2Dir_2, const ARG1 &arg1, const ARG2 &arg2) override;
#ifdef PARDER1PARDER2
    DDRetDArg1DArg2 parDer1ParDer2(const ARG1 &arg1, const ARG2 &arg2) override;
#endif
#ifdef PARDER1
    DRetDArg1 parDer1DirDer2(const ARG2 &arg2Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#endif
    DRetDDir2 dirDer2DirDer1(const ARG2 &arg2Dir, const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#ifdef PARDER2
    DRetDArg2 parDer2DirDer1(const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) override;
#endif
    bool constParDer1() const override;
    bool constParDer2() const override;

  protected:

    Arg1S arg1S;
    Arg2S arg2S;
    Arg1S argDir1S;
    Arg1S argDir1S_2;
    Arg2S argDir2S;
    Arg2S argDir2S_2;
    RetS retS;
#ifdef PARDER1
    typename ReplaceAT<DRetDArg1, SymbolicExpression>::Type pd1;
#endif
    typename ReplaceAT<DRetDDir1, SymbolicExpression>::Type dd1;
#ifdef PARDER2
    typename ReplaceAT<DRetDArg2, SymbolicExpression>::Type pd2;
#endif
    typename ReplaceAT<DRetDDir2, SymbolicExpression>::Type dd2;
#ifdef PARDER1PARDER1
    typename ReplaceAT<DDRetDDArg1, SymbolicExpression>::Type pd1pd1;
#endif
#ifdef PARDER1
    typename ReplaceAT<DRetDArg1, SymbolicExpression>::Type pd1dd1;
#endif
    typename ReplaceAT<DRetDDir1, SymbolicExpression>::Type dd1dd1;
#ifdef PARDER2PARDER2
    typename ReplaceAT<DDRetDDArg2, SymbolicExpression>::Type pd2pd2;
#endif
#ifdef PARDER2
    typename ReplaceAT<DRetDArg2, SymbolicExpression>::Type pd2dd2;
#endif
    typename ReplaceAT<DRetDDir2, SymbolicExpression>::Type dd2dd2;
#ifdef PARDER1PARDER2
    typename ReplaceAT<DDRetDArg1DArg2, SymbolicExpression>::Type pd1pd2;
#endif
#ifdef PARDER1
    typename ReplaceAT<DRetDArg1, SymbolicExpression>::Type pd1dd2;
#endif
    typename ReplaceAT<DRetDDir2, SymbolicExpression>::Type dd2dd1;
#ifdef PARDER2
    typename ReplaceAT<DRetDArg2, SymbolicExpression>::Type pd2dd1;
#endif
    bool isParDer1Const = false;
    bool isParDer2Const = false;
};

template<TEMPLATE>
SymbolicFunction<RET(ARG1, ARG2)>::SymbolicFunction() {
#ifdef PARDER
  isParDer1Const=false;
  isParDer2Const=false;
#endif
}

template<TEMPLATE>
SymbolicFunction<RET(ARG1, ARG2)>::SymbolicFunction(const Arg1S &arg1S_, const Arg2S &arg2S_, const RetS &retS_) :
  arg1S(arg1S_), arg2S(arg2S_), retS(retS_) {
#ifdef PARDER
  isParDer1Const=false;
  isParDer2Const=false;
#endif
  init();
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG1, ARG2)>::setIndependentVariable1(const Arg1S &arg1S_) {
  arg1S<<=arg1S_;
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG1, ARG2)>::setIndependentVariable2(const Arg2S &arg2S_) {
  arg2S<<=arg2S_;
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG1, ARG2)>::setDependentFunction(const RetS &retS_) {
  retS<<=retS_;
}

template<TEMPLATE>
void SymbolicFunction<RET(ARG1, ARG2)>::init() {
  Helper<Arg1S>::initIndep(argDir1S, Helper<Arg1S>::size1(arg1S));
  Helper<Arg1S>::initIndep(argDir1S, Helper<Arg1S>::size1(arg1S));
  Helper<Arg1S>::initIndep(argDir1S_2, Helper<Arg1S>::size1(arg1S));
  Helper<Arg2S>::initIndep(argDir2S, Helper<Arg2S>::size1(arg2S));
  Helper<Arg2S>::initIndep(argDir2S, Helper<Arg2S>::size1(arg2S));
  Helper<Arg2S>::initIndep(argDir2S_2, Helper<Arg2S>::size1(arg2S));
#ifdef PARDER1
  pd1<<=fmatvec::parDer(retS, arg1S);
#endif
  dd1<<=fmatvec::dirDer(retS, argDir1S, arg1S);
#ifdef PARDER2
  pd2<<=fmatvec::parDer(retS, arg2S);
#endif
  dd2<<=fmatvec::dirDer(retS, argDir2S, arg2S);
#ifdef PARDER1PARDER1
  pd1pd1<<=fmatvec::parDer(fmatvec::parDer(retS, arg1S), arg1S);
#endif
#ifdef PARDER1
  pd1dd1<<=fmatvec::dirDer(fmatvec::parDer(retS, arg1S), argDir1S, arg1S);
#endif
  dd1dd1<<=fmatvec::dirDer(fmatvec::dirDer(retS, argDir1S, arg1S), argDir1S_2, arg1S);
#ifdef PARDER2PARDER2
  pd2pd2<<=fmatvec::parDer(fmatvec::parDer(retS, arg2S), arg2S);
#endif
#ifdef PARDER2
  pd2dd2<<=fmatvec::dirDer(fmatvec::parDer(retS, arg2S), argDir2S, arg2S);
#endif
  dd2dd2<<=fmatvec::dirDer(fmatvec::dirDer(retS, argDir2S, arg2S), argDir2S_2, arg2S);
#ifdef PARDER1PARDER2
  pd1pd2<<=fmatvec::parDer(fmatvec::parDer(retS, arg1S), arg2S);
#endif
#ifdef PARDER1
  pd1dd2<<=fmatvec::dirDer(fmatvec::parDer(retS, arg1S), argDir2S, arg2S);
#endif
  dd2dd1<<=fmatvec::dirDer(fmatvec::dirDer(retS, argDir2S, arg2S), argDir1S, arg1S);
#ifdef PARDER2
  pd2dd1<<=fmatvec::dirDer(fmatvec::parDer(retS, arg2S), argDir1S, arg1S);
#endif
}

template<TEMPLATE>
std::pair<int, int> SymbolicFunction<RET(ARG1, ARG2)>::getRetSize() const {
  return std::make_pair(Helper<RetS>::size1(retS), Helper<RetS>::size2(retS));
}

template<TEMPLATE>
int SymbolicFunction<RET(ARG1, ARG2)>::getArg1Size() const{
  return Helper<Arg1S>::size1(arg1S);
}

template<TEMPLATE>
int SymbolicFunction<RET(ARG1, ARG2)>::getArg2Size() const {
  return Helper<Arg2S>::size1(arg2S);
}

template<TEMPLATE>
RET SymbolicFunction<RET(ARG1, ARG2)>::operator()(const ARG1 &arg1, const ARG2 &arg2) {
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(retS);
}

#ifdef PARDER1
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer1(const ARG1 &arg1, const ARG2 &arg2) -> DRetDArg1 {
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd1);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::dirDer1(const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDDir1 {
  argDir1S&=arg1Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(dd1);
}

#ifdef PARDER2
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer2(const ARG1 &arg1, const ARG2 &arg2) -> DRetDArg2 {
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd2);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::dirDer2(const ARG2 &arg2Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDDir2 {
  argDir2S&=arg2Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(dd2);
}

#ifdef PARDER1PARDER1
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer1ParDer1(const ARG1 &arg1, const ARG2 &arg2) -> DDRetDDArg1 {
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd1pd1);
}
#endif

#ifdef PARDER1
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer1DirDer1(const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDArg1 {
  argDir1S&=arg1Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd1dd1);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::dirDer1DirDer1(const ARG1 &arg1Dir_1, const ARG1 &arg1Dir_2, const ARG1 &arg1, const ARG2 &arg2) -> DRetDDir1 {
  argDir1S&=arg1Dir_1;
  argDir1S_2&=arg1Dir_2;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(dd1dd1);
}

#ifdef PARDER2PARDER2
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer2ParDer2(const ARG1 &arg1, const ARG2 &arg2) -> DDRetDDArg2 {
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd2pd2);
}
#endif

#ifdef PARDER2
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer2DirDer2(const ARG2 &arg2Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDArg2 {
  argDir2S&=arg2Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd2dd2);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::dirDer2DirDer2(const ARG2 &arg2Dir_1, const ARG2 &arg2Dir_2, const ARG1 &arg1, const ARG2 &arg2) -> DRetDDir2 {
  argDir2S&=arg2Dir_1;
  argDir2S_2&=arg2Dir_2;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(dd2dd2);
}

#ifdef PARDER1PARDER2
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer1ParDer2(const ARG1 &arg1, const ARG2 &arg2) -> DDRetDArg1DArg2 {
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd1pd2);
}
#endif

#ifdef PARDER1
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer1DirDer2(const ARG2 &arg2Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDArg1 {
  argDir2S&=arg2Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd1dd2);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::dirDer2DirDer1(const ARG2 &arg2Dir, const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDDir2 {
  argDir2S&=arg2Dir;
  argDir1S&=arg1Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(dd2dd1);
}

#ifdef PARDER2
template<TEMPLATE>
auto SymbolicFunction<RET(ARG1, ARG2)>::parDer2DirDer1(const ARG1 &arg1Dir, const ARG1 &arg1, const ARG2 &arg2) -> DRetDArg2 {
  argDir1S&=arg1Dir;
  arg1S&=arg1;
  arg2S&=arg2;
  return eval(pd2dd1);
}
#endif

template<TEMPLATE>
bool SymbolicFunction<RET(ARG1, ARG2)>::constParDer1() const {
  return isParDer1Const;
}

template<TEMPLATE>
bool SymbolicFunction<RET(ARG1, ARG2)>::constParDer2() const {
  return isParDer2Const;
}
