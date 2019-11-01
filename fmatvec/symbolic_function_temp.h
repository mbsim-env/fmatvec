template<TEMPLATE>
class SymbolicFunction<RET(ARG)> : public Function<RET(ARG)> {
  public:
    using DRetDArg = typename Function<RET(ARG)>::DRetDArg;
    using DRetDDir = typename Function<RET(ARG)>::DRetDDir;
    using DDRetDDArg = typename Function<RET(ARG)>::DDRetDDArg;
    using ArgS = typename ReplaceAT<ARG, IndependentVariable>::Type;
    using RetS = typename ReplaceAT<RET, SymbolicExpression>::Type;

    SymbolicFunction(const ArgS &argS_, const RetS &retS_);
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
#ifdef PARDER
    typename ReplaceAT<DRetDArg, SymbolicExpression>::Type pd;
#endif
    typename ReplaceAT<DRetDDir, SymbolicExpression>::Type dd;
#ifdef PARDERPARDER
    typename ReplaceAT<DDRetDDArg, SymbolicExpression>::Type pdpd;
#endif
#ifdef PARDER
    typename ReplaceAT<DRetDArg, SymbolicExpression>::Type pddd;
#endif
    typename ReplaceAT<DRetDDir, SymbolicExpression>::Type dddd;
    bool isParDerConst = false;
};

template<TEMPLATE>
SymbolicFunction<RET(ARG)>::SymbolicFunction(const ArgS &argS_, const RetS &retS_) :
  argS(argS_), argDirS(ARGDIRINIT), argDirS_2(ARGDIRINIT), retS(retS_) {
    
#ifdef PARDER
  pd=fmatvec::parDer(retS, argS);
#endif
  dd=fmatvec::dirDer(retS, argDirS, argS);
#ifdef PARDERPARDER
  pdpd=fmatvec::parDer(fmatvec::parDer(retS, argS), argS);
#endif
#ifdef PARDER
  pddd=fmatvec::dirDer(fmatvec::parDer(retS, argS), argDirS, argS);
#endif
  dddd=fmatvec::dirDer(fmatvec::dirDer(retS, argDirS, argS), argDirS_2, argS);

#ifdef PARDER
  isParDerConst=false;
#endif
}

template<TEMPLATE>
RET SymbolicFunction<RET(ARG)>::operator()(const ARG &arg) {
  argS&=arg;//mfmf use something better here, this disables the caching of ast!?
  return eval(retS);
}

#ifdef PARDER
template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::parDer(const ARG &arg) -> DRetDArg {
  argS&=arg;
  return eval(pd);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::dirDer(const ARG &argDir, const ARG &arg) -> DRetDDir {
  argDirS&=argDir;
  argS&=arg;
  return eval(dd);
}

#ifdef PARDERPARDER
template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::parDerParDer(const ARG &arg) -> DDRetDDArg {
  argS&=arg;
  return eval(pdpd);
}
#endif

#ifdef PARDER
template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::parDerDirDer(const ARG &argDir, const ARG &arg) -> DRetDArg {
  argDirS&=argDir;
  argS&=arg;
  return eval(pddd);
}
#endif

template<TEMPLATE>
auto SymbolicFunction<RET(ARG)>::dirDerDirDer(const ARG &argDir_1, const ARG &argDir_2, const ARG &arg) -> DRetDDir {
  argDirS&=argDir_1;
  argDirS_2&=argDir_2;
  argS&=arg;
  return eval(dddd);
}

template<TEMPLATE>
bool SymbolicFunction<RET(ARG)>::constParDer() const {
  return isParDerConst;
}
