#include <cfenv>
#include <cassert>
#include <iostream>
#include "fmatvec/ast.h"
#include "fmatvec/fmatvec.h"
#include "fmatvec/linear_algebra_complex.h"

using namespace std;
using namespace fmatvec;

void checkDoubleComplex() {
  // check compile of OperatorResult for double / complex<double>
  Matrix<General, Var, Var, complex<double>> M(3,3);
  Vector<Var, double> v(3);
  Vector<Var, complex<double>> r=M*v;
}

void checkComplex() {
  // check operators for mat/vec with double /complex<double>
  typedef Vector<Var, complex<double>> VecComplex;
  typedef Matrix<General, Var, Var, complex<double>> MatComplex;

  VecComplex v2;
  VecComplex v(3);
  v(1)=2.0;
  v(2)=3.0;
  v(2)=complex<double>(2.2, 3.3);
  complex<double> a(3.0);
  VecComplex r1=v*a;
  VecComplex r2=a*v;
  VecComplex r3=v*3.0;
  VecComplex r4=v*3;
  VecComplex r5=3.0*v;
  VecComplex r6=3*v;
  VecComplex r11=v/a;
  VecComplex r31=v/3.0;
  VecComplex r41=v/3;
  VecComplex r12=v+v;
  VecComplex r22=v-v;
  MatComplex m2;
  MatComplex m(3,3);
  m(0,0)=2.0;
  m(0,2)=3.0;
  m(1,0)=4.0;
  m(1,1)=5.0;
  m(1,2)=6.0;
  m(2,0)=7.0;
  m(2,1)=8.0;
  m(2,2)=9.0;
  m(2,2)=complex<double>(2.2, 3.3);
  MatComplex _r1=m*a;
  MatComplex _r2=a*m;
  MatComplex _r3=m*3.0;
  MatComplex _r4=m*3;
  MatComplex _r5=3.0*m;
  MatComplex _r6=3*m;
  MatComplex _r11=m/a;
  MatComplex _r31=m/3.0;
  MatComplex _r41=m/3;
  MatComplex _r12=m+m;
  MatComplex _r22=m-m;
  MatComplex _rmm=m*m;
  VecComplex _rmv=m*v;

  // dump output
  cout<<"complex 0 = "<<_rmv(0)<<endl;
  cout<<"complex 1 = "<<_rmv(1)<<endl;
  cout<<"complex 2 = "<<_rmv(2)<<endl;

  complex<double> v0ref( 6.60,  9.90);
  complex<double> v1ref(23.20, 19.80);
  complex<double> v2ref( 9.95, 14.52);
}

void checkSym(double &pdn0Value,
              string &vSer,
              string &a_Ser,
              string &a5Ser,
              string &a6Ser,
              string &pdn0Ser) {
  {
    // check all vector variants as expression and independent
    RowVector<Fixed<3>, SymbolicExpression> rvf;
    RowVector<Var, SymbolicExpression> rvv(3);
    RowVector<Ref, SymbolicExpression> rvr(3);
    Vector<Fixed<3>, SymbolicExpression> vf;
    Vector<Var, SymbolicExpression> vv(3);
    Vector<Ref, SymbolicExpression> vr(3);

    Vector<Fixed<3>, IndependentVariable> vf_indep;
    Vector<Var, IndependentVariable> vv_indep(3);
    Vector<Ref, IndependentVariable> vr_indep(3);
  }

  typedef Vector<Var, SymbolicExpression> VecSymExpr;
  typedef Vector<Var, IndependentVariable> VecIndep;
  typedef Matrix<General, Var, Var, SymbolicExpression> MatSymExpr;

  {
    VecSymExpr v2;
    VecSymExpr v(3);
  }
  // check operators for mat/vec with IndependentVariable / SymbolicExpression
  VecSymExpr v2;
  VecSymExpr v3(3);
  cout<<"v3 "<<v3<<endl;
  VecIndep v(3, NONINIT);
  SymbolicExpression a(3.0);
  VecSymExpr r1=v*a;
  VecSymExpr r2=a*v;
  VecSymExpr r3=v*3.0;
  VecSymExpr r4=v*3;
  VecSymExpr r5=3.0*v;
  VecSymExpr r6=3*v;
  VecSymExpr r11=v/a;
  VecSymExpr r31=v/3.0;
  VecSymExpr r41=v/3;
  VecSymExpr r12=v+v;
  VecSymExpr r22=v-v;
  MatSymExpr m2;
  MatSymExpr m3(3,3);
  cout<<"m3 "<<m3<<endl;
  MatSymExpr m(3,3);
  IndependentVariable a_=IndependentVariable();
  m(0,1)=a_;
  m(0,0)=2.0;
  m(0,2)=3.0;
  m(1,0)=4.0;
  m(1,1)=5.0;
  m(1,2)=6.0;
  m(2,0)=7.0;
  m(2,1)=8.0;
  m(2,2)=9.0;
  auto a5=IndependentVariable();
  auto a6=IndependentVariable();
  m(2,2)=a5+a6;
  m(2,0)=log(pow(v(0), 3)+v(0));
  MatSymExpr _r1=m*a;
  MatSymExpr _r2=a*m;
  MatSymExpr _r3=m*3.0;
  MatSymExpr _r4=m*3;
  MatSymExpr _r5=3.0*m;
  MatSymExpr _r6=3*m;
  MatSymExpr _r11=m/a;
  MatSymExpr _r31=m/3.0;
  MatSymExpr _r41=m/3;
  MatSymExpr _r12=m+m;
  MatSymExpr _r22=m-m;
  MatSymExpr _rmm=m*m;
  VecSymExpr _rmv=m*v;
  SymbolicExpression norm=nrm2(_rmv);

  // check partial derivatives
  MatSymExpr pd=parDer(_rmv, v);
  RowVector<Var, SymbolicExpression> pdn=parDer(norm, v);

  // dump output
  cout<<"indepVar 0 = "<<v(0)<<endl;
  cout<<"indepVar 1 = "<<v(1)<<endl;
  cout<<"indepVar 2 = "<<v(2)<<endl;
  cout<<"depVar 0 = "<<_rmv(0)<<endl;
  cout<<"depVar 1 = "<<_rmv(1)<<endl;
  cout<<"depVar 2 = "<<_rmv(2)<<endl;
  cout<<"norm = "<<norm<<endl;
  cout<<"depVar0_indepVar 0 = "<<pd(0,0)<<endl;
  cout<<"depVar0_indepVar 1 = "<<pd(0,1)<<endl;
  cout<<"depVar0_indepVar 2 = "<<pd(0,2)<<endl;
  cout<<"depVar1_indepVar 0 = "<<pd(1,0)<<endl;
  cout<<"depVar1_indepVar 1 = "<<pd(1,1)<<endl;
  cout<<"depVar1_indepVar 2 = "<<pd(1,2)<<endl;
  cout<<"depVar2_indepVar 0 = "<<pd(2,0)<<endl;
  cout<<"depVar2_indepVar 1 = "<<pd(2,1)<<endl;
  cout<<"depVar2_indepVar 2 = "<<pd(2,2)<<endl;
  cout<<"norm_indepVar 0 = "<<pdn(0)<<endl;
  cout<<"norm_indepVar 1 = "<<pdn(1)<<endl;
  cout<<"norm_indepVar 2 = "<<pdn(2)<<endl;

  // check evaluation
  Vec3 vdouble3({1.1,2.2,3.3});
  v&=vdouble3;
  a_&=4.4;
  a5&=5.5;
  a6&=6.6;
#ifndef NDEBUG
  SymbolicExpression::evalOperationsCount=0;
#endif
  cout<<"frist eval"<<endl;
  cout<<"pdn.eval = "<<eval(pdn)<<endl;
#ifndef NDEBUG
  cout<<"number of operations evaluated = "<<SymbolicExpression::evalOperationsCount<<endl;
#endif
#ifndef NDEBUG
  SymbolicExpression::evalOperationsCount=0;
#endif
  cout<<"second eval with same values for independent variables"<<endl;
  cout<<"pdn.eval = "<<eval(pdn)<<endl;
#ifndef NDEBUG
  cout<<"number of operations evaluated = "<<SymbolicExpression::evalOperationsCount<<endl;
#endif

  VecV vdoublev({3.1,4.2,5.3});
  v&=vdoublev;
  a_&=6.4;
  a5&=7.5;
  a6&=8.6;
#ifndef NDEBUG
  SymbolicExpression::evalOperationsCount=0;
#endif
  cout<<"third eval with different values for independent variables"<<endl;
  cout<<"pdn.eval = "<<eval(pdn)<<endl;
#ifndef NDEBUG
  cout<<"number of operations evaluated = "<<SymbolicExpression::evalOperationsCount<<endl;
#endif

#ifndef NDEBUG
  SymbolicExpression::evalOperationsCount=0;
#endif
  eval(_rmv(0));
#ifndef NDEBUG
  cout<<"number of operations evaluated = "<<SymbolicExpression::evalOperationsCount<<endl;
#endif

  // check all variants of partial derivatives
  {
    SymbolicExpression dep;
    IndependentVariable indep;
    SymbolicExpression ret=parDer(dep, indep);
  }
  {
    SymbolicExpression dep;
    Vector<Fixed<3>, IndependentVariable> indep;
    RowVector<Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }
  {
    Vector<Fixed<3>, SymbolicExpression> dep;
    IndependentVariable indep;
    Vector<Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }
  {
    Vector<Fixed<3>, SymbolicExpression> dep;
    Vector<Fixed<3>, IndependentVariable> indep;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }
  {
    RowVector<Fixed<3>, SymbolicExpression> dep;
    IndependentVariable indep;
    RowVector<Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }
  {
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> dep;
    IndependentVariable indep;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }
  {
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> dep;
    IndependentVariable indep;
    Vector<Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }
  {
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> dep;
    Vector<Fixed<3>, IndependentVariable> indep;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> ret=parDer(dep, indep);
  }

  // dump all ret* vars to this stream just to avoid set but unused warnings
  {
    stringstream dummy;
    double retd=eval(SymbolicExpression()); dummy<<retd;
    Vector<Var, double> retvv=eval(Vector<Var, SymbolicExpression>(3)); dummy<<retvv;
    Vector<Fixed<3>, double> retvf=eval(Vector<Fixed<3>, SymbolicExpression>()); dummy<<retvf;
    RowVector<Var, double> retrvv=eval(RowVector<Var, SymbolicExpression>(3)); dummy<<retrvv;
    RowVector<Fixed<3>, double> retrvf=eval(RowVector<Fixed<3>, SymbolicExpression>()); dummy<<retrvf;
    Matrix<General, Var, Var, double> retmvv=eval(Matrix<General, Var, Var, SymbolicExpression>(3,3)); dummy<<retmvv;
    Matrix<General, Fixed<3>, Var, double> retmfv=eval(Matrix<General, Fixed<3>, Var, SymbolicExpression>(3)); dummy<<retmfv;
    Matrix<General, Var, Fixed<3>, double> retmvf=eval(Matrix<General, Var, Fixed<3>, SymbolicExpression>(3)); dummy<<retmvf;
    Matrix<General, Fixed<3>, Fixed<3>, double> retmff=eval(Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression>()); dummy<<retmff;
  }

  // write to stream for later reread
  {
    { stringstream ostr; ostr<<v; vSer=ostr.str(); }
    { stringstream ostr; ostr<<a_; a_Ser=ostr.str(); }
    { stringstream ostr; ostr<<a5; a5Ser=ostr.str(); }
    { stringstream ostr; ostr<<a6; a6Ser=ostr.str(); }
    { stringstream ostr; ostr<<pdn(0); pdn0Ser=ostr.str(); }
    pdn0Value=eval(pdn(0));
  }
}

void checkSymReread(double pdn0Value,
                    const string &vSer,
                    const string &a_Ser,
                    const string &a5Ser,
                    const string &a6Ser,
                    const string &pdn0Ser) {
  // read expressions from streams
  Vector<Fixed<3>, IndependentVariable> v; { stringstream istr(vSer); istr>>v; }
  IndependentVariable a_; { stringstream istr(a_Ser); istr>>a_; }
  IndependentVariable a5; { stringstream istr(a5Ser); istr>>a5; }
  IndependentVariable a6; { stringstream istr(a6Ser); istr>>a6; }

  SymbolicExpression pdn0; { stringstream istr(pdn0Ser); istr>>pdn0; }

  // set independent variables
  v&=Vec3({3.1,4.2,5.3});
  a_&=6.4;
  a5&=7.5;
  a6&=8.6;

  // print evaluated value after reading the expressions from stream
  cout<<"reread "<<pdn0Value<<" == "<<eval(pdn0)<<endl;

  // more streamout/in tests
  {
    Vector<Fixed<3>, IndependentVariable> v(NONINIT);
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> M;
    M(0,0)=1;
    M(0,1)=2+v(0);
    M(0,2)=3;
    M(1,0)=4+log(v(1));
    M(1,1)=5;
    M(1,2)=6;
    M(2,0)=7;
    M(2,1)=8;
    M(2,2)=9*v(2);
    Vector<Fixed<3>, SymbolicExpression> y=M*v;
    cout<<"indepvar ser "<<v<<endl;
    cout<<"depvar ser "<<M<<endl;
    cout<<"depvar ser "<<y<<endl;
    { stringstream str; str<<v; Vector<Fixed<3>, IndependentVariable> vrr; str>>vrr; cout<<"reread "<<vrr<<endl; }
    { stringstream str; str<<y; Vector<Fixed<3>, SymbolicExpression>  yrr; str>>yrr; cout<<"reread "<<yrr<<endl; }
  }

  // check for bit-identical write/read to/from stream
  {
    SymbolicExpression x=3;
    SymbolicExpression a=1/x;
    SymbolicExpression b=2/x;
    stringstream ostr;
    ostr<<a<<b;
    stringstream istr(ostr.str());
    SymbolicExpression ar;
    SymbolicExpression br;
    istr>>ar>>br;
    cout<<"reread-precision "<<(a==ar)<<" "<<(eval(a)==eval(ar))<<endl;
    cout<<"reread-precision "<<(b==br)<<" "<<(eval(b)==eval(br))<<endl;
  }
}

void checkSymRereadExistingIndeps(double pdn0Value,
                                  const string &vSer,
                                  const string &a_Ser,
                                  const string &a5Ser,
                                  const string &a6Ser,
                                  const string &pdn0Ser) {
  // read expression and independent variables from stream
  Vector<Fixed<3>, IndependentVariable> v(NONINIT);
  IndependentVariable a_;
  IndependentVariable a5;
  IndependentVariable a6;

  Vector<Fixed<3>, IndependentVariable> vFromSer;  { stringstream istr(vSer);  istr>>vFromSer; }
  IndependentVariable a_FromSer; { stringstream istr(a_Ser); istr>>a_FromSer; }
  IndependentVariable a5FromSer; { stringstream istr(a5Ser); istr>>a5FromSer; }
  IndependentVariable a6FromSer; { stringstream istr(a6Ser); istr>>a6FromSer; }

  // Substitutes the independent variables from the stream with already exising ones from c++
  SymbolicExpression pdn0;
  {
    stringstream istr(pdn0Ser);
    istr>>pdn0;
    pdn0=subst(pdn0, vFromSer, v);
    pdn0=subst(subst(subst(pdn0, a_FromSer, a_), a5FromSer, a5), a6FromSer, a6);
  }

  // set independent variables
  v&=Vec3({3.1,4.2,5.3});
  a_&=6.4;
  a5&=7.5;
  a6&=8.6;

  // dump output
  cout<<"reread "<<pdn0Value<<" == "<<eval(pdn0)<<endl;
}

int main() {
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  checkDoubleComplex();
  checkComplex();
  double pdn0Value;
  string v, a_, a5, a6, pdn0Ser;
  checkSym(pdn0Value, v, a_, a5, a6, pdn0Ser);
  SymbolicExpression::garbageCollect();
  checkSymReread(pdn0Value, v, a_, a5, a6, pdn0Ser);
  checkSymRereadExistingIndeps(pdn0Value, v, a_, a5, a6, pdn0Ser);

  return 0;  
}
