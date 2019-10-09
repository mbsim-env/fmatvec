#include <cfenv>
#include <cassert>
#include <iostream>
#include "fmatvec/ast.h"
#include "fmatvec/fmatvec.h"

//mfmf remove everything with AST

using namespace std;
using namespace fmatvec;

template<class AT>
AT create();

template<>
complex<double> create<complex<double>>() {
  return complex<double>(2.2, 3.3);
}

template<>
Expr create<Expr>() {
  return Expr()+Expr();
}

template<class AT>
int dump(const Vector<Var, AT>& depVar, const Vector<Var, AT>& indepVar);

template<>
int dump<complex<double>>(const Vector<Var, complex<double>>& depVar, const Vector<Var, complex<double>>& indepVar) {
  cout<<"complex 0 = "<<depVar(0)<<endl;
  cout<<"complex 1 = "<<depVar(1)<<endl;
  cout<<"complex 2 = "<<depVar(2)<<endl;
  complex<double> v0ref( 6.60,  9.90);
  complex<double> v1ref(23.20, 19.80);
  complex<double> v2ref( 9.95, 14.52);
  if(abs(v0ref-depVar(0))>1e-12) return 10;
  if(abs(v1ref-depVar(1))>1e-12) return 11;
  if(abs(v2ref-depVar(2))>1e-12) return 12;
  return 0;
}

template<>
int dump<Expr>(const Vector<Var, Expr>& depVar, const Vector<Var, Expr>& indepVar) {
  auto parDer0=depVar(0)->parDer(dynamic_pointer_cast<const AST::Var>(indepVar(0)));

  string indepVar0; { stringstream str; indepVar(0)->writeXMLFile(str); indepVar0=str.str(); }
  string indepVar1; { stringstream str; indepVar(1)->writeXMLFile(str); indepVar1=str.str(); }
  string indepVar2; { stringstream str; indepVar(2)->writeXMLFile(str); indepVar2=str.str(); }
  string depVar0; { stringstream str; depVar(0)->writeXMLFile(str); depVar0=str.str(); }
  string depVar1; { stringstream str; depVar(1)->writeXMLFile(str); depVar1=str.str(); }
  string depVar2; { stringstream str; depVar(2)->writeXMLFile(str); depVar2=str.str(); }
  string depVar_indepVar0; { stringstream str; parDer0->writeXMLFile(str); depVar_indepVar0=str.str(); }

  stringstream str;
  str<<"indepVar 0 = "<<indepVar0<<endl;
  str<<"indepVar 1 = "<<indepVar1<<endl;
  str<<"indepVar 2 = "<<indepVar2<<endl;
  str<<"depVar 0 = "<<depVar0<<endl;
  str<<"depVar 1 = "<<depVar1<<endl;
  str<<"depVar 2 = "<<depVar2<<endl;
  str<<"depVar_indepVar 0 = "<<depVar_indepVar0<<endl;
  cout<<str.str();

  if(indepVar0!="a1") return 100;
  if(indepVar1!="2") return 101;
  if(indepVar2!="(0 a2 a3)") return 102;
  if(depVar0!="(0 (0 (2 2 a1) (2 a4 2)) (2 3 (0 a2 a3)))") return 103;
  if(depVar1!="(0 (0 (2 4 a1) 10) (2 6 (0 a2 a3)))") return 104;
  if(depVar2!="(0 (0 (2 7 a1) 16) (2 (0 a5 a6) (0 a2 a3)))") return 105;
  if(depVar_indepVar0!="2") return 106;

  return 0;
}

template<class AT>
int check() {
  typedef Vector<Var, AT> VecAT;
  typedef Matrix<General, Var, Var, AT> MatAT;

  VecAT v2;
  VecAT v(3);
  v(1)=2.0;
  v(2)=3.0;
  v(2)=create<AT>();
  AT a(3.0);
  auto r1=v*a;
  auto r2=a*v;
  auto r3=v*3.0;
  auto r4=v*3;
  auto r5=3.0*v;
  auto r6=3*v;
  auto r11=v/a;
  auto r31=v/3.0;
  auto r41=v/3;
  auto r12=v+v;
  auto r22=v-v;
  MatAT m2;
  MatAT m(3,3);
  m(0,0)=2.0;
  m(0,2)=3.0;
  m(1,0)=4.0;
  m(1,1)=5.0;
  m(1,2)=6.0;
  m(2,0)=7.0;
  m(2,1)=8.0;
  m(2,2)=9.0;
  m(2,2)=create<AT>();
  auto _r1=m*a;
  auto _r2=a*m;
  auto _r3=m*3.0;
  auto _r4=m*3;
  auto _r5=3.0*m;
  auto _r6=3*m;
  auto _r11=m/a;
  auto _r31=m/3.0;
  auto _r41=m/3;
  auto _r12=m+m;
  auto _r22=m-m;
  auto _rmm=m*m;
  auto _rmv=m*v;

  return dump(_rmv, v);
}

int main() {
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  int retc=check<complex<double>>();
  int rete=check<Expr>();

  if(retc!=0)
    return retc;
  if(rete!=0)
    return rete;

  return 0;  
}
