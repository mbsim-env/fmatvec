#include <cfenv>
#include <cassert>
#include <iostream>
#include "fmatvec/ast.h"
#include "fmatvec/fmatvec.h"

//mfmf remove everything with AST
//mfmf replace auto with a explicit type

using namespace std;
using namespace fmatvec;

int checkComplex() {
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

  cout<<"complex 0 = "<<_rmv(0)<<endl;
  cout<<"complex 1 = "<<_rmv(1)<<endl;
  cout<<"complex 2 = "<<_rmv(2)<<endl;

  complex<double> v0ref( 6.60,  9.90);
  complex<double> v1ref(23.20, 19.80);
  complex<double> v2ref( 9.95, 14.52);
  if(abs(v0ref-_rmv(0))>1e-12) return 10;
  if(abs(v1ref-_rmv(1))>1e-12) return 11;
  if(abs(v2ref-_rmv(2))>1e-12) return 12;

  return 0;
}

int checkSym() {
  typedef Vector<Var, Expr> VecExpr;
  typedef Matrix<General, Var, Var, Expr> MatExpr;

  {
    VecExpr v2;
    VecExpr v(3);
  }
  VecExpr v2;
  VecExpr v(3);
  v(0)=Expr(AST::Var::create());//mfmf must be removed but this does not work for now
  v(1)=Expr(AST::Var::create());//mfmf must be removed but this does not work for now
  v(2)=Expr(AST::Var::create());//mfmf must be removed but this does not work for now
  Expr a(3.0);
  VecExpr r1=v*a;
  VecExpr r2=a*v;
  VecExpr r3=v*3.0;
  VecExpr r4=v*3;
  VecExpr r5=3.0*v;
  VecExpr r6=3*v;
  VecExpr r11=v/a;
  VecExpr r31=v/3.0;
  VecExpr r41=v/3;
  VecExpr r12=v+v;
  VecExpr r22=v-v;
  MatExpr m2;
  MatExpr m(3,3);
  m(0,0)=2.0;
  m(0,2)=3.0;
  m(1,0)=4.0;
  m(1,1)=5.0;
  m(1,2)=6.0;
  m(2,0)=7.0;
  m(2,1)=8.0;
  m(2,2)=9.0;
  m(2,2)=Expr()+Expr();
  MatExpr _r1=m*a;
  MatExpr _r2=a*m;
  MatExpr _r3=m*3.0;
  MatExpr _r4=m*3;
  MatExpr _r5=3.0*m;
  MatExpr _r6=3*m;
  MatExpr _r11=m/a;
  MatExpr _r31=m/3.0;
  MatExpr _r41=m/3;
  MatExpr _r12=m+m;
  MatExpr _r22=m-m;
  MatExpr _rmm=m*m;
  VecExpr _rmv=m*v;

  auto parDer0=_rmv(0)->parDer(dynamic_pointer_cast<const AST::Var>(v(0)));
  auto parDer1=_rmv(0)->parDer(dynamic_pointer_cast<const AST::Var>(v(1)));
  auto parDer2=_rmv(0)->parDer(dynamic_pointer_cast<const AST::Var>(v(2)));

  string indepVar0; { stringstream str; v(0)->writeXMLFile(str); indepVar0=str.str(); }
  string indepVar1; { stringstream str; v(1)->writeXMLFile(str); indepVar1=str.str(); }
  string indepVar2; { stringstream str; v(2)->writeXMLFile(str); indepVar2=str.str(); }
  string depVar0; { stringstream str; _rmv(0)->writeXMLFile(str); depVar0=str.str(); }
  string depVar1; { stringstream str; _rmv(1)->writeXMLFile(str); depVar1=str.str(); }
  string depVar2; { stringstream str; _rmv(2)->writeXMLFile(str); depVar2=str.str(); }
  string depVar_indepVar0; { stringstream str; parDer0->writeXMLFile(str); depVar_indepVar0=str.str(); }
  string depVar_indepVar1; { stringstream str; parDer1->writeXMLFile(str); depVar_indepVar1=str.str(); }
  string depVar_indepVar2; { stringstream str; parDer2->writeXMLFile(str); depVar_indepVar2=str.str(); }

  stringstream str;
  str<<"indepVar 0 = "<<indepVar0<<endl;
  str<<"indepVar 1 = "<<indepVar1<<endl;
  str<<"indepVar 2 = "<<indepVar2<<endl;
  str<<"depVar 0 = "<<depVar0<<endl;
  str<<"depVar 1 = "<<depVar1<<endl;
  str<<"depVar 2 = "<<depVar2<<endl;
  str<<"depVar_indepVar 0 = "<<depVar_indepVar0<<endl;
  str<<"depVar_indepVar 1 = "<<depVar_indepVar1<<endl;
  str<<"depVar_indepVar 2 = "<<depVar_indepVar2<<endl;
  cout<<str.str();

  if(indepVar0!="a1") return 100;
  if(indepVar1!="a2") return 101;
  if(indepVar2!="a3") return 102;
  if(depVar0!="(0 (0 (2 2 a1) (2 a4 a2)) (2 3 a3))") return 103;
  if(depVar1!="(0 (0 (2 4 a1) (2 5 a2)) (2 6 a3))") return 104;
  if(depVar2!="(0 (0 (2 7 a1) (2 8 a2)) (2 (0 a5 a6) a3))") return 105;
  if(depVar_indepVar0!="2") return 106;
  if(depVar_indepVar1!="a4") return 107;
  if(depVar_indepVar2!="3") return 108;

  return 0;
}

int main() {
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  int retc=checkComplex();
  int rete=checkSym();

  if(retc!=0)
    return retc;
  if(rete!=0)
    return rete;

  return 0;  
}
