#include <cfenv>
#include <cassert>
#include <iostream>
#include "fmatvec/ast.h"
#include "fmatvec/fmatvec.h"

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
  typedef Vector<Var, SymbolicExpression> VecSymExpr;
  typedef Matrix<General, Var, Var, SymbolicExpression> MatSymExpr;

  {
    VecSymExpr v2;
    VecSymExpr v(3);
  }
  VecSymExpr v2;
  VecSymExpr v(3);
  v(0)=SymbolicExpression(AST::Symbol::create());//mfmf must be removed but this does not work for now
  v(1)=SymbolicExpression(AST::Symbol::create());//mfmf must be removed but this does not work for now
  v(2)=SymbolicExpression(AST::Symbol::create());//mfmf must be removed but this does not work for now
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
  MatSymExpr m(3,3);
  m(0,0)=2.0;
  m(0,2)=3.0;
  m(1,0)=4.0;
  m(1,1)=5.0;
  m(1,2)=6.0;
  m(2,0)=7.0;
  m(2,1)=8.0;
  m(2,2)=9.0;
  m(2,2)=SymbolicExpression()+SymbolicExpression();
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

  SymbolicExpression pd00=_rmv(0)->parDer(v(0));//mfmf this should be done by MatSymExpr pd=_rmv.parDer(v) 
  SymbolicExpression pd01=_rmv(0)->parDer(v(1));
  SymbolicExpression pd02=_rmv(0)->parDer(v(2));
  SymbolicExpression pd10=_rmv(1)->parDer(v(0));
  SymbolicExpression pd11=_rmv(1)->parDer(v(1));
  SymbolicExpression pd12=_rmv(1)->parDer(v(2));
  SymbolicExpression pd20=_rmv(2)->parDer(v(0));
  SymbolicExpression pd21=_rmv(2)->parDer(v(1));
  SymbolicExpression pd22=_rmv(2)->parDer(v(2));
  SymbolicExpression pdn0=norm->parDer(v(0));//mfmf this should be done by VecSymExpr pdn=norm.parDer(v)
  SymbolicExpression pdn1=norm->parDer(v(1));
  SymbolicExpression pdn2=norm->parDer(v(2));

  string indepVar0; { stringstream str; v(0)->writeXMLFile(str); indepVar0=str.str(); }
  string indepVar1; { stringstream str; v(1)->writeXMLFile(str); indepVar1=str.str(); }
  string indepVar2; { stringstream str; v(2)->writeXMLFile(str); indepVar2=str.str(); }
  string depVar0; { stringstream str; _rmv(0)->writeXMLFile(str); depVar0=str.str(); }
  string depVar1; { stringstream str; _rmv(1)->writeXMLFile(str); depVar1=str.str(); }
  string depVar2; { stringstream str; _rmv(2)->writeXMLFile(str); depVar2=str.str(); }
  string normStr; { stringstream str; norm->writeXMLFile(str); normStr=str.str(); }
  string depVar0_indepVar0; { stringstream str; pd00->writeXMLFile(str); depVar0_indepVar0=str.str(); }
  string depVar0_indepVar1; { stringstream str; pd01->writeXMLFile(str); depVar0_indepVar1=str.str(); }
  string depVar0_indepVar2; { stringstream str; pd02->writeXMLFile(str); depVar0_indepVar2=str.str(); }
  string depVar1_indepVar0; { stringstream str; pd10->writeXMLFile(str); depVar1_indepVar0=str.str(); }
  string depVar1_indepVar1; { stringstream str; pd11->writeXMLFile(str); depVar1_indepVar1=str.str(); }
  string depVar1_indepVar2; { stringstream str; pd12->writeXMLFile(str); depVar1_indepVar2=str.str(); }
  string depVar2_indepVar0; { stringstream str; pd20->writeXMLFile(str); depVar2_indepVar0=str.str(); }
  string depVar2_indepVar1; { stringstream str; pd21->writeXMLFile(str); depVar2_indepVar1=str.str(); }
  string depVar2_indepVar2; { stringstream str; pd22->writeXMLFile(str); depVar2_indepVar2=str.str(); }
  string norm_indepVar0; { stringstream str; pdn0->writeXMLFile(str); norm_indepVar0=str.str(); }
  string norm_indepVar1; { stringstream str; pdn1->writeXMLFile(str); norm_indepVar1=str.str(); }
  string norm_indepVar2; { stringstream str; pdn2->writeXMLFile(str); norm_indepVar2=str.str(); }

  stringstream str;
  str<<"indepVar 0 = "<<indepVar0<<endl;
  str<<"indepVar 1 = "<<indepVar1<<endl;
  str<<"indepVar 2 = "<<indepVar2<<endl;
  str<<"depVar 0 = "<<depVar0<<endl;
  str<<"depVar 1 = "<<depVar1<<endl;
  str<<"depVar 2 = "<<depVar2<<endl;
  str<<"norm = "<<normStr<<endl;
  str<<"depVar0_indepVar 0 = "<<depVar0_indepVar0<<endl;
  str<<"depVar0_indepVar 1 = "<<depVar0_indepVar1<<endl;
  str<<"depVar0_indepVar 2 = "<<depVar0_indepVar2<<endl;
  str<<"depVar1_indepVar 0 = "<<depVar1_indepVar0<<endl;
  str<<"depVar1_indepVar 1 = "<<depVar1_indepVar1<<endl;
  str<<"depVar1_indepVar 2 = "<<depVar1_indepVar2<<endl;
  str<<"depVar2_indepVar 0 = "<<depVar2_indepVar0<<endl;
  str<<"depVar2_indepVar 1 = "<<depVar2_indepVar1<<endl;
  str<<"depVar2_indepVar 2 = "<<depVar2_indepVar2<<endl;
  str<<"norm_indepVar 0 = "<<norm_indepVar0<<endl;
  str<<"norm_indepVar 1 = "<<norm_indepVar1<<endl;
  str<<"norm_indepVar 2 = "<<norm_indepVar2<<endl;
  cout<<str.str();

  if(indepVar0!="a1") return 100;
  if(indepVar1!="a2") return 101;
  if(indepVar2!="a3") return 102;
  if(depVar0!="(0 (0 (2 2 a1) (2 a4 a2)) (2 3 a3))") return 103;
  if(depVar1!="(0 (0 (2 4 a1) (2 5 a2)) (2 6 a3))") return 104;
  if(depVar2!="(0 (0 (2 (5 (0 (4 a1 3) a1)) a1) (2 8 a2)) (2 (0 a5 a6) a3))") return 105;
  //if(normStr!="...") return 200;
  if(depVar0_indepVar0!="2") return 106;
  if(depVar0_indepVar1!="a4") return 107;
  if(depVar0_indepVar2!="3") return 108;
  if(depVar1_indepVar0!="4") return 109;
  if(depVar1_indepVar1!="5") return 110;
  if(depVar1_indepVar2!="6") return 111;
  if(depVar2_indepVar0!="(0 (2 (3 (0 (2 (4 a1 3) (3 3 a1)) 1) (0 (4 a1 3) a1)) a1) (5 (0 (4 a1 3) a1)))") return 112;
  if(depVar2_indepVar1!="8") return 113;
  if(depVar2_indepVar2!="(0 a5 a6)") return 113;
  //if(norm_indepVar0!="...") return 201
  //if(norm_indepVar1!="...") return 202
  //if(norm_indepVar2!="...") return 203

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
