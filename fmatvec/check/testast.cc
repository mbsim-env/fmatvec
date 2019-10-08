#include <cfenv>
#include <cassert>
#include <iostream>
#include "fmatvec/ast.h"
#include "fmatvec/fmatvec.h"

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
int dump(const Vector<Var, AT>& v);

template<>
int dump<complex<double>>(const Vector<Var, complex<double>>& v) {
  cout<<"complex 0 = "<<v(0)<<endl;
  cout<<"complex 1 = "<<v(1)<<endl;
  cout<<"complex 2 = "<<v(2)<<endl;

  complex<double> v0ref(10.60,  9.90);
  complex<double> v1ref(23.20, 19.80);
  complex<double> v2ref( 9.95, 14.52);
  if(abs(v0ref-v(0))>1e-12) return 10;
  if(abs(v1ref-v(1))>1e-12) return 20;
  if(abs(v2ref-v(2))>1e-12) return 30;
  return 0;
}

template<>
int dump<Expr>(const Vector<Var, Expr>& v) {
  cout<<"expr mfmf"<<endl;
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
  m(0,1)=2.0;
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

  return dump(_rmv);
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
