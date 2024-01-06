#include <cfenv>
#include <cassert>
#include "fmatvec/symbolic.h"
#include <fmatvec/symbolic_function.h>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace fmatvec;

void dumpS(const double &s) {
  cout<<s<<endl;
}

template<class V>
void dumpV(const V &v) {
  cout<<"[";
  for(int i=0; i<v.size(); ++i)
    cout<<(i==0?"":", ")<<v(i);
  cout<<"]"<<endl;
}

template<class M>
void dumpM(const M &m) {
  cout<<"[";
  for(int r=0; r<m.rows(); ++r)
    for(int c=0; c<m.cols(); ++c)
      cout<<(r==0?"":((c==0?"; ":", ")))<<m(r,c);
  cout<<"]"<<endl;
}

string check(double a, double b) {
  bool e=true;
  double m=max(abs(a), abs(b));
  if(m<1e-10) {
    if(abs(a-b)>1e-12)
      e=false;
  }
  else {
    if(abs(a-b)/m>1e-12)
      e=false;
  }
  return e ? boost::lexical_cast<string>(a)+" equal" : boost::lexical_cast<string>(a)+" NOT_EQUAL (=a"+", b="+boost::lexical_cast<string>(b)+")";
}

int main() {
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  cout<<setprecision(12);

  {
    IndependentVariable x;
    Vector<Fixed<3>, SymbolicExpression> v; v(0)=x*x; v(1)=x*x; v(2)=x*x;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> M;
    M(0,0)=x*x; M(1,0)=x*x; M(2,0)=x*x;
    M(0,1)=x*x; M(1,1)=x*x; M(2,1)=x*x;
    M(0,2)=x*x; M(1,2)=x*x; M(2,2)=x*x;
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> R;
    R(0,0)=x*x; R(1,0)=x*x; R(2,0)=x*x;
    R(0,1)=x*x; R(1,1)=x*x; R(2,1)=x*x;
    R(0,2)=x*x; R(1,2)=x*x; R(2,2)=x*x;
    SymbolicFunction<double(double)> funcS(x, x*x);
    SymbolicFunction<Vec3(double)> funcV(x, v*x);
    SymbolicFunction<RowVec3(double)> funcRV(x, v.T()*x);
    SymbolicFunction<Mat3x3(double)> funcM(x, M*x);
    SymbolicFunction<RotMat3(double)> funcR(x, R*x);

    double arg=1.2;
    double argd=5.6;
    double argd2=7.8;

    dumpS(funcS(arg));
    dumpV(funcV(arg));
    dumpV(funcRV(arg));
    dumpM(funcM(arg));
    dumpM(funcR(arg));

    dumpS(funcS.parDer(arg));
    dumpV(funcV.parDer(arg));
    dumpV(funcRV.parDer(arg));
    dumpM(funcM.parDer(arg));
    dumpM(funcR.parDer(arg));

    dumpS(funcS.dirDer(argd, arg));
    dumpV(funcV.dirDer(argd, arg));
    dumpV(funcRV.dirDer(argd, arg));
    dumpM(funcM.dirDer(argd, arg));
    dumpM(funcR.dirDer(argd, arg));

    dumpS(funcS.parDerParDer(arg));
    dumpV(funcV.parDerParDer(arg));
    dumpV(funcRV.parDerParDer(arg));
    dumpM(funcM.parDerParDer(arg));
    dumpM(funcR.parDerParDer(arg));

    dumpS(funcS.parDerDirDer(argd, arg));
    dumpV(funcV.parDerDirDer(argd, arg));
    dumpV(funcRV.parDerDirDer(argd, arg));
    dumpM(funcM.parDerDirDer(argd, arg));
    dumpM(funcR.parDerDirDer(argd, arg));

    dumpS(funcS.dirDerDirDer(argd, argd2, arg));
    dumpV(funcV.dirDerDirDer(argd, argd2, arg));
    dumpV(funcRV.dirDerDirDer(argd, argd2, arg));
    dumpM(funcM.dirDerDirDer(argd, argd2, arg));
    dumpM(funcR.dirDerDirDer(argd, argd2, arg));
  }

  {
    Vector<Fixed<3>, IndependentVariable> x(NONINIT);
    Vector<Fixed<3>, SymbolicExpression> v; v(0)=x(0)*x(1); v(1)=x(0)*x(1); v(2)=x(0)*x(1);
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> M;
    M(0,0)=x(0)*x(1); M(1,0)=x(0)*x(1); M(2,0)=x(0)*x(1);
    M(0,1)=x(0)*x(1); M(1,1)=x(0)*x(1); M(2,1)=x(0)*x(1);
    M(0,2)=x(0)*x(1); M(1,2)=x(0)*x(1); M(2,2)=x(0)*x(1);
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> R;
    R(0,0)=x(0)*x(1); R(1,0)=x(0)*x(1); R(2,0)=x(0)*x(1);
    R(0,1)=x(0)*x(1); R(1,1)=x(0)*x(1); R(2,1)=x(0)*x(1);
    R(0,2)=x(0)*x(1); R(1,2)=x(0)*x(1); R(2,2)=x(0)*x(1);
    SymbolicFunction<double(Vec3)> funcS(x, x.T()*v);
    SymbolicFunction<Vec3(Vec3)> funcV(x, M*v);
    SymbolicFunction<RowVec3(Vec3)> funcRV(x, (M*v).T());
    SymbolicFunction<Mat3x3(Vec3)> funcM(x, x*v.T());
    SymbolicFunction<RotMat3(Vec3)> funcR(x, R*x(0)*x(1)*x(2));

    Vec3 arg({1.2, 2.3, 3.4});
    Vec3 argd({5.6, 6.7, 7.8});
    Vec3 argd2({7.8, 8.9, 9.1});

    dumpS(funcS(arg));
    dumpV(funcV(arg));
    dumpV(funcRV(arg));
    dumpM(funcM(arg));
    dumpM(funcR(arg));

    dumpV(funcS.parDer(arg));
    dumpM(funcV.parDer(arg));
//    dumpV(funcRV.parDer(arg));
//    dumpM(funcM.parDer(arg));
    dumpM(funcR.parDer(arg));

    dumpS(funcS.dirDer(argd, arg));
    dumpV(funcV.dirDer(argd, arg));
    dumpV(funcRV.dirDer(argd, arg));
    dumpM(funcM.dirDer(argd, arg));
    dumpM(funcR.dirDer(argd, arg));

//    dumpS(funcS.parDerParDer(arg));
//    dumpV(funcV.parDerParDer(arg));
//    dumpV(funcRV.parDerParDer(arg));
//    dumpM(funcM.parDerParDer(arg));
//    dumpM(funcR.parDerParDer(arg));

    dumpV(funcS.parDerDirDer(argd, arg));
    dumpM(funcV.parDerDirDer(argd, arg));
//    dumpV(funcRV.parDerDirDer(argd, arg));
//    dumpM(funcM.parDerDirDer(argd, arg));
    dumpM(funcR.parDerDirDer(argd, arg));

    dumpS(funcS.dirDerDirDer(argd, argd2, arg));
    dumpV(funcV.dirDerDirDer(argd, argd2, arg));
    dumpV(funcRV.dirDerDirDer(argd, argd2, arg));
    dumpM(funcM.dirDerDirDer(argd, argd2, arg));
    dumpM(funcR.dirDerDirDer(argd, argd2, arg));
  }

  {
    IndependentVariable x1;
    IndependentVariable x2;
    SymbolicExpression scalar;
    Vector<Fixed<3>, SymbolicExpression> vector;
    RowVector<Fixed<3>, SymbolicExpression> rowVector;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> matrix;
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> rotMatrix;
    SymbolicFunction<double(double, double)> funcS(x1, x2, scalar);
    SymbolicFunction<Vec3(double, double)> funcV(x1, x2, vector);
    SymbolicFunction<RowVec3(double, double)> funcRV(x1, x2, rowVector);
    SymbolicFunction<Mat3x3(double, double)> funcM(x1, x2, matrix);
    SymbolicFunction<RotMat3(double, double)> funcR(x1, x2, rotMatrix);

    double arg1(1.2);
    double arg2(2.3);
    double arg1Dir(3.4);
    double arg2Dir(4.5);
    double arg1Dir_2(5.6);
    double arg2Dir_2(6.7);

    funcS (arg1, arg2);
    funcV (arg1, arg2);
    funcRV(arg1, arg2);
    funcM (arg1, arg2);
    funcR (arg1, arg2);

    funcS .parDer1(arg1, arg2);
    funcV .parDer1(arg1, arg2);
    funcRV.parDer1(arg1, arg2);
    funcM .parDer1(arg1, arg2);
    funcR .parDer1(arg1, arg2);

    funcS .dirDer1(arg1Dir, arg1, arg2);
    funcV .dirDer1(arg1Dir, arg1, arg2);
    funcRV.dirDer1(arg1Dir, arg1, arg2);
    funcM .dirDer1(arg1Dir, arg1, arg2);
    funcR .dirDer1(arg1Dir, arg1, arg2);

    funcS .parDer2(arg1, arg2);
    funcV .parDer2(arg1, arg2);
    funcRV.parDer2(arg1, arg2);
    funcM .parDer2(arg1, arg2);
    funcR .parDer2(arg1, arg2);

    funcS .dirDer2(arg2Dir, arg1, arg2);
    funcV .dirDer2(arg2Dir, arg1, arg2);
    funcRV.dirDer2(arg2Dir, arg1, arg2);
    funcM .dirDer2(arg2Dir, arg1, arg2);
    funcR .dirDer2(arg2Dir, arg1, arg2);

    funcS .parDer1ParDer1(arg1, arg2);
    funcV .parDer1ParDer1(arg1, arg2);
    funcRV.parDer1ParDer1(arg1, arg2);
    funcM .parDer1ParDer1(arg1, arg2);
    funcR .parDer1ParDer1(arg1, arg2);

    funcS .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcRV.parDer1DirDer1(arg1Dir, arg1, arg2);
    funcM .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer1DirDer1(arg1Dir, arg1, arg2);

    funcS .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcV .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcRV.dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcM .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcR .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);

    funcS .parDer2ParDer2(arg1, arg2);
    funcV .parDer2ParDer2(arg1, arg2);
    funcRV.parDer2ParDer2(arg1, arg2);
    funcM .parDer2ParDer2(arg1, arg2);
    funcR .parDer2ParDer2(arg1, arg2);

    funcS .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcRV.parDer2DirDer2(arg2Dir, arg1, arg2);
    funcM .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer2DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcV .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcRV.dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcM .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcR .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);

    funcS .parDer1ParDer2(arg1, arg2);
    funcV .parDer1ParDer2(arg1, arg2);
    funcRV.parDer1ParDer2(arg1, arg2);
    funcM .parDer1ParDer2(arg1, arg2);
    funcR .parDer1ParDer2(arg1, arg2);

    funcS .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcRV.parDer1DirDer2(arg2Dir, arg1, arg2);
    funcM .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer1DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcV .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcRV.dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcM .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcR .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);

    funcS .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcRV.parDer2DirDer1(arg1Dir, arg1, arg2);
    funcM .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer2DirDer1(arg1Dir, arg1, arg2);
  }

  {
    Vector<Fixed<3>, IndependentVariable> x1(NONINIT);
    Vector<Fixed<3>, IndependentVariable> x2(NONINIT);
    SymbolicExpression scalar;
    Vector<Fixed<3>, SymbolicExpression> vector;
    RowVector<Fixed<3>, SymbolicExpression> rowVector;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> matrix;
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> rotMatrix;
    SymbolicFunction<double(Vec3, Vec3)> funcS(x1, x2, scalar);
    SymbolicFunction<Vec3(Vec3, Vec3)> funcV(x1, x2, vector);
    SymbolicFunction<RowVec3(Vec3, Vec3)> funcRV(x1, x2, rowVector);
    SymbolicFunction<Mat3x3(Vec3, Vec3)> funcM(x1, x2, matrix);
    SymbolicFunction<RotMat3(Vec3, Vec3)> funcR(x1, x2, rotMatrix);

    Vec3 arg1({1.2,1.2,1.2});
    Vec3 arg2({2.3,2.3,2.3});
    Vec3 arg1Dir({3.4,3.4,3.4});
    Vec3 arg2Dir({4.5,4.5,4.5});
    Vec3 arg1Dir_2({5.6,5.6,5.6});
    Vec3 arg2Dir_2({6.7,6.7,6.7});

    funcS (arg1, arg2);
    funcV (arg1, arg2);
    funcRV(arg1, arg2);
    funcM (arg1, arg2);
    funcR (arg1, arg2);

    funcS .parDer1(arg1, arg2);
    funcV .parDer1(arg1, arg2);
    //funcRV.parDer1(arg1, arg2);
    //funcM .parDer1(arg1, arg2);
    funcR .parDer1(arg1, arg2);

    funcS .dirDer1(arg1Dir, arg1, arg2);
    funcV .dirDer1(arg1Dir, arg1, arg2);
    funcRV.dirDer1(arg1Dir, arg1, arg2);
    funcM .dirDer1(arg1Dir, arg1, arg2);
    funcR .dirDer1(arg1Dir, arg1, arg2);

    funcS .parDer2(arg1, arg2);
    funcV .parDer2(arg1, arg2);
    //funcRV.parDer2(arg1, arg2);
    //funcM .parDer2(arg1, arg2);
    funcR .parDer2(arg1, arg2);

    funcS .dirDer2(arg2Dir, arg1, arg2);
    funcV .dirDer2(arg2Dir, arg1, arg2);
    funcRV.dirDer2(arg2Dir, arg1, arg2);
    funcM .dirDer2(arg2Dir, arg1, arg2);
    funcR .dirDer2(arg2Dir, arg1, arg2);

    //funcS .parDer1ParDer1(arg1, arg2);
    //funcV .parDer1ParDer1(arg1, arg2);
    //funcRV.parDer1ParDer1(arg1, arg2);
    //funcM .parDer1ParDer1(arg1, arg2);
    //funcR .parDer1ParDer1(arg1, arg2);

    funcS .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer1DirDer1(arg1Dir, arg1, arg2);
    //funcRV.parDer1DirDer1(arg1Dir, arg1, arg2);
    //funcM .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer1DirDer1(arg1Dir, arg1, arg2);

    funcS .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcV .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcRV.dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcM .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcR .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);

    //funcS .parDer2ParDer2(arg1, arg2);
    //funcV .parDer2ParDer2(arg1, arg2);
    //funcRV.parDer2ParDer2(arg1, arg2);
    //funcM .parDer2ParDer2(arg1, arg2);
    //funcR .parDer2ParDer2(arg1, arg2);

    funcS .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer2DirDer2(arg2Dir, arg1, arg2);
    //funcRV.parDer2DirDer2(arg2Dir, arg1, arg2);
    //funcM .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer2DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcV .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcRV.dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcM .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcR .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);

    //funcS .parDer1ParDer2(arg1, arg2);
    //funcV .parDer1ParDer2(arg1, arg2);
    //funcRV.parDer1ParDer2(arg1, arg2);
    //funcM .parDer1ParDer2(arg1, arg2);
    //funcR .parDer1ParDer2(arg1, arg2);

    funcS .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer1DirDer2(arg2Dir, arg1, arg2);
    //funcRV.parDer1DirDer2(arg2Dir, arg1, arg2);
    //funcM .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer1DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcV .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcRV.dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcM .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcR .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);

    funcS .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer2DirDer1(arg1Dir, arg1, arg2);
    //funcRV.parDer2DirDer1(arg1Dir, arg1, arg2);
    //funcM .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer2DirDer1(arg1Dir, arg1, arg2);
  }

  {
    IndependentVariable x1;
    Vector<Fixed<3>, IndependentVariable> x2(NONINIT);
    SymbolicExpression scalar;
    Vector<Fixed<3>, SymbolicExpression> vector;
    RowVector<Fixed<3>, SymbolicExpression> rowVector;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> matrix;
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> rotMatrix;
    SymbolicFunction<double(double, Vec3)> funcS(x1, x2, scalar);
    SymbolicFunction<Vec3(double, Vec3)> funcV(x1, x2, vector);
    SymbolicFunction<RowVec3(double, Vec3)> funcRV(x1, x2, rowVector);
    SymbolicFunction<Mat3x3(double, Vec3)> funcM(x1, x2, matrix);
    SymbolicFunction<RotMat3(double, Vec3)> funcR(x1, x2, rotMatrix);

    double arg1(1.2);
    Vec3 arg2({2.3,2.3,2.3});
    double arg1Dir(3.4);
    Vec3 arg2Dir({4.5,4.5,4.5});
    double arg1Dir_2(5.6);
    Vec3 arg2Dir_2({6.7,6.7,6.7});

    funcS (arg1, arg2);
    funcV (arg1, arg2);
    funcRV(arg1, arg2);
    funcM (arg1, arg2);
    funcR (arg1, arg2);

    funcS .parDer1(arg1, arg2);
    funcV .parDer1(arg1, arg2);
    funcRV.parDer1(arg1, arg2);
    funcM .parDer1(arg1, arg2);
    funcR .parDer1(arg1, arg2);

    funcS .dirDer1(arg1Dir, arg1, arg2);
    funcV .dirDer1(arg1Dir, arg1, arg2);
    funcRV.dirDer1(arg1Dir, arg1, arg2);
    funcM .dirDer1(arg1Dir, arg1, arg2);
    funcR .dirDer1(arg1Dir, arg1, arg2);

    funcS .parDer2(arg1, arg2);
    funcV .parDer2(arg1, arg2);
    //funcRV.parDer2(arg1, arg2);
    //funcM .parDer2(arg1, arg2);
    funcR .parDer2(arg1, arg2);

    funcS .dirDer2(arg2Dir, arg1, arg2);
    funcV .dirDer2(arg2Dir, arg1, arg2);
    funcRV.dirDer2(arg2Dir, arg1, arg2);
    funcM .dirDer2(arg2Dir, arg1, arg2);
    funcR .dirDer2(arg2Dir, arg1, arg2);

    funcS .parDer1ParDer1(arg1, arg2);
    funcV .parDer1ParDer1(arg1, arg2);
    funcRV.parDer1ParDer1(arg1, arg2);
    funcM .parDer1ParDer1(arg1, arg2);
    funcR .parDer1ParDer1(arg1, arg2);

    funcS .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcRV.parDer1DirDer1(arg1Dir, arg1, arg2);
    funcM .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer1DirDer1(arg1Dir, arg1, arg2);

    funcS .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcV .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcRV.dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcM .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcR .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);

    //funcS .parDer2ParDer2(arg1, arg2);
    //funcV .parDer2ParDer2(arg1, arg2);
    //funcRV.parDer2ParDer2(arg1, arg2);
    //funcM .parDer2ParDer2(arg1, arg2);
    //funcR .parDer2ParDer2(arg1, arg2);

    funcS .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer2DirDer2(arg2Dir, arg1, arg2);
    //funcRV.parDer2DirDer2(arg2Dir, arg1, arg2);
    //funcM .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer2DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcV .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcRV.dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcM .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcR .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);

    funcS .parDer1ParDer2(arg1, arg2);
    funcV .parDer1ParDer2(arg1, arg2);
    //funcRV.parDer1ParDer2(arg1, arg2);
    //funcM .parDer1ParDer2(arg1, arg2);
    funcR .parDer1ParDer2(arg1, arg2);

    funcS .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcRV.parDer1DirDer2(arg2Dir, arg1, arg2);
    funcM .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer1DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcV .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcRV.dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcM .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcR .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);

    funcS .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer2DirDer1(arg1Dir, arg1, arg2);
    //funcRV.parDer2DirDer1(arg1Dir, arg1, arg2);
    //funcM .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer2DirDer1(arg1Dir, arg1, arg2);
  }

  {
    Vector<Fixed<3>, IndependentVariable> x1(NONINIT);
    IndependentVariable x2;
    SymbolicExpression scalar;
    Vector<Fixed<3>, SymbolicExpression> vector;
    RowVector<Fixed<3>, SymbolicExpression> rowVector;
    Matrix<General, Fixed<3>, Fixed<3>, SymbolicExpression> matrix;
    Matrix<Rotation, Fixed<3>, Fixed<3>, SymbolicExpression> rotMatrix;
    SymbolicFunction<double(Vec3, double)> funcS(x1, x2, scalar);
    SymbolicFunction<Vec3(Vec3, double)> funcV(x1, x2, vector);
    SymbolicFunction<RowVec3(Vec3, double)> funcRV(x1, x2, rowVector);
    SymbolicFunction<Mat3x3(Vec3, double)> funcM(x1, x2, matrix);
    SymbolicFunction<RotMat3(Vec3, double)> funcR(x1, x2, rotMatrix);

    Vec3 arg1({1.2,1.2,1.2});
    double arg2(2.3);
    Vec3 arg1Dir({3.4,3.4,3.4});
    double arg2Dir(4.5);
    Vec3 arg1Dir_2({5.6,5.6,5.6});
    double arg2Dir_2(6.7);

    funcS (arg1, arg2);
    funcV (arg1, arg2);
    funcRV(arg1, arg2);
    funcM (arg1, arg2);
    funcR (arg1, arg2);

    funcS .parDer1(arg1, arg2);
    funcV .parDer1(arg1, arg2);
    //funcRV.parDer1(arg1, arg2);
    //funcM .parDer1(arg1, arg2);
    funcR .parDer1(arg1, arg2);

    funcS .dirDer1(arg1Dir, arg1, arg2);
    funcV .dirDer1(arg1Dir, arg1, arg2);
    funcRV.dirDer1(arg1Dir, arg1, arg2);
    funcM .dirDer1(arg1Dir, arg1, arg2);
    funcR .dirDer1(arg1Dir, arg1, arg2);

    funcS .parDer2(arg1, arg2);
    funcV .parDer2(arg1, arg2);
    funcRV.parDer2(arg1, arg2);
    funcM .parDer2(arg1, arg2);
    funcR .parDer2(arg1, arg2);

    funcS .dirDer2(arg2Dir, arg1, arg2);
    funcV .dirDer2(arg2Dir, arg1, arg2);
    funcRV.dirDer2(arg2Dir, arg1, arg2);
    funcM .dirDer2(arg2Dir, arg1, arg2);
    funcR .dirDer2(arg2Dir, arg1, arg2);

    //funcS .parDer1ParDer1(arg1, arg2);
    //funcV .parDer1ParDer1(arg1, arg2);
    //funcRV.parDer1ParDer1(arg1, arg2);
    //funcM .parDer1ParDer1(arg1, arg2);
    //funcR .parDer1ParDer1(arg1, arg2);

    funcS .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer1DirDer1(arg1Dir, arg1, arg2);
    //funcRV.parDer1DirDer1(arg1Dir, arg1, arg2);
    //funcM .parDer1DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer1DirDer1(arg1Dir, arg1, arg2);

    funcS .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcV .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcRV.dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcM .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);
    funcR .dirDer1DirDer1(arg1Dir, arg1Dir_2, arg1, arg2);

    funcS .parDer2ParDer2(arg1, arg2);
    funcV .parDer2ParDer2(arg1, arg2);
    funcRV.parDer2ParDer2(arg1, arg2);
    funcM .parDer2ParDer2(arg1, arg2);
    funcR .parDer2ParDer2(arg1, arg2);

    funcS .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcRV.parDer2DirDer2(arg2Dir, arg1, arg2);
    funcM .parDer2DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer2DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcV .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcRV.dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcM .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);
    funcR .dirDer2DirDer2(arg2Dir, arg2Dir_2, arg1, arg2);

    funcS .parDer1ParDer2(arg1, arg2);
    funcV .parDer1ParDer2(arg1, arg2);
    //funcRV.parDer1ParDer2(arg1, arg2);
    //funcM .parDer1ParDer2(arg1, arg2);
    funcR .parDer1ParDer2(arg1, arg2);

    funcS .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcV .parDer1DirDer2(arg2Dir, arg1, arg2);
    //funcRV.parDer1DirDer2(arg2Dir, arg1, arg2);
    //funcM .parDer1DirDer2(arg2Dir, arg1, arg2);
    funcR .parDer1DirDer2(arg2Dir, arg1, arg2);

    funcS .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcV .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcRV.dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcM .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);
    funcR .dirDer2DirDer1(arg2Dir, arg1Dir, arg1, arg2);

    funcS .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcV .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcRV.parDer2DirDer1(arg1Dir, arg1, arg2);
    funcM .parDer2DirDer1(arg1Dir, arg1, arg2);
    funcR .parDer2DirDer1(arg1Dir, arg1, arg2);
  }

  // native function with one argument as symbolic function
  {
    IndependentVariable i1, i2;
    auto arg = sin(pow(i1,2))*sinh(pow(i2,3));
    class Func : public Function<double(double)> {
      public:
        double operator()(const double &arg) override {
          return arg*arg;
        }
        double parDer(const double &arg) override {
          return 2*arg;
        }
        // double dirDer(...) override -> use default implemention from fmatvec/function.h
        double parDerParDer(const double &arg) override {
          return 2;
        }
        // double dirDerDirDer(...) override -> use default implemention from fmatvec/function.h
    };
    auto funcN = symbolicFunc<double(double)>(make_shared<Func>(), arg);
    auto funcS = arg*arg;
    auto resN = sin(arg) * pow(funcN,2);
    auto resS = sin(arg) * pow(funcS,2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  // native function with two arguments as symbolic function
  {
    IndependentVariable i1, i2;
    auto arg1 = sin(pow(i1,2))*sinh(pow(i2,3));
    auto arg2 = cos(pow(i2,2))*cosh(pow(i1,3));
    class Func : public Function<double(double,double)> {
      public:
        double operator()(const double &arg1, const double &arg2) override {
          return arg1*arg1+arg2*arg2+arg1*arg1*arg2*arg2;
        }
        double parDer1(const double &arg1, const double &arg2) override {
          return 2*arg1+2*arg1*arg2*arg2;
        }
        double parDer2(const double &arg1, const double &arg2) override {
          return 2*arg2+arg1*arg1*2*arg2;
        }
        double parDer1ParDer2(const double &arg1, const double &arg2) override {
          return 2*arg1*2*arg2;
        }
        double parDer1ParDer1(const double &arg1, const double &arg2) override {
          return 2+2*arg2*arg2;
        }
        double parDer2ParDer2(const double &arg1, const double &arg2) override {
          return 2+arg1*arg1*2;
        }
    };
    auto funcN = symbolicFunc<double(double,double)>(make_shared<Func>(), arg1, arg2);
    auto funcS = arg1*arg1+arg2*arg2+arg1*arg1*arg2*arg2;
    auto resN = cos(arg1) *sin(arg2) * pow(funcN,2);
    auto resS = cos(arg1) *sin(arg2) * pow(funcS,2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  // native function with one vector argument as symbolic function
  {
    IndependentVariable i1, i2;
    Vector<Var, SymbolicExpression> arg({ sin(pow(i1,2))*sinh(pow(i2,3)), cos(pow(i2,2))*cosh(pow(i1,3)) });
    class Func : public Function<double(VecV)> {
      public:
        double operator()(const VecV &arg) override {
          return arg(0)*arg(0)+arg(1)*arg(1)+arg(0)*arg(0)*arg(1)*arg(1);
        }
        RowVecV parDer(const VecV &arg) override {
          return RowVecV({ 2*arg(0)+2*arg(0)*arg(1)*arg(1), 2*arg(1)+arg(0)*arg(0)*2*arg(1) });
        }
        RowVecV parDerDirDer(const VecV &dir, const VecV &arg) override {
          return RowVecV({ (2+2*arg(1)*arg(1))*dir(0) + (2*arg(0)*2*arg(1))*dir(1),
                           (2*arg(0)*2*arg(1))*dir(0) + (2+arg(0)*arg(0)*2)*dir(1) });
        }
    };
    auto funcN = symbolicFunc<double(VecV)>(make_shared<Func>(), arg);
    auto funcS = arg(0)*arg(0)+arg(1)*arg(1)+arg(0)*arg(0)*arg(1)*arg(1);
    auto resN = cos(arg(0)) *sin(arg(1)) * pow(funcN,2);
    auto resS = cos(arg(0)) *sin(arg(1)) * pow(funcS,2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  // native function with two vector arguments as symbolic function
  {
    IndependentVariable i1, i2;
    auto arg1 = sin(pow(i1,2))*sinh(pow(i2,3));
    auto arg2 = Vector<Var,SymbolicExpression>({cos(pow(i2,2))*cosh(pow(i1,3)), log(pow(i2,2))*cosh(pow(i1,5)), cos(2*pow(i2,2))*sqrt(pow(i1,3))});
    class Func : public Function<double(double,VecV)> {
      public:
        double operator()(const double &arg1, const VecV &arg2) override {
          return arg1*arg1+arg2(0)*arg2(1)*arg2(2)+arg1*arg1*arg2(0)*arg2(1)*arg2(2);
        }
        double dirDer1(const double &arg1Dir, const double &arg1, const VecV &arg2) override {
          return (2*arg1+2*arg1*arg2(0)*arg2(1)*arg2(2))*arg1Dir;
        }
        double dirDer2(const VecV &arg2Dir, const double &arg1, const VecV &arg2) override {
          return (arg2(1)*arg2(2)+arg1*arg1*arg2(1)*arg2(2))*arg2Dir(0) +
                 (arg2(0)*arg2(2)+arg1*arg1*arg2(0)*arg2(2))*arg2Dir(1) +
                 (arg2(0)*arg2(1)+arg1*arg1*arg2(0)*arg2(1))*arg2Dir(2);
        }
        double dirDer2DirDer1(const VecV &arg2Dir, const double &arg1Dir, const double &arg1, const VecV &arg2) override {
          return (2*arg1*arg2(1)*arg2(2))*arg2Dir(0)*arg1Dir +
                 (2*arg1*arg2(0)*arg2(2))*arg2Dir(1)*arg1Dir +
                 (2*arg1*arg2(0)*arg2(1))*arg2Dir(2)*arg1Dir;
        }
        double dirDer1DirDer1(const double &arg1Dir_1, const double &arg1Dir_2, const double &arg1, const VecV &arg2) override {
          return (2+2*arg2(0)*arg2(1)*arg2(2))*arg1Dir_1*arg1Dir_2;
        }
        double dirDer2DirDer2(const VecV &arg2Dir_1, const VecV &arg2Dir_2, const double &arg1, const VecV &arg2) override {
          return (arg2(2)+arg1*arg1*arg2(2))*arg2Dir_1(1)*arg2Dir_2(0) +
                 (arg2(1)+arg1*arg1*arg2(1))*arg2Dir_1(2)*arg2Dir_2(0) +
                 (arg2(2)+arg1*arg1*arg2(2))*arg2Dir_1(0)*arg2Dir_2(1) +
                 (arg2(0)+arg1*arg1*arg2(0))*arg2Dir_1(2)*arg2Dir_2(1) +
                 (arg2(1)+arg1*arg1*arg2(1))*arg2Dir_1(0)*arg2Dir_2(2) +
                 (arg2(0)+arg1*arg1*arg2(0))*arg2Dir_1(1)*arg2Dir_2(2);
        }
    };
    auto funcN = symbolicFunc<double(double,VecV)>(make_shared<Func>(), arg1, arg2);
    auto funcS = arg1*arg1+arg2(0)*arg2(1)*arg2(2)+arg1*arg1*arg2(0)*arg2(1)*arg2(2);
    auto resN = cos(arg1+arg1) *sin(arg2(0)+arg2(1)+arg2(2)) * pow(funcN,2);
    auto resS = cos(arg1+arg1) *sin(arg2(0)+arg2(1)+arg2(2)) * pow(funcS,2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  // native function with two vector arguments as symbolic function
  {
    IndependentVariable i1, i2;
    auto arg1 = Vector<Var,SymbolicExpression>({cos(pow(i2,2))*cosh(pow(i1,3)), log(pow(i2,2))*cosh(pow(i1,5)), cos(2*pow(i2,2))*sqrt(pow(i1,3))});
    auto arg2 = sin(pow(i1,2))*sinh(pow(i2,3));
    class Func : public Function<double(VecV,double)> {
      public:
        double operator()(const VecV &arg1, const double &arg2) override {
          return arg2*arg2+arg1(0)*arg1(1)*arg1(2)+arg2*arg2*arg1(0)*arg1(1)*arg1(2);
        }
        double dirDer1(const VecV &arg1Dir, const VecV &arg1, const double &arg2) {
          return (arg1(1)*arg1(2)+arg2*arg2*arg1(1)*arg1(2))*arg1Dir(0) +
                 (arg1(0)*arg1(2)+arg2*arg2*arg1(0)*arg1(2))*arg1Dir(1) +
                 (arg1(0)*arg1(1)+arg2*arg2*arg1(0)*arg1(1))*arg1Dir(2);
        }
        double dirDer2(const double &arg2Dir, const VecV &arg1, const double &arg2) {
          return (2*arg2+2*arg2*arg1(0)*arg1(1)*arg1(2))*arg2Dir;
        }
        double dirDer2DirDer1(const double &arg2Dir, const VecV &arg1Dir, const VecV &arg1, const double &arg2) {
          return (2*arg2*arg1(1)*arg1(2))*arg1Dir(0)*arg2Dir +
                 (2*arg2*arg1(0)*arg1(2))*arg1Dir(1)*arg2Dir +
                 (2*arg2*arg1(0)*arg1(1))*arg1Dir(2)*arg2Dir;
        }
        double dirDer1DirDer1(const VecV &arg1Dir_1, const VecV &arg1Dir_2, const VecV &arg1, const double &arg2) {
          return (arg1(2)+arg2*arg2*arg1(2))*arg1Dir_1(1)*arg1Dir_2(0) +
                 (arg1(1)+arg2*arg2*arg1(1))*arg1Dir_1(2)*arg1Dir_2(0) +
                 (arg1(2)+arg2*arg2*arg1(2))*arg1Dir_1(0)*arg1Dir_2(1) +
                 (arg1(0)+arg2*arg2*arg1(0))*arg1Dir_1(2)*arg1Dir_2(1) +
                 (arg1(1)+arg2*arg2*arg1(1))*arg1Dir_1(0)*arg1Dir_2(2) +
                 (arg1(0)+arg2*arg2*arg1(0))*arg1Dir_1(1)*arg1Dir_2(2);
        }
        double dirDer2DirDer2(const double &arg2Dir_1, const double &arg2Dir_2, const VecV &arg1, const double &arg2) {
          return (2+2*arg1(0)*arg1(1)*arg1(2))*arg2Dir_1*arg2Dir_2;
        }
    };
    auto funcN = symbolicFunc<double(VecV,double)>(make_shared<Func>(), arg1, arg2);
    auto funcS = arg2*arg2+arg1(0)*arg1(1)*arg1(2)+arg2*arg2*arg1(0)*arg1(1)*arg1(2);
    auto resN = cos(arg2+arg2) *sin(arg1(0)+arg1(1)+arg1(2)) * pow(funcN,2);
    auto resS = cos(arg2+arg2) *sin(arg1(0)+arg1(1)+arg1(2)) * pow(funcS,2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  // native function with two vector arguments as symbolic function
  {
    IndependentVariable i1, i2;
    auto arg1 = Vector<Var,SymbolicExpression>({sin(pow(i1,2))*sinh(pow(i2,3)), tan(pow(i1,2))*tanh(pow(i2,3))});
    auto arg2 = Vector<Var,SymbolicExpression>({cos(pow(i2,2))*cosh(pow(i1,3)), log(pow(i2,2))*cosh(pow(i1,5)), cos(2*pow(i2,2))*sqrt(pow(i1,3))});
    class Func : public Function<double(VecV,VecV)> {
      public:
        double operator()(const VecV &arg1, const VecV &arg2) override {
          return arg1(0)*arg1(1)+arg2(0)*arg2(1)*arg2(2)+arg1(0)*arg1(1)*arg2(0)*arg2(1)*arg2(2);
        }
        double dirDer1(const VecV &arg1Dir, const VecV &arg1, const VecV &arg2) {
          return (arg1(1)+arg1(1)*arg2(0)*arg2(1)*arg2(2))*arg1Dir(0) +
                 (arg1(0)+arg1(0)*arg2(0)*arg2(1)*arg2(2))*arg1Dir(1);
        }
        double dirDer2(const VecV &arg2Dir, const VecV &arg1, const VecV &arg2) {
          return (arg2(1)*arg2(2)+arg1(0)*arg1(1)*arg2(1)*arg2(2))*arg2Dir(0) +
                 (arg2(0)*arg2(2)+arg1(0)*arg1(1)*arg2(0)*arg2(2))*arg2Dir(1) +
                 (arg2(0)*arg2(1)+arg1(0)*arg1(1)*arg2(0)*arg2(1))*arg2Dir(2);
        }
        double dirDer2DirDer1(const VecV &arg2Dir, const VecV &arg1Dir, const VecV &arg1, const VecV &arg2) {
          return (arg1(1)*arg2(1)*arg2(2))*arg2Dir(0)*arg1Dir(0) +
                 (arg1(1)*arg2(0)*arg2(2))*arg2Dir(1)*arg1Dir(0) +
                 (arg1(1)*arg2(0)*arg2(1))*arg2Dir(2)*arg1Dir(0) +
                 (arg1(0)*arg2(1)*arg2(2))*arg2Dir(0)*arg1Dir(1) +
                 (arg1(0)*arg2(0)*arg2(2))*arg2Dir(1)*arg1Dir(1) +
                 (arg1(0)*arg2(0)*arg2(1))*arg2Dir(2)*arg1Dir(1);
        }
        double dirDer1DirDer1(const VecV &arg1Dir_1, const VecV &arg1Dir_2, const VecV &arg1, const VecV &arg2) {
          return (1+arg2(0)*arg2(1)*arg2(2))*arg1Dir_1(1)*arg1Dir_2(0) +
                 (1+arg2(0)*arg2(1)*arg2(2))*arg1Dir_1(0)*arg1Dir_2(1);
        }
        double dirDer2DirDer2(const VecV &arg2Dir_1, const VecV &arg2Dir_2, const VecV &arg1, const VecV &arg2) {
          return (arg2(2)+arg1(0)*arg1(1)*arg2(2))*arg2Dir_1(1)*arg2Dir_2(0) +
                 (arg2(1)+arg1(0)*arg1(1)*arg2(1))*arg2Dir_1(2)*arg2Dir_2(0) +
                 (arg2(2)+arg1(0)*arg1(1)*arg2(2))*arg2Dir_1(0)*arg2Dir_2(1) +
                 (arg2(0)+arg1(0)*arg1(1)*arg2(0))*arg2Dir_1(2)*arg2Dir_2(1) +
                 (arg2(1)+arg1(0)*arg1(1)*arg2(1))*arg2Dir_1(0)*arg2Dir_2(2) +
                 (arg2(0)+arg1(0)*arg1(1)*arg2(0))*arg2Dir_1(1)*arg2Dir_2(2);
        }
    };
    auto funcN = symbolicFunc<double(VecV,VecV)>(make_shared<Func>(), arg1, arg2);
    auto funcS = arg1(0)*arg1(1)+arg2(0)*arg2(1)*arg2(2)+arg1(0)*arg1(1)*arg2(0)*arg2(1)*arg2(2);
    auto resN = cos(arg1(0)+arg1(1)) *sin(arg2(0)+arg2(1)+arg2(2)) * pow(funcN,2);
    auto resS = cos(arg1(0)+arg1(1)) *sin(arg2(0)+arg2(1)+arg2(2)) * pow(funcS,2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  {
    IndependentVariable i1, i2;
    auto arg = sin(pow(i1,2))*sinh(pow(i2,3));
    class Func : public Function<VecV(double)> {
      public:
        VecV operator()(const double &arg) override {
          auto x=arg*arg;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        VecV parDer(const double &arg) override {
          auto x=2*arg;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        // double dirDer(...) override -> use default implemention from fmatvec/function.h
        VecV parDerParDer(const double &arg) override {
          auto x=2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        // double dirDerDirDer(...) override -> use default implemention from fmatvec/function.h
    };
    auto funcN = symbolicFunc<VecV(double)>(make_shared<Func>(), arg);
    Vector<Var, SymbolicExpression> funcS(2);
    funcS(0) = arg*arg;
    funcS(1) = 2*funcS(0);
    auto resN = sin(arg) * pow(funcN(0),2)*pow(funcN(1),2);
    auto resS = sin(arg) * pow(funcS(0),2)*pow(funcS(1),2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  {
    IndependentVariable i1, i2;
    auto arg = sin(pow(i1,2))*sinh(pow(i2,3));
    class Func : public Function<MatV(double)> {
      public:
        MatV operator()(const double &arg) override {
          auto x=arg*arg;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        MatV parDer(const double &arg) override {
          auto x=2*arg;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        // double dirDer(...) override -> use default implemention from fmatvec/function.h
        MatV parDerParDer(const double &arg) override {
          auto x=2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        // double dirDerDirDer(...) override -> use default implemention from fmatvec/function.h
    };
    auto funcN = symbolicFunc<MatV(double)>(make_shared<Func>(), arg);
    Matrix<General, Var, Var, SymbolicExpression> funcS(2,2);
    funcS(0,0)=arg*arg; funcS(0,1)=2*funcS(0,0);
    funcS(1,0)=3*arg*arg; funcS(1,1)=4*funcS(0,0);
    auto resN = sin(arg) * pow(funcN(0,0),2)*pow(funcN(0,1),2)*pow(funcN(0,1),2)*pow(funcN(1,1),2);
    auto resS = sin(arg) * pow(funcS(0,0),2)*pow(funcS(0,1),2)*pow(funcS(0,1),2)*pow(funcS(1,1),2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  {
    IndependentVariable i1, i2;
    auto arg1 = sin(pow(i1,2))*sinh(pow(i2,3));
    auto arg2 = cos(pow(i2,2))*cosh(pow(i1,3));
    class Func : public Function<VecV(double,double)> {
      public:
        VecV operator()(const double &arg1, const double &arg2) override {
          auto x=arg1*arg1+arg2*arg2+arg1*arg1*arg2*arg2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        VecV parDer1(const double &arg1, const double &arg2) override {
          auto x=2*arg1+2*arg1*arg2*arg2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        VecV parDer2(const double &arg1, const double &arg2) override {
          auto x=2*arg2+arg1*arg1*2*arg2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        VecV parDer1ParDer2(const double &arg1, const double &arg2) override {
          auto x=2*arg1*2*arg2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        VecV parDer1ParDer1(const double &arg1, const double &arg2) override {
          auto x=2+2*arg2*arg2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
        VecV parDer2ParDer2(const double &arg1, const double &arg2) override {
          auto x=2+arg1*arg1*2;
          VecV r(2);
          r(0)=x;
          r(1)=2*x;
          return r;
        }
    };
    auto funcN = symbolicFunc<VecV(double,double)>(make_shared<Func>(), arg1, arg2);
    Vector<Var, SymbolicExpression> funcS(2);
    funcS(0)= arg1*arg1+arg2*arg2+arg1*arg1*arg2*arg2;
    funcS(1)= 2*funcS(0);
    auto resN = cos(arg1) *sin(arg2) * pow(funcN(0),2)*pow(funcN(1),2);
    auto resS = cos(arg1) *sin(arg2) * pow(funcS(0),2)*pow(funcS(1),2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  {
    IndependentVariable i1, i2;
    auto arg1 = sin(pow(i1,2))*sinh(pow(i2,3));
    auto arg2 = cos(pow(i2,2))*cosh(pow(i1,3));
    class Func : public Function<MatV(double,double)> {
      public:
        MatV operator()(const double &arg1, const double &arg2) override {
          auto x=arg1*arg1+arg2*arg2+arg1*arg1*arg2*arg2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        MatV parDer1(const double &arg1, const double &arg2) override {
          auto x=2*arg1+2*arg1*arg2*arg2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        MatV parDer2(const double &arg1, const double &arg2) override {
          auto x=2*arg2+arg1*arg1*2*arg2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        MatV parDer1ParDer2(const double &arg1, const double &arg2) override {
          auto x=2*arg1*2*arg2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        MatV parDer1ParDer1(const double &arg1, const double &arg2) override {
          auto x=2+2*arg2*arg2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
        MatV parDer2ParDer2(const double &arg1, const double &arg2) override {
          auto x=2+arg1*arg1*2;
          MatV r(2,2);
          r(0,0)=x  ; r(0,1)=2*x;
          r(1,0)=3*x; r(1,1)=4*x;
          return r;
        }
    };
    auto funcN = symbolicFunc<MatV(double,double)>(make_shared<Func>(), arg1, arg2);
    Matrix<General, Var, Var, SymbolicExpression> funcS(2,2);
    funcS(0,0) = arg1*arg1+arg2*arg2+arg1*arg1*arg2*arg2; funcS(0,1) = 2*funcS(0,0);
    funcS(1,0) = 3*funcS(0,0); funcS(1,1) = 4*funcS(0,0);
    auto resN = cos(arg1) *sin(arg2) * pow(funcN(0,0),2)*pow(funcN(0,1),2)*pow(funcN(1,0),2)*pow(funcN(1,1),2);
    auto resS = cos(arg1) *sin(arg2) * pow(funcS(0,0),2)*pow(funcS(0,1),2)*pow(funcS(1,0),2)*pow(funcS(1,1),2);
    i1 ^= 0.9;
    i2 ^= 0.8;
    cout<<check(Eval{resN}(), Eval{resS}())<<endl;
    cout<<check(Eval{parDer(resN,i1)}(), Eval{parDer(resS,i1)}())<<endl;
    cout<<check(Eval{parDer(resN,i2)}(), Eval{parDer(resS,i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i1)}(), Eval{parDer(parDer(resS,i1),i1)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i2)}(), Eval{parDer(parDer(resS,i2),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i1),i2)}(), Eval{parDer(parDer(resS,i1),i2)}())<<endl;
    cout<<check(Eval{parDer(parDer(resN,i2),i1)}(), Eval{parDer(parDer(resS,i2),i1)}())<<endl;
  }

  return 0;  
}
