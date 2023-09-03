#include <cfenv>
#include <cassert>
#include "fmatvec/symbolic.h"
#include <fmatvec/symbolic_function.h>

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

  return 0;  
}
