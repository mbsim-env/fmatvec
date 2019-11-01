#include <cfenv>
#include <cassert>
#include "fmatvec/symbolic.h"
#include <fmatvec/symbolic_function.h>

using namespace std;
using namespace fmatvec;

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

    cout<<funcS(arg)<<endl;
    cout<<funcV(arg)<<endl;
    cout<<funcRV(arg)<<endl;
    cout<<funcM(arg)<<endl;
    cout<<funcR(arg)<<endl;

    cout<<funcS.parDer(arg)<<endl;
    cout<<funcV.parDer(arg)<<endl;
    cout<<funcRV.parDer(arg)<<endl;
    cout<<funcM.parDer(arg)<<endl;
    cout<<funcR.parDer(arg)<<endl;

    cout<<funcS.dirDer(argd, arg)<<endl;
    cout<<funcV.dirDer(argd, arg)<<endl;
    cout<<funcRV.dirDer(argd, arg)<<endl;
    cout<<funcM.dirDer(argd, arg)<<endl;
    cout<<funcR.dirDer(argd, arg)<<endl;

    cout<<funcS.parDerParDer(arg)<<endl;
    cout<<funcV.parDerParDer(arg)<<endl;
    cout<<funcRV.parDerParDer(arg)<<endl;
    cout<<funcM.parDerParDer(arg)<<endl;
    cout<<funcR.parDerParDer(arg)<<endl;

    cout<<funcS.parDerDirDer(argd, arg)<<endl;
    cout<<funcV.parDerDirDer(argd, arg)<<endl;
    cout<<funcRV.parDerDirDer(argd, arg)<<endl;
    cout<<funcM.parDerDirDer(argd, arg)<<endl;
    cout<<funcR.parDerDirDer(argd, arg)<<endl;

    cout<<funcS.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcV.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcRV.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcM.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcR.dirDerDirDer(argd, argd2, arg)<<endl;
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

    cout<<funcS(arg)<<endl;
    cout<<funcV(arg)<<endl;
    cout<<funcRV(arg)<<endl;
    cout<<funcM(arg)<<endl;
    cout<<funcR(arg)<<endl;

    cout<<funcS.parDer(arg)<<endl;
    cout<<funcV.parDer(arg)<<endl;
//    cout<<funcRV.parDer(arg)<<endl;
//    cout<<funcM.parDer(arg)<<endl;
    cout<<funcR.parDer(arg)<<endl;

    cout<<funcS.dirDer(argd, arg)<<endl;
    cout<<funcV.dirDer(argd, arg)<<endl;
    cout<<funcRV.dirDer(argd, arg)<<endl;
    cout<<funcM.dirDer(argd, arg)<<endl;
    cout<<funcR.dirDer(argd, arg)<<endl;

//    cout<<funcS.parDerParDer(arg)<<endl;
//    cout<<funcV.parDerParDer(arg)<<endl;
//    cout<<funcRV.parDerParDer(arg)<<endl;
//    cout<<funcM.parDerParDer(arg)<<endl;
//    cout<<funcR.parDerParDer(arg)<<endl;

    cout<<funcS.parDerDirDer(argd, arg)<<endl;
    cout<<funcV.parDerDirDer(argd, arg)<<endl;
//    cout<<funcRV.parDerDirDer(argd, arg)<<endl;
//    cout<<funcM.parDerDirDer(argd, arg)<<endl;
    cout<<funcR.parDerDirDer(argd, arg)<<endl;

    cout<<funcS.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcV.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcRV.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcM.dirDerDirDer(argd, argd2, arg)<<endl;
    cout<<funcR.dirDerDirDer(argd, argd2, arg)<<endl;
  }

  return 0;  
}
