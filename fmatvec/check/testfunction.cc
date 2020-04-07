#include <config.h>
#include <cassert>
#include <cfenv>
#include <iostream>
#include <fmatvec/fmatvec.h>
#include <fmatvec/function.h>
#include <fmatvec/stream_impl.h>
#include <stdexcept>



using namespace std;
using namespace fmatvec;



// Just a example
// r=r(q(t),t)=[cos(q(0));sin(q(1));0]*(r0-v0*t)
class KinematicSpiral : public Function<Vec3(VecV, double)> {
  public:
    KinematicSpiral(double r0_, double v0_) : r0(r0_), v0(v0_) {}
    int getArg1Size() const {
      return 2;
    }
    Vec3 operator()(const VecV &q, const double &t) {
      // we can also use the user data here (a, b); but not required here
      Vec3 r;
      r(0)=+cos(q(0))*(r0-v0*t);
      r(1)=+sin(q(1))*(r0-v0*t);
      r(2)=0;
      return r;
    }
    Mat3xV parDer1(const VecV &q, const double &t) {
      Mat3xV J(2);
      J(0,0)=-sin(q(0))*(r0-v0*t);  J(0,1)=0;
      J(1,0)=0;                     J(1,1)=+cos(q(1))*(r0-v0*t);
      J(2,0)=0;                     J(2,1)=0;
      return J;
    }
    Vec3 parDer2(const VecV &q, const double &t) {
      Vec3 j;
      j(0)=+cos(q(0))*v0;
      j(1)=+sin(q(1))*v0;
      j(2)=0;
      return j;
    }
    Mat3xV parDer1DirDer1(const VecV &qd, const VecV &q, const double &t) {
      Mat3xV Jq(2);
      Jq(0,0)=-cos(q(0))*(r0-v0*t)*qd(0) + 0;  Jq(0,1)=0 + 0                   *qd(1);
      Jq(1,0)=0                   *qd(0) + 0;  Jq(1,1)=0 + -sin(q(1))*(r0-v0*t)*qd(1);
      Jq(2,0)=0                   *qd(0) + 0;  Jq(2,1)=0 + 0                   *qd(1);
      return Jq;
    }
    Mat3xV parDer1ParDer2(const VecV &q, const double &t) {
      Mat3xV Jt(2);
      Jt(0,0)=-sin(q(0))*v0;  Jt(0,1)=0;
      Jt(1,0)=0;              Jt(1,1)=+cos(q(1))*v0;
      Jt(2,0)=0;              Jt(2,1)=0;
      return Jt;
    }
    // only this function OR the previous is required (see fPdjT later)
    Vec3 parDer2DirDer1(const VecV &qd, const VecV &q, const double &t) {
      // this code equals "return parDer1ParDer2(q, t)*qd;"
      Vec3 jq;
      jq(0)=-sin(q(0))*v0*qd(0) + 0            *qd(1);
      jq(1)=0            *qd(0) + +cos(q(1))*v0*qd(1);
      jq(2)=0            *qd(0) + 0            *qd(1);
      return jq;
    }
    Vec3 parDer2ParDer2(const VecV &q, const double &t) {
      Vec3 jt;
      jt(0)=0;
      jt(1)=0;
      jt(2)=0;
      return jt;
    }
  private:
    double r0, v0;
};



// Just a example
// The transformation matrix A rotates about axis a with with angle alpha=2*q(0)+3*q(1)+3*t
class KinematicRotAboutAxis : public Function<RotMat3(VecV, double)> {
  public:
    KinematicRotAboutAxis(const Vec3 &a_) : a(a_) {
    }
    int getArg1Size() const {
      return 2;
    }
    RotMat3 operator()(const VecV &q, const double &t) {
      RotMat3 A;
      // MISSING
      return A;
    }
    Mat3xV parDer1(const VecV &q, const double &t) {
      Mat3xV J(2);
      // MISSING
      return J;
    }
    Vec3 parDer2(const VecV &q, const double &t) {
      Vec3 j;
      // MISSING
      return j;
    }
    Mat3xV parDer1DirDer1(const VecV &qd, const VecV &q, const double &t) {
      Mat3xV Jq(2);
      // MISSING
      return Jq;
    }
    Mat3xV parDer1ParDer2(const VecV &q, const double &t) {
      Mat3xV Jt(2);
      // MISSING
      return Jt;
    }
    Vec3 parDer2ParDer2(const VecV &q, const double &t) {
      Vec3 jt;
      // MISSING
      return jt;
    }
  private:
    Vec3 a;
};



int main() {
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  // define some arbitary argument for the following function calles
  VecV q(2);
  q(0)=3.5;
  q(1)=2.6;
  VecV qd(2);
  qd(0)=1.2;
  qd(1)=1.7;
  double t=2.2;



  // create the example function object (translation)
  Function<Vec3(VecV, double)> *T=new KinematicSpiral(0.2, 5.3);

  // return some values (required for MBSim)

  cout<<"arg1Size = "<<T->getArg1Size()<<endl;
//cout<<"arg2Size = "<<T->getArg2Size()<<endl; // will give a compile error, since ErrorType is not defined by "cout<<" (Nice)

  cout<<"fPrPK    = "<<(*T)(q, t)<<endl;

  cout<<"fPJT     = "<<T->parDer1(q, t)<<endl;
  cout<<"fPjT     = "<<T->parDer2(q, t)<<endl;

  cout<<"fPdJT    = "<<T->parDer1DirDer1(qd, q, t)+T->parDer1ParDer2(q, t)<<endl;
  // fPdjT can be calculate using parDer2DirDer1 OR parDer1ParDer2 in this special case (see KinematicSpiral above)
  cout<<"fPdjT    = "<<T->parDer2DirDer1(qd, q, t)+T->parDer2ParDer2(q, t)<<endl;
  cout<<"fPdjT    = "<<T->parDer1ParDer2(q, t)*qd+T->parDer2ParDer2(q, t)<<endl;



  // create the example function object (rotation)
  Function<RotMat3(VecV, double)> *R=new KinematicRotAboutAxis(Vec3("[1;0;0]"));

  // return some values (required for MBSim)

  cout<<"arg1Size = "<<R->getArg1Size()<<endl;
//cout<<"arg2Size = "<<R->getArg2Size()<<endl; // will give a compile error, since ErrorType is not defined by "cout<<" (Nice)

  cout<<"fAPK     = "<<(*R)(q, t)<<endl;

  cout<<"fPJR     = "<<R->parDer1(q, t)<<endl;
  cout<<"fPjR     = "<<R->parDer2(q, t)<<endl;

  cout<<"fPdJT    = "<<R->parDer1DirDer1(qd, q, t)+R->parDer1ParDer2(q, t)<<endl;
  cout<<"fPdjR    = "<<R->parDer1ParDer2(q, t)*qd+R->parDer2ParDer2(q, t)<<endl;



  string inputv("[1; 2\n 3 ]");
  { istringstream str(inputv); Vec     v; str>>v; cout<<v<<endl; }
  { istringstream str(inputv); RowVec  v; str>>v; cout<<v<<endl; }
  { istringstream str(inputv); Vec3    v; str>>v; cout<<v<<endl; }
  { istringstream str(inputv); VecV    v; str>>v; cout<<v<<endl; }

  string inputm("[1 2 3; 4, 5 6\n7 8,  9 ]");
  string inputs("[1 2 3; 2, 5 6\n3 6,  9 ]");
  string inputd("[1 0 0; 0, 5 0\n0 0,  9 ]");
  { istringstream str(inputm); Mat       m; str>>m; cout<<m<<endl; }
  { istringstream str(inputm); SqrMat    m; str>>m; cout<<m<<endl; }
  { istringstream str(inputd); DiagMat   m; str>>m; cout<<m<<endl; }
  { istringstream str(inputs); SymMat    m; str>>m; cout<<m<<endl; }
  { istringstream str(inputm); Mat3x3    m; str>>m; cout<<m<<endl; }
  { istringstream str(inputm); SqrMat3   m; str>>m; cout<<m<<endl; }
  { istringstream str(inputs); SymMat3   m; str>>m; cout<<m<<endl; }
  { istringstream str(inputm); MatV      m; str>>m; cout<<m<<endl; }
  { istringstream str(inputm); SqrMatV   m; str>>m; cout<<m<<endl; }
  { istringstream str(inputm); Mat3xV    m; str>>m; cout<<m<<endl; }

  { istringstream str("[3]"); Mat m; str>>m; cout<<m<<endl; }
  { istringstream str("3"  ); Mat m; str>>m; cout<<m<<endl; }
  { istringstream str("[3]"); VecV v; str>>v; cout<<v<<endl; }
  { istringstream str("3"  ); VecV v; str>>v; cout<<v<<endl; }

  return 0;  
}
