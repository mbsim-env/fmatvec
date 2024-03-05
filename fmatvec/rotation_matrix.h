#ifndef _FMATVEC_ROTATION_MATRIX_H_
#define _FMATVEC_ROTATION_MATRIX_H_

#include <fmatvec/fmatvec.h>

namespace fmatvec {

using std::sin;
using std::cos;
using std::asin;
using std::atan2;
using std::min;
using std::max;

template<class AT>
Matrix<Rotation, Fixed<3>, Fixed<3>, AT> al2T(const AT& al) {
  Matrix<Rotation, Fixed<3>, Fixed<3>, AT> T;
  AT sa = sin(al);
  AT ca = cos(al);
  T(0,0)= 1.0; T(0,1)= 0.0; T(0,2)= 0.0;
  T(1,0)= 0.0; T(1,1)= ca ; T(1,2)= -sa;
  T(2,0)= 0.0; T(2,1)= sa ; T(2,2)= ca ;
  return T;
}

template<class AT>
Matrix<Rotation, Fixed<3>, Fixed<3>, AT> be2T(const AT& be) {
  Matrix<Rotation, Fixed<3>, Fixed<3>, AT> T;
  AT sb = sin(be);
  AT cb = cos(be);
  T(0,0)= cb ; T(0,1)= 0.0; T(0,2)= sb ;
  T(1,0)= 0.0; T(1,1)= 1.0; T(1,2)= 0.0;
  T(2,0)= -sb; T(2,1)= 0.0; T(2,2)= cb ;
  return T;
}

template<class AT>
Matrix<Rotation, Fixed<3>, Fixed<3>, AT> ga2T(const AT& ga) {
  Matrix<Rotation, Fixed<3>, Fixed<3>, AT> T;
  AT sg = sin(ga);
  AT cg = cos(ga);
  T(0,0)= cg ; T(0,1)= -sg; T(0,2)= 0.0;
  T(1,0)= sg ; T(1,1)= cg ; T(1,2)= 0.0;
  T(2,0)= 0.0; T(2,1)= 0.0; T(2,2)= 1.0;
  return T;
}

template<class AT>
Matrix<Rotation, Fixed<3>, Fixed<3>, AT> albega2T(const Vector<Fixed<3>, AT> &albega) {
  Matrix<Rotation, Fixed<3>, Fixed<3>, AT> T;
  AT sa = sin(albega(0));
  AT sb = sin(albega(1));
  AT sg = sin(albega(2));
  AT ca = cos(albega(0));
  AT cb = cos(albega(1));
  AT cg = cos(albega(2));
  T(0,0) = cb*cg           ; T(0,1) = -sg*cb           ; T(0,2) = sb    ;
  T(1,0) = sa*sb*cg + sg*ca; T(1,1) = -sa*sb*sg + ca*cg; T(1,2) = -sa*cb;
  T(2,0) = sa*sg - sb*ca*cg; T(2,1) = sa*cg + sb*sg*ca ; T(2,2) = ca*cb ;
  return T;
}

template<class AT>
Matrix<Rotation, Fixed<3>, Fixed<3>, AT> begaal2T(const Vector<Fixed<3>, AT> &begaal) {
  Matrix<Rotation, Fixed<3>, Fixed<3>, AT> T;
  AT sb = sin(begaal(0));
  AT sg = sin(begaal(1));
  AT sa = sin(begaal(2));
  AT cb = cos(begaal(0));
  AT cg = cos(begaal(1));
  AT ca = cos(begaal(2));
  T(0,0) = cb*cg ; T(0,1) = sa*sb - sg*ca*cb; T(0,2) = sa*sg*cb + sb*ca ;
  T(1,0) = sg    ; T(1,1) = ca*cg           ; T(1,2) = -sa*cg           ;
  T(2,0) = -sb*cg; T(2,1) = sa*cb + sb*sg*ca; T(2,2) = -sa*sb*sg + ca*cb;
  return T;
}

template<class AT>
Vector<Fixed<3>, AT> T2albega(const Matrix<Rotation, Fixed<3>, Fixed<3>, AT> &T) {
  Vector<Fixed<3>, AT> albega(NONINIT);
  assert(T(0,2)>=-1-1e-12 && T(0,2)<1+1e-12);
  albega(1)= asin(max(min(T(0,2), 1.0), -1.0));
  double nenner = cos(albega(1));
  if (fabs(nenner)>1e-10) {
    albega(0) = atan2(-T(1,2),T(2,2));
    albega(2) = atan2(-T(0,1),T(0,0));
  } else {
    albega(0)=0;
    albega(2)=atan2(T(1,0),T(1,1));
  }
  return albega;
}

template<class AT>
Vector<Fixed<3>, AT> T2begaal(const Matrix<Rotation, Fixed<3>, Fixed<3>, AT> &T) {
  Vector<Fixed<3>, AT> begaal(NONINIT);
  assert(T(1,0)>=-1-1e-12 && T(1,0)<1+1e-12);
  begaal(1)= asin(max(min(T(1,0), 1.0), -1.0));
  double nenner = cos(begaal(1));
  if (fabs(nenner)>1e-10) {
    begaal(0) = atan2(-T(2,0),T(0,0));
    begaal(2) = atan2(-T(1,2),T(1,1));
  } else {
    begaal(0)=0;
    begaal(2)=atan2(T(2,1),T(2,2));
  }
  return begaal;
}

}

#endif
