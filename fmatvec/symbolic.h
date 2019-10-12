#ifndef _fmatvec_symbolic_h_
#define _fmatvec_symbolic_h_

namespace fmatvec {

template<class DepVec, class AT>
DepVec parDer(const DepVec &dep, const AT &indep) {
  DepVec ret(dep.size());
  for(int i=0; i<dep.size(); ++i)
    ret(i)=parDer(dep(i), indep);
  return ret;
}

template<class Type, class RowShape, class ColShape, class AT>
Matrix<Type, RowShape, ColShape, AT> parDer(const Matrix<Type, RowShape, ColShape, AT> &dep, const AT &indep) {
  Matrix<Type, RowShape, ColShape, AT> ret(dep.rows(), dep.cols());
  for(int r=0; r<dep.rows(); ++r)
    for(int c=0; c<dep.cols(); ++c)
      ret(r,c)=parDer(dep(r,c), indep);
  return ret;
}

template<class DepShape, class IndepShape, class AT>
Matrix<General, DepShape, IndepShape, AT> parDer(const Vector<DepShape, AT> &dep, const Vector<IndepShape, AT> &indep) {
  Matrix<General, DepShape, IndepShape, AT> ret(dep.size(), indep.size());
  for(int r=0; r<dep.size(); ++r)
    for(int c=0; c<indep.size(); ++c)
      ret(r,c)=parDer(dep(r), indep(c));
  return ret;
}

template<class IndepShape, class AT>
RowVector<IndepShape, AT> parDer(const AT &dep, const Vector<IndepShape, AT> &indep) {
  RowVector<IndepShape, AT> ret(indep.size());
  for(int r=0; r<indep.size(); ++r)
    ret(r)=parDer(dep, indep(r));
  return ret;
}

template<class AT>
Vector<Fixed<3>, AT> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, AT> &R, const AT &x) {
  Matrix<General, Fixed<3>, Fixed<3>, AT> Rs;
  for(int r=0; r<3; ++r)
    for(int c=0; c<3; ++c)
      Rs(r,c)=parDer(R(r,c), x);
  auto retTilde=Rs*trans(R);
  Vector<Fixed<3>, AT> ret;
  ret(0)=retTilde(2,1);
  ret(1)=retTilde(0,2);
  ret(2)=retTilde(1,0);
  return ret;
}

template<class AT>
Matrix<General, Fixed<3>, Fixed<3>, AT> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, AT> &R, const Vector<Fixed<3>, AT> &x) {
  Matrix<General, Fixed<3>, Fixed<3>, AT> ret;
  for(int i=0; i<3; ++i)
    ret.set(i, parDer(R, x(i)));
  return ret;
}

}

#endif
