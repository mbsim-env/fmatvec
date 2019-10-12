#ifndef _fmatvec_symbolic_h_
#define _fmatvec_symbolic_h_

namespace fmatvec {

// definitions of matrix/vector operations/functions only usefull for symbolic calculations
// like building partial derivatives, evaluation of symbolic expressions to double values, ...

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

template<class Shape, class AT>
Vector<Shape, double> eval(const Vector<Shape, AT> &x) {
  Vector<Shape, double> ret(x.size());
  for(int i=0; i<x.size(); ++i)
    ret(i)=eval(x(i));
  return ret;
}

template<class Shape, class AT>
RowVector<Shape, double> eval(const RowVector<Shape, AT> &x) {
  RowVector<Shape, double> ret(x.size());
  for(int i=0; i<x.size(); ++i)
    ret(i)=eval(x(i));
  return ret;
}

template<int N, class AT>
Vector<Fixed<N>, double> eval(const Vector<Fixed<N>, AT> &x) {
  Vector<Fixed<N>, double> ret;
  for(int i=0; i<N; ++i)
    ret(i)=eval(x(i));
  return ret;
}

template<int N, class AT>
RowVector<Fixed<N>, double> eval(const RowVector<Fixed<N>, AT> &x) {
  RowVector<Fixed<N>, double> ret;
  for(int i=0; i<N; ++i)
    ret(i)=eval(x(i));
  return ret;
}

template<class Type, class RowShape, class ColShape, class AT>
Matrix<Type, RowShape, ColShape, double> eval(const Matrix<Type, RowShape, ColShape, AT> &x) {
  Matrix<Type, RowShape, ColShape, double> ret(x.rows(), x.cols());
  for(int r=0; r<x.rows(); ++r)
    for(int c=0; c<x.cols(); ++c)
      ret(r,c)=eval(x(r,c));
  return ret;
}

template<class Type, int N, class ColShape, class AT>
Matrix<Type, Fixed<N>, ColShape, double> eval(const Matrix<Type, Fixed<N>, ColShape, AT> &x) {
  Matrix<Type, Fixed<N>, ColShape, double> ret(x.cols());
  for(int r=0; r<N; ++r)
    for(int c=0; c<x.cols(); ++c)
      ret(r,c)=eval(x(r,c));
  return ret;
}

template<class Type, class RowShape, int M, class AT>
Matrix<Type, RowShape, Fixed<M>, double> eval(const Matrix<Type, RowShape, Fixed<M>, AT> &x) {
  Matrix<Type, RowShape, Fixed<M>, double> ret(x.rows());
  for(int r=0; r<x.rows(); ++r)
    for(int c=0; c<M; ++c)
      ret(r,c)=eval(x(r,c));
  return ret;
}

template<class Type, int N, int M, class AT>
Matrix<Type, Fixed<N>, Fixed<M>, double> eval(const Matrix<Type, Fixed<N>, Fixed<M>, AT> &x) {
  Matrix<Type, Fixed<N>, Fixed<M>, double> ret;
  for(int r=0; r<N; ++r)
    for(int c=0; c<M; ++c)
      ret(r,c)=eval(x(r,c));
  return ret;
}

}

#endif
