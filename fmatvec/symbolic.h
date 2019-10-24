#ifndef _fmatvec_symbolic_h_
#define _fmatvec_symbolic_h_

namespace fmatvec {

// definitions of matrix/vector operations/functions only usefull for symbolic calculations
// like building partial derivatives, evaluation of symbolic expressions to double values, ...

template<class DepVec, class ATIndep>
DepVec parDer(const DepVec &dep, const ATIndep &indep) {
  DepVec ret(dep.size());
  for(int i=0; i<dep.size(); ++i)
    ret(i)=parDer(dep(i), indep);
  return ret;
}

template<class Type, class RowShape, class ColShape, class ATDep, class ATIndep>
Matrix<Type, RowShape, ColShape, ATDep> parDer(const Matrix<Type, RowShape, ColShape, ATDep> &dep, const ATIndep &indep) {
  Matrix<Type, RowShape, ColShape, ATDep> ret(dep.rows(), dep.cols());
  for(int r=0; r<dep.rows(); ++r)
    for(int c=0; c<dep.cols(); ++c)
      ret(r,c)=parDer(dep(r,c), indep);
  return ret;
}

template<class DepShape, class IndepShape, class ATDep, class ATIndep>
Matrix<General, DepShape, IndepShape, ATDep> parDer(const Vector<DepShape, ATDep> &dep, const Vector<IndepShape, ATIndep> &indep) {
  Matrix<General, DepShape, IndepShape, ATDep> ret(dep.size(), indep.size());
  for(int r=0; r<dep.size(); ++r)
    for(int c=0; c<indep.size(); ++c)
      ret(r,c)=parDer(dep(r), indep(c));
  return ret;
}

template<class IndepShape, class ATDep, class ATIndep>
RowVector<IndepShape, ATDep> parDer(const ATDep &dep, const Vector<IndepShape, ATIndep> &indep) {
  RowVector<IndepShape, ATDep> ret(indep.size());
  for(int r=0; r<indep.size(); ++r)
    ret(r)=parDer(dep, indep(r));
  return ret;
}

template<class ATDep, class ATIndep>
Vector<Fixed<3>, ATDep> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const ATIndep &x) {
  Matrix<General, Fixed<3>, Fixed<3>, ATDep> Rs;
  for(int r=0; r<3; ++r)
    for(int c=0; c<3; ++c)
      Rs(r,c)=parDer(R(r,c), x);
  auto retTilde=Rs*trans(R);
  Vector<Fixed<3>, ATDep> ret;
  ret(0)=retTilde(2,1);
  ret(1)=retTilde(0,2);
  ret(2)=retTilde(1,0);
  return ret;
}

template<class ATDep, class ATIndep>
Matrix<General, Fixed<3>, Fixed<3>, ATDep> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const Vector<Fixed<3>, ATIndep> &x) {
  Matrix<General, Fixed<3>, Fixed<3>, ATDep> ret;
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

template<class DepR, class ShapeA, class IndepA, class ShapeB, class DepB>
DepR subst(const DepR &se, const Vector<ShapeA, IndepA>& a, const Vector<ShapeB, DepB> &b) {
  if(a.size()!=b.size())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  DepR ret=se;
  for(int i=0; i<a.size(); ++i)
    ret=subst(ret, a(i), b(i));
  return ret;
}

template<class ShapeR, class DepR, class ShapeA, class IndepA, class ShapeB, class DepB>
Vector<ShapeR, DepR> subst(const Vector<ShapeR, DepR> &se, const IndepA& a, const DepB &b) {
  Vector<ShapeR, DepR> ret(se.size());
  for(int i=0; i<se.size(); ++i)
    ret(i)=subst(se(i), a, b);
  return ret;
}

template<class TypeR, class RowShapeR, class ColShapeR, class DepR, class ShapeA, class IndepA, class ShapeB, class DepB>
Matrix<TypeR, RowShapeR, ColShapeR, DepR> subst(const Matrix<TypeR, RowShapeR, ColShapeR, DepR> &se,
                                                const IndepA& a, const DepB &b) {
  Matrix<TypeR, RowShapeR, ColShapeR, DepR> ret(se.rows(), se.cols());
  for(int r=0; r<se.rows(); ++r)
    for(int c=0; c<se.cols(); ++c)
      ret(r,c)=subst(se(r,c), a, b);
  return ret;
}

//mfmf change operator &= to something better
template<class ShapeDst, class ShapeSrc, class ATIndep>
Vector<ShapeDst, ATIndep>& operator&=(Vector<ShapeDst, ATIndep> &dst, const Vector<ShapeSrc, double> &src) {
  assert(dst.size()==src.size());
  for(int i=0; i<dst.size(); ++i)
    dst(i)&=src(i);
  return dst;
}

template<int N, class ATIndep>
Vector<Fixed<N>, ATIndep>& operator&=(Vector<Fixed<N>, ATIndep> &dst, const Vector<Fixed<N>, double> &src) {
  for(int i=0; i<N; ++i)
    dst(i)&=src(i);
  return dst;
}

template<int N, class ShapeSrc, class ATIndep>
Vector<Fixed<N>, ATIndep>& operator&=(Vector<Fixed<N>, ATIndep> &dst, const Vector<ShapeSrc, double> &src) {
  assert(N==src.size());
  for(int i=0; i<N; ++i)
    dst(i)&=src(i);
  return dst;
}

template<class ShapeDst, int N, class ATIndep>
Vector<ShapeDst, ATIndep>& operator&=(Vector<ShapeDst, ATIndep> &dst, const Vector<Fixed<N>, double> &src) {
  assert(N==src.size());
  for(int i=0; i<N; ++i)
    dst(i)&=src(i);
  return dst;
}

}

#endif
