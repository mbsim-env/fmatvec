#ifndef _fmatvec_symbolic_h_
#define _fmatvec_symbolic_h_

#include "ast.h"
#include "function.h"
#include <type_traits>

namespace fmatvec {

// definitions of matrix/vector operations/functions only usefull for symbolic calculations
// like building partial derivatives, evaluation of symbolic expressions to double values, ...

template<class Dep, class ATIndep>
Dep parDer(const Dep &dep, const ATIndep &indep) {
  Dep ret;
  if constexpr (std::is_same_v<Dep, SymbolicExpression>) {
  }
  else if constexpr (Dep::isVector)
    ret.resize(dep.size());
  else
    ret.resize(dep.rows(), dep.cols());
  auto d=dep.begin();
  auto r=ret.begin();
  for(; d!=dep.end(); ++d, ++r)
    *r=parDer(*d, indep);
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

template<class ATDep, class Shape, class ATIndep>
Matrix<General, Fixed<3>, Fixed<3>, ATDep> parDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const Vector<Shape, ATIndep> &x) {
  Matrix<General, Fixed<3>, Shape, ATDep> ret;
  for(int i=0; i<x.size(); ++i)
    ret.set(i, parDer(R, x(i)));
  return ret;
}

template<class Dep, class ATDep, class ATIndep>
Dep dirDer(const Dep &dep, const ATDep &indepdir, const ATIndep &indep) {
  return parDer(dep, indep)*indepdir;
}

template<class ATDep, class ATIndep>
Vector<Fixed<3>, ATDep> dirDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &R, const ATDep &indepdir, const ATIndep &indep) {
  return parDer(R, indep)*indepdir;
}

template<class Dep, class ShapeIndep, class ATDep, class ATIndep>
Dep dirDer(const Dep &dep, const Vector<ShapeIndep, ATDep> &indepdir, const Vector<ShapeIndep, ATIndep> &indep) {
  Dep ret;
  if constexpr (std::is_same_v<Dep, SymbolicExpression>) {
    ret=0.0;
  }
  else if constexpr (Dep::isVector)
    ret=Dep(dep.size(), INIT, 0.0);
  else
    ret=Dep(dep.rows(), dep.cols(), INIT, 0.0);
  for(int i=0; i<indep.size(); ++i)
    ret+=parDer(dep, indep(i))*indepdir(i);
  return ret;
}

template<class ATDep, class ShapeIndep, class ATIndep>
Vector<Fixed<3>, ATDep> dirDer(const Matrix<Rotation, Fixed<3>, Fixed<3>, ATDep> &dep, const Vector<ShapeIndep, ATDep> &indepdir, const Vector<ShapeIndep, ATIndep> &indep) {
  return parDer(dep, indep)*indepdir;
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

template<class Type, class RowShape, class ColShape, class AT>
Matrix<Type, RowShape, ColShape, double> eval(const Matrix<Type, RowShape, ColShape, AT> &x) {
  Matrix<Type, RowShape, ColShape, double> ret(x.rows(), x.cols());
  for(int r=0; r<x.rows(); ++r)
    for(int c=0; c<x.cols(); ++c)
      ret(r,c)=eval(x(r,c));
  return ret;
}

template<class Shape, class AT>
SquareMatrix<Shape, double> eval(const SquareMatrix<Shape, AT> &x) {
  SquareMatrix<Shape, double> ret(x.size());
  for(int r=0; r<x.size(); ++r)
    for(int c=0; c<x.size(); ++c)
      ret(r,c)=eval(x(r,c));
  return ret;
}

// substitute scalar in scalar expression -> defined in ast.h

// substitute vector in scalar expression
template<class DepR, class ShapeA, class IndepA, class ShapeB, class DepB>
DepR subst(const DepR &se, const Vector<ShapeA, IndepA>& a, const Vector<ShapeB, DepB> &b) {
  if(a.size()!=b.size())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  auto ret=se;
  for(int i=0; i<a.size(); ++i)
    ret=subst(ret, a(i), b(i));
  return ret;
}

// substitute matrix in scalar expression
template<class DepR, class TypeA, class RowA, class ColA, class IndepA, class TypeB, class RowB, class ColB, class DepB>
DepR subst(const DepR &se, const Matrix<TypeA, RowA, ColA, IndepA>& A, const Matrix<TypeB, RowB, ColB, DepB> &B) {
  if(A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  auto ret=se;
  for(int r=0; r<A.rows(); ++r)
    for(int c=0; c<A.cols(); ++c)
      ret=subst(ret, A(r,c), B(r,c));
  return ret;
}

// substitute scalar in vector expression
template<class ShapeR, class DepR, class IndepA, class DepB>
Vector<ShapeR, DepR> subst(const Vector<ShapeR, DepR> &se, const IndepA& a, const DepB &b) {
  Vector<ShapeR, DepR> ret(se.size());
  for(int i=0; i<se.size(); ++i)
    ret(i)=subst(se(i), a, b);
  return ret;
}

// substitute vector in vector expression
template<class ShapeR, class DepR, class ShapeA, class IndepA, class ShapeB, class DepB>
Vector<ShapeR, DepR> subst(const Vector<ShapeR, DepR> &se, const Vector<ShapeA, IndepA>& a, const Vector<ShapeB, DepB> &b) {
  if(a.size()!=b.size())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  auto ret=se;
  for(int i=0; i<a.size(); ++i)
    ret=subst(ret, a(i), b(i));
  return ret;
}

// substitute matrix in vector expression
template<class ShapeR, class DepR, class TypeA, class RowA, class ColA, class IndepA, class TypeB, class RowB, class ColB, class DepB>
Vector<ShapeR, DepR> subst(const Vector<ShapeR, DepR> &se, const Matrix<TypeA, RowA, ColA, IndepA>& A, const Matrix<TypeB, RowB, ColB, DepB> &B) {
  if(A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  auto ret=se;
  for(int r=0; r<A.rows(); ++r)
    for(int c=0; c<A.cols(); ++c)
      ret=subst(ret, A(r,c), B(r,c));
  return ret;
}

// substitute scalar in matrix expression
template<class TypeR, class RowR, class ColR, class DepR, class IndepA, class DepB>
Matrix<TypeR, RowR, ColR, DepR> subst(const Matrix<TypeR, RowR, ColR, DepR> &SE, const IndepA& a, const DepB &b) {
  Matrix<TypeR, RowR, ColR, DepR> ret(SE.rows(), SE.cols());
  for(int r=0; r<SE.rows(); ++r)
    for(int c=0; c<SE.cols(); ++c)
      ret(r,c)=subst(SE(r,c), a, b);
  return ret;
}

// substitute vector in matrix expression
template<class TypeR, class RowR, class ColR, class DepR, class ShapeA, class IndepA, class ShapeB, class DepB>
Matrix<TypeR, RowR, ColR, DepR> subst(const Matrix<TypeR, RowR, ColR, DepR> &SE, const Vector<ShapeA, IndepA>& a, const Vector<ShapeB, DepB> &b) {
  if(a.size()!=b.size())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  auto ret=SE;
  for(int i=0; i<a.size(); ++i)
    ret=subst(ret, a(i), b(i));
  return ret;
}

// substitute matrix in matrix expression
template<class TypeR, class RowR, class ColR, class DepR, class TypeA, class RowA, class ColA, class IndepA, class TypeB, class RowB, class ColB, class DepB>
Matrix<TypeR, RowR, ColR, DepR> subst(const Matrix<TypeR, RowR, ColR, DepR> &SE, const Matrix<TypeA, RowA, ColA, IndepA>& A, const Matrix<TypeB, RowB, ColB, DepB> &B) {
  if(A.rows()!=B.rows() || A.cols()!=B.cols())
    throw std::runtime_error("The size of the independent and the dependent substitution variable does not match.");
  auto ret=SE;
  for(int r=0; r<A.rows(); ++r)
    for(int c=0; c<A.cols(); ++c)
      ret=subst(ret, A(r,c), B(r,c));
  return ret;
}



template<class ShapeDst, class ShapeSrc, class ATIndep>
Vector<ShapeDst, ATIndep>& operator^=(Vector<ShapeDst, ATIndep> &dst, const Vector<ShapeSrc, double> &src) {
  FMATVEC_ASSERT(dst.size()==src.size(), ATIndep);
  for(int i=0; i<dst.size(); ++i)
    dst(i)^=src(i);
  return dst;
}

template<int N, class ATIndep>
Vector<Fixed<N>, ATIndep>& operator^=(Vector<Fixed<N>, ATIndep> &dst, const Vector<Fixed<N>, double> &src) {
  for(int i=0; i<N; ++i)
    dst(i)^=src(i);
  return dst;
}

template<int N, class ShapeSrc, class ATIndep>
Vector<Fixed<N>, ATIndep>& operator^=(Vector<Fixed<N>, ATIndep> &dst, const Vector<ShapeSrc, double> &src) {
  FMATVEC_ASSERT(N==src.size(), ATIndep);
  for(int i=0; i<N; ++i)
    dst(i)^=src(i);
  return dst;
}

template<class ShapeDst, int N, class ATIndep>
Vector<ShapeDst, ATIndep>& operator^=(Vector<ShapeDst, ATIndep> &dst, const Vector<Fixed<N>, double> &src) {
  FMATVEC_ASSERT(N==src.size(), ATIndep);
  for(int i=0; i<N; ++i)
    dst(i)^=src(i);
  return dst;
}

template<class TypeDst, class TypeSrc, class ShapeRowDst, class ShapeRowSrc, class ShapeColDst, class ShapeColSrc, class ATIndep>
Matrix<TypeDst, ShapeRowDst, ShapeColDst, ATIndep>& operator^=(Matrix<TypeDst, ShapeRowDst, ShapeColDst, ATIndep> &dst,
                                                               const Matrix<TypeSrc, ShapeRowSrc, ShapeColSrc, double> &src) {
  FMATVEC_ASSERT(dst.rows()==src.rows(), ATIndep);
  FMATVEC_ASSERT(dst.cols()==src.cols(), ATIndep);
  for(int r=0; r<dst.rows(); ++r)
    for(int c=0; c<dst.cols(); ++c)
      dst(r,c)^=src(r,c);
  return dst;
}

}

#endif
