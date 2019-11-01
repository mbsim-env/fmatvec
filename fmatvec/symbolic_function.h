#ifndef _FMATVEC_SYMBOLIC_FUNCTION_H_
#define _FMATVEC_SYMBOLIC_FUNCTION_H_

#include "function.h"
#include "ast.h"
#include "symbolic.h"

namespace fmatvec {

// Replace in a fmatvec scalar, vector or matrix the AT with ATNew.
template<class MatVec, class NewAT>
struct ReplaceAT;

template<class NewAT>
struct ReplaceAT<ErrorType, NewAT> {
  using Type = ErrorType;
};

template<class NewAT>
struct ReplaceAT<double, NewAT> {
  using Type = NewAT;
};

template<class Shape, class NewAT>
struct ReplaceAT<Vector<Shape, double>, NewAT> {
  using Type = Vector<Shape, NewAT>;
};

template<class Shape, class NewAT>
struct ReplaceAT<RowVector<Shape, double>, NewAT> {
  using Type = RowVector<Shape, NewAT>;
};

template<class MatType, class RowShape, class ColShape, class NewAT>
struct ReplaceAT<Matrix<MatType, RowShape, ColShape, double>, NewAT> {
  using Type = Matrix<MatType, RowShape, ColShape, NewAT>;
};

template<typename Sig>
class SymbolicFunction;

#define RET Ret
#define ARG ATArg
#define TEMPLATE typename Ret, typename ATArg
#define PARDER
#define PARDERPARDER
#define ARGDIRINIT 
#include "symbolic_function_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER
#undef ARGDIRINIT

#define RET ATRet
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename ATRet, typename ArgShape, typename ATArg
#define PARDER
//#define PARDERPARDER
#define ARGDIRINIT NONINIT
#include "symbolic_function_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER
#undef ARGDIRINIT

#define RET RowVector<RetShape, ATArg>
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename RetShape, typename ArgShape, typename ATArg
//#define PARDER
//#define PARDERPARDER
#define ARGDIRINIT NONINIT
#include "symbolic_function_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER
#undef ARGDIRINIT

#define RET Matrix<Type, ColShape, RowShape, ATArg>
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename Type, typename ColShape, typename RowShape, typename ArgShape, typename ATArg
//#define PARDER
//#define PARDERPARDER
#define ARGDIRINIT NONINIT
#include "symbolic_function_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER
#undef ARGDIRINIT

#define RET Matrix<Rotation, ColShape, RowShape, ATArg>
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename ColShape, typename RowShape, typename ArgShape, typename ATArg
#define PARDER
//#define PARDERPARDER
#define ARGDIRINIT NONINIT
#include "symbolic_function_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER
#undef ARGDIRINIT

}

#endif
