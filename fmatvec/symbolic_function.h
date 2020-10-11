#ifndef _FMATVEC_SYMBOLIC_FUNCTION_H_
#define _FMATVEC_SYMBOLIC_FUNCTION_H_

#include "function.h"
#include "ast.h"
#include "symbolic.h"

namespace fmatvec {

namespace {
  template<class ATIndep>
  struct Helper {
    static void initIndep(ATIndep &x, int size) { x=ATIndep(); }
    static int size1(const ATIndep &x) { return 1; }
    static int size2(const ATIndep &x) { return 1; }
  };

  template<class Shape, class ATIndep>
  struct Helper<Vector<Shape, ATIndep>> {
    static void initIndep(Vector<Shape, ATIndep> &x, int size) { x <<= Vector<Shape, ATIndep>(size, NONINIT); }
    static int size1(const Vector<Shape, ATIndep> &x) { return x.size(); }
    static int size2(const Vector<Shape, ATIndep> &x) { return 1; }
  };

  template<class Type, class RowShape, class ColShape, class ATIndep>
  struct Helper<Matrix<Type, RowShape, ColShape, ATIndep>> {
    static int size1(const Matrix<Type, RowShape, ColShape, ATIndep> &x) { return x.rows(); }
    static int size2(const Matrix<Type, RowShape, ColShape, ATIndep> &x) { return x.cols(); }
  };
}

template<typename Sig>
class SymbolicFunction;

#define RET Ret
#define ARG ATArg
#define TEMPLATE typename Ret, typename ATArg
#define PARDER
#define PARDERPARDER
#include "symbolic_function1_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER

#define RET Ret
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename Ret, typename ArgShape, typename ATArg
#define PARDER
//#define PARDERPARDER
#include "symbolic_function1_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER

#define RET RowVector<RetShape, ATRet>
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename RetShape, typename ATRet, typename ArgShape, typename ATArg
//#define PARDER
//#define PARDERPARDER
#include "symbolic_function1_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER

#define RET Matrix<Type, ColShape, RowShape, ATRet>
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename Type, typename ColShape, typename RowShape, typename ATRet, typename ArgShape, typename ATArg
//#define PARDER
//#define PARDERPARDER
#include "symbolic_function1_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER

#define RET Matrix<Rotation, Fixed<3>, Fixed<3>, ATRet>
#define ARG Vector<ArgShape, ATArg>
#define TEMPLATE typename ATRet, typename ArgShape, typename ATArg
#define PARDER
//#define PARDERPARDER
#include "symbolic_function1_temp.h"
#undef RET
#undef ARG
#undef TEMPLATE
#undef PARDER
#undef PARDERPARDER



#define RET Ret
#define ARG1 ATArg1
#define ARG2 ATArg2
#define TEMPLATE typename Ret, typename ATArg1, typename ATArg2
#define PARDER1
#define PARDER2
#define PARDER1PARDER1
#define PARDER2PARDER2
#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Ret
#define ARG1 Vector<Shape, ATArg>
#define ARG2 ATArg
#define TEMPLATE typename Ret, typename Shape, typename ATArg
#define PARDER1
#define PARDER2
//#define PARDER1PARDER1
#define PARDER2PARDER2
#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Matrix<Rotation, Fixed<3>, Fixed<3>, ATRet>
#define ARG1 Vector<Shape, ATArg>
#define ARG2 ATArg
#define TEMPLATE typename ATRet, typename Shape, typename ATArg
#define PARDER1
#define PARDER2
//#define PARDER1PARDER1
#define PARDER2PARDER2
#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Ret
#define ARG1 ATArg
#define ARG2 Vector<Shape, ATArg>
#define TEMPLATE typename Ret, typename Shape, typename ATArg
#define PARDER1
#define PARDER2
#define PARDER1PARDER1
//#define PARDER2PARDER2
#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Matrix<Rotation, Fixed<3>, Fixed<3>, ATRet>
#define ARG1 ATArg
#define ARG2 Vector<Shape, ATArg>
#define TEMPLATE typename ATRet, typename Shape, typename ATArg
#define PARDER1
#define PARDER2
#define PARDER1PARDER1
//#define PARDER2PARDER2
#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Ret
#define ARG1 Vector<Shape1, ATArg>
#define ARG2 Vector<Shape2, ATArg>
#define TEMPLATE typename Ret, typename Shape1, typename Shape2, typename ATArg
#define PARDER1
#define PARDER2
//#define PARDER1PARDER1
//#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Matrix<Rotation, Fixed<3>, Fixed<3>, ATRet>
#define ARG1 Vector<Shape1, ATArg>
#define ARG2 Vector<Shape2, ATArg>
#define TEMPLATE typename ATRet, typename Shape1, typename Shape2, typename ATArg
#define PARDER1
#define PARDER2
//#define PARDER1PARDER1
//#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET RowVector<RetShape, ATRet>
#define ARG1 Vector<Shape, ATArg>
#define ARG2 ATArg
#define TEMPLATE typename RetShape, typename ATRet, typename Shape, typename ATArg
//#define PARDER1
#define PARDER2
//#define PARDER1PARDER1
#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Matrix<Type, RowShape, ColShape, ATRet>
#define ARG1 Vector<Shape, ATArg>
#define ARG2 ATArg
#define TEMPLATE typename Type, typename RowShape, typename ColShape, typename ATRet, typename Shape, typename ATArg
//#define PARDER1
#define PARDER2
//#define PARDER1PARDER1
#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET RowVector<RetShape, ATRet>
#define ARG1 ATArg
#define ARG2 Vector<Shape, ATArg>
#define TEMPLATE typename RetShape, typename ATRet, typename Shape, typename ATArg
#define PARDER1
//#define PARDER2
#define PARDER1PARDER1
//#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Matrix<Type, RowShape, ColShape, ATRet>
#define ARG1 ATArg
#define ARG2 Vector<Shape, ATArg>
#define TEMPLATE typename Type, typename RowShape, typename ColShape, typename ATRet, typename Shape, typename ATArg
#define PARDER1
//#define PARDER2
#define PARDER1PARDER1
//#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET RowVector<RetShape, ATRet>
#define ARG1 Vector<Shape1, ATArg>
#define ARG2 Vector<Shape2, ATArg>
#define TEMPLATE typename RetShape, typename ATRet, typename Shape1, typename ATArg, typename Shape2
//#define PARDER1
//#define PARDER2
//#define PARDER1PARDER1
//#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

#define RET Matrix<Type, RowShape, ColShape, ATRet>
#define ARG1 Vector<Shape1, ATArg>
#define ARG2 Vector<Shape2, ATArg>
#define TEMPLATE typename Type, typename RowShape, typename ColShape, typename ATRet, typename Shape1, typename ATArg, typename Shape2
//#define PARDER1
//#define PARDER2
//#define PARDER1PARDER1
//#define PARDER2PARDER2
//#define PARDER1PARDER2
#include "symbolic_function2_temp.h"
#undef RET
#undef ARG1
#undef ARG2
#undef TEMPLATE
#undef PARDER1
#undef PARDER2
#undef PARDER1PARDER1
#undef PARDER2PARDER2
#undef PARDER1PARDER2

}

#endif





// 1: definition
// 2: value
// 3: parDer1
// 4: dirDer1
// 5: parDer2
// 6: dirDer2
// 7: parDer1ParDer1
// 8: parDer1DirDer1
// 9: dirDer1DirDer1
// 10: parDer2ParDer2
// 11: parDer2DirDer2
// 12: dirDer2DirDer2
// 13: parDer1ParDer2
// 14: parDer1DirDer2
// 15: dirDer2DirDer1
// 16: parDer2DirDer1
// 
// 1   2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
// 
// sss s s s s s s s s s  s  s  s  s  s  s        
// vss v v v v v v v v v  v  v  v  v  v  v        
// rss r r r r r r r r r  r  r  r  r  r  r        
// Mss M M M M M M M M M  M  M  M  M  M  M        
// Rss R v v v v v v v v  v  v  v  v  v  v        
// svv s r s r s E r s E  r  s  E  r  s  r        
// vvv v M v M v E M v E  M  v  E  M  v  M        
// rvv r E r E r E E r E  E  r  E  E  r  E        
// Mvv M E M E M E E M E  E  M  E  E  M  E        
// Rvv R M v M v E M v E  M  v  E  M  v  M        
// ssv s s s r s s s s E  r  s  r  s  s  r        
// vsv v v v M v v v v E  M  v  M  v  v  M        
// rsv r r r E r r r r E  E  r  E  r  r  E        
// Msv M M M E M M M M E  E  M  E  M  M  E        
// Rsv R v v M v v v v E  M  v  M  v  v  M        
// svs s r s s s E r s s  s  s  r  r  s  s        
// vvs v M v v v E M v v  v  v  M  M  v  v        
// rvs r E r r r E E r r  r  r  E  E  r  r        
// Mvs M E M M M E E M M  M  M  E  E  M  M        
// Rvs R M v v v E M v v  v  v  M  M  v  v        
// 
// 1   2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
// 
// 
// 
// 1   2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
// 
// aaa s s s s s s s s s  s  s  s  s  s  s        
// vaa v v v v v v v v v  v  v  v  v  v  v        
// raa r r r r r r r r r  r  r  r  r  r  r        
// Maa M M M M M M M M M  M  M  M  M  M  M        
// Raa R v v v v v v v v  v  v  v  v  v  v        
// 
// ava s r s s s E r s s  s  s  r  r  s  s        
// vva v M v v v E M v v  v  v  M  M  v  v        
// Rva R M v v v E M v v  v  v  M  M  v  v        
// 
// aav s s s r s s s s E  r  s  r  s  s  r        
// vav v v v M v v v v E  M  v  M  v  v  M        
// Rav R v v M v v v v E  M  v  M  v  v  M        
// 
// avv s r s r s E r s E  r  s  E  r  s  r        
// vvv v M v M v E M v E  M  v  E  M  v  M        
// Rvv R M v M v E M v E  M  v  E  M  v  M        
// 
// rva r E r r r E E r r  r  r  E  E  r  r        
// Mva M E M M M E E M M  M  M  E  E  M  M        
// 
// rav r r r E r r r r E  E  r  E  r  r  E        
// Mav M M M E M M M M E  E  M  E  M  M  E        
// 
// rvv r E r E r E E r E  E  r  E  E  r  E        
// Mvv M E M E M E E M E  E  M  E  E  M  E        
// 
// 1   2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
