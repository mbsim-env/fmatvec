import sympy
import numpy

def dump(x):
  def s(x):
    return str(round(x, 12))
  xx=numpy.array(x).astype(numpy.float64)
  if xx.shape==(3,1):
    print("["+s(xx[0,0])+"; "+s(xx[1,0])+"; "+s(xx[2,0])+"]")
  elif xx.shape==(1,3):
    print("["+s(xx[0,0])+", "+s(xx[0,1])+", "+s(xx[0,2])+"]")
  elif xx.shape==(3,3):
    print("["+s(xx[0][0])+", "+s(xx[0][1])+", "+s(xx[0][2])+"; "+\
              s(xx[1][0])+", "+s(xx[1][1])+", "+s(xx[1][2])+"; "+\
              s(xx[2][0])+", "+s(xx[2][1])+", "+s(xx[2][2])+"]")
  else:
    print(round(x,12))

def diffR(R, x):
  if type(x)==sympy.Symbol:
    A=sympy.diff(R, x)*R.T
    return sympy.Matrix([A[2,1],A[0,2],A[1,0]])
  else:
    r0=diffR(R, x[0])
    r1=diffR(R, x[1])
    r2=diffR(R, x[2])
    return r2.col_insert(0, r1).col_insert(0, r0)

def diffV(y, x):
  matrix=False
  try:
    y.shape
    matrix=True
  except:
    pass
  if not matrix:
    r0=sympy.diff(y, x[0])
    r1=sympy.diff(y, x[1])
    r2=sympy.diff(y, x[2])
    return sympy.Matrix([[r0,r1,r2]])
  else:
    r0=sympy.diff(y, x[0])
    r1=sympy.diff(y, x[1])
    r2=sympy.diff(y, x[2])
    return r2.col_insert(0, r1).col_insert(0, r0)

def diffVdir(y, x, xd):
  return sympy.diff(y, x[0])*xd[0]+sympy.diff(y, x[1])*xd[1]+sympy.diff(y, x[2])*xd[2]

def scalar():
  x=sympy.Symbol('x')
  v=sympy.Matrix([x*x, x*x, x*x])
  M=sympy.Matrix([
    [x*x, x*x, x*x],
    [x*x, x*x, x*x],
    [x*x, x*x, x*x]
  ])
  R=sympy.Matrix([
    [x*x, x*x, x*x],
    [x*x, x*x, x*x],
    [x*x, x*x, x*x]
  ])
  funcS=x*x
  funcV=v*x
  funcRV=(v.T*x)
  funcM=M*x
  funcR=R*x

  arg=1.2
  argd=5.6
  argd2=7.8

  dump(funcS.subs(x, arg))
  dump(funcV.subs(x, arg))
  dump(funcRV.subs(x, arg))
  dump(funcM.subs(x, arg))
  dump(funcR.subs(x, arg))

  dump(sympy.diff(funcS, x).subs(x, arg))
  dump(sympy.diff(funcV, x).subs(x, arg))
  dump(sympy.diff(funcRV, x).subs(x, arg))
  dump(sympy.diff(funcM, x).subs(x, arg))
  dump(diffR(funcR, x).subs(x, arg))

  dump((sympy.diff(funcS, x)*argd).subs(x, arg))
  dump((sympy.diff(funcV, x)*argd).subs(x, arg))
  dump((sympy.diff(funcRV, x)*argd).subs(x, arg))
  dump((sympy.diff(funcM, x)*argd).subs(x, arg))
  dump((diffR(funcR, x)*argd).subs(x, arg))

  dump(sympy.diff(sympy.diff(funcS, x), x).subs(x, arg))
  dump(sympy.diff(sympy.diff(funcV, x), x).subs(x, arg))
  dump(sympy.diff(sympy.diff(funcRV, x), x).subs(x, arg))
  dump(sympy.diff(sympy.diff(funcM, x), x).subs(x, arg))
  dump(sympy.diff(diffR(funcR, x), x).subs(x, arg))

  dump((sympy.diff(sympy.diff(funcS, x), x)*argd).subs(x, arg))
  dump((sympy.diff(sympy.diff(funcV, x), x)*argd).subs(x, arg))
  dump((sympy.diff(sympy.diff(funcRV, x), x)*argd).subs(x, arg))
  dump((sympy.diff(sympy.diff(funcM, x), x)*argd).subs(x, arg))
  dump((sympy.diff(diffR(funcR, x), x)*argd).subs(x, arg))

  dump((sympy.diff((sympy.diff(funcS, x)*argd), x)*argd2).subs(x, arg))
  dump((sympy.diff((sympy.diff(funcV, x)*argd), x)*argd2).subs(x, arg))
  dump((sympy.diff((sympy.diff(funcRV, x)*argd), x)*argd2).subs(x, arg))
  dump((sympy.diff((sympy.diff(funcM, x)*argd), x)*argd2).subs(x, arg))
  dump((sympy.diff((diffR(funcR, x)*argd), x)*argd2).subs(x, arg))

scalar()



def vector():
  x0=sympy.Symbol('x0')
  x1=sympy.Symbol('x1')
  x2=sympy.Symbol('x2')
  x=sympy.Matrix([x0, x1, x2])
  v=sympy.Matrix([x[0]*x[1], x[0]*x[1], x[0]*x[1]])
  M=sympy.Matrix([
    [x[0]*x[1], x[0]*x[1], x[0]*x[1]],
    [x[0]*x[1], x[0]*x[1], x[0]*x[1]],
    [x[0]*x[1], x[0]*x[1], x[0]*x[1]],
  ])
  R=sympy.Matrix([
    [x[0]*x[1], x[0]*x[1], x[0]*x[1]],
    [x[0]*x[1], x[0]*x[1], x[0]*x[1]],
    [x[0]*x[1], x[0]*x[1], x[0]*x[1]],
  ])
  funcS=(x.T*v)[0]
  funcV=M*v
  funcRV=(M*v).T
  funcM=x*v.T
  funcR=R*x[0]*x[1]*x[2]

  arg=sympy.Matrix([1.2, 2.3, 3.4]);
  argd=sympy.Matrix([5.6, 6.7, 7.8]);
  argd2=sympy.Matrix([7.8, 8.9, 9.1]);

  subst=((x0, arg[0]), (x1, arg[1]), (x2, arg[2]))

  dump(funcS.subs(subst))
  dump(funcV.subs(subst))
  dump(funcRV.subs(subst))
  dump(funcM.subs(subst))
  dump(funcR.subs(subst))

  dump(diffV(funcS, x).subs(subst))
  dump(diffV(funcV, x).subs(subst))
  #dump(diffV(funcRV, x).subs(subst))
  #dump(diffV(funcM, x).subs(subst))
  dump(diffR(funcR, x).subs(subst))

  dump(diffVdir(funcS, x,argd).subs(subst))
  dump(diffVdir(funcV, x,argd).subs(subst))
  dump(diffVdir(funcRV, x,argd).subs(subst))
  dump(diffVdir(funcM, x,argd).subs(subst))
  dump((diffR(funcR, x)*argd).subs(subst))

  #dump(diffV(diffV(funcS, x), x).subs(subst))
  #dump(diffV(diffV(funcV, x), x).subs(subst))
  #dump(diffV(diffV(funcRV, x), x).subs(subst))
  #dump(diffV(diffV(funcM, x), x).subs(subst))
  #dump(diffV(diffR(funcR, x), x).subs(subst))

  dump(diffVdir(diffV(funcS, x), x,argd).subs(subst))
  dump(diffVdir(diffV(funcV, x), x,argd).subs(subst))
  #dump(diffVdir(diffV(funcRV, x), x,argd).subs(subst))
  #dump(diffVdir(diffV(funcM, x), x,argd).subs(subst))
  dump(diffVdir(diffR(funcR, x), x,argd).subs(subst))

  dump(diffVdir(diffVdir(funcS, x,argd), x,argd2).subs(subst))
  dump(diffVdir(diffVdir(funcV, x,argd), x,argd2).subs(subst))
  dump(diffVdir(diffVdir(funcRV, x,argd), x,argd2).subs(subst))
  dump(diffVdir(diffVdir(funcM, x,argd), x,argd2).subs(subst))
  dump(diffVdir((diffR(funcR, x)*argd), x,argd2).subs(subst))

vector()
