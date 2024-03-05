# This is a python helper module for mathematic code generation.
# Its mainly useful during C++ development.

import sympy
import numpy

def symbol(name, rows=None, cols=None):
  if rows is None and cols is None:
    return sympy.Symbol(name)
  if cols is None:
    return numpy.array(list(map(lambda r: sympy.Symbol(name+"("+str(r)+")"), range(0,rows))))
  m=[]
  for r in range(0,rows):
    m.append(list(map(lambda c: sympy.Symbol(name+"("+str(r)+","+str(c)+")"), range(0,cols))))
  return numpy.array(m)

def ccode(x, retName="ret", oneOutputPerLine=False):
  def dumpScalar(x):
    x = sympy.trigsimp(x)
    x = sympy.simplify(x)
    x = x.evalf()
    return sympy.ccode(x)
  if not hasattr(x, "shape"):
    return retName + " = " + dumpScalar(x)+";"
  if len(x.shape)==1:
    lr=len(str(x.shape[0]-1))
    ret="{\n"
    for r in range(0, x.shape[0]):
      ret+=f"{retName}({r:{lr}}) = "+dumpScalar(x[r])+";\n"
    ret+="}"
    return ret
  if len(x.shape)==2:
    lr=len(str(x.shape[0]-1))
    lc=len(str(x.shape[1]-1))
    l=numpy.zeros((x.shape[1],))
    mat=[]
    for r in range(0, x.shape[0]):
      row=[]
      for c in range(0, x.shape[1]):
        code=dumpScalar(x[r,c])
        row.append(code)
        l[c]=max(l[c], len(code))
      mat.append(row)
    ret="{\n"
    for r in range(0, x.shape[0]):
      for c in range(0, x.shape[1]):
        if oneOutputPerLine:
          ret+=f"{retName}({r:{lr}},{c:{lc}}) = {mat[r][c]};\n"
        else:
          ret+=f"{retName}({r:{lr}},{c:{lc}}) = {mat[r][c].ljust(int(l[c]))};"+(" " if c<x.shape[1]-1 else "")
      if not oneOutputPerLine:
        ret+="\n"
    ret+="}"
    return ret

def __getattr__(name):
  def wrap(x, *args,  **kwargs):
    return numpy.reshape(getattr(sympy, name)(sympy.Matrix(x), *args, **kwargs), x.shape)
  return wrap
