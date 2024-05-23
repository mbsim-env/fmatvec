# This is a python helper module for mathematic code generation.
# Its mainly useful during C++ development.

import sympy
import numpy

# Create a scalar (rows=None, cols=None), vector (cols=None) or matrix symbolic variable with basename name.
# For vector/matrix the element access in in fmatvec notation.
def symbol(name, rows=None, cols=None):
  if rows is None and cols is None:
    return sympy.Symbol(name, real=True)
  if cols is None:
    return numpy.array(list(map(lambda r: sympy.Symbol(name+"("+str(r)+")", real=True), range(0,rows))))
  m=[]
  for r in range(0,rows):
    m.append(list(map(lambda c: sympy.Symbol(name+"("+str(r)+","+str(c)+")", real=True), range(0,cols))))
  return numpy.array(m)

# Generate c-code for all expressions *xx.
# If only one expression is provided use retName as output variable name.
# If more then one expression is provided and retName is a string use "{retName}[{idx}]" as output variable name.
# If more then one expression is provided and retName is a list of string use "{retName[idx]] as output variable name.
# If common subexpressions are extracted (cse=True) then subsName is used as prefix for cse variables.
# If oneOutputPerLine=True matrix expressions a printed one line per element, else matrices are printed as a matrix.
def ccode(*xx, retName="ret", subsName="x", oneOutputPerLine=False, cse=True):
  if cse:
    # convert matrix/vector to sympy
    sympy_xx = []
    for x in xx:
      if type(x)==numpy.ndarray:
        sympy_xx.append(sympy.Matrix(x))
      else:
        sympy_xx.append(x)
    # extract common subexpressions
    (subs,sympy_xx_cse)=sympy.cse(sympy_xx, symbols=sympy.numbered_symbols(subsName))
    # convert sympy_xx_cse back to numpy (while keeping the original scalar/vector/matrix form)
    numpy_xx_cse = []
    for (x,xorg) in zip(sympy_xx_cse, xx):
      if x.__class__.__name__=="MutableDenseMatrix" or x.__class__.__name__=="ImmutableDenseMatrix":
        numpy_xx_cse.append(numpy.reshape(x, xorg.shape))
      else:
        numpy_xx_cse.append(x)
  else:
    # do not extract cse (no subs and numpy_xx_cse=xx)
    subs=[]
    numpy_xx_cse=xx

  def dumpScalar(x):
    x = x.evalf()
    return sympy.ccode(x)

  fullRet=""
  
  # dump common subexpressions
  for x in subs:
    fullRet += f"double {x[0]} = {dumpScalar(x[1])};\n"

  # dump the expressions
  exprIdx=-1
  for x in numpy_xx_cse:
    exprIdx+=1
    # expression output variable name
    if len(xx)>1:
      if type(retName) == list:
        retNameFull=f"{retName[exprIdx]}"
      else:
        retNameFull=f"{retName}[{exprIdx}]"
    else:
      retNameFull=f"{retName}"

    if not hasattr(x, "shape"):
      # dump scalar
      fullRet += f"{retNameFull} = {dumpScalar(x)};\n"
    elif len(x.shape)==1:
      # dump vector
      lr=len(str(x.shape[0]-1))
      ret=""
      for r in range(0, x.shape[0]):
        ret+=f"{retNameFull}({r:{lr}}) = {dumpScalar(x[r])};\n"
      fullRet += ret
    elif len(x.shape)==2:
      # dump matrix
      lr=len(str(x.shape[0]-1))
      lc=len(str(x.shape[1]-1))
      l=numpy.zeros((x.shape[1],))
      mat=[]
      # calculate column size
      for r in range(0, x.shape[0]):
        row=[]
        for c in range(0, x.shape[1]):
          code=dumpScalar(x[r,c])
          row.append(code)
          l[c]=max(l[c], len(code))
        mat.append(row)
      ret=""
      # dump
      for r in range(0, x.shape[0]):
        for c in range(0, x.shape[1]):
          if oneOutputPerLine:
            ret+=f"{retNameFull}({r:{lr}},{c:{lc}}) = {mat[r][c]};\n"
          else:
            ret+=f"{retNameFull}({r:{lr}},{c:{lc}}) = {mat[r][c].ljust(int(l[c]))};"+(" " if c<x.shape[1]-1 else "")
        if not oneOutputPerLine:
          ret+="\n"
      fullRet += ret
  return fullRet

# Convert from numpy to sympy scalar/vector/matrix.
# Note that a sympy matrix of size Nx1 is converted to a numpy vector of shape (N,) if asVectorIfPossible=True
# (there is no sympy vector)
def C(x, asVectorIfPossible=True):
  if type(x)==numpy.ndarray:
    # convert numpy array to sympy matrix
    return sympy.Matrix(x)
  if x.__class__.__name__=="MutableDenseMatrix" or x.__class__.__name__=="ImmutableDenseMatrix" or x.__class__.__name__=="ImmutableDenseNDimArray":
    # convert sympy matrix to numpy array (as 1D or 2D)
    if asVectorIfPossible and x.shape[1]==1:
      return numpy.reshape(x, (x.shape[0],))
    return numpy.array(x)
  # no conversion needed for scalars
  return x
