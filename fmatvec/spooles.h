/* Copyright (C) 2003-2022  Martin Förg

 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact:
 *   martin.o.foerg@googlemail.com
 *
 */

#ifndef spooles_h
#define spooles_h

#include "sparse_matrix.h"
#include "symmetric_sparse_matrix.h"

typedef struct _IV   IV;
typedef struct _InpMtx InpMtx;
typedef struct _FrontMtx FrontMtx;
typedef struct _SubMtxManager SubMtxManager;
typedef struct _ETree ETree;

namespace fmatvec {

  class Spooles {
    public:
      Spooles(const Matrix<SymmetricSparse, Ref, Ref, double> &A);
      Spooles(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Matrix<SymmetricSparse, Ref, Ref, double> &M, double sigma);
      virtual ~Spooles();
      void factorize();
      Matrix<General, Ref, Ref, double> solve(const Matrix<General, Ref, Ref, double> &B);
      Vector<Ref, double> solve(const Vector<Ref, double> &b);
    private:
      int neq{0};
      InpMtx *mtxA{nullptr};
      double cpus[11];
      IV *newToOldIV{nullptr};
      IV *oldToNewIV{nullptr};
      FrontMtx *frontmtx{nullptr};
      SubMtxManager *mtxmanager{nullptr};
      ETree *frontETree{nullptr};
  };

}

#endif
