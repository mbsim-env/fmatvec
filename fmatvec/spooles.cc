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

#include "config.h"
#include "spooles.h"
#include "vector.h"
extern "C" {
#include <FrontMtx.h>
#include <misc.h>
#include <SymbFac.h>
}

namespace fmatvec {

  Spooles::Spooles(const Matrix<SymmetricSparse, Ref, Ref, double> &A) : neq(A.rows()) {
    mtxA = InpMtx_new();
    InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, A.nonZeroElements(), neq);
    for(int i=0; i<neq; i++) {
      for (int j = A.Ip()[i]; j < A.Ip()[i+1]; j++)
	InpMtx_inputRealEntry(mtxA, i, A.Jp()[j], A()[j]);
    }
    InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
  }

  Spooles::Spooles(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Matrix<SymmetricSparse, Ref, Ref, double> &M, double sigma) : neq(A.rows()) {
    mtxA = InpMtx_new();
    InpMtx_init(mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, A.nonZeroElements(), neq);
    for(int i=0; i<neq; i++) {
      for (int j = A.Ip()[i]; j < A.Ip()[i+1]; j++)
	InpMtx_inputRealEntry(mtxA, i, A.Jp()[j], A()[j]+sigma*M()[j]);
    }
    InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
  }

  Spooles::~Spooles() {
    FrontMtx_free(frontmtx);
    IV_free(newToOldIV);
    IV_free(oldToNewIV);
    SubMtxManager_free(mtxmanager);
    ETree_free(frontETree);
  }

  void Spooles::factorize() {
    DVfill(10, cpus, 0.0);

    Graph *graph = Graph_new();
    IVL *adjIVL = InpMtx_fullAdjacency(mtxA);
    int nedges = IVL_tsize(adjIVL);
    Graph_init2(graph, 0, neq, 0, nedges, neq, nedges, adjIVL, nullptr, nullptr);

    frontETree = orderViaBestOfNDandMS(graph, 800, 1000, 64, 7892713, 0, nullptr);

    oldToNewIV = ETree_oldToNewVtxPerm(frontETree);
    int *oldToNew = IV_entries(oldToNewIV);
    newToOldIV = ETree_newToOldVtxPerm(frontETree);
    ETree_permuteVertices(frontETree, oldToNewIV);
    InpMtx_permute(mtxA, oldToNew, oldToNew);
    InpMtx_mapToUpperTriangle(mtxA);
    InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS);
    InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
    IVL *symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA);

    frontmtx = FrontMtx_new();
    mtxmanager = SubMtxManager_new();
    SubMtxManager_init(mtxmanager, NO_LOCK, 0);
    FrontMtx_init(frontmtx, frontETree, symbfacIVL, SPOOLES_REAL, 0, FRONTMTX_DENSE_FRONTS, SPOOLES_PIVOTING, NO_LOCK, 0, nullptr, mtxmanager, 0, nullptr);

    int stats[20];
    int error;

    ChvManager *chvmanager = ChvManager_new();
    ChvManager_init(chvmanager, NO_LOCK, 1);
    IVfill(20, stats, 0);
    Chv *rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, 100, 0, chvmanager, &error, cpus, stats, 0, nullptr);
    ChvManager_free(chvmanager);
    if(rootchv != nullptr)
      throw std::runtime_error("Singular matrix");
    if(error >= 0)
      throw std::runtime_error("Error in front");

    FrontMtx_postProcess(frontmtx, 0, nullptr);
    IVL_free(symbfacIVL);
    InpMtx_free(mtxA);
    Graph_free(graph);
  }

  Vector<Ref, double> Spooles::solve(const Vector<Ref, double> &b) {
    DenseMtx *mtxB = DenseMtx_new();
    DenseMtx_init(mtxB, SPOOLES_REAL, 0, 0, neq, 1, 1, neq);
    DenseMtx_zero(mtxB);
    for(int i = 0; i < neq; i++)
      DenseMtx_setRealEntry(mtxB, i, 0, b(i));

    DenseMtx_permuteRows(mtxB, oldToNewIV);
    DenseMtx *mtxX = DenseMtx_new();
    DenseMtx_init(mtxX, SPOOLES_REAL, 0, 0, neq, 1, 1, neq);
    DenseMtx_zero(mtxX);
    FrontMtx_solve(frontmtx, mtxX, mtxB, mtxmanager, cpus, 0, nullptr);
    DenseMtx_permuteRows(mtxX, newToOldIV);
    DenseMtx_free(mtxB);

    Vector<Ref, double> x(b.size(),NONINIT);
    for(int i = 0; i < neq; i++)
      x(i) = DenseMtx_entries(mtxX)[i];
    DenseMtx_free(mtxX);

    return x;
  }

  Matrix<General, Ref, Ref, double> Spooles::solve(const Matrix<General, Ref, Ref, double> &B) {
    if(B.rows() != neq)
      throw std::runtime_error("Rows of B must be neq");
    DenseMtx *mtxB = DenseMtx_new();
    DenseMtx_init(mtxB, SPOOLES_REAL, 0, 0, B.rows(), B.cols(), 1, B.ldim());
    DenseMtx_zero(mtxB);
    for(int i=0; i<B.rows(); i++) {
      for(int j=0; j<B.cols(); j++)
	DenseMtx_setRealEntry(mtxB, i, j, B(i,j));
    }

    DenseMtx_permuteRows(mtxB, oldToNewIV);
    DenseMtx *mtxX = DenseMtx_new();
    DenseMtx_init(mtxX, SPOOLES_REAL, 0, 0, B.rows(), B.cols(), 1, B.ldim());
    DenseMtx_zero(mtxX);
    FrontMtx_solve(frontmtx, mtxX, mtxB, mtxmanager, cpus, 0, nullptr);
    DenseMtx_permuteRows(mtxX, newToOldIV);
    DenseMtx_free(mtxB);

    Matrix<General, Ref, Ref, double> X(B.rows(),B.cols(),NONINIT);
    for(int i=0; i<B.rows(); i++) {
      for(int j=0; j<B.cols(); j++)
	X(i,j) = DenseMtx_entries(mtxX)[i+j*neq];
    }
    DenseMtx_free(mtxX);

    return X;
  }

}
