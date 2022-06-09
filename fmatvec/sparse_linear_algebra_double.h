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

#ifndef sparse_linear_algebra_double_h
#define sparse_linear_algebra_double_h

#include "vector.h"

namespace fmatvec {

  /*! \brief Systems of linear equations
   *
   * This function solves systems of sparse linear equations
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LU decompostion.
   * \param A A symmetric sparse matrix.
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, double> slvLU(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Vector<Ref, double> &x);

  /*! \brief System of linear equations
   *
   * This function solves a system of sparse linear equations
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A A symmetric sparse matrix.
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvLU(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &X);

  /*! \brief Eigenvectors and Eigenvalues
   *
   * This function computes a limited number of eigenvectors and eigenvalues of the general eigenvalue equation A*x=(lambda)*B*x.
   * \param A A symmetric sparse matrix.
   * \param B A symmetric sparse, positive definite matrix.
   * \param nev The number of eigenvalues und eigenvectors to be calculated.
   * \param sigma A scalar that determines which eigenvalues are returned.
   * \param eigenvector A matrix, containing the nev eigenvectors et the end of the function.
   * \param eigenvalues A vector, containing the nev eigenvalues closest to sigma at the end of the function.
   * \return void
   * */
  FMATVEC_EXPORT int eigvec(const Matrix<SymmetricSparse, Ref, Ref, double> &A, const Matrix<SymmetricSparse, Ref, Ref, double> &M, int nev, double sigma, Matrix<General, Ref, Ref, double> &eigenvectors, Vector<Ref, double> &eigenvalues);
}

#endif
