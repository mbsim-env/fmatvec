/* Copyright (C) 2003-2005  Jan Clauberg

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
 *   jan.p.clauberg@googlemail.com
 *
 */

#ifndef TESTBENCH
#define TESTBENCH

#include "fmatvec.h"
#include "fmatvecTestbench/matgenerator.h"
#include <iostream>

using namespace std;

namespace fmatvec {

  class Testbench {
  
    private:
      int Dim;
      int NumRuns;
  
      int Range;
      int DecimalPlaces;

      Matgenerator matgenerator;

      Vec Vec1, Vec2, Vec3, Vec4, Vec5;
      RowVec RowVec1, RowVec2, RowVec3, RowVec4, RowVec5;
      Mat GenMat1, GenMat2, GenMat3, GenMat4, GenMat5;
      SymMat SymMat1, SymMat2, SymMat3, SymMat4, SymMat5, Hposdef1, hposdefSymFacLLMat1;
      SqrMat SqrMat1, SqrMat2, SqrMat3, SqrMat4, SqrMat5;
      DiagMat DiagMat1, DiagMat2, DiagMat3, DiagMat4, DiagMat5;
  
    public:
      /*! \brief Standard constructor
      *
      * Constructs a Matgenerator. 
      * */      
      Testbench() {}
      
      /*! \brief Standard destructor
       *
       * Destructs a Randomgenerator. 
       * */      
      ~Testbench() {}
  
      /*! \brief Setter: number of runs
       *
       * Sets the number of runs for a individual operation
       * \param NumRuns_ numbers of runs.
       * */      
      void setNumRuns(int NumRuns_) {NumRuns=NumRuns_;}
      
      /*! \brief Setter: dimension of operation
       *
       * Sets the dimension of the operations, e.g. Mat(Dim,Dim)*Vec(Dim)
       * \param Dim_ dimension of operation
       * */ 
      void setDim(int Dim_) {Dim=Dim_;}
      
      /*! \brief Setter: range of random numbers
       *
       * Sets the range from ]-Range/2;Range/2[
       * \param Range_ range for random numbers.
       * */
      void setRange(int Range_) {Range=Range_;}
      
      /*! \brief Setter: number of decimal places
       *
       * Sets the number of decimal places
       * \param DecimalPlaces_ number of decimal places.
       * */      
      void setDecimalPlaces(int DecimalPlaces_) {DecimalPlaces=DecimalPlaces_;}
  
      /*! \brief Initializing the Matgenerator
       *
       * */      
      void init_Matgenerator();
      
      /*! \brief Initializing the used matrices and vectors
       *
       * */
      void init_MatVec();
      
      /*! \brief Calls init_Matgenerator() and init_MatVec()
       *
       * */
      void init();

      /*! \brief Operation: general matrix operations
       *
       * Does some operations with general matrices
       * */      
      void operate_GenMatGenMat();

      /*! \brief Operation: operations similar not MultiBody Dynamics
       *
       * Does some operations similar to Multibody Dynamics
       * */  
      void operate_MultibodyDynamics();
      
      /*! \brief Operation: symmetric matrix operations
       *
       * Does some operations with symmetric matrices
       * */      
      void operate_SymMatSymMat();

      /*! \brief Operation: square matrix operations
       *
       * Does some operations with square matrices
       * */       
      void operate_SqrMatSqrMat();

      /*! \brief Operation: diagonal matrix operations
       *
       * Does some operations with diagonal matrices
       * */       
      void operate_DiagMatDiagMat();

      /*! \brief Operation: vector operations
       *
       * Does some operations with vectors
       * */       
      void operate_VecVec();

      /*! \brief Operation: inversion of a matrix 
       *
       * Does inversions of matrices
       * */      
      void operate_inv();

      /*! \brief Operation: calculation of eigenvalues
       *
       * Does the calculation of eigenvalues
       * */      
      void operate_eigval();

      /*! \brief Operation: factorization of a matrix
       *
       * Does the factorization of a matrix
       * */      
      void operate_facLL();

      /*! \brief Operation: solution of a system of linear equations (A*x=b)
       *
       * Does the solution of a system of linear equations (matrix already factorized)
       * */      
      void operate_slvLLFac();
      
      /*! \brief Operation: solution of a system of linear equations (A*X=b)
       *
       * Does the solution of a system of linear equations (matrix already factorized)
       * */
      void operate_slvLLFac_MatMat();
      
      /*! \brief Operation: operation similar to DynamicSystemSolver::updateG
       *
       * Does operations similar to DynamicSystemSolver::update
       * */
      void operate_G();
      
      /*! \brief Operation: solution of a system of linear equations (A*x=b)
       *
       * Does the solution of a system of linear equations (by LU decomposition)
       * */
      void operate_slvLU();
      
      /*! \brief Operation: solution of a system of linear equations (A*X=b)
       *
       * Does the solution of a system of linear equations (by LU decomposition)
       * */
      void operate_slvLU_MatMat();
  };
}

#endif /* TESTBENCH */
