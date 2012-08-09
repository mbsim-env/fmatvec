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

#ifndef MATGENERATOR
#define MATGENERATOR

#include "fmatvec.h"
#include "fmatvecTestbench/randomgenerator.h"

using namespace std;

namespace fmatvec {

  class Matgenerator : public Randomgenerator {
  
    private:

    public:
      /*! \brief Standard constructor
      *
      * Constructs a Matgenerator. 
      * */
      Matgenerator() {}
      
      /*! \brief Standard destructor
       *
       * Destructs a Randomgenerator. 
       * */
      ~Matgenerator() {}
      
      /*! \brief Initializing the RandomGenerator
       *
       * */ 
      void init();

      /*! \brief Generates a random general matrix
       *
       * Generates a general matrix
       * \param Rows_ number of rows.
       * \param Cols_ number of cols.
       * \param norm_ only positiv entries
       * \return A general Matrix
       * */ 
      Mat getMat(int Rows_, int Cols_, bool norm_=false);

      /*! \brief Generates a random symmetric matrix
       *
       * Generates a symmetric matrix
       * \param Rows_ number of rows/cols.
       * \param norm_ only positiv entries
       * \return A symmetric Matrix
       * */       
      SymMat getSymMat(int Rows_, bool norm_=false);

      /*! \brief Generates a random square matrix
       *
       * Generates a square matrix
       * \param Rows_ number of rows/cols.
       * \param norm_ only positiv entries
       * \return A square Matrix
       * */       
      SqrMat getSqrMat(int Rows_, bool norm_=false);

      /*! \brief Generates a random diagonal matrix
       *
       * Generates a diagonal matrix
       * \param Rows_ number of rows/cols.
       * \param norm_ only positiv entries
       * \return A diagonal Matrix
       * */       
      DiagMat getDiagMat(int Rows_, bool norm_=false);

      /*! \brief Generates a random vector
       *
       * Generates a vector
       * \param Rows_ number of Rows.
       * \param norm_ only positiv entries
       * \return A vector
       * */       
      Vec getVec(int Rows_, bool norm_=false);

      /*! \brief Generates a random rowvector
       *
       * Generates a rowvector
       * \param Cols_ number of cols.
       * \param norm_ only positiv entries
       * \return A general rowvector
       * */       
      RowVec getRowVec(int Cols_, bool norm_=false);

      /*! \brief Generates a positiv definit symmetric random matrix
       *
       * Generates a positiv definit symmetric matrix
       * \param Rows_ number of rows/cols.
       * \return A positiv definit symmetric matrix
       * */       
      SymMat hposdefSymMat(int Rows_);
      
      /*! \brief Generates a decomposited (Cholesky) symmetric random matrix
       *
       * Generates a decomposited (Cholesky) symmetric matrix
       * \param Rows_ number of rows/cols.
       * \return A decomposited (Cholesky) symmetric matrix
       * */
      SymMat hposdefSymFacLLMat(int Rows_);
  };
}

#endif /* MATGENERATOR */
