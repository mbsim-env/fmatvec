/* Copyright (C) 2003-2005  Martin Förg

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

#ifndef types_h
#define types_h

#ifndef HAVE_LIBMKL_INTEL_LP64
#ifndef CBLAS_ENUM_DEFINED_H
   #define CBLAS_ENUM_DEFINED_H
   enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
   enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
                         AtlasConj=114};
   enum CBLAS_UPLO  {CblasUpper=121, CblasLower=122};
   enum CBLAS_DIAG  {CblasNonUnit=131, CblasUnit=132};
   enum CBLAS_SIDE  {CblasLeft=141, CblasRight=142};
#endif
#define CBLAS_INDEX int

#ifndef ATLAS_ENUM_H
   #define ATLAS_ENUM_H
   #define ATLAS_ORDER CBLAS_ORDER
      #define AtlasRowMajor CblasRowMajor
      #define AtlasColMajor CblasColMajor
   #define ATLAS_TRANS CBLAS_TRANSPOSE
      #define AtlasNoTrans CblasNoTrans
      #define AtlasTrans CblasTrans
      #define AtlasConjTrans CblasConjTrans
   #define ATLAS_UPLO CBLAS_UPLO
      #define AtlasUpper CblasUpper
      #define AtlasLower CblasLower
   #define ATLAS_DIAG CBLAS_DIAG
      #define AtlasNonUnit CblasNonUnit
      #define AtlasUnit CblasUnit
   #define ATLAS_SIDE CBLAS_SIDE
      #define AtlasLeft  CblasLeft
      #define AtlasRight CblasRight
#endif
#else
#include "mkl_cblas.h"
#endif

namespace fmatvec {

  /*! 
   *  \brief Basic shape class for matrices.
   *
   * Class BasicType is the basic shape type.
   * */
  class BasicType {
  };

  /*! 
   *  \brief Shape class for general matrices.
   *
   * Class General is a shape class for general matrices.
   * */
  class General : public BasicType {
  };

  /*! 
   *  \brief Shape class for general band matrices.
   *
   * Class GeneralBand is a shape class for general band matrices.
   * */
  class GeneralBand : public BasicType {
  };

  /*! 
   *  \brief Shape class for symmetric matrices.
   *
   * Class Symmetric is a shape class for symmetric matrices.
   * */
  class Symmetric : public BasicType {
  };

  /*! 
   *  \brief Shape class for diagonal matrices.
   *
   * Class Diagonal is a shape class for diagonal matrices.
   * */
  class Diagonal : public BasicType {
  };

  /*! 
   *  \brief Shape class for sparse matrices.
   *
   * Class Sparse is a shape class for sparse matrices.
   * */
  class Sparse : public BasicType {
  };

  template<int M, int N>
  class FixedGeneral: public BasicType {
  };

  template<int M>
  class FixedSymmetric: public BasicType {
  };

  class VarGeneral: public BasicType {
  };

  template<int M>
  class FixedVarGeneral: public BasicType {
  };

}

#endif
