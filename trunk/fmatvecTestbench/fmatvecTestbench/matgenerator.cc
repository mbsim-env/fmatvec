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

#include <matgenerator.h>
#include <config.h>
#include <iostream>
#include "math.h"

namespace fmatvec {

  void Matgenerator::init() {
    Randomgenerator::init();
  }
 
  Mat Matgenerator::getMat(int Rows_, int Cols_, bool norm_) {
    Mat A(Rows_,Cols_,INIT,0.0);

    for(int i=0; i<Rows_; i++) {
      for(int j=0; j<Cols_; j++) {
        if(norm_) A(i,j)=(*this)(true);
        else A(i,j)=(*this)();
      }
    }
  return A;
  }
  
  SymMat Matgenerator::getSymMat(int Rows_, bool norm_) {
    SymMat A(Rows_,INIT,0.0);

    for(int i=0; i<Rows_; i++) {
      for(int j=0; j<=i; j++) {
        if(norm_) A(i,j)=(*this)(true);
        else A(i,j)=(*this)();
      }
    }
  return A;
  }
  
  SqrMat Matgenerator::getSqrMat(int Rows_, bool norm_) {
    SqrMat A(Rows_,INIT,0.0);

    for(int i=0; i<Rows_; i++) {
      for(int j=0; j<Rows_; j++) {
        if(norm_) A(i,j)=(*this)(true);
        else A(i,j)=(*this)();        
      }
    }
  return A;
  }
  
  DiagMat Matgenerator::getDiagMat(int Rows_, bool norm_) {
    DiagMat A(Rows_,INIT,0.0);

    for(int i=0; i<Rows_; i++) {
        if(norm_) A(i)=(*this)(true);
        else A(i)=(*this)();
      }
  return A;
  }
  
  Vec Matgenerator::getVec(int Rows_, bool norm_) {
    Vec a(Rows_,INIT,0.0);

    for(int i=0; i<Rows_; i++) {
        if(norm_) a(i)=(*this)(true);
        else a(i)=(*this)();
      }
  return a;  
  }
  
  RowVec Matgenerator::getRowVec(int Cols_, bool norm_) {
    RowVec a(Cols_,INIT,0.0);

    for(int i=0; i<Cols_; i++) {
        if(norm_) a(i)=(*this)(true);
        else a(i)=(*this)();
      }
  return a;  
  }

  SymMat Matgenerator::hposdefSymMat(int Rows_) {
    Mat tmp1, tmp2;
    SymMat hposdefMat(Rows_,INIT,0.0);
    
    tmp1=Matgenerator::getMat(Rows_,Rows_,false);
    tmp2=trans(tmp1)*tmp1;
    bool notposdef=false;

    do {
      for (int i=0; i<tmp2.rows();i++) {
        for (int j=0; j<tmp2.cols();j++) {
          if (j>=i){
          hposdefMat(i,j)=tmp2(i,j);
          }  
        }
      }
    
      Vec eig=eigval(hposdefMat);
      for(int i=0; i<hposdefMat.rows(); i++) {
        if(eig(i)<=0) notposdef=true;
      }
    } while (notposdef);
  
  return hposdefMat;
  }

  SymMat Matgenerator::hposdefSymFacLLMat(int Rows_) {
    SymMat A(Rows_,INIT,0.0);
    A=Matgenerator::hposdefSymMat(Rows_);
  return facLL(A);
  }

}
