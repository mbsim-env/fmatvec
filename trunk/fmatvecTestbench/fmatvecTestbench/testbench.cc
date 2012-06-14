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

#include <testbench.h>
#include <config.h>

namespace fmatvec {

  void Testbench::init_Matgenerator() {
  
    // initializing Matgenerator
    matgenerator.setRange(Range);
    matgenerator.setDecimalPlaces(DecimalPlaces);
    matgenerator.init();
  }

  void Testbench::init_MatVec() {
  
    // creating Mat
    GenMat1.resize(Dim,Dim); GenMat1<<matgenerator.getMat(Dim,Dim);
    GenMat2.resize(Dim,Dim); GenMat2<<matgenerator.getMat(Dim,Dim);
    //GenMat3.resize(Dim,Dim); GenMat3<<matgenerator.getMat(Dim,Dim);
    //GenMat4.resize(Dim,Dim); GenMat4<<matgenerator.getMat(Dim,Dim);
    //GenMat5.resize(Dim,Dim); GenMat5<<matgenerator.getMat(Dim,Dim);

    // creating SymMat
    //SymMat1.resize(Dim); SymMat1<<matgenerator.getSymMat(Dim);
    //SymMat2.resize(Dim); SymMat2<<matgenerator.getSymMat(Dim);
    //SymMat3.resize(Dim); SymMat3<<matgenerator.getSymMat(Dim);
    //SymMat4.resize(Dim); SymMat4<<matgenerator.getSymMat(Dim);
    //SymMat5.resize(Dim); SymMat5<<matgenerator.getSymMat(Dim);   
    
    // creating SqrMat
    SqrMat1.resize(Dim); SqrMat1<<matgenerator.getSqrMat(Dim);
    //SqrMat2.resize(Dim); SqrMat2<<matgenerator.getSqrMat(Dim);
    //SqrMat3.resize(Dim); SqrMat3<<matgenerator.getSqrMat(Dim);
    //SqrMat4.resize(Dim); SqrMat4<<matgenerator.getSqrMat(Dim);
    //SqrMat5.resize(Dim); SqrMat5<<matgenerator.getSqrMat(Dim);    
    
    // creating DiagMat
    //DiagMat1.resize(Dim); DiagMat1<<matgenerator.getDiagMat(Dim);
    //DiagMat2.resize(Dim); DiagMat2<<matgenerator.getDiagMat(Dim);
    //DiagMat3.resize(Dim); DiagMat3<<matgenerator.getDiagMat(Dim);
    //DiagMat4.resize(Dim); DiagMat4<<matgenerator.getDiagMat(Dim);
    //DiagMat5.resize(Dim); DiagMat5<<matgenerator.getDiagMat(Dim);    
    
    // creating Vec
    //Vec1.resize(Dim); Vec1<<matgenerator.getVec(Dim);
    //Vec2.resize(Dim); Vec2<<matgenerator.getVec(Dim);
    //Vec3.resize(Dim); Vec3<<matgenerator.getVec(Dim);
    //Vec4.resize(Dim); Vec4<<matgenerator.getVec(Dim);
    //Vec5.resize(Dim); Vec5<<matgenerator.getVec(Dim);    
    
    // creating RowVec
    //RowVec1.resize(Dim); RowVec1<<matgenerator.getRowVec(Dim);
    //RowVec2.resize(Dim); RowVec2<<matgenerator.getRowVec(Dim);
    //RowVec3.resize(Dim); RowVec3<<matgenerator.getRowVec(Dim);
    //RowVec4.resize(Dim); RowVec4<<matgenerator.getRowVec(Dim);
    //RowVec5.resize(Dim); RowVec5<<matgenerator.getRowVec(Dim);

    // creating positiv definit matrix
    //Hposdef1.resize(Dim); Hposdef1<<matgenerator.hposdefSymMat(Dim);

    // creating LL decomposited matrix
    hposdefSymFacLLMat1.resize(Dim); hposdefSymFacLLMat1<<matgenerator.hposdefSymFacLLMat(Dim);

  }


  void Testbench::operate_GenMatGenMat() {
    for(int i=0; i<NumRuns; i++){
      GenMat1=GenMat5-GenMat2*GenMat3-GenMat4*GenMat1-GenMat5*GenMat3;
    }
  }
  
  void Testbench::operate_MultibodyDynamics() {
    for(int i=0; i<NumRuns; i++){
      Vec1 = Vec2 + GenMat1*Vec3 + GenMat2*Vec4;
    }
  }
  
  void Testbench::operate_SymMatSymMat() {
    for(int i=0; i<NumRuns; i++){
      GenMat1=SymMat5-SymMat2*SymMat3-SymMat4*SymMat1-SymMat5*SymMat3;
    }
  }

  void Testbench::operate_SqrMatSqrMat() {
    for(int i=0; i<NumRuns; i++){
      GenMat1=SqrMat5-SqrMat2*SqrMat3-SqrMat4*SqrMat1-SqrMat5*SqrMat3;
    }
  }

  void Testbench::operate_DiagMatDiagMat() {
    for(int i=0; i<NumRuns; i++){
      GenMat1=DiagMat5-DiagMat2*DiagMat3-DiagMat4*DiagMat1-DiagMat5*DiagMat3;
    }
  }

  void Testbench::operate_VecVec() {
    for(int i=0; i<NumRuns; i++){
      GenMat1=Vec2*trans(Vec3)-Vec4*trans(Vec1)-Vec5*trans(Vec3);
    }
  }

  void Testbench::operate_inv() {
    for(int i=0; i<NumRuns; i++){
      GenMat1=inv(SqrMat1)*inv(SqrMat2);
    }
  }

  void Testbench::operate_eigval() {
    for(int i=0; i<NumRuns; i++){
      eigval(SqrMat1);
      eigval(SqrMat2);
    }
  } 

  void Testbench::operate_facLL() {
    for(int i=0; i<NumRuns; i++){
      SymMat1 = facLL(Hposdef1);
    }
  }  
  
  void Testbench::operate_slvLLFac() {
    for(int i=0; i<NumRuns; i++){
      Vec1 = slvLLFac(hposdefSymFacLLMat1,Vec2);
    }
  }
  
  void Testbench::operate_slvLLFac_MatMat() {
    for(int i=0; i<NumRuns; i++){
      slvLLFac(hposdefSymFacLLMat1,GenMat2);
    }
  }
  
  void Testbench::operate_G() {
    for(int i=0; i<NumRuns; i++){
      SqrMat1 = SqrMat(GenMat1.T()*slvLLFac(hposdefSymFacLLMat1,GenMat2));
    }
  }

  void Testbench::operate_slvLU() {
    for(int i=0; i<NumRuns; i++){
      Vec1 = slvLU(SqrMat1,Vec2);
    }
  }
  
  void Testbench::operate_slvLU_MatMat() {
    for(int i=0; i<NumRuns; i++){
      GenMat1 = slvLU(SqrMat1,GenMat2);
    }
  }



  void Testbench::init() {
    Testbench::init_Matgenerator();
    Testbench::init_MatVec();
  }
}
