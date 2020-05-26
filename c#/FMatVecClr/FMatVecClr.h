#pragma once

#include <complex>

#include "ManagedObject.h"
#include "fmatvec/square_matrix.h"
#include "fmatvec/vector.h"

using namespace System;
using namespace System::Numerics;

namespace FMatVecClr {

	public ref class LinearAlgebraDouble
	{
        public:
            static array<double>^ eigval(int size, array<double>^ in);
	};

    public ref class ComplexVec : public ManagedObject<fmatvec::Vector<fmatvec::Ref, std::complex<double>>>
    {
    public:
        ComplexVec(); // cannot be created from .Net
        ComplexVec(const fmatvec::Vector<fmatvec::Ref, std::complex<double>>& vec);
        Complex^ e(int i);
        property int Size
        {
        public:
            int get()
            {
                return GetInstance()->size();
            }
        }
    private:

    };

    public ref class SquareMatrix : public ManagedObject<fmatvec::SquareMatrix<fmatvec::Ref, double>>
    {
    public:
        SquareMatrix(int size, array<double>^ in);
        array<double>^ eigval();
        ComplexVec^ eigvalComplexVec();
        property int Size
        {
        public:
            int get()
            {
                return GetInstance()->size();
            }
        }
    private:

    };



    


}
