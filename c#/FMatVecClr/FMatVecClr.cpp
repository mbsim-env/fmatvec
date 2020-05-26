#include "FMatVecClr.h"
#include <iostream>

#include "fmatvec/linear_algebra_double.h"

namespace FMatVecClr {
    
    
    // *************************************************************************************
    // Static function exchanging CLR arrays
    // *************************************************************************************

    array<double>^ LinearAlgebraDouble::eigval(int size, array<double>^ in)
    {
        // a more elegant way may exist
        double* inC{ new double[size] };
        for (int i = 0; i < (size*size); ++i) 
            inC[i] = in[i];

        const fmatvec::SquareMatrix<fmatvec::Ref, double> A(size, inC);
        fmatvec::Vector<fmatvec::Ref, std::complex<double>> eigenValues = fmatvec::eigval(A);
        delete[] inC;

        array< double >^ result = gcnew array< double >(size);
        int cnt = 0;
        for (const std::complex<double> ev : eigenValues)
        {
            result[cnt] = ev.real();
            ++cnt;
        }
        return result;
    }

    // *************************************************************************************
    // Class SquareMatrix
    // *************************************************************************************

    SquareMatrix::SquareMatrix(int size, array<double>^ in) 
        : ManagedObject< fmatvec::SquareMatrix<fmatvec::Ref, double> >(new fmatvec::SquareMatrix<fmatvec::Ref, double>(size))
    {
        // a more elegant way may exist
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                GetInstance()->e(i, j) = in[i * size + j];
        
    }
    array<double>^ SquareMatrix::eigval()
    {
        fmatvec::Vector<fmatvec::Ref, std::complex<double>> eigenValues = fmatvec::eigval(*GetInstance());

        array< double >^ result = gcnew array< double >(GetInstance()->size());
        int cnt = 0;
        for (const std::complex<double> ev : eigenValues)
        {
            result[cnt] = ev.real();
            ++cnt;
        }
        return result;
    }

    ComplexVec^ SquareMatrix::eigvalComplexVec()
    {
        fmatvec::Vector<fmatvec::Ref, std::complex<double>> eigenValues = fmatvec::eigval(*GetInstance());
        eigenValues.size();
        return gcnew ComplexVec(eigenValues);
    }


    // *************************************************************************************
    // Class ComplexVec
    // *************************************************************************************
    ComplexVec::ComplexVec()
        : ManagedObject<fmatvec::Vector<fmatvec::Ref, std::complex<double>>>(new fmatvec::Vector<fmatvec::Ref, std::complex<double>>())
    {}

    ComplexVec::ComplexVec(const fmatvec::Vector<fmatvec::Ref, std::complex<double>>& vec)
        : ManagedObject<fmatvec::Vector<fmatvec::Ref, std::complex<double>>>(new fmatvec::Vector<fmatvec::Ref, std::complex<double>>(vec))
    {}

    Complex^ ComplexVec::e(int i)
    {
        return gcnew Complex( GetInstance()->e(i).real(), GetInstance()->e(i).imag());
    }
}

