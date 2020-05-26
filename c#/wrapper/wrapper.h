#ifndef WRAPPER_HPP
#define WRAPPER_HPP


#include <limits>

extern "C"
{
    void __declspec(dllexport) eigval(std::size_t s, double* in, double* out);
}

#endif