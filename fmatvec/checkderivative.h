#ifndef _FMATVEC_CHECKDERIVATIVE_H_
#define _FMATVEC_CHECKDERIVATIVE_H_

#include <map>
#include <functional>
#include <string>
#include <iostream>
#include <fmatvec/linear_algebra.h>

namespace fmatvec {

  #define FMATVEC_CHECKDERIVATIVE_MEMBERFUNC(func, anaDiff, value, /*indep,*/ ...) \
    fmatvec::checkDerivative(this, std::string(__FILE__)+":"+std::to_string(__LINE__)+":"+#anaDiff, func, anaDiff, value, /*indep,*/ __VA_ARGS__)
  #define FMATVEC_CHECKDERIVATIVE_GLOBALFUNC(func, anaDiff, value, /*indep,*/ ...) \
    fmatvec::checkDerivative(nullptr, std::string(__FILE__)+":"+std::to_string(__LINE__)+":"+#anaDiff, func, anaDiff, value, /*indep,*/ __VA_ARGS__)

  // Check a analytical derivative against the finite difference.
  // Any local variable for which a analytical derivative is given can be checked.
  // - idPtr and idName must, together, be a unique identifier for this call. Usually idPtr is "this" and idName is the name of the variable to be checked.
  // - func is the callers function with the independent variable as input.
  // - anaDiff is the analytical derivative which should be checked.
  // - value is the value for which the analytical derivative should be checked.
  // - indep is the independent variable with respect to which the derivative is calculated.
  // - eps is the tolerance. If anaDiff differs more than this value a error is printed.
  // - delta is the finite difference for the finite derivative calculation.
  template<class Value>
  void checkDerivative(const void *idPtr, const std::string &idName, const std::function<void(double)> &func, const Value &anaDiff, const Value &value, double indep, double eps=1e-6, double delta=1e-8) {
    #ifdef NDEBUG
      #error "fmatvec::checkDerivative should not be active for release builds"
    #endif

    static int distributedCount = 0;
    static std::map<std::pair<const void*, std::string>, std::pair<bool, Value>> ele;

    auto &[disturbed, disturbedValue] = ele.emplace(std::make_pair(idPtr, idName), std::make_pair(false, Value())).first->second;

    if(disturbed) { // this is a distributed call (for the key in "ele") ..
      // ... save the distributed value (in "ele" map)
      disturbedValue = value;
    }
    else if(distributedCount == 0) { // this is a normal call (not distributed) ...
      std::array<Value, 2> disturbedValueRightLeft;

      // .. make a distributed call for the current key in "ele"
      disturbed = true;
      distributedCount++;

      func(indep + delta); // right disturbed
      // save the distributed value (scalar or fmatvec vec/mat)
      if constexpr (std::is_same_v<Value, double>)
        disturbedValueRightLeft[0] = disturbedValue;
      else
        disturbedValueRightLeft[0] <<= disturbedValue;

      func(indep - delta); // left disturbed
      // save the distributed value (scalar or fmatvec vec/mat)
      if constexpr (std::is_same_v<Value, double>)
        disturbedValueRightLeft[1] = disturbedValue;
      else
        disturbedValueRightLeft[1] <<= disturbedValue;

      func(indep); // undisturbed again to reset everything

      disturbed = false;
      distributedCount--;

      // calculate the right and left finite difference and its norm
      std::array<double, 2> dist;
      std::array<Value, 2> finiteDiffRightLeft;
      for(int i=0; i<2; ++i) {
        finiteDiffRightLeft[i] = (i==0?1.0:-1.0) * (disturbedValueRightLeft[i] - value)/delta;

        if constexpr (std::is_same_v<Value, double>)
          dist[i] = std::abs(finiteDiffRightLeft[i] - anaDiff);
        else if constexpr (Value::isVector)
          dist[i] = fmatvec::nrmInf(finiteDiffRightLeft[i] - anaDiff);
        else {
          dist[i] = 0;
          for(int c=0; c<anaDiff.cols(); c++)
            dist[i] = std::max(dist[i], fmatvec::nrm2(finiteDiffRightLeft[i].col(c) - anaDiff.col(c)));
        }
      }

      // if both, the right and the left finite difference, differs from analytically derivative value print a error message
      if(dist[0] > eps && dist[1] > eps) {
        std::stringstream str;
        str<<"checkDerivative failed in "<<idPtr<<": ID = "<<idName<<":"<<std::endl
           <<"- anaDiff         = "<<anaDiff<<std::endl
           <<"- finiteDiffRight = "<<finiteDiffRightLeft[0]<<std::endl
           <<"- finiteDiffLeft  = "<<finiteDiffRightLeft[1]<<std::endl
           <<"- indep = "<<indep<<std::endl
           <<"- nrmInf (right/left) = "<<dist[0]<<" "<<dist[1];
        //throw std::runtime_error(str.str());
        std::cerr<<str.str()<<std::endl;
        assert(0);
        throw std::runtime_error(str.str());
      }
    }
  }

}

#endif
