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

#include "stream_impl.h"
#include "general_matrix.h"
#include "var_general_matrix.h"
#include "fixed_var_general_matrix.h"
#include "var_fixed_general_matrix.h"
#include "fixed_general_matrix.h"
#include "symmetric_matrix.h"
#include "var_symmetric_matrix.h"
#include "fixed_symmetric_matrix.h"
#include "diagonal_matrix.h"
#include "ast.h"

namespace fmatvec {

  template<>
  boost::spirit::qi::rule<boost::spirit::istream_iterator, double()>& getBoostSpiritQiRule<double>() {
    static boost::spirit::qi::rule<boost::spirit::istream_iterator, double()> ret=boost::spirit::qi::double_;
    return ret;
  }

  template<>
  boost::spirit::qi::rule<boost::spirit::istream_iterator, int()>& getBoostSpiritQiRule<int>() {
    static boost::spirit::qi::rule<boost::spirit::istream_iterator, int()> ret=boost::spirit::qi::int_;
    return ret;
  }
  template<>
  boost::spirit::qi::rule<boost::spirit::istream_iterator, long()>& getBoostSpiritQiRule<long>() {
    static boost::spirit::qi::rule<boost::spirit::istream_iterator, long()> ret=boost::spirit::qi::long_;
    return ret;
  }

  template<>
  boost::spirit::qi::rule<boost::spirit::istream_iterator, std::complex<double>()>& getBoostSpiritQiRule<std::complex<double>>() {
    static boost::spirit::qi::rule<boost::spirit::istream_iterator, std::complex<double>()> ret =
      '(' >> boost::spirit::qi::double_[boost::spirit::qi::_val=boost::spirit::qi::_1] >>
      '+' >> boost::spirit::qi::double_[boost::spirit::qi::_val=boost::spirit::qi::_val+boost::spirit::qi::_1*std::complex<double>(0,1)] >> "*i)";
    return ret;
  }

  template<>
  boost::spirit::karma::rule<std::ostream_iterator<char>, double()>& getBoostSpiritKarmaRule<double>() {
    struct DoublePolicy : boost::spirit::karma::real_policies<double> {
      static constexpr int floatfield(double n) { return fmtflags::scientific; }
      static constexpr unsigned precision(double n) { return std::numeric_limits<double>::digits10+1; }
      // BEGIN: added bugfix from https://github.com/boostorg/spirit/issues/529
      // computes the number of digits left of the decimal in base 10
      static unsigned int intDigits(double n) {
          if (n < 0)
              n = std::abs(n);
  
          // need to do this in 64-bit to be able to represent large enough numbers
          const uint64_t limit = UINT64_MAX / 10;
          uint64_t num = floor(n);
          for (uint64_t x = 10, i = 1;; x *= 10, i++) {
              if (num < x)
                  return i;
              if (x > limit)
                  return i + 1;
          }
      }
  
      // Output the fractional part of the number
      //
      // By default, if the fractional part is zero we do nothing, unless the policy
      // was explicitly told to do so
      #if !defined(BOOST_MPL_CFG_NO_ADL_BARRIER_NAMESPACE)
        #define BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE mpl_
      #else
        #define BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE boost::mpl
      #endif
      bool fraction_part(boost::spirit::karma::detail::output_iterator<std::ostream_iterator<char>, BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::int_<15>, boost::spirit::unused_type>& sink, double n, unsigned int precision_, unsigned int precision) const
      {
          // Inline code for karma::real_policies<double>::fraction_part() and change it to use
          // a method to calculate number of digits for n that doesn't overflow when n is very large.
          unsigned int digits = boost::spirit::traits::test_zero(n) ? 1 : intDigits(n);
          
          bool r = true;
          for (/**/; r && digits < precision_; digits = digits + 1)
              r = boost::spirit::karma::char_inserter<>::call(sink, '0');
          if (precision && r)
              r = boost::spirit::karma::int_inserter<10>::call(sink, n);
          return r;
      }
      // END: added bugfix from https://github.com/boostorg/spirit/issues/529
    };
    static const boost::spirit::karma::real_generator<double, DoublePolicy> doubleBitIdentical;
    static boost::spirit::karma::rule<std::ostream_iterator<char>, double()> mydouble=doubleBitIdentical;
    return mydouble;
  }

  template<>
  boost::spirit::karma::rule<std::ostream_iterator<char>, int()>& getBoostSpiritKarmaRule<int>() {
    static boost::spirit::karma::rule<std::ostream_iterator<char>, int()> myint=boost::spirit::karma::int_;
    return myint;
  }
  template<>
  boost::spirit::karma::rule<std::ostream_iterator<char>, long()>& getBoostSpiritKarmaRule<long>() {
    static boost::spirit::karma::rule<std::ostream_iterator<char>, long()> myint=boost::spirit::karma::long_;
    return myint;
  }

  namespace {
    // std::real/std::image cannot be used as a function pointer in boost::phoenix!??????
    double myReal(const std::complex<double> &c) { return std::real(c); }
    double myImag(const std::complex<double> &c) { return std::imag(c); }
  }
  template<>
  boost::spirit::karma::rule<std::ostream_iterator<char>, std::complex<double>()>& getBoostSpiritKarmaRule<std::complex<double>>() {
    auto &mydouble = getBoostSpiritKarmaRule<double>();
    static boost::spirit::karma::rule<std::ostream_iterator<char>, std::complex<double>()> mycomplex =
      '(' << mydouble[boost::spirit::karma::_1=boost::phoenix::bind(&myReal, boost::spirit::karma::_val)] <<
      '+' << mydouble[boost::spirit::karma::_1=boost::phoenix::bind(&myImag, boost::spirit::karma::_val)] << "*i)";
    return mycomplex;
  }

}
