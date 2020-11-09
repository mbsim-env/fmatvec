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
