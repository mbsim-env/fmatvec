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

  //!!! we explicitly instantiate a lot of mat/vec stream operators to avoid compiling in each included file boost::spirit !!!
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,std::complex<double>> &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,std::complex<double>> &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,int                 > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,int                 > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Var     ,int                 > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Var     ,int                 > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,int                 > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,int                 > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Rotation ,Fixed<3>,Fixed<3>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Rotation ,Fixed<3>,Fixed<3>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Ref     ,Ref     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Ref     ,Ref     ,double              > &A);
//template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Var     ,Var     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Var     ,Var     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Fixed<2>,Fixed<2>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<2>,Fixed<2>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Fixed<3>,Fixed<3>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<3>,Fixed<3>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Var     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Var     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<2>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<2>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<3>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<3>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Var     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Var     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Var     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Var     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Var     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Var     ,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<2>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<2>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<3>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<3>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<6>,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<6>,double              > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Rotation ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Rotation ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Ref     ,Ref     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Ref     ,Ref     ,SymbolicExpression  > &A);
//template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Var     ,Var     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Var     ,Var     ,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Symmetric,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Var     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Var     ,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<2>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<2>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<3>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<3>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Var     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Var     ,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Var     ,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Var     ,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<6>,SymbolicExpression  > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<6>,SymbolicExpression  > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Rotation ,Fixed<3>,Fixed<3>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Rotation ,Fixed<3>,Fixed<3>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Diagonal ,Ref     ,Ref     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Diagonal ,Ref     ,Ref     ,double              > &A);

}
