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

#ifndef stream_h
#define stream_h

#include <vector>
#include "matrix.h"
#include <iterator>
#include <boost/spirit/home/qi/nonterminal/nonterminal_fwd.hpp>
#include <boost/spirit/home/karma/nonterminal/nonterminal_fwd.hpp>

namespace boost {
  namespace spirit {
    template<typename Elem, typename Traits> class basic_istream_iterator;
    typedef basic_istream_iterator<char, std::char_traits<char>> istream_iterator;
  }
}

namespace fmatvec {

  class SymbolicExpression;
  class IndependentVariable;

  template<class AT>
  boost::spirit::qi::rule<boost::spirit::istream_iterator, AT()>& getBoostSpiritQiRule();

  template<class AT>
  boost::spirit::karma::rule<std::ostream_iterator<char>, AT()>& getBoostSpiritKarmaRule();

  template<> boost::spirit::qi::rule<boost::spirit::istream_iterator, double()>& getBoostSpiritQiRule<double>();
  template<> boost::spirit::qi::rule<boost::spirit::istream_iterator, int()>& getBoostSpiritQiRule<int>();
  template<> boost::spirit::qi::rule<boost::spirit::istream_iterator, std::complex<double>()>& getBoostSpiritQiRule<std::complex<double>>();
  template<> boost::spirit::karma::rule<std::ostream_iterator<char>, double()>& getBoostSpiritKarmaRule<double>();
  template<> boost::spirit::karma::rule<std::ostream_iterator<char>, int()>& getBoostSpiritKarmaRule<int>();
  template<> boost::spirit::karma::rule<std::ostream_iterator<char>, std::complex<double>()>& getBoostSpiritKarmaRule<std::complex<double>>();

  /*! \brief Matrix input
   *
   * This function loads a matrix from a stream.
   * \param is An input stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the input stream.
   * */
  template<class Type, class Row, class Col, class AT> std::istream& operator>>(std::istream &s, Matrix<Type,Row,Col,AT> &A);

  /*! \brief Matrix output 
   *
   * This function writes a matrix into a stream.
   * \param os An output stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the output stream.
   * */
  template<class Type, class Row, class Col, class AT> std::ostream& operator<<(std::ostream &s, const Matrix<Type,Row,Col,AT> &A);

  //!!! we explicitly instantiate a lot of mat/vec stream operators to avoid compiling in each included file boost::spirit !!!
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,std::complex<double>> &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,std::complex<double>> &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,int                 > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,int                 > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Var     ,int                 > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Var     ,int                 > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,int                 > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,int                 > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Rotation ,Fixed<3>,Fixed<3>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Rotation ,Fixed<3>,Fixed<3>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Ref     ,Ref     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Ref     ,Ref     ,double              > &A);
//extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Var     ,Var     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Var     ,Var     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Fixed<2>,Fixed<2>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<2>,Fixed<2>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Fixed<3>,Fixed<3>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<3>,Fixed<3>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Var     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Var     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<2>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<2>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<3>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<3>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Var     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Var     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Var     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Var     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Var     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Var     ,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<2>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<2>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<3>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<3>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<6>,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<6>,double              > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Rotation ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Rotation ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Ref     ,Ref     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Ref     ,Ref     ,SymbolicExpression  > &A);
//extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Var     ,Var     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Var     ,Var     ,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Symmetric,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Symmetric,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Var     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Var     ,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<2>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<2>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<3>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<3>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Var     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Var     ,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Var     ,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Var     ,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<2>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<3>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<6>,SymbolicExpression  > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<6>,SymbolicExpression  > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Rotation ,Fixed<3>,Fixed<3>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Rotation ,Fixed<3>,Fixed<3>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,IndependentVariable > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,IndependentVariable > &A);
  extern template std::istream FMATVEC_EXPORT & operator>>(std::istream &s,       Matrix<Diagonal ,Ref     ,Ref     ,double              > &A);
  extern template std::ostream FMATVEC_EXPORT & operator<<(std::ostream &s, const Matrix<Diagonal ,Ref     ,Ref     ,double              > &A);

}

#endif
