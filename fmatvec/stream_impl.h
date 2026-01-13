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

#ifndef stream_impl_h
#define stream_impl_h

#include "stream.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <boost/scope_exit.hpp>
#include "range.h"
#include "types.h"
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/bind/bind_function.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/spirit/include/karma_real.hpp>
#include <boost/spirit/include/karma.hpp>

namespace {
  template<class AT>
  std::vector<std::vector<AT>> scalarToVecVec(const AT& x) {
    return std::vector<std::vector<AT>>(1, std::vector<AT>(1, x));
  }
}

namespace fmatvec {

  /*! \brief Matrix input
   *
   * This function loads a matrix from a stream.
   * \param is An input stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the input stream.
   * */
  template <class Type, class Row, class Col, class AT> std::istream& operator>>(std::istream &s, Matrix<Type,Row,Col,AT> &A) {
    namespace qi = boost::spirit::qi;
    namespace phx = boost::phoenix;
    using It = boost::spirit::istream_iterator;

    qi::rule<It, std::vector<std::vector<AT>>()> scalar;
    qi::rule<It, std::vector<AT>()> row;
    qi::rule<It, std::vector<std::vector<AT>>()> matrix;
    qi::rule<It, std::vector<std::vector<AT>>()> scalarOrMatrix;

    qi::rule<It, AT()> &atomicType = getBoostSpiritQiRule<AT>();
    scalar = *qi::space >> atomicType[qi::_val=phx::bind(&scalarToVecVec<AT>, qi::_1)];
    row = atomicType % (+qi::blank | (*qi::blank >> ',' >> *qi::blank));
    matrix = *qi::space >> '[' >> *qi::space >>
               -(row % (*qi::blank >> (';' | qi::eol) >> *qi::blank))[qi::_val=qi::_1] >>
             *qi::space >> ']';
    scalarOrMatrix = scalar | matrix;

    auto savedFlags=s.flags();
    s.unsetf(std::ios::skipws);
    BOOST_SCOPE_EXIT_TPL(&s, &savedFlags) {
      s.flags(savedFlags);
    } BOOST_SCOPE_EXIT_END

    std::vector<std::vector<AT>> Avecvec;
    if(!qi::parse(It(s), It(), scalarOrMatrix, Avecvec))
      throw std::runtime_error("The stream does not contain a valid scalar, vector or matrix expression. Not parsed content of stream:\n"+
                               std::string(std::istreambuf_iterator<char>(s), std::istreambuf_iterator<char>()));
    A<<=Matrix<Type,Row,Col,AT>(Avecvec);
    return s;
  }

  /*! \brief Matrix output 
   *
   * This function writes a matrix into a stream.
   * \param os An output stream.
   * \param A A matrix of any shape and type.
   * \return A reference to the output stream.
   * */
  template <class Type, class Row, class Col, class AT> std::ostream& operator<<(std::ostream &os, const Matrix<Type,Row,Col,AT> &A) {
    namespace karma = boost::spirit::karma;
    using It = std::ostream_iterator<char>;

    static boost::spirit::karma::rule<It, std::vector<std::vector<AT>>()> matrix;
    static bool init=false;
    if(!init) {
      static karma::rule<It, std::vector<AT>()> row;

      auto &atomicType=getBoostSpiritKarmaRule<AT>();

      row = -(atomicType % ", ")[karma::_1=karma::_val];
      matrix = '[' << -(row % "; ")[karma::_1=karma::_val] << ']';
    }

    auto Avecvec=static_cast<std::vector<std::vector<AT>>>(A);
    if(!karma::generate(It(os), matrix, Avecvec))
      throw std::runtime_error("Failed to write matrix to stream");
    return os;
  }

}

#endif
