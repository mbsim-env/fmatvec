#include <config.h>
#include "fmatvec/atom.h"
#include <boost/preprocessor/iteration/local.hpp>

using namespace std;

namespace fmatvec {

// a shared null stream object
shared_ptr<ostream> Atom::_nullStream = make_shared<ostream>(static_cast<streambuf*>(nullptr));

// initialize the static streams with cout/cerr and corresponding active flags
array<shared_ptr<bool>   , Atom::SIZE> Atom::_msgActStatic = {{
  // initial active flag of the message stream and/or envvar postfix being used as initial value
  make_shared<bool>(true ), // Info
  make_shared<bool>(true ), // Warn
  make_shared<bool>(false), // Debug
  make_shared<bool>(true ), // Error
  make_shared<bool>(true ), // Deprecated
  make_shared<bool>(true )  // Status
  // NEW TYPES HERE: initial active flag
}};
array<shared_ptr<ostream>, Atom::SIZE> Atom::_msgSavedStatic = {{
  // stream to use for the message stream
  make_shared<PrePostfixedStream>("Info: ", ""  , cout), // Info
  make_shared<PrePostfixedStream>("Warn: ", ""  , cerr), // Warn
  make_shared<PrePostfixedStream>("Debug:", ""  , cout), // Debug
  make_shared<PrePostfixedStream>(""      , ""  , cerr), // Error
  make_shared<PrePostfixedStream>("Depr: ", ""  , cerr), // Deprecated
  make_shared<PrePostfixedStream>(""      , "\r", cout)  // Status
  // NEW TYPES HERE: stream
}};

#define BOOST_PP_LOCAL_LIMITS (0, FMATVEC_ATOM_MSGTYPE_SIZE-1)
#define BOOST_PP_LOCAL_MACRO(type) *_msgActStatic[type] ? _msgSavedStatic[type] : _nullStream,
array<shared_ptr<ostream>, Atom::SIZE> Atom::_msgStatic = {{
  #include BOOST_PP_LOCAL_ITERATE()
}};

Atom::Atom() :
  _msgAct(_msgActStatic),
  _msgSaved(_msgSavedStatic),
  _msg(_msgStatic) {
}

Atom::Atom(const Atom &) :
  _msgAct(_msgActStatic),
  _msgSaved(_msgSavedStatic),
  _msg(_msgStatic) {
}

Atom::~Atom() = default;

Atom& Atom::operator=(const Atom &) {
  return *this;
}

void Atom::setCurrentMessageStream(MsgType type, const shared_ptr<bool> &a, const shared_ptr<ostream> &s) {
  _msgActStatic[type]=a;
  _msgSavedStatic[type]=s;
  _msgStatic[type] = *_msgActStatic[type] ? _msgSavedStatic[type] : _nullStream;
}

void Atom::setMessageStreamActive(MsgType type, bool activeFlag) {
  *_msgAct[type] = activeFlag;
  _msg[type] = *_msgAct[type] ? _msgSaved[type] : _nullStream;
}

void Atom::getMessageStream(MsgType type,
       shared_ptr<bool> &a,
       shared_ptr<ostream> &s) const {
  a=_msgAct[type];
  s=_msgSaved[type];
}

void Atom::adoptMessageStreams(const Atom *src) {
  if(src) {
    _msgAct  =src->_msgAct;
    _msgSaved=src->_msgSaved;
    _msg     =src->_msg;
  }
  else {
    _msgAct  =_msgActStatic;
    _msgSaved=_msgSavedStatic;
    _msg     =_msgStatic;
  }
}

PrePostfixedStream::PrePostfixedStream(const string &prefix, const string &postfix,
                                       const function<void(const string&)> &outputFunc,
                                       const function<void(string&)> &escapingFunc) :
  ostream(nullptr), buffer(*this, prefix, postfix, outputFunc, escapingFunc) {
  unsetf(skipws);
  rdbuf(&buffer);
}

PrePostfixedStream::PrePostfixedStream(const string &prefix, const string &postfix, ostream &outstr,
                                       const function<void(string&)> &escapingFunc) :
  ostream(nullptr), buffer(*this, prefix, postfix, [&outstr](const string &s){
    outstr<<s<<std::flush;
  }, escapingFunc) {
  unsetf(skipws);
  rdbuf(&buffer);
}

PrePostfixedStream::StringBuf::StringBuf(const PrePostfixedStream &stream_, const string prefix_, string postfix_,
                                         const function<void(const string&)> &outputFunc_,
                                         const function<void(string&)> &escapingFunc_) :
  stream(stream_), stringbuf(ios_base::out), prefix(std::move(prefix_)), postfix(std::move(postfix_)),
  outputFunc(outputFunc_), escapingFunc(escapingFunc_) {}

int PrePostfixedStream::StringBuf::sync() {
  auto msg=str();
  str("");
  if(msg.empty())
    return 0;
  if(msg=="\n") {
    outputFunc("\n");
    return 0;
  }
  if(!escapingFunc || (stream.flags() & skipws))
    outputFunc(prefix+msg+postfix);
  else {
    escapingFunc(msg);
    outputFunc(prefix+msg+postfix);
  }
  return 0;
}

AdoptCurrentMessageStreamsUntilScopeExit::AdoptCurrentMessageStreamsUntilScopeExit(const Atom* src) {
  for(int t=0; t<Atom::SIZE; ++t) {
    auto type=static_cast<Atom::MsgType>(t);
    savedStreams[type].first =Atom::_msgActStatic[type];
    savedStreams[type].second=Atom::_msgSavedStatic[type];
    Atom::setCurrentMessageStream(type, src->_msgAct[type], src->_msgSaved[type]);
  }
}

AdoptCurrentMessageStreamsUntilScopeExit::~AdoptCurrentMessageStreamsUntilScopeExit() {
  for(int t=0; t<Atom::SIZE; ++t) {
    auto type=static_cast<Atom::MsgType>(t);
    Atom::setCurrentMessageStream(type, savedStreams[type].first, savedStreams[type].second);
  }
}

}
