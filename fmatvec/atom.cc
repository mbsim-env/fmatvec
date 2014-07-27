#include "fmatvec/atom.h"
#include <boost/preprocessor/iteration/local.hpp>

using namespace std;
using namespace boost;

namespace {
  bool msgTypeActive(bool defaultValue, const string &envVarPostfix) {
    char *env=getenv((string("FMATVEC_MSG_")+envVarPostfix).c_str());
    if(env) {
      if(string(env)=="0") return false;
      if(string(env)=="1") return true;
    }
    return defaultValue;
  }
}

namespace fmatvec {

// a shared null stream object
shared_ptr<ostream> Atom::_nullStream = make_shared<ostream>(static_cast<streambuf*>(NULL));

// initialize the static streams with cout/cerr and corresponding active flags
array<shared_ptr<bool>   , Atom::SIZE> Atom::_msgActStatic = {{
  // initial active flag of the message stream and/or envvar postfix being used as initial value
  make_shared<bool>(msgTypeActive(true , "Info" )), // Info
  make_shared<bool>(msgTypeActive(true , "Warn" )), // Warn
  make_shared<bool>(msgTypeActive(false, "Debug"))  // Debug
  // NEW TYPES HERE: initial active flag
}};
array<shared_ptr<ostream>, Atom::SIZE> Atom::_msgSavedStatic = {{
  // stream to use for the message stream
  make_shared<ostream>(cout.rdbuf()), // Info
  make_shared<ostream>(cerr.rdbuf()), // Warn
  make_shared<ostream>(cout.rdbuf())  // Debug
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

Atom::~Atom() {
}

Atom& Atom::operator=(const Atom &) {
  return *this;
}

void Atom::setCurrentMessageStream(MsgType type, const boost::shared_ptr<std::ostream> &s) {
  _msgActStatic[type]=boost::make_shared<bool>(true);
  _msgSavedStatic[type]=s;
  _msgStatic[type] = _msgSavedStatic[type];
}

void Atom::setCurrentMessageStreamActive(MsgType type, bool activeFlag) {
  *_msgActStatic[type] = activeFlag;
  _msgStatic[type] = *_msgActStatic[type] ? _msgSavedStatic[type] : _nullStream;
}

void Atom::setMessageStreamActive(MsgType type, bool activeFlag) {
  *_msgAct[type] = activeFlag;
  _msg[type] = *_msgAct[type] ? _msgSaved[type] : _nullStream;
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

}
