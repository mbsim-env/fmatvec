#ifndef _FMATVEC_ATOM_H
#define _FMATVEC_ATOM_H

#include <ostream>
#include <iostream>
#include <functional>
#include <memory>
#include <array>
#include <sstream>
#include <chrono>
#include "stream.h"

namespace fmatvec {
FMATVEC_MSVC_DISABLEW4251_BEGIN

class AdoptCurrentMessageStreamsUntilScopeExit;

class DisableEscaping {};
class EnableEscaping {};
FMATVEC_EXPORT extern DisableEscaping disableEscaping;
FMATVEC_EXPORT extern EnableEscaping enableEscaping;

// Write synchronized to str_.
// The synchronization is done globally (not seperately for each str_)
// The content is buffered and written to str_ at destruction time of this object.
// The two special output manipulators disableEscaping and enableEscaping can be used (always as a pair)
// to enforce a none escaping mode of str_ if str_ supports this.
// If so, disableEscaping is translated to str_<<flush<<skipws and
// enableEscaping to std_<<flush<<noskipws.
// What str_ actually does on skipws/noskipws is defined by the used class for str_.
class FMATVEC_EXPORT osyncstream
#ifndef SWIG // swig can not parse this however it is not needed for swig
  : public std::ostream
#endif
{
  public:
    osyncstream(std::ostream &str_);
    osyncstream(osyncstream &) = delete;
    osyncstream(osyncstream &&) noexcept;
    osyncstream& operator=(osyncstream &) = delete;
    osyncstream& operator=(osyncstream &&) = delete;
    ~osyncstream() override;
    osyncstream& operator<<(DisableEscaping&);
    osyncstream& operator<<(EnableEscaping&);
    osyncstream& operator<<(bool value)                             { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(long value)                             { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(unsigned long value)                    { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(long long value)                        { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(unsigned long long value)               { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(double value)                           { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(long double value)                      { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(const void* value)                      { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(std::nullptr_t)                         { std::ostream::operator<<(nullptr); return *this; }
    osyncstream& operator<<(short value)                            { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(int value)                              { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(unsigned short value)                   { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(unsigned int value)                     { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(float value)                            { std::ostream::operator<<(value  ); return *this; }
    osyncstream& operator<<(std::streambuf* sb)                     { std::ostream::operator<<(sb     ); return *this; }
    osyncstream& operator<<(std::ios_base& (*func)(std::ios_base&)) { std::ostream::operator<<(func   ); return *this; }
    osyncstream& operator<<(std::ostream& (*func)(std::ostream&))   { std::ostream::operator<<(func   ); return *this; }

    // Use this function with care!!! It's mainly not indented to be used.
    std::ostream &getOStream() const { return str; }
  private:
    std::stringbuf buf;
    std::ostream &str;
    bool moved { false };
};

inline osyncstream& operator<<(osyncstream& os, char ch)                { operator<<(static_cast<std::ostream&>(os), ch); return os; }
inline osyncstream& operator<<(osyncstream& os, signed char ch)         { operator<<(static_cast<std::ostream&>(os), ch); return os; }
inline osyncstream& operator<<(osyncstream& os, unsigned char ch)       { operator<<(static_cast<std::ostream&>(os), ch); return os; }
inline osyncstream& operator<<(osyncstream& os, const char* s)          { operator<<(static_cast<std::ostream&>(os), s ); return os; }
inline osyncstream& operator<<(osyncstream& os, const signed char* s)   { operator<<(static_cast<std::ostream&>(os), s ); return os; }
inline osyncstream& operator<<(osyncstream& os, const unsigned char* s) { operator<<(static_cast<std::ostream&>(os), s ); return os; }
inline osyncstream& operator<<(osyncstream& os, const std::string &s)   { operator<<(static_cast<std::ostream&>(os), s ); return os; }

/*! Top level class.
 * This is the top level class which is used for (at least) all classes
 * which may be created using a object factory.
 * This class contains only totally basic functionallity like streams for
 * printing messages. No mathematical or other "none" basic content should
 * be added here.
 */
class FMATVEC_EXPORT Atom {
  friend class AdoptCurrentMessageStreamsUntilScopeExit;
  public:

    //! Messages can be printed to different message types named here.
    // When adding new message type a stream and a initial active flag must be provided in atom.cc (see NEW TYPES HERE)
    enum MsgType {
      Info,       // Informational messages
      Warn,       // Warning messages
      Debug,      // Debugging messages
      Error,      // Error messages
      Deprecated, // Error messages
      Status,     // Status messages (only the last message is relevant)
      SIZE        // Must be the last enum in this list
    };
    #define FMATVEC_ATOM_MSGTYPE_SIZE 6
#ifndef SWIG // swig can not parse this however it is not needed for swig
    static_assert(SIZE==FMATVEC_ATOM_MSGTYPE_SIZE, "The proprocessor define FMATVEC_ATOM_MSGTYPE_SIZE must be equal Atom::SIZE.");
#endif

  protected:
    //! When a Atom is default constructed use the current statically set message streams.
    Atom();
    //! When a Atom is copy constructed use the current statically set message streams, not the message streams from src.
    Atom(const Atom &src);
  public:
    //! dtor.
    virtual ~Atom();
#ifndef SWIG // no assignment operator for swig
    //! When a Atom is assinged do not change the messsage streams since we always use the message streams being active at ctor time.
    Atom& operator=(const Atom &);
#endif

    //! Set the current message stream used by all subsequent created objects.
    //! type defines the message type which should be set using this call.
    //! If s is not defined, the message of type type prints to cout.
    //! If a is not defined, a new shared bool flag set to true is used.
    //! Be aware of data races in streams if objects of type Atom print messages in threads.
    static void setCurrentMessageStream(MsgType type,
                  const std::shared_ptr<bool> &a=std::make_shared<bool>(true),
                  const std::shared_ptr<std::ostream> &s=std::make_shared<std::ostream>(std::cout.rdbuf()));

    //! Set the active flag of this object and all objects which were created using the same message stream as this object.
    void setMessageStreamActive(MsgType type, bool activeFlag);

    //! Get the shared message stream active flag and the shared message stream of this object
    void getMessageStream(MsgType type,
           std::shared_ptr<bool> &a,
           std::shared_ptr<std::ostream> &s) const;

    //! Adopt the message streams from src to this object.
    //! If src is NULL adopt the current (static) message streams.
    //! Normally always the streams at ctor time are used. But in some special cassed this function is usefull.
    void adoptMessageStreams(const Atom *src=nullptr);

    //! Return the message stream of type type.
    //! Node: If the code is performance critical you should check first whether this stream is really
    //! printed using msgAct(type). If this return false just skip the complete message.
    osyncstream msg(MsgType type) const {
      return *_msg[type];
    }
    //! Return true if the the message of type type is currently active.
    //! Note: If the code is not performance critical their is no need to check this flag. You can
    //! just print using msg(type)<<"Hello world"<<endl; and it is not really printed.
    bool msgAct(MsgType type) const {
      return *_msgAct[type];
    }

    //! Same as msg(type).
    //! Use this function only if not object is available. This should normally not be the case.
    static osyncstream msgStatic(MsgType type) {
      return *_msgStatic[type];
    }
    //! Same as msgAct(type).
    //! Use this function only if not object is available. This should normally not be the case.
    static bool msgActStatic(MsgType type) {
      return *_msgActStatic[type];
    }

  private:

    // A stream which prints to null.
    static std::shared_ptr<std::ostream> _nullStream;

    // Static pointer arrays of streams and active flags which were used for newly created objects.
    // These can be changed using setCurrentMessageStream(...)
    static std::array<std::shared_ptr<bool        >, SIZE> _msgActStatic;
    static std::array<std::shared_ptr<std::ostream>, SIZE> _msgSavedStatic;
    static std::array<std::shared_ptr<std::ostream>, SIZE> _msgStatic;

    // Pointer arrays to streams and active flags this object uses.
    // (these have a life-time at least as long as the object itself, ensured by reference counting)
    std::array<std::shared_ptr<bool        >, SIZE> _msgAct;
    std::array<std::shared_ptr<std::ostream>, SIZE> _msgSaved;
    std::array<std::shared_ptr<std::ostream>, SIZE> _msg;
};

#ifndef SWIG
//! A ostream object which prefix/postfix every messasge.
//! Before the output happens the message is converted using escapingFunc, if given.
//! Then result is then passed to outputFunc of printed to outstr_ (dependent on the used ctor).
//! The conversion using escapingFunc is skipped if the stream has set the skipws formatting flag, which is not set at the beginning.
//! However, you have ALWAYS to call flush before skipws and to call flush and noskipws afterwards!!!
//! But note that normally this is done by using disableEscaping and enableEscaping on the osyncstream of Atom::msg/Atom::msgStatic.
//! (this flag is used since it is not relevant for a ostream and can hence be "misused" here).
class FMATVEC_EXPORT PrePostfixedStream : public std::ostream {
  public:
    //! Call f with every message of the stream, prefixed/postfixed.
    PrePostfixedStream(const std::string &prefix_, const std::string &postfix_,
                       const std::function<void(const std::string &)> &outputFunc,
                       const std::function<void(std::string&)> &escapingFunc=nullptr);
    //! Convinence function to print to outstr_.
    PrePostfixedStream(const std::string &prefix_, const std::string &postfix_, std::ostream &outstr_,
                       const std::function<void(std::string&)> &escapingFunc=nullptr);
  private:
    class StringBuf : public std::streambuf {
      public:
        StringBuf(const PrePostfixedStream &stream_, std::string prefix_, std::string postfix_,
                  const std::function<void(const std::string &)> &outputFunc_,
                  const std::function<void(std::string&)> &escapingFunc_);
      protected:
        int overflow(int c) override;
        std::streamsize xsputn(const char_type* s, std::streamsize count) override;
        int sync() override;
      private:
        void flushBuffer();
        std::string buffer;
        std::function<void(const std::string &)> outputFunc;
        const PrePostfixedStream &stream;
        const std::string prefix;
        const std::string postfix;
        const std::function<void(std::string&)> escapingFunc;
    };
    StringBuf buffer;
};
#endif

//! Adopt the current message streams for the livetime of this object.
//! This can be used to set the current message streams to an existing object.
class FMATVEC_EXPORT AdoptCurrentMessageStreamsUntilScopeExit {
  public:
    //! Save the current message streams and then set it to be equal to src.
    //! The current message streams are reset to the saved ones when this object goes out of scope.
    AdoptCurrentMessageStreamsUntilScopeExit(const Atom* src);
    ~AdoptCurrentMessageStreamsUntilScopeExit();
    AdoptCurrentMessageStreamsUntilScopeExit(const AdoptCurrentMessageStreamsUntilScopeExit& other) = delete;
    AdoptCurrentMessageStreamsUntilScopeExit(AdoptCurrentMessageStreamsUntilScopeExit&& other) noexcept = delete;
    AdoptCurrentMessageStreamsUntilScopeExit& operator=(const AdoptCurrentMessageStreamsUntilScopeExit& other) = delete;
    AdoptCurrentMessageStreamsUntilScopeExit& operator=(AdoptCurrentMessageStreamsUntilScopeExit&& other) noexcept = delete;
  private:
    std::array<std::pair<std::shared_ptr<bool>, std::shared_ptr<std::ostream>>, Atom::SIZE> savedStreams;
};

// BEGIN: measure time (just used for performance debugging)
class DebugMeasureAccumulatedTime {
  public:
    DebugMeasureAccumulatedTime(std::chrono::duration<double> &sum_, size_t &count_) : sum(sum_), count(count_) {
      start = std::chrono::steady_clock::now();
    }
    ~DebugMeasureAccumulatedTime() {
      auto end = std::chrono::steady_clock::now();
      auto delta = std::chrono::duration<double>(end - start);
      sum += delta;
      count++;
    }
  private:
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::duration<double> &sum;
    size_t &count;
};
class DebugMeasureAccumulatedTimeOnExit {
  public:
    DebugMeasureAccumulatedTimeOnExit(const std::string &id_, std::chrono::duration<double> &sum_, size_t &count_) : id(id_), sum(sum_), count(count_) {
    }
    ~DebugMeasureAccumulatedTimeOnExit() {
      std::cout<<"fmatvec_measure_id,"<<id<<",totalsec,"<<sum.count()<<",count,"<<count<<",avgsec,"<<sum.count()/count<<std::endl;
    }
  private:
    std::string id;
    std::chrono::duration<double> &sum;
    size_t &count;
};
// measure the accumulated time and number of calls of a block
#define FMATVEC_MEASURE_BLOCK(id) \
  static std::chrono::duration<double> fmatvec_measure_sum_##id {}; \
  static size_t fmatvec_measure_count_##id {}; \
  static fmatvec::DebugMeasureAccumulatedTimeOnExit fmatvec_measure_onexit_##id(#id, fmatvec_measure_sum_##id, fmatvec_measure_count_##id); \
  fmatvec::DebugMeasureAccumulatedTime fmatvec_measure_dummy_##id(fmatvec_measure_sum_##id, fmatvec_measure_count_##id);
// measure the accumulated time and number of calls of all blocks which contain FMATVEC_MEASURE_RUN(id)
#define FMATVEC_MEASURE_BLOCKINIT(id) \
  static std::chrono::duration<double> fmatvec_measure_sum_##id {}; \
  static size_t fmatvec_measure_count_##id {}; \
  static fmatvec::DebugMeasureAccumulatedTimeOnExit fmatvec_measure_onexit_##id(#id, fmatvec_measure_sum_##id, fmatvec_measure_count_##id);
// see, the define above
#define FMATVEC_MEASURE_BLOCKEXE(id) \
  fmatvec::DebugMeasureAccumulatedTime fmatvec_measure_dummy_##id(fmatvec_measure_sum_##id, fmatvec_measure_count_##id);
// start time measure for id
#define FMATVEC_MEASURE_START(id) \
  auto fmatvec_measure_start_##id = std::chrono::steady_clock::now();
// dump elapsed time since the start of id (see the define above)
#define FMATVEC_MEASURE_DUMP(id) \
  { \
    auto end = std::chrono::steady_clock::now(); \
    auto delta = std::chrono::duration<double>(end - fmatvec_measure_start_##id); \
    std::cout<<"fmatvec_measure_start_end_id,"<<#id<<",sec,"<<delta.count()<<std::endl; \
  }
// BEGIN: measure time (just used for performance debugging)

FMATVEC_MSVC_DISABLEW4251_END
}

#endif
