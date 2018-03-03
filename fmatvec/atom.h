#ifndef _FMATVEC_ATOM_H
#define _FMATVEC_ATOM_H

#include <ostream>
#include <iostream>
#include <memory>
#include <array>

namespace fmatvec {

/*! Top level class.
 * This is the top level class which is used for (at least) all classes
 * which may be created using a object factory.
 * This class contains only totally basic functionallity like streams for
 * printing messages. No mathematical or other "none" basic content should
 * be added here.
 */
class Atom {
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
           std::shared_ptr<std::ostream> &s);

    //! Adopt the message streams from src to this object.
    //! If src is NULL adopt the current (static) message streams.
    //! Normally always the streams at ctor time are used. But in some special cassed this function is usefull.
    void adoptMessageStreams(const Atom *src=nullptr);

    //! Return the message stream of type type.
    //! Node: If the code is performance critical you should check first whether this stream is really
    //! printed using msgAct(type). If this return false just skip the complete message.
    std::ostream &msg(MsgType type) {
      return *_msg[type];
    }
    //! Return true if the the message of type type is currently active.
    //! Note: If the code is not performance critical their is no need to check this flag. You can
    //! just print using msg(type)<<"Hello world"<<endl; and it is not really printed.
    bool msgAct(MsgType type) {
      return *_msgAct[type];
    }

    //! Same as msg(type).
    //! Use this function only if not object is available. This should normally not be the case.
    static std::ostream &msgStatic(MsgType type) {
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

}

#endif
