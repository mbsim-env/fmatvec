#ifndef _FMATVEC_ATOM_H
#define _FMATVEC_ATOM_H

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
    virtual ~Atom();
};

}

#endif
