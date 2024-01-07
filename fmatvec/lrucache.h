#ifndef _FMATVEC_LRUCACHE_H_
#define _FMATVEC_LRUCACHE_H_

#include <list>
#include <cmath>
#include <unordered_map>
#include <functional>
#include <boost/container_hash/hash.hpp>
#include <fmatvec/fmatvec.h>
#include <fmatvec/atom.h>

namespace boost {

  template<class Row, class AT>
  struct hash<fmatvec::Vector<Row,AT>> {
    size_t operator()(const fmatvec::Vector<Row,AT> &v) const {
      size_t h=0;
      for(int i=0; i<v.size(); ++i)
        hash_combine(h, v(i));
      return h;
    };
  };

  template<class Col, class AT>
  struct hash<fmatvec::RowVector<Col,AT>> {
    size_t operator()(const fmatvec::RowVector<Col,AT> &v) const {
      size_t h=0;
      for(int i=0; i<v.size(); ++i)
        hash_combine(h, v(i));
      return h;
    };
  };

  template<class Type, class Row, class Col, class AT>
  struct hash<fmatvec::Matrix<Type,Row,Col,AT>> {
    size_t operator()(const fmatvec::Matrix<Type,Row,Col,AT> &v) const {
      size_t h=0;
      for(int r=0; r<v.rows(); ++r)
        for(int c=0; c<v.cols(); ++c)
          hash_combine(h, v(r,c));
      return h;
    };
  };

}

namespace fmatvec {

  /*! A generic fast (O(1)) LRU cache */
  template<class In,
           class Out,
           class InHash = boost::hash<In>,
           class InKeyEqual = std::equal_to<In>,
           class UnorderedMapAllocator = std::allocator<std::pair<const In, typename std::list<std::pair<const In&,Out>>::iterator>>,
           class ListAllocator = std::allocator<std::pair<const In&, Out>> >
  class LRUCache {
    public:
      //! Create a LRU cache of at most size items.
      LRUCache(size_t size_, double bucketSizeFactor = 1.5);
      //! Prints the cache hit rate in debug builds to the stream Atom::Debug
      ~LRUCache();
      /*! Add the key "in" to the cache.
       * Ensures that no more then size entries exist in the cache.
       * Returns a reference to the cache value and a bool which indicates if a new entry was created.
       * If a new entry was created, than the reference is a uninizialized new value which must be set.
       * If a existing entry is returned, than the reference is a reference to the existing cache entry.
       */
      std::pair<Out&, bool> operator()(const In &in);
      double getCacheHitRate() const;
    private:
      std::list<std::pair<const In&,Out>, ListAllocator> itemList;
      std::unordered_map<In, decltype(itemList.begin()), InHash, InKeyEqual, UnorderedMapAllocator> itemMap;
      size_t size;
      size_t allCalls { 0 }, cacheHits { 0 };
  };
  
  template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
  LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::LRUCache(size_t size_, double bucketSizeFactor) {
    // we can override the LRU cache size with FMATVEC_LRUCACHE_SIZE
    static auto envSize = static_cast<size_t>(-1);
    if(envSize == static_cast<size_t>(-1)) {
      const char *env = getenv("FMATVEC_LRUCACHE_SIZE");
      if(env)
        envSize = atoi(env);
      else
        envSize = 0;
    }
    if(envSize != 0)
      size_ = envSize;

    itemMap.rehash(floor(bucketSizeFactor*size_));
    size = size_;
  }

  template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
  LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::~LRUCache() {
    if(Atom::msgActStatic(Atom::Debug))
      Atom::msgStatic(Atom::Debug)<<"Cache hit rate "<<getCacheHitRate()*100<<"%"<<std::endl;
  }
  
  template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
  std::pair<Out&, bool> LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::operator()(const In &in) {
    allCalls++;
    auto [it, created] = itemMap.emplace(in, itemList.begin());
  
    if(created) {
      // new itemMap{in,it_dummy} was created by the above call (as "it")
      itemList.emplace_front(it->first, Out()); // create a new first itemList with a default ctor out value
      // remove old item(s) if too many items exist
      if(itemMap.size() > size) {
        auto itLast = --itemList.end();
        itemMap.erase(itLast->first);
        itemList.pop_back();
      }
    }
    else {
      cacheHits++;
      // a existing itemMap{in,it_to_itemList} was selected by the above call (as "it")
      itemList.emplace_front(it->first, std::move(it->second->second)); // create a new first itemList by move-ctor out from the existing
      itemList.erase(it->second); // remove the old itemList
    }
    it->second = itemList.begin(); // fix second value of itemMap (it was either created with a dummy iterator or has changed)
    return {it->second->second,created};
  }

  template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
  double LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::getCacheHitRate() const {
    return allCalls==0 ? -1 : static_cast<double>(cacheHits)/allCalls;
  }

}

#endif
