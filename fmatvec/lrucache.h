#ifndef _FMATVEC_LRUCACHE_H_
#define _FMATVEC_LRUCACHE_H_

#include <list>
#include <cmath>
#include <unordered_map>
#include <functional>
#include <boost/container_hash/hash.hpp>
#include <fmatvec/fmatvec.h>

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
      LRUCache(size_t size_, size_t bucketSize = 0);
      /*! Add the key "in" to the cache.
       * Ensures that no more then size entries exist in the cache.
       * Returns a reference to the cache value and a bool which indicates if a new entry was created.
       * If a new entry was created, than the reference is a uninizialized new value which must be set.
       * If a existing entry is returned, than the reference is a reference to the existing cache entry.
       */
      std::pair<Out&, bool> operator()(const In &in);
    private:
      std::list<std::pair<const In&,Out>, ListAllocator> itemList;
      std::unordered_map<In, decltype(itemList.begin()), InHash, InKeyEqual, UnorderedMapAllocator> itemMap;
      size_t size;
  };
  
  template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
  LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::LRUCache(size_t size_, size_t bucketSize) :
    itemMap(bucketSize==0 ? floor(1.5*size_) : bucketSize), size(size_) {
  }
  
  template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
  std::pair<Out&, bool> LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::operator()(const In &in) {
    auto [it, created] = itemMap.emplace(in, itemList.begin());
  
    // item found -> save out and remove it from itemList and itemMap to be readded a last recently used
    if(!created) {
      itemList.emplace_front(it->first, std::move(it->second->second)); // move the exiting out value a new first list (Mark *)
      itemList.erase(it->second);
      itemMap.erase(it);
    }
    else
      itemList.emplace_front(it->first, Out()); // create a new first list item the list with a default ctor out value (Mark *)
  
    // add a new last recently used element
    // Mark *: to use only the default ctor and move ctor of Out the creation of the new item in itemList
    //         is done in the above if/else construct. Doing it here will be more readable but requires also a move assignment operator.
    auto newIt = itemList.begin();
    if(created)
      it->second=newIt;
    
    // remove old item(s) if too many items exist
    while(itemMap.size() > size) {
      auto it = itemList.end();
      --it;
      itemMap.erase(it->first);
      itemList.pop_back();
    }
  
    return {newIt->second,created};
  }

}

#endif
