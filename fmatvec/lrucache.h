#ifndef _FMATVEC_LRUCACHE_H_
#define _FMATVEC_LRUCACHE_H_

#include <list>
#include <unordered_map>

namespace fmatvec {

/*! A generic fast (O(1)) LRU cache */
template<class In,
         class Out,
         class InHash = std::hash<In>,
         class InKeyEqual = std::equal_to<In>,
         class UnorderedMapAllocator = std::allocator<std::pair<const In, typename std::list<std::pair<In,Out>>::iterator>>,
         class ListAllocator = std::allocator<std::pair<In, Out>> >
class LRUCache {
  public:
    //! Create a LRU cache of at most size items.
    LRUCache(size_t size_, size_t bucketSize = 0) : itemMap(bucketSize==0?size_:bucketSize), size(size_) {};
    /*! Add the key "in" to the cache.
     * Ensures that no more then size entries exist in the cache.
     * Returns a reference to the cache value and a bool which indicates if a new entry was created.
     * If a new entry was created, than the reference is a uninizialized new value which must be set.
     * If a existing entry is returned, than the reference is a reference to the existing cache entry.
     */
    std::pair<Out&, bool> operator()(const In &in);
  private:
    std::list<std::pair<In,Out>, ListAllocator> itemList;
    std::unordered_map<In, decltype(itemList.begin()), InHash, InKeyEqual, UnorderedMapAllocator> itemMap;
    size_t size;
};

template<class In, class Out, class InHash, class InKeyEqual, class UnorderedMapAllocator, class ListAllocator>
std::pair<Out&, bool> LRUCache<In,Out,InHash,InKeyEqual,UnorderedMapAllocator,ListAllocator>::operator()(const In &in) {
  auto [it, created] = itemMap.emplace(in, itemList.begin());
  Out out;

  // item found -> save out and remove it from itemList and itemMap to be readded a last reasondly used
  if(!created) {
    out = it->second->second;
    itemList.erase(it->second);
    itemMap.erase(it);
  }

  // add a new last reaseondly used element
  itemList.emplace_front(in, out);
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
