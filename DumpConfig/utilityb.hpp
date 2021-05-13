//
//  utilityb.hpp
//  IsingTMat
//
//  Created by Yi Hu on 8/5/20.
//

#ifndef utilityb_h
#define utilityb_h

#include <unordered_map>
#include <list>

namespace MyUtility{
/// \brief Set random generator seed with current time
  void setseed(void);
/// \brief Set random generator seed
  void setseed(unsigned int seed);
/// \brief get random number
  double rand();
};

/// \brief LRU cache
template <typename K, typename C>
class Cache
{
private:
  std::list< std::pair<K,C> > item_list;
  std::unordered_map<K, decltype(item_list.begin()) > item_map;
  int cache_size;
public:
  Cache(int size): cache_size(size){}
  void put(const K &key, const C &val)
  {
    auto it = item_map.find(key);
    if(it != item_map.end() && it->second != item_list.begin() ){
      item_list.splice(item_list.begin(), item_list, it->second);
      it->second = item_list.begin();
    }
    else
    {
      if(item_map.size()>cache_size)
      {
        auto last_it = std::prev(item_list.end());
        item_map.erase(last_it->first);
        item_list.pop_back();
      }
      item_list.push_front({key,val});
      item_map.insert({key, item_list.begin()});
    }
  }
  const C *get(K key)
  {
    auto it = item_map.find(key);
    if(it == item_map.end()) return nullptr;
    else return &(it->second->second);
  }
};

#endif /* utilityb_h */
