#ifndef InstantStore_hpp
#define InstantStore_hpp

#include <map>
#include <vector>
#include <unordered_map>
#include <list>
#include "ExaaltDef.hpp"
#include "StoreEntry.hpp"

class InstantStore {
public:
    InstantStore() = default;

    ~InstantStore() = default;

    int initialize(const std::string &path) {
        return 0;
    };

    int initialize(size_t capacity_) {
        this->capacity = capacity_;
        return 0;
    };

    bool contains(unsigned int dbKey, int64 key) {
        return contains(std::make_pair(dbKey, key));
    };

    bool contains(KeyPair keys) {
        return locations.find(keys) != locations.end();
    };

    int put(const Entry &entry) {
        if (!contains(entry.getKeys())) {
            data.emplace_front(entry);
            locations[entry.getKeys()] = data.begin();
        }
        return 0;
    };

    int put(unsigned int dbKey, int64 &key, const RawDataVector &rawData) {
        return put(Entry(std::make_pair(dbKey, key), rawData));
    };

    int get(unsigned int dbKey, int64 &key, Entry &entry) {
        auto k = std::make_pair(dbKey, key);
        if (contains(dbKey, key)) {
            data.splice(data.begin(), data, locations[k]);
            entry = *locations[k];
            return 0;
        } else {
            return KEY_NOTFOUND;
        }
    }

    int get(unsigned int dbKey, int64 &key, RawDataVector &data) {
        Entry entry;
        auto result = get(dbKey, key, entry);
        if (result == 0) {
            data = *entry.getData();
        }
        return result;
    };

    Entry push(const Entry &entry) {
        Entry delEntry;
        if (!contains(entry.getKeys())) {
            if (locations.size() >= capacity) {
                delEntry = pop();
            }
            data.emplace_front(entry);
            locations[entry.getKeys()] = data.begin();
        }
        return delEntry;
    };

    Entry pop() {
        auto delEntry = data.back();
        data.pop_back();
        locations.erase(delEntry.getKeys());
        return delEntry;
    }

    std::list<Entry> &getAll() {
        return data;
    }

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet) { return 0; }

    int sync() { return 0; }

    int erase(unsigned int dbKey, int64 &key) {
        auto k = std::make_pair(dbKey, key);
        data.erase(locations[k]);
        locations.erase(k);
        return 0;
    }

    size_t size() {
        return locations.size();
    }

    void print() {
        std::cout << "InstantStore num of blocks = " << locations.size() << std::endl;
    }

private:
    std::list<Entry> data;
    std::unordered_map<KeyPair, std::list<Entry>::iterator, pair_hash> locations;
    size_t capacity = 10000;

};


#endif /* InstantStore_hpp */
