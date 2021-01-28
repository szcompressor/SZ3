//
//  MultiLevelStore.hpp
//  ParSplice-refact
//
//  Created on 9/22/20.
//  Copyright © 2020 dp. All rights reserved.
//

#ifndef MultiLevelStoreL1_hpp
#define MultiLevelStoreL1_hpp

#include "ExaaltDef.hpp"
#include "InstantStore.hpp"
#include "PersistentStore.hpp"
#include "CompressedStore.hpp"

template<class T>
class MultiLevelStoreL1 {
public:
    MultiLevelStoreL1() = default;

    ~MultiLevelStoreL1() {
        for (auto &e: instantStore.getAll()) {
            persistentStore.put(e);
        }
    };

    int initialize(const std::string &path, size_t l1_capacity, size_t l2_capacity,
                   size_t blockSize_, double errorBound_, size_t batchSize_) {
        instantStore.initialize(l1_capacity);
        persistentStore.initialize(path + "/exaalt_multilevel_nocompress");
        return 0;
    }

    void print() {
        instantStore.print();
        persistentStore.print();
    }


    int put(Entry &entry) {
        auto l1_out = instantStore.push(entry);
        if (l1_out.isValid() && l1_out.isDirty()) {
            persistentStore.put(l1_out);
        }
        return 0;
    }

    int put(unsigned int dbKey, int64 &key, RawDataVector &data) {
        Entry entry(std::make_pair(dbKey, key), data);
        return put(entry);
    }

    int get(unsigned int dbKey, int64 &key, RawDataVector &data) {
        auto l1_code = instantStore.get(dbKey, key, data);
        if (l1_code == KEY_NOTFOUND) {
            Entry entry;
            auto disk_code = persistentStore.get(dbKey, key, entry);
            if (disk_code == KEY_NOTFOUND) {
                return KEY_NOTFOUND;
            } else {
                data = *entry.getData();
                put(entry);
            }
        };
        return 0;
    }

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet);

    int sync();

    int erase(unsigned int dbKey, int64 &key);

private:

    InstantStore instantStore;
    PersistentStore persistentStore;

};


#endif /* MultiLevelCache_hpp */
