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
#include "utils/Timer.hpp"

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
        printf("Get Count: instant = %lu , persistent = %lu\n", l1_count, l2_count);
        printf("Get Time: persistent = %.2f\n", l2_time);

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
            SZ::Timer timer;
            timer.start();
            Entry entry;
            auto disk_code = persistentStore.get(dbKey, key, entry);
            l2_time += timer.stop();
            if (disk_code == KEY_NOTFOUND) {
                return KEY_NOTFOUND;
            } else {
                l2_count++;
                data = *entry.getData();
                put(entry);
            }
        } else {
            l1_count++;
        };
        return 0;
    }

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet);

    int sync();

    int erase(unsigned int dbKey, int64 &key);

private:

    InstantStore instantStore;
    PersistentStore persistentStore;
    size_t l1_count = 0;
    size_t l2_count = 0;
    double l2_time = 0;
};


#endif /* MultiLevelCache_hpp */
