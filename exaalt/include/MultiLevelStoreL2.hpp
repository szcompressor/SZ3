//
//  MultiLevelStore.hpp
//  ParSplice-refact
//
//  Created on 9/22/20.
//  Copyright © 2020 dp. All rights reserved.
//

#ifndef MultiLevelStoreL2_hpp
#define MultiLevelStoreL2_hpp

#include "ExaaltDef.hpp"
#include "InstantStore.hpp"
#include "PersistentStore.hpp"
#include "CompressedStore.hpp"
#include "utils/Timer.hpp"

template<class T>
class MultiLevelStoreL2 {
public:
    MultiLevelStoreL2() = default;

    ~MultiLevelStoreL2() {
        for (auto &e: instantStore.getAll()) {
            compressedStore.put(e);
        }
        compressedStore.compressBuffer();
        for (auto &e : compressedStore.getAllData()) {
            persistentDataStore.put(e);
        }
        for (auto &e : compressedStore.getAllHuffman()) {
            persistentHuffmanStore.put(e);
        }
    };

    int initialize(const std::string &path, size_t l1_capacity, size_t l2_capacity,
                   size_t blockSize_, double errorBound_, size_t batchSize_) {
        instantStore.initialize(l1_capacity);
        compressedStore.initialize("", l2_capacity, blockSize_, errorBound_, batchSize_);
        persistentDataStore.initialize(path + "/exaalt_multilevel_data");
        persistentHuffmanStore.initialize(path + "/exaalt_multilevel_huffman");
        return 0;
    }

    void print() {
        instantStore.print();
        compressedStore.print();
        persistentDataStore.print();
        printf("Get Count: instant = %lu , compressed = %lu , persistent = %lu\n", l1_count, l2_count, l3_count);
        printf("Get Time: Compressed = %.2f , persistent = %.2f\n", l2_time, l3_time);
    }

    int put(Entry &entry) {
        auto l1_out = instantStore.push(entry);
        if (l1_out.isValid() && l1_out.isDirty()) {
            auto l2_out = compressedStore.push(l1_out);
            if (l2_out.first.isValid() && l2_out.first.isDirty()) {
                persistentDataStore.put(l2_out.first);
            }
            if (l2_out.second.isValid() && l2_out.second.isDirty()) {
                persistentHuffmanStore.put(l2_out.second);
            }
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
            SZ::Timer timer;
            timer.start();
            auto l2_code = compressedStore.get(dbKey, key, entry);
            l2_time += timer.stop();
            if (l2_code == KEY_NOTFOUND) {
                l3_count++;
                timer.start();
                Entry compressedEntry;
                auto disk_code = persistentDataStore.get(dbKey, key, compressedEntry);
                l3_time += timer.stop();
                if (disk_code == KEY_NOTFOUND) {
                    return KEY_NOTFOUND;
                } else {
                    timer.start();
                    std::pair<Entry, Entry> l2_out;
                    if (!compressedStore.containsHuffman(compressedEntry.getHuffmanIndex())) {
                        Entry huffmanEntry;
                        int64 huffmanIndex = compressedEntry.getHuffmanIndex();
                        persistentHuffmanStore.get(0, huffmanIndex, huffmanEntry);
                        l2_out = compressedStore.push_compressed(compressedEntry, huffmanEntry);
                    } else {
                        l2_out = compressedStore.push_compressed(compressedEntry);
                    }

                    if (l2_out.first.isValid() && l2_out.first.isDirty()) {
                        persistentDataStore.put(l2_out.first);
                    }
                    if (l2_out.second.isValid() && l2_out.second.isDirty()) {
                        persistentHuffmanStore.put(l2_out.second);
                    }

                    compressedStore.get(dbKey, key, entry);
                    l2_time += timer.stop();
                }
            }
            l2_count++;
            data = *entry.getData();
            put(entry);
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
    CompressedStore<T, InstantStore> compressedStore;
    PersistentStore persistentDataStore;
    PersistentStore persistentHuffmanStore;
    size_t l1_count = 0;
    size_t l2_count = 0;
    size_t l3_count = 0;
    double l2_time = 0;
    double l3_time = 0;
};


#endif /* MultiLevelCache_hpp */
