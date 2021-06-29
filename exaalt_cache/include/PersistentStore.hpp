//
//  LocalStore.hpp
//  ParSplice-refact
//
//  Created by Danny Perez on 1/9/17.
//  Copyright © 2017 dp. All rights reserved.
//

#ifndef PersistentStore_hpp
#define PersistentStore_hpp

#include <iostream>
#include <filesystem>
#include <unordered_map>
#include "ExaaltDef.hpp"
#include "StoreEntry.hpp"

namespace fs = std::filesystem;

class PersistentStore {
public:
//    PersistentStore() = default;

    ~PersistentStore() {
        offsets.flush();
        io.flush();
        offsets.close();
        io.close();
    };

    int initialize(const std::string &path) {

        std::string dbName = path + ".db";
        std::string offsetsName = path + ".offsets";
        fs::create_directories(fs::path(dbName).parent_path().string());

//        bool restore = fs::exists(dbName) && fs::exists(offsetsName);
        bool restore = false;
        if (!restore) {
            fs::remove(fs::path(dbName));
            fs::remove(fs::path(offsetsName));
        }
        offsets.open(offsetsName, std::fstream::in | std::fstream::out | std::fstream::app);
        io.open(dbName, std::fstream::in | std::fstream::out | std::fstream::app | std::fstream::binary);

//        std::cout << "STORE INITIALIZE: " << offsets.is_open() << " " << io.is_open() << std::endl;

        if (restore) {
            std::cout << "RESTORING DATABASE" << std::endl;
            //read in the offsets
            int dbk;
            int64 k;
            long int off;
            long int s;
            offsets.seekg(0, offsets.beg);
            while (offsets >> dbk >> k >> off >> s) {
                //std::cout<<"RESTORING DB: "<<dbk<<" "<<k<<" "<<off<<" "<<s<<std::endl;
                locations[std::make_pair(dbk, k)] = std::pair<std::streampos, std::streamsize>(off, s);
            }
            offsets.close();
            offsets.open(offsetsName, std::fstream::in | std::fstream::out | std::fstream::app);
        }
        return 0;
    };

    bool contains(unsigned int dbKey, int64 &key) {
        return contains(std::make_pair(dbKey, key));
    };

    bool contains(KeyPair keys) {
        return locations.count(keys) > 0;
    };

    int put(const Entry &entry) {
        auto k = entry.getKeys();
        if (locations.count(k) == 0) {
            //insert the item

            //move to the end of the stream
            io.seekp(0, io.end);

            locations[k] = std::pair<std::streampos, std::streamsize>(io.tellp(), entry.size());
            offsets.seekp(0, offsets.end);
            offsets << k.first << " " << k.second << " " << locations[k].first << " " << locations[k].second << " " << std::flush;

            //std::cout<<"PUT-IO "<<data.size()<<std::endl;
            auto data = entry.getData();
            io.write(data->data(), data->size());
            auto index = entry.getHuffmanIndex();
            io.write((char *) &index, sizeof(index));

            io.flush();
            //std::cout<<"PUT-sk"<<std::endl;

        }
        return 0;
    };

    int put(unsigned int dbKey, int64 &key, const RawDataVector &data) {
        return put(Entry(std::make_pair(dbKey, key), data));
    };

    int get(unsigned int dbKey, int64 &key, RawDataVector &data) {
        auto k = std::make_pair(dbKey, key);
        if (locations.count(k) != 0) {
            auto loc = locations[k];
            io.seekg(loc.first, io.beg);
            int64 index;
            data.resize(loc.second);
            io.read(data.data(), loc.second - sizeof(index));
            io.read((char *) &index, sizeof(index));
            io.sync();
        } else {
            return KEY_NOTFOUND;
        }
        return 0;
    }

    int get(unsigned int dbKey, int64 &key, Entry &entry) {
        auto k = std::make_pair(dbKey, key);
        if (locations.count(k) != 0) {
            auto loc = locations[k];
            io.seekg(loc.first, io.beg);
            auto data = std::make_shared<RawDataVector>(loc.second);
            int64 index;
            io.read(data->data(), loc.second - sizeof(index));
            io.read((char *) &index, sizeof(index));
            entry.setHuffmanIndex(index);
            entry.setData(std::move(data));
            entry.setDirtyBit(false);
            entry.setKeys(std::make_pair(dbKey, key));
            io.sync();
        } else {
            return KEY_NOTFOUND;
        }
        return 0;
    };

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet) { return 0; }

    int sync() {
        offsets.flush();
        io.flush();
        //std::cout<<"SYNC"<<std::endl;
        return 0;
    };

    virtual int erase(unsigned int dbKey, int64 &key) { return 0; }

    void print() {
        std::cout << "PersistentStore num of blocks = " << locations.size() << std::endl;
    }

private:
    std::unordered_map<KeyPair, std::pair<std::streampos, std::streamsize>, pair_hash> locations;
    std::fstream io;
    std::fstream offsets;

};

#endif /* PersistentStore_hpp */
