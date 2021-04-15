//
//  LocalStore.hpp
//  ParSplice-refact
//
//  Created by Danny Perez on 1/9/17.
//  Copyright © 2017 dp. All rights reserved.
//

#ifndef LocalStore_hpp
#define LocalStore_hpp

#include <pthread.h>
#include <unistd.h>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <list>
#include "compressor/SZZoneCompressor.hpp"
#include "predictor/LorenzoPredictor.hpp"


#define DBKEY_NOTFOUND 1
#define KEY_NOTFOUND 2

typedef uint64_t Label;
typedef char RawData;
typedef std::vector<RawData> RawDataVector;

class AbstractDataStore {
public:
    virtual int put(unsigned int dbKey, Label &key, RawDataVector &data) = 0;

    virtual int get(unsigned int dbKey, Label &key, RawDataVector &data) = 0;

    virtual unsigned int count(unsigned int dbKey, Label &key) = 0;

    virtual int erase(unsigned int dbKey, Label &key) = 0;
};

//abstract base class for the local data store
class AbstractLocalDataStore : public AbstractDataStore {
public:
    AbstractLocalDataStore();

    virtual ~AbstractLocalDataStore() {
    };

    virtual int put(unsigned int dbKey, Label &key, RawDataVector &data) = 0;

    virtual int get(unsigned int dbKey, Label &key, RawDataVector &data) = 0;

    virtual int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet) = 0;

    virtual int initialize(std::string homeDir, std::string baseName) = 0;

    virtual int sync() = 0;

    virtual int erase(unsigned int dbKey, Label &key) = 0;

    void setMaximumSize(unsigned long maxSize_);

    virtual unsigned int count(unsigned int dbKey, Label &key);

    virtual std::set<Label> availableKeys(unsigned int dbKey);

    virtual std::set<unsigned int> availableDbKeys();

    virtual std::unordered_map<unsigned int, std::set<Label> > getStoredKeys();

//virtual std::multimap< std::pair<unsigned int, Label>, RawDataVector> purge();
    void purge();

    virtual void touch(unsigned int dbKey, Label &key);

protected:

    std::unordered_map<unsigned int, std::set<Label> > storedKeys;
//MRU< std::pair<unsigned int, Label> > mru;
    unsigned long maxSize;
    unsigned long currentSize;

    std::unordered_map<unsigned int, std::pair<bool, bool> > dbAttributes;
};


class PersistentLocalStore : public AbstractLocalDataStore {
public:
    PersistentLocalStore();

    ~PersistentLocalStore();

    virtual unsigned int count(unsigned int dbKey, Label &key);


    int initialize(std::string homeDir, std::string baseName);

    virtual int put(unsigned int dbKey, uint64_t &key, RawDataVector &data);

    virtual int get(unsigned int dbKey, uint64_t &key, RawDataVector &data);

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet);

    int sync();

    virtual int erase(unsigned int dbKey, uint64_t &key);

private:
    std::map<std::pair<int, uint64_t>, std::pair<std::streampos, std::streamsize> > locations;
    std::fstream io;
    std::fstream offsets;

};

class CompressedPersistentLocalStore : public AbstractLocalDataStore {
public:
    CompressedPersistentLocalStore();

    ~CompressedPersistentLocalStore();

    int initialize(std::string homeDir, std::string baseName);

    int initialize(std::string homeDir, std::string baseName, size_t block_size, double realPrecision, size_t maxBufferSize);

    virtual int put(unsigned int dbKey, uint64_t &key, RawDataVector &data);

    virtual int get(unsigned int dbKey, uint64_t &key, RawDataVector &data);

    int put_compressible(unsigned int dbKey, uint64_t &key, std::vector<double> &data);

    int get_compressible(unsigned int dbKey, uint64_t &key, std::vector<double> &data);

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet);

    int sync();

    virtual int erase(unsigned int dbKey, uint64_t &key);

    void compressBuffer();

private:

    std::map<std::pair<unsigned int, uint64_t>, std::vector<double> > buffer;
    double realPrecision;
    size_t block_size;
    size_t maxBufferSize;

    PersistentLocalStore compressedStore;
    PersistentLocalStore huffmanStore;
    std::map<std::pair<unsigned int, uint64_t>, Label> huffmanIndexStore;
    Label treeIndex = 0;
    std::map<Label, std::unique_ptr<SZ::SZ_Zone_Compressor<double, 1, SZ::LorenzoPredictor<double, 1, 1>,
            SZ::LinearQuantizer<double>, SZ::HuffmanEncoder<int>>>> szmap;

};

#endif /* LocalStore_hpp */
