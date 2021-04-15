//
// Created by Kai Zhao on 10/28/20.
//
#ifndef SZ_COMPRESSEDSTORE_HPP
#define SZ_COMPRESSEDSTORE_HPP

#include <unordered_map>
#include <list>
#include <vector>
#include <string>
#include "ExaaltDef.hpp"
#include "StoreEntry.hpp"
#include "compressor/SZZoneCompressor.hpp"
#include "predictor/LorenzoPredictor.hpp"

template<class T, class Store>
class CompressedStore {
    typedef SZ::SZZoneCompressor<T, 1, SZ::LorenzoPredictor<T, 1, 1>,
            SZ::LinearQuantizer<T>, SZ::HuffmanEncoder<int>, SZ::Lossless_zstd> SZ_Compressor;
public:
//    CompressedStore() = default;

//    ~CompressedStore() = default;

    bool contains(unsigned int dbKey, int64 &key) {
        return contains(std::make_pair(dbKey, key));
    }

    bool contains(KeyPair keys) {
        return (dataStore.contains(keys) || buffer.count(keys) > 0);
    }

    bool containsHuffman(int64 index) {
        return huffmanCacheIters.count(index) > 0 || huffmanStore.contains(0, index);
    }

    int
    initialize(const std::string &path, size_t capacity_, size_t blockSize_, double errorBound_, size_t batchSize_) {
        dataStore.initialize(path + "/exaalt_compressed_data");
        huffmanStore.initialize(path + "/exaalt_compressed_huffman");
        this->batchSize = batchSize_;
        this->blockSize = blockSize_;
        this->capacity = capacity_;
        this->errorBound = errorBound_;
        return 0;
    }

    std::list<Entry> &getAllData() {
        return dataStore.getAll();
    }

    std::list<Entry> &getAllHuffman() {
        return huffmanStore.getAll();
    }

    std::pair<Entry, Entry> push_compressed(const Entry &entry) {
        if (!dataStore.contains(entry.getKeys())) {
            dataStore.put(entry);
        }
        huffmanCount[entry.getHuffmanIndex()]++;
        return (dataStore.size() + buffer.size() >= capacity) ? pop() : std::make_pair(Entry(), Entry());
    }

    std::pair<Entry, Entry> push_compressed(const Entry &entry, const Entry &huffmanEntry) {
        if (!huffmanStore.contains(huffmanEntry.getKeys())) {
            huffmanStore.put(huffmanEntry);
            huffmanCount[entry.getHuffmanIndex()] = 0;
        }
        return push_compressed(entry);
    }

    std::pair<Entry, Entry> push(const Entry &entry) {
        if (!contains(entry.getKeys())) {
            buffer[entry.getKeys()] = entry;
        }
        if (buffer.size() >= batchSize) {
            compressBuffer();
        }
        return (dataStore.size() >= capacity) ? pop() : std::make_pair(Entry(), Entry());
    }

    std::pair<Entry, Entry> pop() {
        Entry delDataEntry, delHuffmanEntry;
        delDataEntry = dataStore.pop();
        auto huffman = delDataEntry.getHuffmanIndex();
        auto it = huffmanCount.find(huffman);
        auto removeHuffman = (it == huffmanCount.end());
        if (it != huffmanCount.end()) {
            if (--(it->second) == 0) {
                huffmanCount.erase(it);
                removeHuffman = true;
            }
        }
        if (removeHuffman && huffmanStore.contains(0, huffman)) {
            huffmanStore.get(0, huffman, delHuffmanEntry);
            huffmanStore.erase(0, huffman);
        }
        return std::make_pair(delDataEntry, delHuffmanEntry);
    }

    int put(unsigned int dbKey, int64 &key, RawDataVector &data) {
        return put(Entry(std::make_pair(dbKey, key), data));
    }

    int put(const Entry &entry) {
        auto k = entry.getKeys();
        if (!contains(k)) {
            buffer[k] = entry;
        }
        if (buffer.size() >= batchSize) {
            compressBuffer();
        }
        return 0;
    }

    int get(unsigned int dbKey, int64 &key, RawDataVector &data) {
        Entry entry;
        auto result = get(dbKey, key, entry);
        if (result == 0) {
            data = *entry.getData();
        }
        return result;
    }

    int get(unsigned int dbKey, int64 &key, Entry &entry) {

        auto k = std::make_pair(dbKey, key);
        if (buffer.count(k) > 0) {
            entry = buffer[k];
            return 0;
        }

        Entry compressedEntry;
        int check = dataStore.get(dbKey, key, compressedEntry);
        if (check == 0) {
            //found the key
            int64 huffman = compressedEntry.getHuffmanIndex();

            std::shared_ptr<SZ_Compressor> sz;
            auto it = huffmanCacheIters.find(huffman);
            if (it == huffmanCacheIters.end()) {
                auto conf = SZ::Config<T, 1>(errorBound, std::array<size_t, 1>{blockSize});
                sz = std::make_shared<SZ_Compressor>(conf, SZ::LorenzoPredictor<T, 1, 1>(errorBound),
                                                     SZ::LinearQuantizer<T>(errorBound), SZ::HuffmanEncoder<int>(),
                                                     SZ::Lossless_zstd());
                RawDataVector tdata;
                huffmanStore.get(0, huffman, tdata);
                sz->decompress_encoder((unsigned char *) tdata.data(), tdata.size());
                huffmanCache.emplace_front(huffman, sz);
                huffmanCacheIters[huffman] = huffmanCache.begin();

                if (huffmanCacheIters.size() > huffmanCacheCapacity) {
                    auto delIndex = huffmanCache.back().first;
                    huffmanCacheIters.erase(delIndex);
                    huffmanCache.pop_back();
                }
            } else {
                sz = it->second->second;
                huffmanCache.splice(huffmanCache.begin(), huffmanCache, it->second);
            }
            auto compressed = compressedEntry.getData();
            auto decompressed = entry.initData(blockSize * sizeof(T));

            sz->decompress_zone((unsigned char *) compressed->data(), compressed->size(), (T *) decompressed->data());

            //std::cout << "compression ratio = " << (3137*8.0/ddata->size()) << "\n";
            return 0;
        } else {
            return KEY_NOTFOUND;
        }
    }

    int createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet) { return 0; }

    int sync() {
        dataStore.sync();
        huffmanStore.sync();
        return 0;
    }

    virtual int erase(unsigned int dbKey, int64 &key) { return 0; };

    int getBufferSize() { return buffer.size(); };

    void compressBuffer() {
        std::vector<T> rdata(buffer.size() * blockSize); //the double data to be compressed

        //construct compression buffer
        int i = 0;
        for (auto it = buffer.begin(); it != buffer.end(); it++, i++) {
            memcpy(&rdata[i * blockSize], it->second.getData()->data(), blockSize * sizeof(T));
        }

        std::vector<size_t> compressed_size(buffer.size() + 1, 0); //to store the size of each zone

        auto conf = SZ::Config<T, 1>(errorBound, std::array<size_t, 1>{buffer.size() * blockSize});
        auto sz = SZ::make_sz_zone_compressor(conf, SZ::LorenzoPredictor<T, 1, 1>(errorBound),
                                              SZ::LinearQuantizer<T>(errorBound, 128), SZ::HuffmanEncoder<int>());
//    std::cout << rdata.size() << "," << conf.num<<std::endl;
        auto szresult = sz.compress(rdata.data(), compressed_size, buffer.size());
        //save the tree
        RawDataVector tree(szresult[buffer.size()], szresult[buffer.size()] + compressed_size[buffer.size()]);
        huffmanStore.put(0, huffmanIndex, tree);
        huffmanCount[huffmanIndex] = 0;

        i = 0;
        for (auto it = buffer.begin(); it != buffer.end(); it++, i++) {
            Entry entry(it->first, szresult[i], szresult[i] + compressed_size[i]);
            entry.setHuffmanIndex(huffmanIndex);
            dataStore.put(entry);
            huffmanCount[huffmanIndex]++;
        }
        for (i = 0; i < buffer.size() + 1; i++) {
            delete[] szresult[i];
        }

        huffmanIndex++;
        buffer.clear();
    }

    void print() {
        std::cout << "CompressedStore num of blocks = " << dataStore.size()
                  << " buffer size = " << buffer.size()
                  << " huffman size = " << huffmanStore.size()
                  << std::endl;
    }

private:
    double errorBound;
    size_t blockSize;
    size_t batchSize;
    size_t capacity = 10000;

    Store dataStore; //compressed data
    Store huffmanStore; // compressed huffman
    std::unordered_map<int64, int64> huffmanCount;
    int64 huffmanIndex = 0;

    size_t huffmanCacheCapacity = 10000;
    typedef std::list<std::pair<int64, std::shared_ptr<SZ_Compressor>>> HuffmanCacheType;
    std::unordered_map<int64, typename HuffmanCacheType::iterator> huffmanCacheIters;
    HuffmanCacheType huffmanCache;

    std::unordered_map<KeyPair, Entry, pair_hash> buffer;

};

#endif //SZ_COMPRESSEDSTORE_HPP
