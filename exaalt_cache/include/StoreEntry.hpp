#ifndef INSTANT_STORE_ENTRY
#define INSTANT_STORE_ENTRY

#include <utility>


#include "ExaaltDef.hpp"

class Entry {
public:
    Entry() = default;

    Entry(const KeyPair keys_, std::shared_ptr<RawDataVector> data_) : keys(keys_), data(std::move(data_)) {
    }

    Entry(const KeyPair keys_, u_char *begin, u_char *end) : keys(keys_), data(std::make_shared<RawDataVector>(begin, end)) {
    }

    Entry(const KeyPair keys_, RawDataVector data_) : keys(keys_), data(std::make_shared<RawDataVector>(data_)) {
    }

    ~Entry() = default;


    void setData(std::shared_ptr<RawDataVector> data_) {
        this->data = std::move(data_);
    }

    std::shared_ptr<RawDataVector> initData(size_t size) {
        data = std::make_shared<RawDataVector>(size);
        return data;
    }

    std::shared_ptr<RawDataVector> getData() const {
        return data;
    }

    void setHuffmanIndex(int64 huffman) {
        this->huffmanIndex = huffman;
    }

    int64 getHuffmanIndex() const { return huffmanIndex; }

    bool isDirty() const { return dirtyBit; };

    void setDirtyBit(bool dirty) {
        dirtyBit = dirty;
    }

    size_t size() const {
        return data->size() + sizeof(huffmanIndex);
    }

    void setKeys(KeyPair key_) {
        keys = key_;
    }

    KeyPair getKeys() const {
        return keys;
    }

    bool isValid() {
        return data != nullptr;
    }

private:
    int64 huffmanIndex = 0;
    bool dirtyBit = true;
    std::shared_ptr<RawDataVector> data;
    KeyPair keys;
};


#endif
