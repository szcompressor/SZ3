// Decompresses chunks written by the SZ3 3.1.x HDF5 filter, which hdf5plugin 7.1.0 and earlier ship, with the code of
// SZ3 3.1.7 that H5Z_SZ3_READ_V3_1 fetches. This is the only file that includes it: 3.1.7 lives in namespace SZ and its
// global SZ_* functions take SZ::Config, so nothing here collides with SZ3.
#include <cstdint>  // 3.1.7's Config.hpp uses uint8_t without including it
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>

#include "SZ3/api/sz.hpp"

template <class T>
static void* decompress(char* cmpData, size_t cmpSize, size_t* outBytes) {
    // 3.1.x chunks end with the Config and its size, which give the element count to allocate
    int confSize = 0;
    memcpy(&confSize, cmpData + cmpSize - sizeof(int), sizeof(int));
    const SZ::uchar* confPos = reinterpret_cast<SZ::uchar*>(cmpData) + cmpSize - sizeof(int) - confSize;
    SZ::Config conf;
    conf.load(confPos);
    // HDF5 frees the chunk with free()
    std::unique_ptr<T, decltype(&free)> out(static_cast<T*>(malloc(conf.num * sizeof(T))), &free);
    if (!out) throw std::bad_alloc();
    T* decData = out.get();
    SZ_decompress<T>(conf, cmpData, cmpSize, decData);
    *outBytes = conf.num * sizeof(T);
    return out.release();
}

// dataType is the 3.1.x filter's code: 0 float, 1 double, 2 uint8, 3 int8, 4 uint16, 5 int16, 6 uint32, 7 int32,
// 8 uint64, 9 int64
void* H5Z_SZ3_decompress_v3_1(int dataType, char* cmpData, size_t cmpSize, size_t* outBytes) {
    switch (dataType) {
        case 0:
            return decompress<float>(cmpData, cmpSize, outBytes);
        case 1:
            return decompress<double>(cmpData, cmpSize, outBytes);
        case 2:
            return decompress<uint8_t>(cmpData, cmpSize, outBytes);
        case 3:
            return decompress<int8_t>(cmpData, cmpSize, outBytes);
        case 4:
            return decompress<uint16_t>(cmpData, cmpSize, outBytes);
        case 5:
            return decompress<int16_t>(cmpData, cmpSize, outBytes);
        case 6:
            return decompress<uint32_t>(cmpData, cmpSize, outBytes);
        case 7:
            return decompress<int32_t>(cmpData, cmpSize, outBytes);
        case 8:
            return decompress<uint64_t>(cmpData, cmpSize, outBytes);
        case 9:
            return decompress<int64_t>(cmpData, cmpSize, outBytes);
        default:
            return nullptr;
    }
}
