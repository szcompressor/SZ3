#ifndef EXAALT_DEF
#define EXAALT_DEF

#include "compressor/SZZoneCompressor.hpp"
#include "predictor/LorenzoPredictor.hpp"

#define DBKEY_NOTFOUND 1
#define KEY_NOTFOUND 2

#define ERROR_BOUND 0.0001
#define BLOCK_SIZE 3137
#define LEVEL1_CAPACITY 10000
#define LEVEL2_CAPACITY 20000
#define LEVEL2_BATCH_SIZE 400

typedef uint64_t int64;
typedef std::vector<char> RawDataVector;
typedef std::pair<unsigned int, int64> KeyPair;

// from boost (functional/hash):
// see http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html template
template <typename T>
inline void hash_combine(std::size_t &seed, const T &val) {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
// auxiliary generic functions to create a hash value using a seed
template <typename T> inline void hash_val(std::size_t &seed, const T &val) {
    hash_combine(seed, val);
}
template <typename T, typename... Types>
inline void hash_val(std::size_t &seed, const T &val, const Types &... args) {
    hash_combine(seed, val);
    hash_val(seed, args...);
}

template <typename... Types>
inline std::size_t hash_val(const Types &... args) {
    std::size_t seed = 0;
    hash_val(seed, args...);
    return seed;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        return hash_val(p.first, p.second);
    }
};

#endif
