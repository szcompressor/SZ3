#ifndef _SZ_HUFFMAN_ENCODER_V2_HPP
#define _SZ_HUFFMAN_ENCODER_V2_HPP

#include <cassert>
#include <iostream>
#include <map>
#include <queue>
#include <stack>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/ska_hash/unordered_map.hpp"

namespace SZ3 {

template <class T>
class HuffmanEncoderV2 : public concepts::EncoderInterface<T> {
   private:
    class Node {
       public:
        Node(T c_ = 0, Node *lp = nullptr, Node *rp = nullptr) {
            c = c_;
            p[0] = lp;
            p[1] = rp;
        }

        T c;
        Node *p[2];

        inline uchar isLeaf() { return p[0] == nullptr; }
    };

    class HuffmanTree {
       private:
        uchar _constructed = 0;

        uchar len = 0;
        int vec = 0;

        class cmp {
           public:
            bool operator()(const std::pair<int, size_t> &u, const std::pair<int, size_t> &v) {
                return u.second == v.second ? u.first > v.first : u.second > v.second;
            }
        };

       public:
        void dfs_mp(Node *u) {
            if (u->isLeaf()) {
                mplen[u->c] = len;
                mpcode[u->c] = vec;

                limit = std::max(limit, len);

                return;
            }

            ++len;
            dfs_mp(u->p[0]);
            --len;

            vec ^= 1 << len++;
            dfs_mp(u->p[1]);
            vec ^= 1 << --len;
        }

        void dfs_vec(Node *u) {
            if (u->isLeaf()) {
                veclen[u->c] = len;
                veccode[u->c] = vec;

                limit = std::max(limit, len);

                return;
            }

            ++len;
            dfs_vec(u->p[0]);
            --len;

            vec ^= 1 << len++;
            dfs_vec(u->p[1]);
            vec ^= 1 << --len;
        }

        uchar usemp;
        // 0 : vec
        // 1 : mp

        std::vector<uchar> veclen;
        std::vector<int> veccode;
        ska::unordered_map<size_t, uchar> mplen;
        ska::unordered_map<size_t, int> mpcode;
        //            std::unordered_map<size_t,uchar> mplen;
        //            std::unordered_map<size_t,int> mpcode;
        //            std::map<size_t,uchar> mplen;
        //            std::map<size_t,int> mpcode;

        T offset;
        // minimum bits for T
        uchar mbft;
        uchar limit;

        void init() {
            _constructed = 0;
            ht.clear();
            veclen.clear();
            veccode.clear();
            mplen.clear();
            mpcode.clear();
            vecfreq.clear();
            mpfreq.clear();

            offset = 0;
            mbft = 0;
            root = 0;
            n = 0;
            maxval = 0;
            limit = 0;
        }

        HuffmanTree() { init(); }

        int root;
        int n;
        T maxval;
        std::vector<Node> ht;
        std::vector<size_t> vecfreq;
        ska::unordered_map<T, size_t> mpfreq;
        //            std::unordered_map<T,size_t> mpfreq;
        //            std::map<T,size_t> mpfreq;

        void addElementInMap(T c, size_t freqc) {
            assert(!_constructed);

            ht.push_back(Node(c));
            mpfreq[c] = freqc;
            ++n;
        }

        void addElementInVector(T c, size_t freqc) {
            assert(!_constructed);

            ht.push_back(Node(c));
            vecfreq[c] = freqc;
            ++n;
        }

        void constructHuffmanTree() {
            assert(!_constructed);

            if (n == 1 || maxval == 1) {
                mbft = 1;
                maxval = 1;
                offset += ht[0].c;
                ht[0].c = 0;
                ht.push_back(Node(0, &ht[0], nullptr));
                if (usemp) {
                    mplen[0] = 1;
                    mpcode[0] = 0;
                } else {
                    veclen[0] = 1;
                    veccode[0] = 0;
                }
                limit = 1;
                setConstructed();
                return;
            }

            // Timer timer(true);

            mbft = 1;
            while ((1llu << mbft) < maxval) ++mbft;

            std::priority_queue<std::pair<int, size_t>, std::vector<std::pair<int, size_t>>, cmp> q;

            if (usemp) {
                for (int i = 0; i < ht.size(); i++) {
                    q.push({i, mpfreq[ht[i].c]});
                }
            } else {
                for (int i = 0; i < ht.size(); i++) {
                    q.push({i, vecfreq[ht[i].c]});
                }
            }

            while (q.size() > 1) {
                int u = q.top().first;
                size_t freq_u = q.top().second;
                q.pop();
                int v = q.top().first;
                size_t freq_v = q.top().second;
                q.pop();

                ht.push_back(Node(0, &ht[u], &ht[v]));

                q.push({ht.size() - 1, freq_u + freq_v});
            }

            root = ht.size() - 1;

            if (usemp)
                dfs_mp(&ht[root]);
            else
                dfs_vec(&ht[root]);

            setConstructed();

            // timer.stop("construct huffman tree");
        }

        uchar isConstructed() { return _constructed; }

        void setConstructed() { _constructed = 1; }
    };

    HuffmanTree tree;

   public:
    void preprocess_encode(const T *const bins, size_t num_bin, int stateNum, uchar flag = 0x00) {
        // Timer timer(true);

        tree.init();

        T __minval, __maxval;

        if (stateNum == 0) {
            __minval = *bins;
            __maxval = *bins;
            for (int i = 1; i < num_bin; i++) {
                __minval = std::min(__minval, *(bins + i));
                __maxval = std::max(__maxval, *(bins + i));
            }
        } else {
            __minval = 0;
            __maxval = stateNum - 1;
        }

        tree.offset = __minval;
        tree.maxval = __maxval - __minval + 1;

        switch ((flag & 0xc0) >> 6) {
            case 0: {
                tree.usemp = tree.maxval >= (1 << 12) && num_bin < 2 * __maxval || tree.maxval >= (1 << 28) ? 1 : 0;
                break;
            }
            case 1: {
                tree.usemp = 1;
                break;
            }
            case 2: {
                tree.usemp = 0;
                break;
            }
            case 3: {
                tree.usemp = 0;
                break;
            }
        }

        if ((flag & 0x01) == 0x01) {
            tree.mbft = 1;
            while ((1llu << tree.mbft) < tree.maxval) ++tree.mbft;
            tree.n = 0;
            tree.setConstructed();
            return;
        }

        if (tree.usemp) {
            //                tree.mpfreq.reserve(num_bin);
            //                tree.mplen.reserve(num_bin);
            //                tree.mpcode.reserve(num_bin);

            //                ska::unordered_map<T,size_t> freq;
            //                std::unordered_map<T,size_t> freq;
            std::map<T, size_t> freq;
            //                freq.reserve(num_bin);

            if (tree.offset == 0) {
                for (int i = 0; i < num_bin; i++) {
                    ++freq[bins[i]];
                }
            } else {
                for (int i = 0; i < num_bin; i++) {
                    ++freq[bins[i] - tree.offset];
                }
            }

            tree.ht.reserve(freq.size() << 1);

            for (auto it : freq) {
                tree.addElementInMap(it.first, it.second);
            }
        } else {
            tree.vecfreq.resize(tree.maxval);
            tree.veclen.resize(tree.maxval);
            tree.veccode.resize(tree.maxval);

            std::vector<size_t> freq(tree.maxval);

            if (tree.offset == 0) {
                for (int i = 0; i < num_bin; i++) {
                    ++freq[bins[i]];
                }
            } else {
                for (int i = 0; i < num_bin; i++) {
                    ++freq[bins[i] - tree.offset];
                }
            }

            tree.ht.reserve(freq.size() << 1);

            for (int i = 0; i < tree.maxval; i++) {
                if (freq[i]) tree.addElementInVector(i, freq[i]);
            }
        }

        // printf("begins to construct huffman tree\n");

        tree.constructHuffmanTree();

        // timer.stop("preprocess_encode");
    }

    void preprocess_encode(const std::vector<T> &bins, int stateNum) override {
        preprocess_encode(bins.data(), bins.size(), stateNum);
    }

    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        return encode(bins.data(), bins.size(), bytes);
    }

    size_t encode(const T *bins, size_t num_bin, uchar *&bytes) {
        if (tree.maxval == 1) {
            int64ToBytes_bigEndian(bytes, num_bin ^ 0x1234abcd);
            bytes += 8;
            return 8;
        }

        // Timer timer(true);

        assert(tree.isConstructed());

        uchar *head = bytes;
        bytes += 8;

        size_t len = 0;

        uchar mask = 0;
        uchar index = 0;

        if (tree.n == 0) {
            if (tree.offset == 0) {
                for (size_t i = 0; i < num_bin; i++) {
                    writeBytes(bytes, bins[i], tree.mbft, mask, index);
                }
            } else {
                for (size_t i = 0; i < num_bin; i++) {
                    writeBytes(bytes, bins[i] - tree.offset, tree.mbft, mask, index);
                }
            }
            writeBytesClearMask(bytes, mask, index);
            len = tree.mbft * num_bin;
            int64ToBytes_bigEndian(head, len ^ 0x1234abcd);
            return bytes - head;
        }

        if (tree.offset == 0) {
            if (tree.usemp) {
                for (size_t i = 0; i < num_bin; i++) {
                    const T &it = bins[i];
                    const uchar &len_i = tree.mplen[it];
                    const int &code_i = tree.mpcode[it];
                    len += len_i;
                    writeBytes(bytes, code_i, len_i, mask, index);
                }
            } else {
                for (size_t i = 0; i < num_bin; i++) {
                    const T &it = bins[i];
                    const uchar &len_i = tree.veclen[it];
                    const int &code_i = tree.veccode[it];
                    len += len_i;
                    writeBytes(bytes, code_i, len_i, mask, index);
                }
            }
        } else {
            if (tree.usemp) {
                for (size_t i = 0; i < num_bin; i++) {
                    const T &it = bins[i];
                    const uchar &len_i = tree.mplen[it - tree.offset];
                    const int &code_i = tree.mpcode[it - tree.offset];
                    len += len_i;
                    writeBytes(bytes, code_i, len_i, mask, index);
                }
            } else {
                for (size_t i = 0; i < num_bin; i++) {
                    const T &it = bins[i];
                    const uchar &len_i = tree.veclen[it - tree.offset];
                    const int &code_i = tree.veccode[it - tree.offset];
                    len += len_i;
                    writeBytes(bytes, code_i, len_i, mask, index);
                }
            }
        }

        writeBytesClearMask(bytes, mask, index);

        int64ToBytes_bigEndian(head, len ^ 0x1234abcd);

        // timer.stop("encode");

        // printf("code size = %d\n",(int)(bytes-head));

        // Lossless_zstd zstd;
        // size_t compressed_code_size;

        // delete[] zstd.compress(head,bytes-head,compressed_code_size);

        // printf("compressed code size = %d\n",(int)compressed_code_size);

        return bytes - head;
    }

    void postprocess_encode() override {}

    void preprocess_decode() override {}

    std::vector<T> decode(const uchar *&bytes, size_t targetLength) override {
        if (tree.maxval == 1) {
            size_t len = bytesToInt64_bigEndian(bytes) ^ 0x1234abcd;
            bytes += 8;
            //                assert(len==targetLength);

            return std::vector<T>(len, tree.offset);
        }

        // Timer timer(true);

        assert(tree.isConstructed());

        // Node *u = &tree.ht[tree.root];

        size_t len = bytesToInt64_bigEndian(bytes) ^ 0x1234abcd;
        bytes += 8;
        std::vector<T> out(targetLength);
        size_t outLen = 0;

        // For fixed length encoding
        if (tree.n == 0) {
            size_t byteIndex = 0;
            size_t i = 0;
            size_t b;

            T c = 0;
            int j = 0;

            if (tree.offset == 0) {
                for (; i + 8 < len; i += 8, byteIndex++) {
                    b = bytes[byteIndex];

                    c |= ((b) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 1) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 2) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 3) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 4) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 5) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 6) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                    c |= ((b >> 7) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c, c = j = 0;
                }
            } else {
                for (; i + 8 < len; i += 8, byteIndex++) {
                    b = bytes[byteIndex];

                    c |= ((b) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 1) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 2) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 3) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 4) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 5) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 6) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                    c |= ((b >> 7) & 1) << j, ++j;
                    if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
                }
            }

            b = bytes[byteIndex];

            for (size_t k = 0; k < len - i; k++) {
                c |= ((b >> k) & 1) << j, ++j;
                if (j == tree.mbft) out[outLen++] = c + tree.offset, c = j = 0;
            }

            bytes += (len + 7) >> 3;

            return out;
        }

        if (tree.limit > 16) {
            // if huffman tree is large, a cache of huffman codebook is used to increase the performance
            // Reference paper: Xiangyu Zou, Tao Lu, Wen Xia, Xuan Wang, Weizhe Zhang, Haijun Zhang, Sheng Di, Dingwen
            // Tao, and Franck Cappello, "Performance Optimization for Relative-Error-Bounded Lossy Compression on
            // Scientific Data", IEEE Transactions on Parallel and Distributed Systems (IEEE TPDS), 2020. Reference
            // code:
            // https://github.com/szcompressor/SZ/blob/a92658e785c072de1061f549c6cbc6d42d0f7f22/sz/src/Huffman.c#L345

            int maxBits = 16;
            size_t count = 0;
            Node *t = &tree.ht[tree.root];
            Node *n = t;

            size_t tableSize = 1 << maxBits;
            std::vector<T> valueTable(tableSize);
            std::vector<uint8_t> lengthTable(tableSize);
            std::vector<Node *> nodeTable(tableSize);
            size_t j;
            for (size_t i = 0; i < tableSize; i++) {
                n = t;
                j = 0;
                size_t res = i;
                while (!n->isLeaf() && j < maxBits) {
                    n = n->p[res & 0x00000001];
                    res >>= 1;
                    j++;
                }
                if (!n->isLeaf()) {
                    nodeTable[i] = n;
                    valueTable[i] = -1;
                    lengthTable[i] = maxBits;
                } else {
                    valueTable[i] = n->c + tree.offset;
                    lengthTable[i] = j;
                }
            }

            size_t leftBits = 0;
            T currentValue = 0;
            size_t i = 0;

            while (count < targetLength) {
                while (leftBits < maxBits) {
                    currentValue += (bytes[i] << leftBits);
                    leftBits += 8;
                    i++;
                }

                size_t index = currentValue & ((1 << maxBits) - 1);
                T value = valueTable[index];
                if (value != -1) {
                    out[count] = value;
                    int bitLength = lengthTable[index];
                    leftBits -= bitLength;
                    currentValue >>= bitLength;
                    count++;
                } else {
                    int bitLength = lengthTable[index];
                    leftBits -= bitLength;
                    currentValue >>= bitLength;
                    n = nodeTable[index];
                    while (!n->isLeaf()) {
                        if (!leftBits) {
                            currentValue += (bytes[i] << leftBits);
                            leftBits += 8;
                            i++;
                        }
                        n = n->p[(currentValue & 0x01)];
                        leftBits--;
                        currentValue >>= 1;
                    }
                    out[count] = n->c + tree.offset;
                    count++;
                }
            }
        } else {
            // for small huffman tree, use loop unrolling to increase the performance
            // for(int i=0;i<len;){
            //     u=u->p[readBit(bytes,i++)];
            //     if(u->isLeaf()){
            //         out[outLen++]=u->c+tree.offset;
            //         u=&tree.ht[tree.root];
            //     }
            // }

            int byteIndex = 0;
            int i = 0;
            uchar b;
            Node *u = &tree.ht[tree.root];
            auto offset = tree.offset;

            for (; i + 8 < len; i += 8, byteIndex++) {
                b = bytes[byteIndex];

                u = u->p[b & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 1) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 2) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 3) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 4) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 5) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 6) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
                u = u->p[(b >> 7) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + offset;
                    u = &tree.ht[tree.root];
                }
            }

            b = bytes[byteIndex];

            for (int j = 0; j < len - i; j++) {
                u = u->p[(b >> j) & 1];
                if (u->isLeaf()) {
                    out[outLen++] = u->c + tree.offset;
                    u = &tree.ht[tree.root];
                }
            }
        }
        bytes += (len + 7) >> 3;

        // timer.stop("decode");

        return out;
    }

    void postprocess_decode() override {}

    void save(uchar *&c) override {
        //            saveAsCode(c);
        saveAsDFSOrder(c);
    }

    void load(const uchar *&c, size_t &remaining_length) override {
        //            loadAsCode(c,remaining_length);
        loadAsDFSOrder(c, remaining_length);
    }

   private:
    static void writeBytesBit(uchar *&c, uchar val, uchar &mask, uchar &index) {
        assert(val == 0 || val == 1);

        mask |= val << index++;
        if (index == 8) {
            *c++ = mask;
            mask = index = 0;
        }
    }

    static void writeBytes(uchar *&c, T val, uchar len, uchar &mask, uchar &index) {
        assert(len >= 1 && len <= sizeof(T) * 8);

        if (len + index >= 8) {
            mask |= (val & ((1 << (8 - index)) - 1)) << index;
            val >>= 8 - index;
            len -= 8 - index;
            *c++ = mask;
            mask = index = 0;

            while (len >= 8) {
                *c++ = val & (1 << 8) - 1;
                val >>= 8;
                len -= 8;
            }
        }

        mask |= (val & (1 << len) - 1) << index;
        index += len;

        // for(int i=0;i<len;i++){
        //     mask|=(val&1)<<index++;
        //     if(index==8){
        //         *c++=mask;
        //         mask=index=0;
        //     }
        //     val>>=1;
        // }
    }

    static void writeBytesByte(uchar *&c, uchar val) { *c++ = val; }

    static void writeBytesClearMask(uchar *&c, uchar &mask, uchar &index) {
        if (index > 0) {
            *c++ = mask;
            // mask=i=0;
        }
    }

    static uchar readBit(const uchar *const &c, int i) { return ((*(c + (i >> 3))) >> (i & 7)) & 1; }

    void saveAsCode(uchar *&c) {
        // Timer timer(true);

        // uchar *head = c;
        // whether the tree is full binary tree

        uchar &limit = tree.limit;

        std::vector<std::deque<T>> mp(limit + 1);

        if (tree.usemp) {
            for (auto it : tree.mplen) {
                mp[it.second].push_back(it.first);
            }
        } else {
            for (int i = 0; i < tree.maxval; i++) {
                mp[tree.veclen[i]].push_back(i);
            }
        }

        uchar mask = 0;
        uchar index = 0;

        assert(sizeof(T) <= 8);

        if (mp[limit].size() == tree.n) {
            // 00 XXXXXX (mbft)
            if (tree.maxval > 1)
                writeBytesByte(c, tree.mbft);
            else
                writeBytesByte(c, 0x80 | tree.mbft);

            writeBytesByte(c, ((sizeof(T) - 1) << 5) | (limit - 1));

            writeBytes(c, tree.offset, sizeof(T) << 3, mask, index);

            int32ToBytes_bigEndian(c, tree.n);
            c += 4;

            int cnt = mp[limit].size();

            uchar logcnt = 0;
            while (logcnt < 32 && (1 << logcnt) != cnt) ++logcnt;
            assert(logcnt != 32);

            if (tree.n > 1) {
                for (T it : mp[limit]) {
                    writeBytes(c, it, tree.mbft, mask, index);

                    const int code = tree.usemp ? tree.mpcode[it] : tree.veccode[it];

                    writeBytes(c, code, logcnt, mask, index);
                }

                writeBytesClearMask(c, mask, index);
            }

            return;
        }

        writeBytesByte(c, 0x40 | tree.mbft);

        writeBytesByte(c, ((sizeof(T) - 1) << 5) | (limit - 1));

        writeBytes(c, tree.offset, sizeof(T) << 3, mask, index);

        int32ToBytes_bigEndian(c, tree.maxval);
        c += 4;

        for (uchar len = 1; len <= limit; len++) {
            int cnt = mp[len].size();

            writeBytes(c, cnt, len, mask, index);

            if (cnt) {
                for (const T &it : mp[len]) {
                    writeBytes(c, it, tree.mbft, mask, index);

                    const int code = tree.usemp ? tree.mpcode[it] : tree.veccode[it];

                    writeBytes(c, code, len, mask, index);
                }
            }
        }

        writeBytesClearMask(c, mask, index);

        // timer.stop("saveAsCode");

        //            printf("n = %d\n",tree.n);
        //
        //            printf("huffman tree size = %d\n",(int)(c-head));
        //
        //            Lossless_zstd zstd;
        //            size_t compressed_tree_size;
        //
        //            // uchar *compressed_tree = zstd.compress(head,c-head,compressed_tree_size);
        //            delete[] zstd.compress(head,c-head,compressed_tree_size);
        //
        //            printf("compressed huffman tree size = %d\n",(int)compressed_tree_size);

        // return;
    }

    void saveAsDFSOrder(uchar *&c) {
        // uchar *head = c;

        uchar mask = 0, index = 0;

        writeBytesByte(c, (tree.usemp << 7) | ((tree.n == 1) << 6) | tree.mbft);
        writeBytes(c, tree.offset, sizeof(T) << 3, mask, index);
        int64ToBytes_bigEndian(c, tree.n);
        c += sizeof(size_t);
        if (tree.usemp == 0x00) {
            int64ToBytes_bigEndian(c, tree.maxval);
            c += sizeof(size_t);
        }

        if (tree.n == 0 || tree.n == 1) {
            writeBytesClearMask(c, mask, index);
            return;
        }

        std::stack<Node *> stk;

        stk.push(&tree.ht[tree.root]);

        while (!stk.empty()) {
            Node *u = stk.top();
            stk.pop();
            if (u->isLeaf()) {
                writeBytesBit(c, 0x01, mask, index);
                writeBytes(c, u->c, tree.mbft, mask, index);
            } else {
                writeBytesBit(c, 0x00, mask, index);
                stk.push(u->p[1]);
                stk.push(u->p[0]);
            }
        }

        writeBytesClearMask(c, mask, index);

        //            printf("n = %d\n",tree.n);
        //
        //            printf("huffman tree size = %d\n",(int)(c-head));
        //
        //            Lossless_zstd zstd;
        //            size_t compressed_tree_size;
        //
        //            // uchar *compressed_tree = zstd.compress(head,c-head,compressed_tree_size);
        //            delete[] zstd.compress(head,c-head,compressed_tree_size);
        //
        //            printf("compressed huffman tree size = %d\n",(int)compressed_tree_size);
    }

    void loadAsCode(const uchar *&bytes, size_t &remaining_length) {
        // Timer timer(true);

        tree.init();

        uchar feature = (*bytes) >> 6;
        tree.mbft = (*bytes) & 0x3f;
        ++bytes;

        // uchar szT = ((*bytes) >> 5) + 1;
        tree.limit = ((*bytes) & 0x1f) + 1;
        ++bytes;

        // assert(szT == sizeof(T));

        for (int i = 0; i < sizeof(T); i++) {
            tree.offset |= static_cast<T>(*bytes) << (i << 3);
            ++bytes;
        }

        tree.maxval = bytesToInt32_bigEndian(bytes);
        bytes += 4;

        tree.usemp = tree.maxval >= (1 << 12) && (1 << (tree.limit - 1)) < tree.maxval ? 1 : 0;

        if (tree.usemp) {
            tree.ht.reserve(2 << tree.limit);
            //                tree.mpfreq.reserve(1<<tree.limit);
            //                tree.mplen.reserve(1<<tree.limit);
            //                tree.mpcode.reserve(1<<tree.limit);
        } else {
            tree.ht.reserve(tree.maxval << 1);
            tree.vecfreq.resize(tree.maxval);
            tree.veclen.resize(tree.maxval);
            tree.veccode.resize(tree.maxval);
        }

        tree.ht.push_back(Node());

        if (feature == 0x00 || feature == 0x02) {
            int i = 0;
            tree.n = 1 << tree.limit;
            if (feature == 0x02) {
                tree.n = 1;
                tree.ht.resize(2);
                tree.root = 0;
                tree.ht[0] = Node(0, &tree.ht[1]);
                tree.ht[1] = Node(0);
                tree.mplen[0] = 1;
                tree.mpcode[0] = 0;
                tree.setConstructed();
                return;
            }

            for (int j = 0; j < tree.n; j++) {
                T c = 0;

                for (uchar k = 0; k < tree.mbft; k++) {
                    c |= static_cast<T>(readBit(bytes, i++)) << k;
                }

                Node *u = &tree.ht[tree.root];
                int vec = 0;

                for (uchar k = 0; k < tree.limit; k++) {
                    int e = readBit(bytes, i++);
                    vec |= e << k;

                    if (u->p[e] == nullptr) {
                        tree.ht.push_back(Node());
                        u->p[e] = &tree.ht[tree.ht.size() - 1];
                    }

                    u = u->p[e];
                }

                u->c = c;
                if (tree.usemp) {
                    tree.mplen[c] = tree.limit;
                    tree.mpcode[c] = vec;
                } else {
                    tree.veclen[c] = tree.limit;
                    tree.veccode[c] = vec;
                }
            }

            bytes += (i + 7) >> 3;

            tree.setConstructed();

            return;
        }

        tree.n = 0;

        int i = 0;

        for (uchar len = 1; len <= tree.limit; len++) {
            int cnt = 0;

            for (uchar j = 0; j < len; j++) {
                cnt |= static_cast<int>(readBit(bytes, i++)) << j;
            }

            for (int j = 0; j < cnt; j++) {
                T c = 0;

                for (uchar k = 0; k < tree.mbft; k++) {
                    c |= static_cast<T>(readBit(bytes, i++)) << k;
                }

                Node *u = &tree.ht[0];
                int vec = 0;

                for (int k = 0; k < len; k++) {
                    int e = readBit(bytes, i++);
                    vec |= e << k;
                    if (u->p[e] == nullptr) {
                        tree.ht.push_back(Node());
                        u->p[e] = &tree.ht[tree.ht.size() - 1];
                    }

                    u = u->p[e];
                }

                u->c = c;
                ++tree.n;
                if (tree.usemp) {
                    tree.mplen[c] = len;
                    tree.mpcode[c] = vec;
                } else {
                    tree.veclen[c] = len;
                    tree.veccode[c] = vec;
                }
            }
        }

        bytes += (i + 7) >> 3;

        tree.setConstructed();

        // timer.stop("loadAsCode");
    }

    void loadAsDFSOrder(const uchar *&bytes, size_t &remaining_length) {
        tree.init();

        tree.usemp = (*bytes) >> 7;
        tree.mbft = (*bytes) & 0x3f;
        ++bytes;

        for (int i = 0; i < sizeof(T); i++) {
            tree.offset |= static_cast<T>(*bytes) << (i << 3);
            ++bytes;
        }

        tree.n = bytesToInt64_bigEndian(bytes);
        bytes += sizeof(size_t);
        tree.ht.reserve(tree.n << 1);

        if (tree.usemp == 0x00) {
            tree.maxval = bytesToInt64_bigEndian(bytes);
            bytes += sizeof(size_t);
            if (tree.n > 0) {
                tree.veccode.resize(tree.maxval);
                tree.veclen.resize(tree.maxval);
            }
        }

        if (tree.n == 0) {
            tree.setConstructed();
            return;
        }

        if (tree.n == 1) {
            tree.ht.resize(2);
            tree.root = 0;
            tree.ht[0] = Node(0, &tree.ht[1]);
            tree.ht[1] = Node(0);
            if (tree.usemp) {
                tree.mplen[0] = 1;
                tree.mpcode[0] = 0;
            } else {
                tree.veclen[0] = 1;
                tree.veccode[0] = 0;
            }
            tree.setConstructed();
            return;
        }

        tree.root = 0;
        tree.ht.push_back(Node());
        std::stack<Node *> stk;
        stk.push(&tree.ht[0]);
        size_t i = 1;

        while (!stk.empty()) {
            Node *u = stk.top();

            if (readBit(bytes, i++) == 0x00) {
                tree.ht.push_back(Node());
                if (u->p[0] == nullptr) {
                    u->p[0] = &tree.ht[tree.ht.size() - 1];
                } else {
                    u->p[1] = &tree.ht[tree.ht.size() - 1];
                }
                stk.push(&tree.ht[tree.ht.size() - 1]);
            } else {
                T c = 0;
                for (int j = 0; j < tree.mbft; j++) c |= static_cast<T>(readBit(bytes, i++)) << j;
                tree.ht.push_back(Node(c));
                if (u->p[0] == nullptr)
                    u->p[0] = &tree.ht[tree.ht.size() - 1];
                else
                    u->p[1] = &tree.ht[tree.ht.size() - 1];
                while (!stk.empty() && stk.top()->p[1] != nullptr) {
                    //                        Node *tem=stk.top();
                    stk.pop();
                    //                        if(!stk.empty()){
                    //                            if(stk.top()->p[0]==nullptr) stk.top()->p[0]=tem;
                    //                            else stk.top()->p[1]=tem;
                    //                        }
                }
            }
        }

        bytes += (i + 7) >> 3;

        if (tree.usemp) {
            tree.dfs_mp(&tree.ht[tree.root]);
        } else {
            tree.dfs_vec(&tree.ht[tree.root]);
        }

        tree.setConstructed();
    }
};
}  // namespace SZ3

#endif
