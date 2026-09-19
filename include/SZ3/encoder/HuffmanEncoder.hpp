#ifndef SZ3_HUFFMAN_ENCODER_HPP
#define SZ3_HUFFMAN_ENCODER_HPP

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_set>

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"

namespace SZ3 {

template <class T>
class HuffmanEncoder : public concepts::EncoderInterface<T> {
   public:
    typedef struct node_t {
        struct node_t *left, *right;
        size_t freq;
        char t;  // in_node:0; otherwise:1
        T c;
    } *node;

    typedef struct HuffmanTree {
        unsigned int stateNum;
        unsigned int allNodes;
        struct node_t *pool;
        node *qqq, *qq;  // the root node of the HuffmanTree is qq[1]
        int n_nodes;     // n_nodes is for compression
        int qend;
        uint64_t **code;
        unsigned char *cout;
        int n_inode;  // n_inode is for decompression
        int maxBitCount;
    } HuffmanTree;

    HuffmanEncoder() {
        int x = 1;
        char *y = reinterpret_cast<char *>(&x);
        if (*y == 1)
            sysEndianType = 0;
        else  //=0
            sysEndianType = 1;
    }

    ~HuffmanEncoder() override { SZ_FreeHuffman(); }

    // Declared because the destructor above suppresses the implicit ones. Copying is shallow and
    // has to stay -- the SZAlgo* wirings pass encoders by value -- so do not copy one that has
    // built its tree: both copies would free huffmanTree. Today every copy happens before that.
    HuffmanEncoder(const HuffmanEncoder &) = default;
    HuffmanEncoder &operator=(const HuffmanEncoder &) = default;
    HuffmanEncoder(HuffmanEncoder &&) noexcept = default;
    HuffmanEncoder &operator=(HuffmanEncoder &&) noexcept = default;

    // build huffman tree
    HuffmanTree *createHuffmanTree(int stateNum) {
        HuffmanTree *tree = static_cast<HuffmanTree *>(malloc(sizeof(HuffmanTree)));
        memset(tree, 0, sizeof(HuffmanTree));
        tree->stateNum = stateNum;
        tree->allNodes = 2 * stateNum;

        tree->pool = static_cast<struct node_t *>(malloc(tree->allNodes * 2 * sizeof(struct node_t)));
        tree->qqq = static_cast<node *>(malloc(tree->allNodes * 2 * sizeof(node)));
        tree->code = static_cast<uint64_t **>(malloc(tree->stateNum * sizeof(uint64_t *)));
        tree->cout = static_cast<unsigned char *>(malloc(tree->stateNum * sizeof(unsigned char)));

        memset(tree->pool, 0, tree->allNodes * 2 * sizeof(struct node_t));
        memset(tree->qqq, 0, tree->allNodes * 2 * sizeof(node));
        memset(tree->code, 0, tree->stateNum * sizeof(uint64_t *));
        memset(tree->cout, 0, tree->stateNum * sizeof(unsigned char));
        tree->qq = tree->qqq - 1;
        tree->n_nodes = 0;
        tree->n_inode = 0;
        tree->qend = 1;

        return tree;
    }

    /**
     * build huffman tree using bins
     * @param stateNum ignored. The interface offers it as a promise that the bins fall in [0, stateNum); this encoder
     * takes its range from the bins themselves, so the promise buys it nothing.
     */
    void preprocess_encode(const std::vector<T> &bins, int stateNum) override {
        preprocess_encode(bins.data(), bins.size(), stateNum);
    }

    /**
     * build huffman tree using bins
     * @param num_bin how many bins `bins` points at; a raw pointer carries no length of its own
     */
    void preprocess_encode(const T *bins, size_t num_bin, int /*stateNum*/) {
        nodeCount = 0;
        if (num_bin == 0) {
            throw std::invalid_argument("Huffman bins should not be empty");
        }
        init(bins, num_bin);
        for (unsigned int i = 0; i < huffmanTree->stateNum; i++)
            if (huffmanTree->code[i]) nodeCount++;
        nodeCount = nodeCount * 2 - 1;
    }

    // save the huffman Tree in the compressed data
    void save(uchar *&c) override {
        // auto cc = c;
        write(offset, c);
        int32ToBytes_bigEndian(c, nodeCount);
        c += sizeof(int);
        int32ToBytes_bigEndian(c, huffmanTree->stateNum / 2);
        c += sizeof(int);
        uint totalSize = 0;  // = convert_HuffTree_to_bytes_anyStates(nodeCount, c);
        // std::cout << "nodeCount = " << nodeCount << std::endl;
        if (nodeCount <= 256)
            totalSize = convert_HuffTree_to_bytes_anyStates<unsigned char>(nodeCount, c);
        else if (nodeCount <= 65536)
            totalSize = convert_HuffTree_to_bytes_anyStates<unsigned short>(nodeCount, c);
        else
            totalSize = convert_HuffTree_to_bytes_anyStates<unsigned int>(nodeCount, c);
        c += totalSize;
        //            return c - cc;
    }

    size_t size_est() override {
        size_t b = (nodeCount <= 256) ? sizeof(unsigned char)
                                      : ((nodeCount <= 65536) ? sizeof(unsigned short) : sizeof(unsigned int));
        return 1 + 2 * nodeCount * b + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T) + sizeof(int) +
               sizeof(int) + sizeof(T);
    }

    // perform encoding
    size_t encode(const std::vector<T> &bins, uchar *&bytes) override {
        return encode(bins.data(), bins.size(), bytes);
    }

    // perform encoding
    size_t encode(const T *bins, size_t num_bin, uchar *&bytes) {
        size_t outSize = 0;
        size_t i = 0;
        unsigned char bitSize = 0, byteSize, byteSizep;
        int state;
        uchar *p = bytes + sizeof(size_t);
        int lackBits = 0;
        // int64_t totalBitSize = 0, maxBitSize = 0, bitSize21 = 0, bitSize32 = 0;
        for (i = 0; i < num_bin; i++) {
            state = bins[i] - offset;
            bitSize = huffmanTree->cout[state];

            if (lackBits == 0) {
                byteSize = bitSize % 8 == 0
                               ? bitSize / 8
                               : bitSize / 8 + 1;  // it's equal to the number of bytes involved (for *outSize)
                byteSizep = bitSize / 8;           // it's used to move the pointer p for next data
                if (byteSize <= 8) {
                    int64ToBytes_bigEndian(p, (huffmanTree->code[state])[0]);
                    p += byteSizep;
                } else  // byteSize>8
                {
                    int64ToBytes_bigEndian(p, (huffmanTree->code[state])[0]);
                    p += 8;
                    int64ToBytes_bigEndian(p, (huffmanTree->code[state])[1]);
                    p += (byteSizep - 8);
                }
                outSize += byteSize;
                lackBits = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
            } else {
                *p = (*p) | static_cast<unsigned char>((huffmanTree->code[state])[0] >> (64 - lackBits));
                if (lackBits < bitSize) {
                    p++;

                    int64_t newCode = (huffmanTree->code[state])[0] << lackBits;
                    int64ToBytes_bigEndian(p, newCode);

                    if (bitSize <= 64) {
                        bitSize -= lackBits;
                        byteSize = bitSize % 8 == 0 ? bitSize / 8 : bitSize / 8 + 1;
                        byteSizep = bitSize / 8;
                        p += byteSizep;
                        outSize += byteSize;
                        lackBits = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
                    } else  // bitSize > 64
                    {
                        byteSizep = 7;  // must be 7 bytes, because lackBits!=0
                        p += byteSizep;
                        outSize += byteSize;

                        bitSize -= 64;
                        if (lackBits < bitSize) {
                            *p = (*p) | static_cast<unsigned char>((huffmanTree->code[state])[0] >> (64 - lackBits));
                            p++;
                            newCode = (huffmanTree->code[state])[1] << lackBits;
                            int64ToBytes_bigEndian(p, newCode);
                            bitSize -= lackBits;
                            byteSize = bitSize % 8 == 0 ? bitSize / 8 : bitSize / 8 + 1;
                            byteSizep = bitSize / 8;
                            p += byteSizep;
                            outSize += byteSize;
                            lackBits = bitSize % 8 == 0 ? 0 : 8 - bitSize % 8;
                        } else  // lackBits >= bitSize
                        {
                            *p = (*p) | static_cast<unsigned char>((huffmanTree->code[state])[0] >> (64 - bitSize));
                            lackBits -= bitSize;
                        }
                    }
                } else  // lackBits >= bitSize
                {
                    lackBits -= bitSize;
                    if (lackBits == 0) p++;
                }
            }
        }
        write(outSize, bytes);
        bytes += outSize;  // move pointer to end of encoded array
        return outSize;
    }

    void postprocess_encode() override { SZ_FreeHuffman(); }

    void preprocess_decode() override {}

    // perform decoding
    std::vector<T> decode(const uchar *&bytes, size_t targetLength, size_t &remaining_length) override {
        node t = treeRoot;
        std::vector<T> out(targetLength);
        size_t i = 0, byteIndex = 0, count = 0;
        int r;
        node n = treeRoot;
        size_t encodedLength = 0;
        read(encodedLength, bytes, remaining_length);
        if (n->t)  // root->t==1 means that all state values are the same (constant)
        {
            for (count = 0; count < targetLength; count++) out[count] = n->c + offset;
            return out;
        }

        if (encodedLength > remaining_length)
            throw std::out_of_range("SZ3 Huffman: encoded length exceeds compressed buffer");

        // Walk at most the bits the stream holds, and stop once targetLength symbols are out.
        const size_t maxBits = targetLength > 0 ? encodedLength * 8 : 0;
        for (i = 0; i < maxBits; i++) {
            byteIndex = i >> 3;  // i/8
            r = i % 8;
            if (((bytes[byteIndex] >> (7 - r)) & 0x01) == 0)
                n = n->left;
            else
                n = n->right;

            if (n->t) {
                out[count] = n->c + offset;
                n = t;
                if (++count == targetLength) break;
            }
        }
        if (count < targetLength) throw std::out_of_range("SZ3 Huffman: corrupted encoded stream");
        bytes += encodedLength;
        remaining_length -= encodedLength;
        return out;
    }

    // empty function
    void postprocess_decode() override { SZ_FreeHuffman(); }

    // load Huffman tree
    void load(const uchar *&c, size_t &remaining_length) override {
        read(offset, c, remaining_length);
        if (remaining_length < 2 * sizeof(int)) throw std::out_of_range("SZ3 Huffman: truncated tree header");
        nodeCount = bytesToInt32_bigEndian(c);
        // The stored state count is skipped: it sizes the encode-side code tables, which decoding never
        // touches. nodeCount is bounded before it sizes anything, or the encodeStartIndex arithmetic
        // below overflows and the tree is read past the buffer.
        if (nodeCount <= 0 || static_cast<size_t>(nodeCount) > remaining_length)
            throw std::out_of_range("SZ3 Huffman: invalid node count");
        size_t encodeStartIndex;
        if (nodeCount <= 256)
            encodeStartIndex = 1 + 3 * nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        else if (nodeCount <= 65536)
            encodeStartIndex =
                1 + 2 * nodeCount * sizeof(unsigned short) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);
        else
            encodeStartIndex =
                1 + 2 * nodeCount * sizeof(unsigned int) + nodeCount * sizeof(unsigned char) + nodeCount * sizeof(T);

        size_t tree_bytes = sizeof(int) + sizeof(int) + encodeStartIndex;
        if (tree_bytes > remaining_length) throw std::out_of_range("SZ3 Huffman: tree exceeds compressed buffer");

        // The pool is 4x what is asked for, and unpad_tree builds each of the nodeCount nodes once.
        huffmanTree = createHuffmanTree(nodeCount);
        treeRoot = reconstruct_HuffTree_from_bytes_anyStates(c + sizeof(int) + sizeof(int), nodeCount);
        c += tree_bytes;
        remaining_length -= tree_bytes;
        loaded = true;
    }

    bool isLoaded() const { return loaded; }

   private:
    HuffmanTree *huffmanTree = nullptr;
    node treeRoot;
    unsigned int nodeCount = 0;
    uchar sysEndianType;  // 0: little endian, 1: big endian
    bool loaded = false;
    T offset;

    node reconstruct_HuffTree_from_bytes_anyStates(const unsigned char *bytes, uint nodeCount_) {
        if (nodeCount_ <= 256) {
            std::vector<unsigned char> L(nodeCount_);
            std::vector<unsigned char> R(nodeCount_);
            std::vector<T> C(nodeCount_);
            std::vector<unsigned char> t(nodeCount_);
            // TODO: Endian type
            // unsigned char cmpSysEndianType = bytes[0];
            // if(cmpSysEndianType!=(unsigned char)sysEndianType)
            // {
            // 	unsigned char* p = (unsigned char*)(bytes+1+2*nodeCount_*sizeof(unsigned char));
            // 	size_t i = 0, size = nodeCount_*sizeof(unsigned int);
            // 	while(1)
            // 	{
            // 		symTransform_4bytes(p);
            // 		i+=sizeof(unsigned int);
            // 		if(i<size)
            // 			p+=sizeof(unsigned int);
            // 		else
            // 			break;
            // 	}
            // }
            memcpy(L.data(), bytes + 1, nodeCount_ * sizeof(unsigned char));
            memcpy(R.data(), bytes + 1 + nodeCount_ * sizeof(unsigned char), nodeCount_ * sizeof(unsigned char));
            memcpy(C.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned char), nodeCount_ * sizeof(T));
            memcpy(t.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned char) + nodeCount_ * sizeof(T),
                   nodeCount_ * sizeof(unsigned char));
            std::vector<bool> seen(nodeCount_, false);
            seen[0] = true;
            node root = this->new_node2(C[0], t[0]);
            this->unpad_tree<uchar>(L.data(), R.data(), C.data(), t.data(), 0, root, nodeCount_, seen);
            return root;
        } else if (nodeCount_ <= 65536) {
            std::vector<unsigned short> L(nodeCount_);
            std::vector<unsigned short> R(nodeCount_);
            std::vector<T> C(nodeCount_);
            std::vector<unsigned char> t(nodeCount_);

            // TODO: Endian type
            // unsigned char cmpSysEndianType = bytes[0];
            // if(cmpSysEndianType!=(unsigned char)sysEndianType)
            // {
            // 	unsigned char* p = (unsigned char*)(bytes+1);
            // 	size_t i = 0, size = 3*nodeCount_*sizeof(unsigned int);
            // 	while(1)
            // 	{
            // 		symTransform_4bytes(p);
            // 		i+=sizeof(unsigned int);
            // 		if(i<size)
            // 			p+=sizeof(unsigned int);
            // 		else
            // 			break;
            // 	}
            // }

            memcpy(L.data(), bytes + 1, nodeCount_ * sizeof(unsigned short));
            memcpy(R.data(), bytes + 1 + nodeCount_ * sizeof(unsigned short), nodeCount_ * sizeof(unsigned short));
            memcpy(C.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned short), nodeCount_ * sizeof(T));

            memcpy(t.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned short) + nodeCount_ * sizeof(T),
                   nodeCount_ * sizeof(unsigned char));

            std::vector<bool> seen(nodeCount_, false);
            seen[0] = true;
            node root = this->new_node2(0, 0);
            this->unpad_tree<unsigned short>(L.data(), R.data(), C.data(), t.data(), 0, root, nodeCount_, seen);
            return root;
        } else  // nodeCount_>65536
        {
            std::vector<unsigned int> L(nodeCount_);
            std::vector<unsigned int> R(nodeCount_);
            std::vector<T> C(nodeCount_);
            std::vector<unsigned char> t(nodeCount_);
            // TODO: Endian type
            // unsigned char cmpSysEndianType = bytes[0];
            // if(cmpSysEndianType!=(unsigned char)sysEndianType)
            // {
            // 	unsigned char* p = (unsigned char*)(bytes+1);
            // 	size_t i = 0, size = 3*nodeCount_*sizeof(unsigned int);
            // 	while(1)
            // 	{
            // 		symTransform_4bytes(p);
            // 		i+=sizeof(unsigned int);
            // 		if(i<size)
            // 			p+=sizeof(unsigned int);
            // 		else
            // 			break;
            // 	}
            // }

            memcpy(L.data(), bytes + 1, nodeCount_ * sizeof(unsigned int));
            memcpy(R.data(), bytes + 1 + nodeCount_ * sizeof(unsigned int), nodeCount_ * sizeof(unsigned int));
            memcpy(C.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned int), nodeCount_ * sizeof(T));

            memcpy(t.data(), bytes + 1 + 2 * nodeCount_ * sizeof(unsigned int) + nodeCount_ * sizeof(T),
                   nodeCount_ * sizeof(unsigned char));

            std::vector<bool> seen(nodeCount_, false);
            seen[0] = true;
            node root = this->new_node2(0, 0);
            this->unpad_tree<unsigned int>(L.data(), R.data(), C.data(), t.data(), 0, root, nodeCount_, seen);
            return root;
        }
    }

    node new_node(size_t freq, T c, node a, node b) {
        node n = huffmanTree->pool + huffmanTree->n_nodes++;
        if (freq) {
            n->c = c;
            n->freq = freq;
            n->t = 1;
            // printf("new_node: c = %d, freq = %zu, t = %d \n", n->c, n->freq, n->t);
        } else {
            n->left = a;
            n->right = b;
            n->freq = a->freq + b->freq;
            n->t = 0;
            // printf("new_node: c = %d, freq = %zu, t = %d, left = %d, right = %d \n", n->c, n->freq, n->t, n->left->c,
            // n->right->c);
            // n->c = 0;
        }
        return n;
    }

    node new_node2(T c, unsigned char t) {
        huffmanTree->pool[huffmanTree->n_nodes].c = c;
        huffmanTree->pool[huffmanTree->n_nodes].t = t;
        return huffmanTree->pool + huffmanTree->n_nodes++;
    }

    /* priority queue */
    void qinsert(node n) {
        int j, i = huffmanTree->qend++;
        while ((j = (i >> 1)))  // j=i/2
        {
            if (huffmanTree->qq[j]->freq <= n->freq) break;
            huffmanTree->qq[i] = huffmanTree->qq[j], i = j;
        }
        huffmanTree->qq[i] = n;
    }

    node qremove() {
        int i = 1, l;
        node n = huffmanTree->qq[i = 1];
        node p;
        if (huffmanTree->qend < 2) return nullptr;
        huffmanTree->qend--;
        huffmanTree->qq[i] = huffmanTree->qq[huffmanTree->qend];

        while ((l = (i << 1)) < huffmanTree->qend) {  // l=(i*2)
            if (l + 1 < huffmanTree->qend && huffmanTree->qq[l + 1]->freq < huffmanTree->qq[l]->freq) l++;
            if (huffmanTree->qq[i]->freq > huffmanTree->qq[l]->freq) {
                p = huffmanTree->qq[i];
                huffmanTree->qq[i] = huffmanTree->qq[l];
                huffmanTree->qq[l] = p;
                i = l;
            } else {
                break;
            }
        }
        return n;
    }

    /* walk the tree and put 0s and 1s */
    /**
     * @out1 should be set to 0.
     * @out2 should be 0 as well.
     * @index: the index of the byte
     * */
    void build_code(node n, int len, uint64_t out1, uint64_t out2) {
        if (n->t) {
            huffmanTree->code[n->c] = static_cast<uint64_t *>(malloc(2 * sizeof(uint64_t)));
            if (len <= 64) {
                // A single-symbol tree gives the root a zero-length code, and shifting by 64 is undefined.
                (huffmanTree->code[n->c])[0] = (len == 0) ? 0 : (out1 << (64 - len));
                (huffmanTree->code[n->c])[1] = out2;
            } else {
                (huffmanTree->code[n->c])[0] = out1;
                // len >= 128 would shift by >= 64, and such a code does not fit in 128 bits anyway.
                (huffmanTree->code[n->c])[1] = (len >= 128) ? out2 : (out2 << (128 - len));
            }
            huffmanTree->cout[n->c] = static_cast<unsigned char>(len);
            // std::cout << "build_code: c = " << n->c << ", len = " << len << ", out1 = " << out1 << ", out2 = " << out2
            //           << ", code0 = " << (huffmanTree->code[n->c])[0] << ", code1 = " << (huffmanTree->code[n->c])[1]
            //           << std::endl;
            return;
        }
        int index = len >> 6;  //=len/64
        if (index == 0) {
            out1 = out1 << 1;
            out1 = out1 | 0;
            build_code(n->left, len + 1, out1, 0);
            out1 = out1 | 1;
            build_code(n->right, len + 1, out1, 0);
        } else {
            if (len % 64 != 0) out2 = out2 << 1;
            out2 = out2 | 0;
            build_code(n->left, len + 1, out1, out2);
            out2 = out2 | 1;
            build_code(n->right, len + 1, out1, out2);
        }
    }

    /**
     * Compute the frequency of the data and build the Huffman tree
     * @param s the bins to measure
     * @param length how many bins `s` points at; a raw pointer carries no length of its own
     */
    void init(const T *s, size_t length) {
        // Locals, not `offset` itself: a store to a member of type T may alias the T array being read,
        // and that is enough to stop this reduction vectorising.
        T max = s[0];
        T min = s[0];

        // The range decides how the counts are stored, so it has to be known before anything is sized by it.
        for (size_t i = 0; i < length; i++) {
            if (s[i] > max) {
                max = s[i];
            }
            if (s[i] < min) {
                min = s[i];
            }
        }
        offset = min;  // offset is min

        // The state table is sized by the bin range rather than the distinct count, so a sparse wide-range
        // stream overflows this narrowing. Checked before anything is allocated from that range.
        if (static_cast<double>(max) - static_cast<double>(min) > 2e9) {
            throw std::invalid_argument("HuffmanEncoder: bin range too wide; use HuffmanEncoderV2");
        }
        int stateNum = max - min + 2;
        huffmanTree = createHuffmanTree(stateNum);

        // Counted straight into a dense array indexed by bin: the tree is then built in index order, so it
        // does not depend on the iteration order of any hash container, on linux or on win.
        //
        // Real bins repeat -- long runs land in the same bin -- and `count[bin]++` on one table turns that
        // into a chain of same-address store-to-load forwards, which costs more than the hash map did.
        // Four tables counted in parallel break the chain and are worth 3x on this loop. They are used
        // only while all four fit in kLaneBudget, so a sparse wide range neither pays for three extra
        // tables nor scatters them past the last level of cache, where they would stop helping anyway.
        constexpr size_t kLaneBudget = 8u << 20;
        const size_t lanes = (static_cast<size_t>(stateNum) * sizeof(size_t) * 4 <= kLaneBudget) ? 4 : 1;
        std::vector<size_t> frequencyList(static_cast<size_t>(stateNum) * lanes, 0);
        size_t *freq = frequencyList.data();
        if (lanes == 1) {
            for (size_t i = 0; i < length; i++) {
                freq[s[i] - min] += 1;
            }
        } else {
            size_t *f1 = freq + stateNum, *f2 = f1 + stateNum, *f3 = f2 + stateNum;
            size_t i = 0;
            for (; i + 4 <= length; i += 4) {
                freq[s[i] - min] += 1;
                f1[s[i + 1] - min] += 1;
                f2[s[i + 2] - min] += 1;
                f3[s[i + 3] - min] += 1;
            }
            for (; i < length; i++) {
                freq[s[i] - min] += 1;
            }
            for (int j = 0; j < stateNum; j++) {
                freq[j] += f1[j] + f2[j] + f3[j];
            }
        }
        for (int i = 0; i < stateNum; i++) {
            if (frequencyList[i] != 0) {
                qinsert(new_node(frequencyList[i], i, nullptr, nullptr));
            }
        }

        while (huffmanTree->qend > 2) {
            auto left = qremove();
            auto right = qremove();
            qinsert(new_node(0, 0, left, right));
        }

        build_code(huffmanTree->qq[1], 0, 0, 0);
        treeRoot = huffmanTree->qq[1];
    }

    template <class T1>
    void pad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root) {
        C[i] = root->c;
        t[i] = root->t;
        node lroot = root->left;
        if (lroot != nullptr) {
            huffmanTree->n_inode++;
            L[i] = huffmanTree->n_inode;
            pad_tree(L, R, C, t, huffmanTree->n_inode, lroot);
        }
        node rroot = root->right;
        if (rroot != nullptr) {
            huffmanTree->n_inode++;
            R[i] = huffmanTree->n_inode;
            pad_tree(L, R, C, t, huffmanTree->n_inode, rroot);
        }
    }

    template <class T1>
    void unpad_tree(T1 *L, T1 *R, T *C, unsigned char *t, unsigned int i, node root, unsigned int nodeCount_,
                    std::vector<bool> &seen) {
        // root->c = C[i];
        if (root->t == 0) {
            T1 l, r;
            l = L[i];
            if (l != 0) {
                // pad_tree gives a child a higher index than its parent, so a valid index satisfies i < l < nodeCount_.
                // Enforcing it keeps L/R/C/t reads inside the pool and rules out a cycle.
                if (l <= i || l >= nodeCount_) throw std::out_of_range("SZ3 Huffman: invalid left child index in tree");
                // Increasing indices rule out a cycle but not two parents naming one child, which would
                // expand the tree exponentially instead of building nodeCount_ nodes.
                if (seen[l]) throw std::out_of_range("SZ3 Huffman: tree node reached twice");
                seen[l] = true;
                node lroot = new_node2(C[l], t[l]);
                root->left = lroot;
                unpad_tree(L, R, C, t, l, lroot, nodeCount_, seen);
            }
            r = R[i];
            if (r != 0) {
                if (r <= i || r >= nodeCount_)
                    throw std::out_of_range("SZ3 Huffman: invalid right child index in tree");
                if (seen[r]) throw std::out_of_range("SZ3 Huffman: tree node reached twice");
                seen[r] = true;
                node rroot = new_node2(C[r], t[r]);
                root->right = rroot;
                unpad_tree(L, R, C, t, r, rroot, nodeCount_, seen);
            }
            if (root->left == nullptr || root->right == nullptr) {
                throw std::out_of_range("SZ3 Huffman: internal tree node is missing a child");
            }
        }
    }

    template <class T1>
    unsigned int convert_HuffTree_to_bytes_anyStates(unsigned int nodeCount_, unsigned char *out) {
        T1 *L = static_cast<T1 *>(malloc(nodeCount_ * sizeof(T1)));
        memset(L, 0, nodeCount_ * sizeof(T1));
        T1 *R = static_cast<T1 *>(malloc(nodeCount_ * sizeof(T1)));
        memset(R, 0, nodeCount_ * sizeof(T1));
        T *C = static_cast<T *>(malloc(nodeCount_ * sizeof(T)));
        memset(C, 0, nodeCount_ * sizeof(T));
        unsigned char *t = static_cast<unsigned char *>(malloc(nodeCount_ * sizeof(unsigned char)));
        memset(t, 0, nodeCount_ * sizeof(unsigned char));

        pad_tree(L, R, C, t, 0, huffmanTree->qq[1]);

        unsigned int totalSize =
            1 + 2 * nodeCount_ * sizeof(T1) + nodeCount_ * sizeof(unsigned char) + nodeCount_ * sizeof(T);
        //*out = (unsigned char*)malloc(totalSize);
        out[0] = sysEndianType;
        memcpy(out + 1, L, nodeCount_ * sizeof(T1));
        memcpy(out + 1 + nodeCount_ * sizeof(T1), R, nodeCount_ * sizeof(T1));
        memcpy(out + 1 + 2 * nodeCount_ * sizeof(T1), C, nodeCount_ * sizeof(T));
        memcpy(out + 1 + 2 * nodeCount_ * sizeof(T1) + nodeCount_ * sizeof(T), t, nodeCount_ * sizeof(unsigned char));

        free(L);
        free(R);
        free(C);
        free(t);
        return totalSize;
    }

    void SZ_FreeHuffman() {
        if (huffmanTree != nullptr) {
            size_t i;
            free(huffmanTree->pool);
            huffmanTree->pool = nullptr;
            free(huffmanTree->qqq);
            huffmanTree->qqq = nullptr;
            for (i = 0; i < huffmanTree->stateNum; i++) {
                if (huffmanTree->code[i] != nullptr) free(huffmanTree->code[i]);
            }
            free(huffmanTree->code);
            huffmanTree->code = nullptr;
            free(huffmanTree->cout);
            huffmanTree->cout = nullptr;
            free(huffmanTree);
            huffmanTree = nullptr;
        }
    }
};
}  // namespace SZ3

#endif
