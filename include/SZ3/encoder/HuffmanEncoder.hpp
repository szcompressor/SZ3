#ifndef _SZ_HUFFMAN_ENCODER_LZ_HPP
#define _SZ_HUFFMAN_ENCODER_LZ_HPP

#include "SZ3/def.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/utils/ska_hash/unordered_map.hpp"
#include <iostream>
#include <assert.h>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

namespace SZ{

    template<class T>
    class HuffmanEncoder:public concepts::EncoderInterface<T>{

    private:

        class Node{

        public:

            Node(T c_=0,Node *lp=nullptr,Node *rp=nullptr){

                c=c_;
                p[0]=lp;
                p[1]=rp;
            }

            T c;
            Node *p[2];

            inline uchar isLeaf(){

                return p[0]==nullptr;
            }
        };

        class HuffmanTree{

        private:

            uchar _constructed=0;

            uchar len=0;
            int vec=0;

            void dfs(Node* u){

                if(u->isLeaf()){

                    mplen[u->c]=len;
                    mpcode[u->c]=vec;

                    limit=std::max(limit,len);

                    return;
                }

                ++len;
                dfs(u->p[0]);
                --len;

                vec^=1<<len++;
                dfs(u->p[1]);
                vec^=1<<--len;
            }

            class cmp{
            public:
                bool operator()(const std::pair<int,size_t>& u, const std::pair<int,size_t>& v) {
                    return u.second==v.second?u.first>v.first:u.second>v.second;
                }
            };

        public:

            std::vector<uchar> mplen;
            std::vector<int> mpcode;

            T offset;
            // minimum bits for T
            uchar mbft;
            uchar limit;

            void init(){

                _constructed=0;
                ht.clear();
                mplen.clear();
                mpcode.clear();
                freq.clear();

                offset=0;
                mbft=0;
                root=0;
                n=0;
                maxval=0;
                limit=0;
            }

            HuffmanTree(){

                init();
            }

            int root;
            int n;
            int maxval;
            std::vector<Node> ht;
            std::vector<size_t> freq;

            void addElement(T c,size_t freqc){

                assert(!_constructed);

                ht.push_back(Node(c));
                freq[c]=freqc;
                ++n;
            }

            void constructHuffmanTree(){

                assert(!_constructed);
                assert(ht.size()>1);

                if(maxval==1){

                    mbft=1;
                    ht.push_back(Node(0,&ht[0],nullptr));
                    mplen[0]=1;
                    mpcode[0]=0;
                    limit=1;
                    setConstructed();
                    return;
                }

                Timer timer(true);

                mbft=1;
                while((1<<mbft)<maxval) ++mbft;

                std::priority_queue<std::pair<int,size_t>,std::vector<std::pair<int,size_t>>,cmp> q;

                for(int i=0;i<ht.size();i++){

                    q.push({i,freq[ht[i].c]});
                }

                while(q.size()>1){

                    int u=q.top().first;
                    size_t freq_u=q.top().second;
                    q.pop();
                    int v=q.top().first;
                    size_t freq_v=q.top().second;
                    q.pop();

                    ht.push_back(Node(0,&ht[u],&ht[v]));

                    q.push({ht.size()-1,freq_u+freq_v});
                }

                root=ht.size()-1;

                dfs(&ht[root]);

                setConstructed();

                timer.stop("construct huffman tree");
            }

            uchar isConstructed(){

                return _constructed;
            }

            void setConstructed(){

                _constructed=1;
            }
        };

        HuffmanTree tree;

    public:

        void preprocess_encode(const T *const bins,size_t num_bin,int stateNum){

            Timer timer(true);

            tree.init();

            T __minval,__maxval;

            if(stateNum==0){

                __minval=*bins;
                __maxval=*bins;
                for(int i=1;i<num_bin;i++){
                    __minval=std::min(__minval,*(bins+i));
                    __maxval=std::max(__maxval,*(bins+i));
                }
            }
            else{

                __minval=0;
                __maxval=stateNum-1;
            }

            tree.offset=__minval;
            tree.maxval=__maxval-__minval+1;
            tree.freq.resize(tree.maxval);
            tree.mplen.resize(tree.maxval);
            tree.mpcode.resize(tree.maxval);

            // printf("begins to cal freq\n");

            // ska::unordered_map<T,size_t> freq;
            std::vector<size_t> freq(tree.maxval);
            // freq.reserve(4*stateNum);

            for(int i=0;i<num_bin;i++){

                const T& it=*(bins+i);
                ++freq[it-tree.offset];
            }

            tree.ht.reserve(freq.size()<<1);

            for(int i=0;i<tree.maxval;i++){

                if(freq[i]) tree.addElement(i,freq[i]);
            }

            // printf("begins to construct huffman tree\n");

            tree.constructHuffmanTree();

            timer.stop("preprocess_encode");
        }

        void preprocess_encode(const std::vector<T> &bins,int stateNum){

            preprocess_encode(bins.data(),bins.size(),stateNum);
        }

        void saveAsCode(uchar *&c){

            Timer timer(true);

            uchar *head=c;

            // whether the tree is full binary tree

            uchar& limit=tree.limit;

            std::vector<std::deque<T>> mp(limit+1);

            for(int i=0;i<tree.maxval;i++){

                mp[tree.mplen[i]].push_back(i);
            }

            uchar mask=0;
            uchar index=0;

            assert(sizeof(T)<=8);

            if(mp[limit].size()==tree.n){

                // 00 XXXXXX (mbft)
                if(tree.maxval>1) writeBytesByte(c,tree.mbft);
                else writeBytesByte(c,0x80|tree.mbft);

                writeBytesByte(c,((sizeof(T)-1)<<5)|(limit-1));

                writeBytes(c,tree.offset,sizeof(T)<<3,mask,index);

                int32ToBytes_bigEndian(c,tree.n);
                c+=4;

                int cnt=mp[limit].size();

                uchar logcnt=0;
                while(logcnt<32&&(1<<logcnt)!=cnt) ++logcnt;
                assert(logcnt!=32);

                for(T it:mp[limit]){

                    writeBytes(c,it,tree.mbft,mask,index);

                    const int code=tree.mpcode[it];

                    writeBytes(c,code,logcnt,mask,index);
                }

                writeBytesClearMask(c,mask,index);

                return;
            }

            writeBytesByte(c,0x40|tree.mbft);

            writeBytesByte(c,((sizeof(T)-1)<<5)|(limit-1));

            writeBytes(c,tree.offset,sizeof(T)<<3,mask,index);

            int32ToBytes_bigEndian(c,tree.maxval);
            c+=4;

            for(uchar len=1;len<=limit;len++){

                int cnt=mp[len].size();

                writeBytes(c,cnt,len,mask,index);

                if(cnt){

                    for(const T& it:mp[len]){

                        writeBytes(c,it,tree.mbft,mask,index);

                        const int code=tree.mpcode[it];

                        writeBytes(c,code,len,mask,index);
                    }
                }
            }

            writeBytesClearMask(c,mask,index);

            timer.stop("saveAsCode");

            // printf("huffman tree size = %d\n",(int)(c-head));

            // Lossless_zstd zstd;
            // size_t compressed_tree_size;

            // // uchar *compressed_tree = zstd.compress(head,c-head,compressed_tree_size);
            // delete[] zstd.compress(head,c-head,compressed_tree_size);

            // printf("compressed huffman tree size = %d\n",(int)compressed_tree_size);

            return;
        }

        void loadAsCode(const uchar *&bytes,size_t &remaining_length){

            Timer timer(true);

            tree.init();

            uchar feature=(*bytes)>>6;
            tree.mbft=(*bytes)&0x3f;
            ++bytes;

            uchar szT=((*bytes)>>5)+1;
            tree.limit=((*bytes)&0x1f)+1;
            ++bytes;

            assert(szT==sizeof(T));

            for(int i=0;i<sizeof(T);i++){

                tree.offset|=(T)(*bytes)<<(i<<3);
                ++bytes;
            }

            tree.maxval=bytesToInt32_bigEndian(bytes);
            bytes+=4;

            tree.ht.reserve(tree.maxval<<1);
            tree.freq.resize(tree.maxval);
            tree.mplen.resize(tree.maxval);
            tree.mpcode.resize(tree.maxval);

            tree.ht.push_back(Node());

            if(feature==0x00||feature==0x02){

                int i=0;
                tree.n=1<<tree.limit;
                if(feature==0x02) tree.n=1;

                for(int j=0;j<tree.n;j++){

                    T c=0;

                    for(uchar k=0;k<tree.mbft;k++){

                        c|=(T)readBit(bytes,i++)<<k;
                    }

                    Node* u=&tree.ht[tree.root];
                    int vec=0;

                    for(uchar k=0;k<tree.limit;k++){

                        int e=readBit(bytes,i++);
                        vec|=e<<k;

                        if(u->p[e]==nullptr){

                            tree.ht.push_back(Node());
                            u->p[e]=&tree.ht[tree.ht.size()-1];
                        }

                        u=u->p[e];
                    }

                    u->c=c;
                    assert(c<tree.mplen.size());
                    tree.mplen[c]=tree.limit;
                    tree.mpcode[c]=vec;
                }

                bytes+=(i+7)>>3;

                return;
            }

            tree.n=0;

            int i=0;

            for(uchar len=1;len<=tree.limit;len++){

                int cnt=0;

                for(uchar j=0;j<len;j++){

                    cnt|=(int)(readBit(bytes,i++))<<j;
                }

                for(int j=0;j<cnt;j++){

                    T c=0;

                    for(uchar k=0;k<tree.mbft;k++){

                        c|=(T)readBit(bytes,i++)<<k;
                    }

                    Node* u=&tree.ht[0];
                    int vec=0;

                    for(int k=0;k<len;k++){

                        int e=readBit(bytes,i++);
                        vec|=e<<k;
                        if(u->p[e]==nullptr){

                            tree.ht.push_back(Node());
                            u->p[e]=&tree.ht[tree.ht.size()-1];
                        }

                        u=u->p[e];
                    }

                    u->c=c;
                    ++tree.n;
                    tree.mplen[c]=len;
                    tree.mpcode[c]=vec;
                }
            }

            bytes+=(i+7)>>3;

            timer.stop("loadAsCode");

            tree.setConstructed();
        }

        size_t encode(const T *bins, size_t num_bin, uchar *&bytes){

            if(tree.maxval==1){

                int32ToBytes_bigEndian(bytes,num_bin^0x1234abcd);
                bytes+=4;
                return 4;
            }

            Timer timer(true);

            assert(tree.isConstructed());

            uchar *head=bytes;
            bytes+=4;

            int len=0;

            uchar mask=0;
            uchar index=0;

            for(int i=0;i<num_bin;i++){

                const T& it=*(bins+i);

                const uchar len_i=tree.mplen[it-tree.offset];
                const int code_i=tree.mpcode[it-tree.offset];

                len+=len_i;

                writeBytes(bytes,code_i,len_i,mask,index);
            }

            writeBytesClearMask(bytes,mask,index);

            int32ToBytes_bigEndian(head,len^0x1234abcd);

            timer.stop("encode");

            // printf("code size = %d\n",(int)(bytes-head));

            // Lossless_zstd zstd;
            // size_t compressed_code_size;

            // delete[] zstd.compress(head,bytes-head,compressed_code_size);

            // printf("compressed code size = %d\n",(int)compressed_code_size);

            return bytes-head;
        }

        size_t encode(const std::vector<T> &bins, uchar *&bytes){

            return encode(bins.data(),bins.size(),bytes);
        }

        void postprocess_encode(){

        }

        void preprocess_decode(){

        }

        std::vector<T> decode(const uchar *&bytes, size_t targetLength){

            if(tree.maxval==1){

                int len=bytesToInt32_bigEndian(bytes)^0x1234abcd;
                bytes+=4;
                assert(len==targetLength);
                return std::vector<T>(len,tree.offset);
            }

            Timer timer(true);

            assert(tree.isConstructed());

            assert(targetLength>4);

            Node *u=&tree.ht[tree.root];

            int len=bytesToInt32_bigEndian(bytes)^0x1234abcd;
            bytes+=4;

            std::vector<T> a(targetLength);
            int sza=0;
            // a.reserve(targetLength);

            // for(int i=0;i<len;){

            //     u=u->p[readBit(bytes,i++)];

            //     if(u->isLeaf()){

            //         a[sza++]=u->c+tree.offset;
            //         u=&tree.ht[tree.root];
            //     }
            // }

            // use unroll loops to optimize the above code

            int byteIndex=0;
            int i=0;
            uchar b;
            for(;i+8<len;i+=8,byteIndex++){

                b=bytes[byteIndex];

                u=u->p[b&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>1)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>2)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>3)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>4)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>5)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>6)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
                u=u->p[(b>>7)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
            }

            b=bytes[byteIndex];

            for(int j=0;j<len-i;j++){

                u=u->p[(b>>j)&1];
                if(u->isLeaf()){
                    a[sza++]=u->c+tree.offset;
                    u=&tree.ht[tree.root];
                }
            }

            bytes+=(len+7)>>3;

            timer.stop("decode");

            return a;
        }

        void postprocess_decode(){

        }

        void save(uchar *&c){

            saveAsCode(c);
        }

        void load(const uchar *&c,size_t &remaining_length){

            loadAsCode(c,remaining_length);
        }

    };

}

#endif