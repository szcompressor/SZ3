//
//  LocalStore.cpp
//  ParSplice-refact
//
//  Created by Danny Perez on 1/9/17.
//  Copyright © 2017 dp. All rights reserved.
//
#include <iostream>
#include <filesystem>
#include <utils/Verification.hpp>
#include "LocalStore.hpp"

namespace fs = std::filesystem;

/*
 * Definition of AbstractLocalDataStore
 */
AbstractLocalDataStore::AbstractLocalDataStore() {
    currentSize = 0;
    maxSize = 0;
};


unsigned int AbstractLocalDataStore::count(unsigned int dbKey, Label &key) {
    return storedKeys[dbKey].count(key);
};

void AbstractLocalDataStore::touch(unsigned int dbKey, Label &key) {
//    if (count(dbKey, key) > 0) {
//        mru.touch(std::make_pair(dbKey, key));
//    }
};

std::set<Label> AbstractLocalDataStore::availableKeys(unsigned int dbKey) {
    if (storedKeys.count(dbKey) == 0) {
        return std::set<Label>();
    }
    return storedKeys[dbKey];
};


std::set<unsigned int> AbstractLocalDataStore::availableDbKeys() {

    std::set<unsigned int> dbKeys;
    for (auto it = storedKeys.begin(); it != storedKeys.end(); it++) {
        dbKeys.insert(it->first);
    }
    return dbKeys;
};

std::unordered_map<unsigned int, std::set<Label> > AbstractLocalDataStore::getStoredKeys() {
    return storedKeys;
};


#if 0
std::multimap< std::pair<unsigned int, Label>, RawDataVector> AbstractLocalDataStore::purge(){
    std::multimap< std::pair<unsigned int, Label>, RawDataVector> deletedItems;

    while(maxSize>0 and currentSize>maxSize) {

        auto p=mru.oldest();

        if(not p.first) {
            std::cout<<"AbstractLocalDataStore::purge ERROR: "<<currentSize<<" "<<maxSize<< " "<<p.second.first<<" "<<p.second.second<<std::endl;
        }
        RawDataVector data;
        get(p.second.first,p.second.second,data);
        deletedItems.insert(std::make_pair(p.second,data));
        //currentSize-=data.size();
        std::cout<<"AbstractLocalDataStore::purge b "<<currentSize<<" "<<maxSize<< " "<<p.second.first<<" "<<p.second.second<<std::endl;
        int sd=erase(p.second.first,p.second.second);
        std::cout<<"AbstractLocalDataStore::purge a "<<currentSize<<" "<<maxSize<< " "<<p.second.first<<" "<<p.second.second<<" "<<sd<<std::endl;
        //mru.pop_oldest();

    }
    return deletedItems;
};
#endif

void AbstractLocalDataStore::purge() {

//    while (maxSize > 0 and currentSize > maxSize) {
//        auto p = mru.oldest();
//        if (not p.first) {
//            std::cout << "AbstractLocalDataStore::purge ERROR: " << currentSize << " " << maxSize << " " << p.second.first << " "
//                      << p.second.second << std::endl;
//            break;
//        }
//        erase(p.second.first, p.second.second);
//    }

};


void AbstractLocalDataStore::setMaximumSize(unsigned long maxSize_) {
    maxSize = maxSize_;
};

unsigned int PersistentLocalStore::count(unsigned int dbKey, Label &key) {
    //return 0;
    return storedKeys[dbKey].count(key);
};

PersistentLocalStore::PersistentLocalStore() {
};

PersistentLocalStore::~PersistentLocalStore() {
    offsets.flush();
    io.flush();
    offsets.close();
    io.close();
};

int PersistentLocalStore::initialize(std::string homeDir, std::string baseName) {

    std::string dbName = homeDir + "/" + baseName + ".db";
    std::string offsetsName = homeDir + "/" + baseName + ".offsets";
    bool restore = fs::exists(dbName) && fs::exists(offsetsName);
//    bool restore = false;

    fs::path p(homeDir);
//    std::cout<<homeDir<<" "<<baseName<<" "<<p.parent_path().string()<<std::endl;
    fs::create_directories(p.parent_path().string());
    fs::create_directories(homeDir);

    offsets.open(offsetsName, std::fstream::in | std::fstream::out | std::fstream::app);
    io.open(dbName, std::fstream::in | std::fstream::out | std::fstream::app | std::fstream::binary);

    std::cout << "STORE INITIALIZE: " << offsets.is_open() << " " << io.is_open() << std::endl;

    if (restore) {
        std::cout << "RESTORING DATABASE" << std::endl;
        //read in the offsets
        int dbk;
        uint64_t k;
        long int off;
        long int s;
        offsets.seekg(0, offsets.beg);
        while (offsets >> dbk >> k >> off >> s) {
            //std::cout<<"RESTORING DB: "<<dbk<<" "<<k<<" "<<off<<" "<<s<<std::endl;
            locations[std::make_pair(dbk, k)] = std::pair<std::streampos, std::streamsize>(off, s);
            storedKeys[dbk].insert(k);
        }
        offsets.close();
        offsets.open(offsetsName, std::fstream::in | std::fstream::out | std::fstream::app);
    }
    return 0;
};

int PersistentLocalStore::put(unsigned int dbKey, uint64_t &key, RawDataVector &data) {
    auto k = std::make_pair(dbKey, key);
    if (locations.count(k) == 0) {
        //insert the item

        //move to the end of the stream
        io.seekp(0, io.end);
        locations[k] = std::pair<std::streampos, std::streamsize>(io.tellp(), data.size());
        offsets.seekp(0, offsets.end);
        offsets << k.first << " " << k.second << " " << locations[k].first << " " << locations[k].second << " " << std::flush;

        //std::cout<<"PUT-IO "<<data.size()<<std::endl;
        io.write(&data[0], data.size());
        io.flush();

        //std::cout<<"PUT-sk"<<std::endl;
        storedKeys[dbKey].insert(key);

        //std::cout<<"PUT-E: "<<k.first<<" "<<k.second<<" "<<locations[k].first<<" "<<locations[k].second<<" "<<offsets.tellp()<<" "<<dbKey<<" "<<key<<" "<<storedKeys[dbKey].size()<<" "<<storedKeys.size()<<std::endl;

        /*
           {
                //TEST THE ROUND TRIP
                RawDataVector dd;
                get(dbKey,key,dd);
                std::cout<<"TEST: "<<data.size()<<" "<<dd.size()<<" "<<(data==dd)<<std::endl;
           }
         */

    }
    return 0;
};


int PersistentLocalStore::get(unsigned int dbKey, uint64_t &key, RawDataVector &data) {
    data.clear();
    auto k = std::make_pair(dbKey, key);
    if (locations.count(k) != 0) {
        auto loc = locations[k];
        io.seekg(loc.first, io.beg);
        data = RawDataVector(loc.second);
        io.read(&data[0], loc.second);
        io.sync();
        //std::cout<<"GET: "<<k.first<<" "<<k.second<<" "<<loc.first<<" "<<loc.second<<" "<<data.size()<<" "<<io.gcount()<<std::endl;
        //std::cout<<io.good()<<" "<<io.eof()<<" "<<io.fail()<<" "<<io.bad()<<std::endl;
    } else {
        return KEY_NOTFOUND;
    }
    return 0;
};

int PersistentLocalStore::createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet) {
    return 0;
};

int PersistentLocalStore::sync() {
    offsets.flush();
    io.flush();
    //std::cout<<"SYNC"<<std::endl;
    return 0;
};

int PersistentLocalStore::erase(unsigned int dbKey, uint64_t &key) {
    return 0;
};


CompressedPersistentLocalStore::CompressedPersistentLocalStore() {
//    SZ_Init(NULL);
};

CompressedPersistentLocalStore::~CompressedPersistentLocalStore() {

};

int CompressedPersistentLocalStore::initialize(std::string homeDir, std::string baseName) {
    return initialize(homeDir, baseName, 3137, 0.0001, 100);
}

int CompressedPersistentLocalStore::initialize(std::string homeDir, std::string baseName, size_t block_size,
                                               double realPrecision, size_t maxBufferSize) {

    compressedStore.initialize(homeDir, baseName + ".c");
    huffmanStore.initialize(homeDir, baseName + ".tree");
    storedKeys = compressedStore.getStoredKeys();
    this->maxBufferSize = maxBufferSize;
    this->block_size = block_size;
    this->realPrecision = realPrecision;
    return 0;
};

int CompressedPersistentLocalStore::put_compressible(unsigned int dbKey, uint64_t &key, std::vector<double> &data) {
    //std::cout<<"============ PUT " <<dbKey <<" "<< key  << "=============="<<std::endl;
    if (count(dbKey, key) == 0) {
        storedKeys[dbKey].insert(key);

        //put in the buffer
        auto k = std::make_pair(dbKey, key);
        buffer[k] = data;

    }

    //if buffer is full, compress and add to datastore
    if (buffer.size() >= maxBufferSize) {
        compressBuffer();
    }

    return 0;
};

int CompressedPersistentLocalStore::put(unsigned int dbKey, uint64_t &key, RawDataVector &data) {
    //std::cout<<"============ PUT " <<dbKey <<" "<< key  << "=============="<<std::endl;
    if (count(dbKey, key) == 0) {
        //put in the buffer

        //unpack
//        SystemType s;
//        unpack(data,s);

        //extract compressible data
        std::vector<double> d;
//        d=s.extractCompressibleData();

        auto k = std::make_pair(dbKey, key);
        buffer[k] = d;

        //pack and store the rest of the object
//        uncompressedStore.put(dbKey, key, data);

        storedKeys[dbKey].insert(key);
    }

    //if buffer is full, compress and add to datastore
    if (buffer.size() >= maxBufferSize) {
        compressBuffer();
    }

    //std::cout<<"============ PUT " <<dbKey <<" "<< key  << " DONE =============="<<std::endl;
    /*
       SystemType s;
       unpack(data,s);

       //test round trip
       RawDataVector dout;
       get(dbKey,key,dout);

       std::cout<<"ROUND TRIP: "<<dout.size()<<std::endl;
       SystemType sout;
       unpack(dout,sout);

       double maxd=0;
       for(int i=0; i<sout.x.size(); i++) {
            maxd=(fabs(sout.x[i]-s.x[i])  > maxd ? fabs(sout.x[i]-s.x[i])  : maxd);
       }
       std::cout<<"MAX COMPRESSION ERROR "<<maxd<<std::endl;
     */

    return 0;
};

void CompressedPersistentLocalStore::compressBuffer() {
    if (buffer.size() == 0) {
        return;
    }
    //std::cout<<"============  COMPRESS "<< buffer.size()<< "=============="<<std::endl;
//    SZ_Init("sz.config");

    std::vector<double> rdata;

    //construct compression buffer
    for (auto it = buffer.begin(); it != buffer.end(); it++) {
        rdata.insert(rdata.end(), it->second.begin(), it->second.end());
    }

    std::vector<size_t> compressed_size(buffer.size() + 1, 0);

    auto conf = SZ::Config<double, 1>(realPrecision, std::array<size_t, 1>{buffer.size() * block_size});
    auto sz = SZ::make_sz_zone_compressor(conf, SZ::LorenzoPredictor<double, 1, 1>(realPrecision),
                                          SZ::LinearQuantizer<double>(realPrecision), SZ::HuffmanEncoder<int>());
//    std::cout << rdata.size() << "," << conf.num<<std::endl;
    auto szresult = sz.compress(&(rdata[0]), compressed_size, buffer.size());

    //save the tree
    RawDataVector tree(szresult[buffer.size()], szresult[buffer.size()] + compressed_size[buffer.size()]);
    huffmanStore.put(0, treeIndex, tree);
//    std::cout << treeIndex << std::endl;
    //std::cout<<"PUT T: "<<treeIndex<<" "<<tree.size()<<" "<<treeByteSize<<" "<<ddd<<std::endl;


    int ii = 0;
    for (auto it = buffer.begin(); it != buffer.end(); it++, ii++) {
//        auto k = std::make_pair(it->first.first, it->first.second);
        RawDataVector d(szresult[ii], szresult[ii] + compressed_size[ii]);
        auto key = it->first.second;
        compressedStore.put(it->first.first, key, d);
        huffmanIndexStore[it->first] = treeIndex;
        //std::cout<<"PUT C: "<<cmprSize<<" "<<treeIndex<<" "<<blockSize<<" "<<cmprBuffer.size()<<" "<<d.size()<<std::endl;

    }

    treeIndex++;
    buffer.clear();
//    SZ_ReleaseHuffman();
//    SZ_Finalize();
    //std::cout<<"============  COMPRESS "<< buffer.size()<< " DONE =============="<<std::endl;
}

int CompressedPersistentLocalStore::get_compressible(unsigned int dbKey, uint64_t &key, std::vector<double> &data) {
    //std::cout<<"============ GET " <<dbKey <<" "<< key  << "=============="<<std::endl;
    auto k = std::make_pair(dbKey, key);
    if (storedKeys[dbKey].count(key) == 0) {
        return KEY_NOTFOUND;
    }

    if (buffer.count(k) > 0) {
        data = buffer[k];
    } else {

        if (szmap.count(huffmanIndexStore[k]) == 0) {
            auto conf = SZ::Config<double, 1>(realPrecision, std::array<size_t, 1>{block_size});
//            auto sz = SZ::make_sz_zone_compressor(conf, SZ::LorenzoPredictor<double, 1, 1>(realPrecision),
//                                                  SZ::LinearQuantizer<double>(realPrecision), SZ::HuffmanEncoder<int>());

            auto sz = new SZ::SZ_Zone_Compressor<double, 1, SZ::LorenzoPredictor<double, 1, 1>,
                    SZ::LinearQuantizer<double>, SZ::HuffmanEncoder<int>>
                    (conf, SZ::LorenzoPredictor<double, 1, 1>(realPrecision),
                     SZ::LinearQuantizer<double>(realPrecision), SZ::HuffmanEncoder<int>());

            RawDataVector tdata;
            huffmanStore.get(0, huffmanIndexStore[k], tdata);
            //std::cout<<"GET T: "<<treeLabel<<" "<<tdata.size()<<std::endl;
            sz->decompress_encoder((unsigned char *) &(tdata[0]), tdata.size());
            szmap[huffmanIndexStore[k]].reset(sz);
        }
        RawDataVector ddata;
        compressedStore.get(dbKey, key, ddata);
        szmap[huffmanIndexStore[k]]->decompress_zone((unsigned char *) &(ddata[0]), ddata.size(), data.data());


//        struct timespec start, end;
//        clock_gettime(CLOCK_REALTIME, &start);
//        clock_gettime(CLOCK_REALTIME, &end);
//        std::cout << "get2 time = "
//                  << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
//                  << "s" << std::endl;

    }

    //std::cout<<dd.size()<<" "<<s.getNAtoms()<<std::endl;

    //std::cout<<"============ GET " <<dbKey <<" "<< key  << " DONE =============="<<std::endl;

    return 0;
};


int CompressedPersistentLocalStore::get(unsigned int dbKey, uint64_t &key, RawDataVector &data) {
    //std::cout<<"============ GET " <<dbKey <<" "<< key  << "=============="<<std::endl;
    auto k = std::make_pair(dbKey, key);
//    if (uncompressedStore.count(dbKey, key) == 0) {
//        return KEY_NOTFOUND;
//    }

    RawDataVector d;
//    uncompressedStore.get(dbKey, key, d);
    //std::cout<<"GET UC: "<<d.size()<<std::endl;

//    SystemType s;
//    unpack(d,s);
    std::vector<double> dd;
    std::vector<unsigned char> dc;
    if (buffer.count(k) > 0) {
        dd = buffer[k];
    } else {
//        SZ_Init("sz.config");
        RawDataVector dt;
        Label treeLabel;
        int blockSize;
        double realPrecision2 = 0;
        compressedStore.get(dbKey, key, dt);
        //std::cout<<"GET C: "<<dt.size()<<std::endl;
//        unpack(dt,treeLabel,blockSize,dc);
        //std::cout<<"GET C: "<<dt.size()<<" "<<treeLabel<<" "<<blockSize<<std::endl;
        dd.resize(blockSize);


        RawDataVector tdata;
        int ddd = 0;
        huffmanStore.get(ddd, treeLabel, tdata);
        //std::cout<<"GET T: "<<treeLabel<<" "<<tdata.size()<<std::endl;

        unsigned char *p = (unsigned char *) &(tdata[0]);
//        node huffmanTreeRootNode = SZ_reconstruct_HuffmanEncoder(p, tdata.size(), &realPrecision2);


        size_t cmprBlockSize;
        //decompress
        p = (unsigned char *) &(dc[0]);
//        SZ_decompress_double_1D_exaalt_block(huffmanTreeRootNode, realPrecision2, p, dc.size(), &(dd[0]), blockSize);
//        SZ_ReleaseHuffman();
//        SZ_Finalize();
    }

    //std::cout<<dd.size()<<" "<<s.getNAtoms()<<std::endl;
    //reassemble item
//    s.restoreCompressibleData(dd);
//    pack(data,s);

    //std::cout<<"============ GET " <<dbKey <<" "<< key  << " DONE =============="<<std::endl;

    return 0;
};

int CompressedPersistentLocalStore::createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet) {
    return 0;
};

int CompressedPersistentLocalStore::sync() {
    compressedStore.sync();
//    uncompressedStore.sync();
    huffmanStore.sync();
    return 0;
};

int CompressedPersistentLocalStore::erase(unsigned int dbKey, uint64_t &key) {
    return 0;
}

int main(int argc, char **argv) {
    CompressedPersistentLocalStore store;
    size_t block_size = 3137;

    std::string dbdir(argv[1]);
    std::string datafile(argv[2]);
    int buffersize = 100;
    double eb = 0.0001;
    if (argc > 3) {
        buffersize = atoi(argv[3]);
        eb = atof(argv[4]);
    }

//    std::filesystem::remove_all(dbdir);
    store.initialize(dbdir, "exaalt_sz", block_size, eb, buffersize);
    size_t input_num;
//    auto input = SZ::readfile<float>(std::string(homedir + "/data/exaalt-2869440/xx.dat2").c_str(), input_num);
    auto input = SZ::readfile<float>(datafile.c_str(), input_num);

    srand(time(0));
    int testcase = input_num / block_size;
    std::vector<std::vector<double>> raw_double(testcase);
    std::vector<RawDataVector> raw_char(testcase);
    std::vector<uint64_t> keys(testcase);
    for (int i = 0; i < testcase; i++) {
        raw_double[i] = std::vector<double>(input.get() + i * block_size, input.get() + (i + 1) * block_size);
        raw_char[i] = std::vector<char>((char *) raw_double[i].data(),
                                        (char *) (raw_double[i].data() + raw_double[i].size()));
        keys[i] = rand() % testcase;
//        keys[i] = i;
    }

    std::cout << "data size: " << input_num << ", block size: " << block_size
              << ", number of blocks: " << testcase
              << ", buffer size: " << buffersize
              << ", error bound: " << eb
              << std::endl;

    std::cout << "raw_double vector size: " << raw_double[0].size() << std::endl;
    std::cout << "raw_char vector size: " << raw_char[0].size() << std::endl;
//    std::cout << input[0] << "," << input[1] << std::endl;
//    std::cout << raw_double[0][0] << "," << raw_double[0][1] << std::endl;
    std::vector<double> ori(testcase * block_size), dec(testcase * block_size);

    std::cout << "==============buffer with compression===============" << std::endl;

    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    for (int i = 0; i < testcase; i++) {
        store.put_compressible(1, keys[i], raw_double[i]);
    }
    store.compressBuffer();
    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "put time = "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
              << "s" << std::endl;


    clock_gettime(CLOCK_REALTIME, &start);
    std::vector<double> tmp_double(block_size);
    for (int i = 0; i < testcase; i++) {
        store.get_compressible(1, keys[i], tmp_double);
        memcpy(&dec[i * block_size], tmp_double.data(), tmp_double.size() * sizeof(double));
    }

    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "get time = "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
              << "s" << std::endl;

    std::cout << "==============buffer without compression===============" << std::endl;

    PersistentLocalStore pstore;
    pstore.initialize(dbdir, "exaalt");
    clock_gettime(CLOCK_REALTIME, &start);

    for (int i = 0; i < testcase; i++) {
        pstore.put(1, keys[i], raw_char[i]);
    }
    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "put time = "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
              << "s" << std::endl;


    clock_gettime(CLOCK_REALTIME, &start);
    RawDataVector tmp_char(raw_char[0].size());
    for (int i = 0; i < testcase; i++) {
        pstore.get(1, keys[i], tmp_char);
        memcpy(&ori[i * block_size], tmp_char.data(), tmp_char.size());
    }

    clock_gettime(CLOCK_REALTIME, &end);
    std::cout << "get time = "
              << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
              << "s" << std::endl;

    std::cout << "==============verification===============" << std::endl;
    SZ::verify(ori.data(), dec.data(), ori.size());
}