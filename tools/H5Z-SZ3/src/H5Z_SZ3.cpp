#include "H5Z_SZ3.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

// MSG is a printf format string, so anything that is not a literal goes through "%s".
#define H5Z_SZ_PUSH_AND_GOTO(MAJ, MIN, RET, ...)                                                  \
    do {                                                                                          \
        H5Epush(H5E_DEFAULT, __FILE__, _funcname_, __LINE__, H5E_ERR_CLS, MAJ, MIN, __VA_ARGS__); \
        return RET;                                                                               \
    } while (0)

// The filter's own callbacks. HDF5 reaches them through H5Z_SZ3 below, so they stay out of the
// public header, where every consumer of the filter would see them declared and never defined.
extern "C" {
static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id);

static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes,
                             size_t* buf_size, void** buf);
}

const H5Z_class2_t H5Z_SZ3[1] = {{
    H5Z_CLASS_T_VERS,                                       /* H5Z_class_t version */
    H5Z_FILTER_SZ3,                                         /* Filter id number */
    1,                                                      /* encoder_present flag (set to true) */
    1,                                                      /* decoder_present flag (set to true) */
    // Written into every dataset that uses the filter, and quoted back by HDF5 when the filter
    // is missing, so it has to say where to get it.
    "H5Z-SZ3-" SZ3_VER " (data format " SZ3_DATA_VER "); see "
    "https://github.com/szcompressor/SZ3/tree/master/tools/H5Z-SZ3",
    NULL,                                                   /* The "can apply" callback */
    H5Z_sz3_set_local,                                      /* The "set local" callback */
    static_cast<H5Z_func_t>(H5Z_filter_sz3),                /* The actual filter function */
}};

HDF5SZ3_EXPORT H5PL_type_t H5PLget_plugin_type(void) { return H5PL_TYPE_FILTER; }

HDF5SZ3_EXPORT const void* H5PLget_plugin_info(void) { return H5Z_SZ3; }

static std::vector<unsigned int> save_cd_values(const SZ3::Config& conf) {
    std::vector<unsigned char> bytes(conf.size_est());
    auto pos = bytes.data();
    size_t len = conf.save(pos);
    std::vector<unsigned int> cd(1 + (len + sizeof(unsigned int) - 1) / sizeof(unsigned int), 0);
    cd[0] = versionInt(SZ3_DATA_VER);
    memcpy(cd.data() + 1, bytes.data(), len);
    return cd;
}

static void load_cd_values(const unsigned int* cd, size_t cd_nelmts, SZ3::Config& conf) {
    auto bytes = reinterpret_cast<const unsigned char*>(cd);
    size_t len = cd_nelmts * sizeof(unsigned int);
    // 3.3.2 stored the Config without the version in front. A Config starts with its length, which
    // is never 0, and versionInt() always leaves the low byte 0: that byte tells the two apart.
    if (cd_nelmts > 0 && (cd[0] & 0xFFu) == 0) {
        if (versionStr(cd[0]) != SZ3_DATA_VER) {
            throw std::invalid_argument("SZ3 " SZ3_VER " reads data version " SZ3_DATA_VER
                                        ", but this dataset's SZ3 settings are version " +
                                        versionStr(cd[0]) + ". Use SZ3 v" + versionStr(cd[0]) + ".");
        }
        bytes += sizeof(unsigned int);
        len -= sizeof(unsigned int);
    }
    conf.load(bytes, len);
}

// Do not use H5Zfilter_avail() here: it answers for the library, not for this property list.
static bool sz3_filter_on_plist(const hid_t propertyList) {
    const int nfilters = H5Pget_nfilters(propertyList);
    for (int i = 0; i < nfilters; i++) {
        unsigned int flags = 0;
        size_t cd_nelmts = 0;
        unsigned int cd_values[1] = {0};
        if (H5Z_FILTER_SZ3 ==
            H5Pget_filter2(propertyList, static_cast<unsigned>(i), &flags, &cd_nelmts, cd_values, 0, NULL, NULL)) {
            return true;
        }
    }
    return false;
}

herr_t set_SZ3_conf_to_H5(const hid_t propertyList, SZ3::Config& conf) {
    static char const* _funcname_ = "set_SZ3_conf_to_H5";

    std::vector<unsigned int> cd_values = save_cd_values(conf);
    size_t cd_nelmts = cd_values.size();

    if (sz3_filter_on_plist(propertyList)) {
        if (0 > H5Pmodify_filter(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data())) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
        }
    } else {
        // filter not set, set filter. Notice that calling H5Pset_filter twice with the same filter id
        // will cause unexpected errors for decompression
        if (0 > H5Pset_filter(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values.data())) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
        }
    }

    return 1;
}

static_assert(H5Z_SZ3_EB_ABS == SZ3::EB_ABS && H5Z_SZ3_EB_ABS_OR_REL == SZ3::EB_ABS_OR_REL,
              "H5Z_SZ3.hpp error-bound modes must match SZ3::EB");
static_assert(H5Z_SZ3_ALGO_LORENZO_REG == SZ3::ALGO_LORENZO_REG && H5Z_SZ3_ALGO_BIOMDXTC == SZ3::ALGO_BIOMDXTC,
              "H5Z_SZ3.hpp algorithms must match SZ3::ALGO");

herr_t H5Pset_sz3(hid_t plist, int algo, int eb_mode, double abs_bound, double rel_bound, double psnr_bound,
                  double l2norm_bound) {
    static char const* _funcname_ = "H5Pset_sz3";
    if (algo < H5Z_SZ3_ALGO_LORENZO_REG || algo > H5Z_SZ3_ALGO_BIOMDXTC) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "unknown SZ3 algorithm %d", algo);
    }
    if (eb_mode < H5Z_SZ3_EB_ABS || eb_mode > H5Z_SZ3_EB_ABS_OR_REL) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "unknown SZ3 error-bound mode %d", eb_mode);
    }
    SZ3::Config conf;
    conf.cmprAlgo = static_cast<uint8_t>(algo);
    conf.errorBoundMode = static_cast<uint8_t>(eb_mode);
    conf.absErrorBound = abs_bound;
    conf.relErrorBound = rel_bound;
    conf.psnrErrorBound = psnr_bound;
    conf.l2normErrorBound = l2norm_bound;
    return set_SZ3_conf_to_H5(plist, conf) > 0 ? 1 : -1;
}

herr_t get_SZ3_conf_from_H5(const hid_t propertyList, SZ3::Config& conf) {
    static char const* _funcname_ = "get_SZ3_conf_from_H5";

    // H5Pget_filter_by_id fails when the list carries no SZ3 filter, which is not an error to report.
    if (!sz3_filter_on_plist(propertyList)) {
        return 0;
    }
    // Ask for the count before the values: a stored config carries the dataset's dimensions, so its
    // length is not the caller's to guess, and HDF5 reports the real count even when it filled less.
    size_t cd_nelmts = 0;
    unsigned int probe[1] = {0};
    if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, &cd_nelmts, probe, 0, NULL, NULL)) {
        return -1;
    }
    if (cd_nelmts > 0) {
        std::vector<unsigned int> cd_values(cd_nelmts, 0);
        if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, &cd_nelmts, cd_values.data(), 0,
                                    NULL, NULL)) {
            return -1;
        }
        try {
            SZ3::Config loaded;
            load_cd_values(cd_values.data(), cd_values.size(), loaded);
            conf = loaded;
        } catch (const std::exception& e) {
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, -1, "%s", e.what());
        }
    }
    return 1;
}

static herr_t H5Z_sz3_set_local_impl(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    // printf("start H5Z_sz3_set_local\n");

    // printf("start in H5Z_sz3_set_local, dcpl_id = %d\n", dcpl_id);
    static char const* _funcname_ = "H5Z_sz3_set_local";

    // herr_t ret = H5Zregister(H5Z_SZ3);

    SZ3::Config conf;
    if (get_SZ3_conf_from_H5(dcpl_id, conf) < 0) return -1;

    // read datatype and dims from HDF5
    H5T_class_t dclass;
    if (0 > (dclass = H5Tget_class(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

    size_t dsize;
    if (0 == (dsize = H5Tget_size(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

    int ndims;
    hsize_t dims_all[H5S_MAX_RANK];
    if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims_all, NULL)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");
    std::vector<size_t> dims(dims_all, dims_all + ndims);
    // update conf with datatype
    conf.dataType = SZ_FLOAT;
    if (dclass == H5T_FLOAT)
        conf.dataType = dsize == 4 ? SZ_FLOAT : SZ_DOUBLE;
    else if (dclass == H5T_INTEGER) {
        H5T_sign_t dsign;
        if (0 > (dsign = H5Tget_sign(type_id)))
            H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");
        if (dsign == H5T_SGN_NONE) // unsigned
        {
            switch (dsize) {
                case 1:
                    conf.dataType = SZ_UINT8;
                    break;
                case 2:
                    conf.dataType = SZ_UINT16;
                    break;
                case 4:
                    conf.dataType = SZ_UINT32;
                    break;
                case 8:
                    conf.dataType = SZ_UINT64;
                    break;
            }
        } else {
            switch (dsize) {
                case 1:
                    conf.dataType = SZ_INT8;
                    break;
                case 2:
                    conf.dataType = SZ_INT16;
                    break;
                case 4:
                    conf.dataType = SZ_INT32;
                    break;
                case 8:
                    conf.dataType = SZ_INT64;
                    break;
            }
        }
    } else {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
    }
    // update conf with dims
    conf.setDims(std::begin(dims), std::end(dims));
    //  need to update magic number and data version,
    //  as the config may be from cd_values passed by users
    conf.sz3MagicNumber = SZ3_MAGIC_NUMBER;
    conf.sz3DataVer = versionInt(SZ3_DATA_VER);

    set_SZ3_conf_to_H5(dcpl_id, conf);
    return 1;
}

template <typename T>
void process_data(SZ3::Config& conf, void** buf, size_t* buf_size, size_t nbytes, bool is_decompress) {
    if (is_decompress) {
        T* processedData = static_cast<T*>(malloc(conf.num * sizeof(T)));
        // HDF5 frees what this returns, so it has to come from malloc. On null SZ_decompress would
        // allocate with new[] instead, and that pairing is undefined.
        if (processedData == nullptr) throw std::bad_alloc();
        SZ_decompress(conf, static_cast<char*>(*buf), nbytes, processedData);
        free(*buf);
        *buf = processedData;
        *buf_size = conf.num * sizeof(T);
    } else {
        // The bound assumes the payload fits in the raw size, so leave headroom on top of it for
        // algorithms whose output can reach or exceed that.
        size_t cmpCap = std::max(SZ3::SZ_compress_size_bound<T>(conf), sizeof(T) * conf.num * 2);
        char* cmpData = static_cast<char*>(malloc(cmpCap));
        *buf_size = SZ_compress(conf, static_cast<T*>(*buf), cmpData, cmpCap);
        free(*buf);
        *buf = cmpData;
    }
}

/**
 * https://docs.hdfgroup.org/hdf5/v1_14/_f_i_l_t_e_r.html
 * The flags, cd_nelmts, and cd_values are the same as for the H5Pset_filter() function with the additional flag
 * H5Z_FLAG_REVERSE which is set when the filter is called as part of the input pipeline. The input buffer is pointed to
 * by *buf and has a total size of *buf_size bytes but only nbytes are valid data. The filter should perform the
 * transformation in place if possible and return the number of valid bytes or zero for failure. If the transformation
 * cannot be done in place then the filter should allocate a new buffer with malloc() and assign it to *buf, assigning
 * the allocated size of that buffer to *buf_size. The old buffer should be freed by calling free().
 */
static size_t H5Z_filter_sz3_impl(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[],
                                  size_t nbytes, size_t* buf_size, void** buf) {
    // printf("start H5Z_filter_sz3\n");

    if (cd_nelmts == 0) // this is special data such as string, which should not be treated as values.
        return nbytes;

    bool is_decompress = flags & H5Z_FLAG_REVERSE;
    SZ3::Config conf;

    // Ahead of conf.load: every chunk this filter writes carries an SZ3 header, and a file from
    // another version wrote cd_values in a layout this build would misread.
    if (is_decompress) {
        uint32_t magic = 0, dataVer = 0;
        if (nbytes >= 8) {
            auto header = reinterpret_cast<const unsigned char*>(*buf);
            SZ3::read(magic, header);
            SZ3::read(dataVer, header);
        }
        if (magic != SZ3_MAGIC_NUMBER) {
            // v3.2.0 through v3.3.2 wrote and read a chunk raw when cd_values held fewer than 20
            // elements, so those carry no header.
            SZ3::Config legacy;
            load_cd_values(cd_values, cd_nelmts, legacy);
            if (legacy.num > 0 && legacy.num < 20) return nbytes;
            throw std::invalid_argument("SZ3 HDF5 filter: chunk was not written by SZ3");
        }
        if (versionStr(dataVer) != SZ3_DATA_VER)
            throw std::invalid_argument("SZ3 HDF5 filter: data is in SZ3 data format v" + versionStr(dataVer) +
                                        ", this build reads v" SZ3_DATA_VER);
    }

    load_cd_values(cd_values, cd_nelmts, conf);

    switch (conf.dataType) {
        case SZ_FLOAT:
            process_data<float>(conf, buf, buf_size, nbytes, is_decompress);
            break;
#if (!SZ3_DEBUG_TIMINGS)
        case SZ_DOUBLE:
            process_data<double>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT8:
            process_data<int8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT8:
            process_data<uint8_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT16:
            process_data<int16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT16:
            process_data<uint16_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT32:
            process_data<int32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT32:
            process_data<uint32_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_INT64:
            process_data<int64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
        case SZ_UINT64:
            process_data<uint64_t>(conf, buf, buf_size, nbytes, is_decompress);
            break;
#endif
        default:
            throw std::invalid_argument("SZ3 HDF5 filter: unknown datatype in cd_values");
    }
    return *buf_size;
}

// HDF5 calls these from C, where an exception that escapes ends the application. Zero is how a
// filter reports failure; set_local uses a negative return.
static size_t H5Z_filter_sz3(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes,
                             size_t* buf_size, void** buf) {
    static char const* _funcname_ = "H5Z_filter_sz3";
    try {
        return H5Z_filter_sz3_impl(flags, cd_nelmts, cd_values, nbytes, buf_size, buf);
    } catch (const std::exception& e) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, 0, "%s", e.what());
    } catch (...) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, 0, "unknown error");
    }
}

static herr_t H5Z_sz3_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    static char const* _funcname_ = "H5Z_sz3_set_local";
    try {
        return H5Z_sz3_set_local_impl(dcpl_id, type_id, chunk_space_id);
    } catch (const std::exception& e) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, -1, "%s", e.what());
    } catch (...) {
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CALLBACK, -1, "unknown error");
    }
}
