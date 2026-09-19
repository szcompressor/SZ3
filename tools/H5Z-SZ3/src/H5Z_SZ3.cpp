#include "H5Z_SZ3.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <memory>


const H5Z_class2_t H5Z_SZ3[1] = {{
    H5Z_CLASS_T_VERS,                                       /* H5Z_class_t version */
    H5Z_FILTER_SZ3,                                         /* Filter id number */
    1,                                                      /* encoder_present flag (set to true) */
    1,                                                      /* decoder_present flag (set to true) */
    "H5Z-SZ3-" SZ3_VER,                                     /* Filter name for debugging */
    NULL,                                                   /* The "can apply" callback */
    H5Z_sz3_set_local,                                      /* The "set local" callback */
    static_cast<H5Z_func_t>(H5Z_filter_sz3),                /* The actual filter function */
}};

HDF5SZ3_EXPORT H5PL_type_t H5PLget_plugin_type(void) { return H5PL_TYPE_FILTER; }

HDF5SZ3_EXPORT const void* H5PLget_plugin_info(void) { return H5Z_SZ3; }

// Whether this library registered the filter. Anything else's registration is not ours to undo.
static int h5z_sz3_was_registered = 0;

herr_t H5Z_SZ3_initialize(void) {
    // H5Zfilter_avail registers on a miss, from HDF5_PLUGIN_PATH; it is not a passive query.
    if (H5Zfilter_avail(H5Z_FILTER_SZ3) > 0) {
        return 0;
    }
    if (0 > H5Zregister(H5Z_SZ3)) {
        return -1;
    }
    h5z_sz3_was_registered = 1;
    return 1;
}

herr_t H5Z_SZ3_finalize(void) {
    herr_t ret = 0;
    if (h5z_sz3_was_registered) {
        ret = H5Zunregister(H5Z_FILTER_SZ3);
    }
    h5z_sz3_was_registered = 0;
    return ret < 0 ? -1 : 1;
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

    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts, 0);
    auto buffer = reinterpret_cast<unsigned char*>(cd_values.data());

    // conf.save(buffer);
    auto confSizeReal = conf.save(buffer);
    cd_nelmts = std::ceil(confSizeReal / 1.0 / sizeof(int));

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

herr_t get_SZ3_conf_from_H5(const hid_t propertyList, SZ3::Config& conf) {
    // static char const* _funcname_ = "get_SZ3_conf_from_H5";

    size_t cd_nelmts = std::ceil(conf.size_est() / 1.0 / sizeof(int));
    std::vector<unsigned int> cd_values(cd_nelmts, 0);

    // Check both returns: on failure cd_nelmts keeps its input value and conf would load zeros.
    if (!sz3_filter_on_plist(propertyList)) {
        return 0;
    }
    if (0 > H5Pget_filter_by_id(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, &cd_nelmts, cd_values.data(), 0, NULL,
                                NULL)) {
        return -1;
    }
    if (cd_nelmts > 0) {
        auto buffer = reinterpret_cast<const unsigned char*>(cd_values.data());
        size_t cd_bytes = cd_nelmts * sizeof(unsigned int);
        conf.load(buffer, cd_bytes);
    }
    return 1;
}

static herr_t H5Z_sz3_set_local_impl(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id) {
    // printf("start H5Z_sz3_set_local\n");

    // printf("start in H5Z_sz3_set_local, dcpl_id = %d\n", dcpl_id);
    static char const* _funcname_ = "H5Z_sz3_set_local";

    // herr_t ret = H5Zregister(H5Z_SZ3);

    SZ3::Config conf;
    get_SZ3_conf_from_H5(dcpl_id, conf);

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

    auto buffer = reinterpret_cast<const unsigned char*>(cd_values);
    size_t cd_bytes = cd_nelmts * sizeof(unsigned int);
    // Ahead of conf.load: every chunk this filter writes carries an SZ3 header, and a file from
    // another version wrote cd_values in a layout this build would misread.
    if (is_decompress) {
        if (nbytes < 8) throw std::invalid_argument("SZ3 HDF5 filter: chunk is smaller than an SZ3 header");
        auto header = reinterpret_cast<const unsigned char*>(*buf);
        uint32_t magic = 0, dataVer = 0;
        SZ3::read(magic, header);
        SZ3::read(dataVer, header);
        if (magic != SZ3_MAGIC_NUMBER)
            throw std::invalid_argument("SZ3 HDF5 filter: chunk was not written by SZ3");
        if (versionStr(dataVer) != SZ3_DATA_VER)
            throw std::invalid_argument("SZ3 HDF5 filter: data is in SZ3 data format v" + versionStr(dataVer) +
                                        ", this build reads v" SZ3_DATA_VER);
    }

    conf.load(buffer, cd_bytes);

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
