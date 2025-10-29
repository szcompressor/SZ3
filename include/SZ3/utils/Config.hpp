/**
 * @file Config.hpp
 * @brief This file defines the configuration class and related utilities for SZ3.
 *
 * @note Do not delete any existing fields for backward compatibility.
 *       These fields are critical for ensuring that serialized configurations
 *       from older versions of SZ3 remain compatible with newer versions.
 *       Removing or altering fields could break deserialization or cause
 *       unexpected behavior in applications relying on SZ3.
 */

#ifndef SZ3_Config_HPP
#define SZ3_Config_HPP

#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>

#include "SZ3/def.hpp"
#include "SZ3/utils/ByteUtil.hpp"
#include "SZ3/utils/MemoryUtil.hpp"
#include "SZ3/version.hpp"

#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9


namespace SZ3 {

/**
 * @enum EB
 * @brief Enumeration for error bound modes.
 *
 * - EB_ABS: Absolute error bound.
 * - EB_REL: Relative error bound.
 * - EB_PSNR: Peak Signal-to-Noise Ratio error bound.
 * - EB_L2NORM: L2 norm error bound.
 * - EB_ABS_AND_REL: Combined absolute and relative error bound.
 * - EB_ABS_OR_REL: Either absolute or relative error bound.
 */
enum EB { EB_ABS, EB_REL, EB_PSNR, EB_L2NORM, EB_ABS_AND_REL, EB_ABS_OR_REL };

/**
 * @enum ALGO
 * @brief Enumeration for compression algorithms.
 *
 * - ALGO_LORENZO_REG: SZ2 algorithm (Lorenzo with regression).
 * - ALGO_INTERP_LORENZO: SZ3 algorithm (Interpolation with Lorenzo, with tuning).
 * - ALGO_INTERP: Interpolation-only (without tuning).
 * - ALGO_NOPRED: No prediction.
 * - ALGO_LOSSLESS: Lossless compression.
 */
enum ALGO { ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_NOPRED, ALGO_LOSSLESS };

/**
 * @enum INTERP_ALGO
 * @brief Enumeration for interpolation algorithms.
 *
 * - INTERP_ALGO_LINEAR: Linear interpolation.
 * - INTERP_ALGO_CUBIC: Cubic interpolation.
 */
enum INTERP_ALGO { INTERP_ALGO_LINEAR, INTERP_ALGO_CUBIC };

const std::map<std::string, ALGO> ALGO_MAP = {
    {"ALGO_LORENZO_REG", ALGO_LORENZO_REG}, {"ALGO_INTERP_LORENZO", ALGO_INTERP_LORENZO},
    {"ALGO_INTERP", ALGO_INTERP}, {"ALGO_NOPRED", ALGO_NOPRED},
    {"ALGO_LOSSLESS", ALGO_LOSSLESS},
};

const std::map<std::string, EB> EB_MAP = {
    {"ABS", EB_ABS},
    {"REL", EB_REL},
    {"PSNR", EB_PSNR},
    {"NORM", EB_L2NORM},
    {"ABS_AND_REL", EB_ABS_AND_REL},
    {"ABS_OR_REL", EB_ABS_OR_REL},
};

const std::map<std::string, INTERP_ALGO> INTERP_ALGO_MAP = {
    {"INTERP_ALGO_LINEAR", INTERP_ALGO_LINEAR},
    {"INTERP_ALGO_CUBIC", INTERP_ALGO_CUBIC},
};

ALWAYS_INLINE std::string to_lower(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

template <typename EnumType>
ALWAYS_INLINE void match_enum(const std::string& input, const std::map<std::string, EnumType>& table, uint8_t& out) {
    std::string input_lc = to_lower(input);
    for (const auto& [key, val] : table) {
        if (to_lower(key) == input_lc) {
            out = static_cast<int>(val);
        }
    }
}

template <typename EnumType>
ALWAYS_INLINE std::string enum_to_string(EnumType value, const std::map<std::string, EnumType>& forward_map) {
    for (const auto& [str, val] : forward_map) {
        if (val == value) {
            return str;
        }
    }
    return "";
}

/**
 * @class Config
 * @brief Configuration class for SZ3.
 *
 * The Config class stores various parameters and settings used by SZ3, such as:
 * - Dimensions of the data.
 * - Compression algorithm and error bound settings.
 * - Module-specific parameters like interpolation settings.
 *
 * It provides methods to load/save configurations from/to INI files, serialize/deserialize configurations,
 * and estimate the size of the configuration.
 */
class Config {
public:
    /**
     * @brief Constructor to initialize the configuration with dimensions.
     *
     * @tparam Dims Variadic template for dimension arguments.
     * @param args Dimensions of the data.
     */
    template <class... Dims>
    Config(Dims... args) {
        dims = std::vector<size_t>{static_cast<size_t>(std::forward<Dims>(args))...};
        setDims(dims.begin(), dims.end());
    }

    /**
     * @brief Set the dimensions of the data.
     *
     * @tparam Iter Iterator type for the dimension range.
     * @param begin Iterator to the beginning of the dimensions.
     * @param end Iterator to the end of the dimensions.
     * @return Total number of elements in the data.
     */
    template <class Iter>
    size_t setDims(Iter begin, Iter end) {
        auto dims_ = std::vector<size_t>(begin, end);
        dims.clear();
        for (auto dim : dims_) {
            if (dim > 1) {
                dims.push_back(dim);
            }
        }
        if (dims.empty()) {
            dims = {1};
        }
        N = dims.size();
        num = std::accumulate(dims.begin(), dims.end(), static_cast<size_t>(1), std::multiplies<size_t>());
        predDim = N;
        blockSize = (N == 1 ? 128 : (N == 2 ? 16 : 6));
        return num;
    }

    /**
     * @brief Load configuration from an INI file.
     *
     * @param ini_file_path Path to the INI file.
     * @throws std::runtime_error If the file cannot be opened.
     */
    void loadcfg(const std::string& ini_file_path) {
        std::ifstream file(ini_file_path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open INI file: " + ini_file_path);
        }
        std::ostringstream buffer;
        buffer << file.rdbuf();
        load_ini(buffer.str());
    }

    /**
     * @brief Load configuration from an INI content string.
     *
     * @param ini_content Content of the INI file as a string.
     */
    void load_ini(const std::string& ini_content) {
        std::istringstream ss(ini_content);
        std::string line, section;

        auto trim = [](std::string& s) {
            s.erase(0, s.find_first_not_of(" \t\r\n"));
            s.erase(s.find_last_not_of(" \t\r\n") + 1);
        };

        auto eq = [&](const std::string& a, const std::string& b) { return to_lower(a) == to_lower(b); };

        auto parse_bool = [&](const std::string& s) {
            auto ls = to_lower(s);
            return ls == "true" || ls == "1" || ls == "yes" || ls == "on";
        };

        while (std::getline(ss, line)) {
            trim(line);
            if (line.empty() || line[0] == '#') continue;
            if (line.front() == '[') {
                section = line.substr(1, line.find(']') - 1);
                continue;
            }

            auto sep = line.find('=');
            if (sep == std::string::npos) continue;

            std::string key = line.substr(0, sep);
            std::string value = line.substr(sep + 1);
            trim(key);
            trim(value);

            if (eq(section, "GlobalSettings")) {
                if (eq(key, "CmprAlgo"))
                    match_enum(value, ALGO_MAP, cmprAlgo);
                else if (eq(key, "ErrorBoundMode"))
                    match_enum(value, EB_MAP, errorBoundMode);
                else if (eq(key, "AbsErrorBound"))
                    absErrorBound = std::stod(value);
                else if (eq(key, "RelErrorBound"))
                    relErrorBound = std::stod(value);
                else if (eq(key, "PSNRErrorBound"))
                    psnrErrorBound = std::stod(value);
                else if (eq(key, "L2NormErrorBound"))
                    l2normErrorBound = std::stod(value);
                else if (eq(key, "OpenMP"))
                    openmp = parse_bool(value);
            } else if (eq(section, "AlgoSettings")) {
                if (eq(key, "Lorenzo"))
                    lorenzo = parse_bool(value);
                else if (eq(key, "Lorenzo2ndOrder"))
                    lorenzo2 = parse_bool(value);
                else if (eq(key, "Regression"))
                    regression = parse_bool(value);
                else if (eq(key, "Regression2ndOrder"))
                    regression2 = parse_bool(value);
                else if (eq(key, "InterpolationAlgo"))
                    match_enum(value, INTERP_ALGO_MAP, interpAlgo);
                else if (eq(key, "InterpolationDirection"))
                    interpDirection = std::stoi(value);
                else if (eq(key, "BlockSize"))
                    blockSize = std::stoi(value);
                else if (eq(key, "QuantizationBinTotal"))
                    quantbinCnt = std::stoi(value);
                else if (eq(key, "InterpolationAnchorStride"))
                    interpAnchorStride = std::stoi(value);
                else if (eq(key, "InterpolationAlpha"))
                    interpAlpha = std::stod(value);
                else if (eq(key, "InterpolationBeta"))
                    interpBeta = std::stod(value);
            }
        }
    }

    /**
     * @brief Save the current configuration to an INI format string.
     *
     * @return INI format string representing the current configuration.
     */
    std::string save_ini() const {
        std::ostringstream ss;

        ss << "[GlobalSettings]\n";
        ss << "CmprAlgo = " << enum_to_string(static_cast<ALGO>(cmprAlgo), ALGO_MAP) << "\n";
        ss << "ErrorBoundMode = " << enum_to_string(static_cast<EB>(errorBoundMode), EB_MAP) << "\n";
        ss << "AbsErrorBound = " << absErrorBound << "\n";
        ss << "RelErrorBound = " << relErrorBound << "\n";
        ss << "PSNRErrorBound = " << psnrErrorBound << "\n";
        ss << "L2NormErrorBound = " << l2normErrorBound << "\n";
        ss << "OpenMP = " << (openmp ? "true" : "false") << "\n";

        ss << "\n[AlgoSettings]\n";
        ss << "Lorenzo = " << (lorenzo ? "true" : "false") << "\n";
        ss << "Lorenzo2ndOrder = " << (lorenzo2 ? "true" : "false") << "\n";
        ss << "Regression = " << (regression ? "true" : "false") << "\n";
        ss << "Regression2ndOrder = " << (regression2 ? "true" : "false") << "\n";
        ss << "BlockSize = " << blockSize << "\n";
        ss << "QuantizationBinTotal = " << quantbinCnt << "\n";
        ss << "InterpolationAlgo = " << enum_to_string(static_cast<INTERP_ALGO>(interpAlgo), INTERP_ALGO_MAP) << "\n";
        ss << "InterpolationDirection = " << static_cast<int>(interpDirection) << "\n";
        ss << "InterpolationAnchorStride = " << interpAnchorStride << "\n";
        ss << "InterpolationAlpha = " << interpAlpha << "\n";
        ss << "InterpolationBeta = " << interpBeta << "\n";
        return ss.str();
    }

    /**
     * @brief Serialize the configuration to a byte array.
     *
     * @param c Pointer to the byte array.
     * @return Size of the serialized configuration.
     */
    size_t save(unsigned char*& c) const {
        auto c0 = c;
        c += sizeof(uchar); //reserve space for conf size

        write(sz3MagicNumber, c);
        write(sz3DataVer, c);
        write(N, c);
        // write(dims.data(), dims.size(), c);
        auto bitWidth = vector_bit_width(dims);
        write(bitWidth, c);
        vector2bytes(dims, bitWidth, c);

        write(num, c);
        write(cmprAlgo, c);

        write(errorBoundMode, c);
        if (errorBoundMode == EB_ABS) {
            write(absErrorBound, c);
        } else if (errorBoundMode == EB_REL) {
            write(relErrorBound, c);
        } else if (errorBoundMode == EB_PSNR) {
            write(psnrErrorBound, c);
        } else if (errorBoundMode == EB_L2NORM) {
            write(l2normErrorBound, c);
        } else if (errorBoundMode == EB_ABS_OR_REL) {
            write(absErrorBound, c);
            write(relErrorBound, c);
        } else if (errorBoundMode == EB_ABS_AND_REL) {
            write(absErrorBound, c);
            write(relErrorBound, c);
        }

        uint8_t boolvals = (lorenzo & 1) << 7 | (lorenzo2 & 1) << 6 | (regression & 1) << 5 | (regression2 & 1) << 4 |
                           (openmp & 1) << 3;
        write(boolvals, c);

        write(dataType, c);
        write(quantbinCnt, c);
        write(blockSize, c);
        write(predDim, c);

        auto confSize = static_cast<uchar>(c - c0);
        write(confSize, c0); //write conf size at reserved space
        return confSize;
    }

    /**
     * @brief Deserialize the configuration from a byte array.
     *
     * @param c Pointer to the byte array.
     */
    void load(const unsigned char*& c) {
        uchar confSize = 0;
        read(confSize, c);
        auto c1 = c + confSize;
        read(sz3MagicNumber, c);
        read(sz3DataVer, c);

        read(N, c);
        uint8_t bitWidth;
        read(bitWidth, c);
        dims = bytes2vector<size_t>(c, bitWidth, N);
        // dims.resize(N);
        // read(dims.data(), N, c);
        read(num, c);
        read(cmprAlgo, c);

        read(errorBoundMode, c);
        if (errorBoundMode == EB_ABS) {
            read(absErrorBound, c);
        } else if (errorBoundMode == EB_REL) {
            read(relErrorBound, c);
        } else if (errorBoundMode == EB_PSNR) {
            read(psnrErrorBound, c);
        } else if (errorBoundMode == EB_L2NORM) {
            read(l2normErrorBound, c);
        } else if (errorBoundMode == EB_ABS_OR_REL) {
            read(absErrorBound, c);
            read(relErrorBound, c);
        } else if (errorBoundMode == EB_ABS_AND_REL) {
            read(absErrorBound, c);
            read(relErrorBound, c);
        }

        if (c < c1) {
            uint8_t boolvals;
            read(boolvals, c);
            lorenzo = (boolvals >> 7) & 1;
            lorenzo2 = (boolvals >> 6) & 1;
            regression = (boolvals >> 5) & 1;
            regression2 = (boolvals >> 4) & 1;
            openmp = (boolvals >> 3) & 1;
        }
        if (c < c1) {
            read(dataType, c);
        }
        if (c < c1) {
            read(quantbinCnt, c);
        }
        if (c < c1) {
            read(blockSize, c);
        }
        if (c < c1) {
            read(predDim, c);
        }
    }

    /**
     * @brief Print the current configuration to the console.
     */
    void print() {
        // printf("===================== Begin SZ3 Configuration =====================\n");
        printf("\nsz3MagicNumber = %u\n", sz3MagicNumber);
        printf("sz3DataVer = %s\n", versionStr(sz3DataVer).data());
        std::cout << "Dimensions =";
        for (auto d : dims) {
            std::cout << " " << d;
        }
        std::cout << "\n";
        std::cout << save_ini();
    }

    /**
     * @brief Estimate the size of the serialized configuration.
     *
     * @return Estimated size in bytes.
     */
    size_t size_est() const {
        std::vector<uchar> buffer(sizeof(Config) + 1024);
        auto buffer_pos = buffer.data();
        return save(buffer_pos);
    }

    uint32_t sz3MagicNumber = SZ3_MAGIC_NUMBER;
    uint32_t sz3DataVer = versionInt(SZ3_DATA_VER);
    char N = 0;
    std::vector<size_t> dims;
    size_t num = 0;
    uint8_t cmprAlgo = ALGO_INTERP_LORENZO;
    uint8_t errorBoundMode = EB_ABS;
    double absErrorBound = 1e-3;
    double relErrorBound = 0.0;
    double psnrErrorBound = 0.0;
    double l2normErrorBound = 0.0;
    bool lorenzo = true;
    bool lorenzo2 = false;
    bool regression = true;
    bool regression2 = false;
    bool openmp = false;
    int quantbinCnt = 65536;
    int blockSize = 0;
    uint8_t predDim = 0;         // not used now
    uint8_t dataType = SZ_FLOAT; // dataType is only used in HDF5 filter

    /**
     * The following parameters are only used by specific modules.
     * They are saved and loaded in their module implementation, not in Config.
     */
    uint8_t interpAlgo = INTERP_ALGO_CUBIC;
    uint8_t interpDirection = 0;
    int interpAnchorStride = -1; // -1: using dynamic default setting
    double interpAlpha = 1.25;   // level-wise eb reduction rate
    double interpBeta = 2.0;     // maximum eb reduction rate
};

} // namespace SZ3

#endif  // SZ_CONFIG_HPP

