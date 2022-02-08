//
// Created by Kai Zhao on 1/28/21.
//

#ifndef SZ3_BYTEUTIL_HPP
#define SZ3_BYTEUTIL_HPP

#include "SZ3/def.hpp"
#include <cstring>

namespace SZ {

    typedef union lint16 {
        unsigned short usvalue;
        short svalue;
        unsigned char byte[2];
    } lint16;

    typedef union lint32 {
        int ivalue;
        unsigned int uivalue;
        unsigned char byte[4];
    } lint32;

    typedef union lint64 {
        long lvalue;
        unsigned long ulvalue;
        unsigned char byte[8];
    } lint64;

    typedef union ldouble {
        double value;
        unsigned long lvalue;
        unsigned char byte[8];
    } ldouble;

    typedef union lfloat {
        float value;
        unsigned int ivalue;
        unsigned char byte[4];
        uint16_t int16[2];
    } lfloat;

    inline void symTransform_4bytes(uchar data[4]) {
        unsigned char tmp = data[0];
        data[0] = data[3];
        data[3] = tmp;

        tmp = data[1];
        data[1] = data[2];
        data[2] = tmp;
    }

    inline unsigned short bytesToUInt16_bigEndian(const uchar *bytes) {
        int temp = 0;
        unsigned short res = 0;

        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        return res;
    }

    inline unsigned int bytesToUInt32_bigEndian(const uchar *bytes) {
        unsigned int temp = 0;
        unsigned int res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline unsigned long bytesToUInt64_bigEndian(const uchar *b) {
        unsigned long temp = 0;
        unsigned long res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }

    inline short bytesToInt16_bigEndian(unsigned char *bytes) {
        int temp = 0;
        short res = 0;

        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        return res;
    }

    inline int bytesToInt32_bigEndian(unsigned char *bytes) {
        int temp = 0;
        int res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline long bytesToInt64_bigEndian(const unsigned char *b) {
        long temp = 0;
        long res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }

    inline int bytesToInt_bigEndian(const unsigned char *bytes) {
        int temp = 0;
        int res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline long bytesToLong_bigEndian(const uchar *b) {
        long temp = 0;
        long res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }

    inline void longToBytes_bigEndian(unsigned char *b, unsigned long num) {
        b[0] = (unsigned char) (num >> 56);
        b[1] = (unsigned char) (num >> 48);
        b[2] = (unsigned char) (num >> 40);
        b[3] = (unsigned char) (num >> 32);
        b[4] = (unsigned char) (num >> 24);
        b[5] = (unsigned char) (num >> 16);
        b[6] = (unsigned char) (num >> 8);
        b[7] = (unsigned char) (num);
    }

    inline void int16ToBytes_bigEndian(unsigned char *b, uint16_t num) {
        b[0] = (unsigned char) (num >> 8);
        b[1] = (unsigned char) (num);
    }

    inline void int32ToBytes_bigEndian(unsigned char *b, uint32_t num) {
        b[0] = (unsigned char) (num >> 24);
        b[1] = (unsigned char) (num >> 16);
        b[2] = (unsigned char) (num >> 8);
        b[3] = (unsigned char) (num);
    }

    inline void intToBytes_bigEndian(unsigned char *b, unsigned int num) {
        b[0] = (unsigned char) (num >> 24);
        b[1] = (unsigned char) (num >> 16);
        b[2] = (unsigned char) (num >> 8);
        b[3] = (unsigned char) (num);
    }

    inline void int64ToBytes_bigEndian(unsigned char *b, uint64_t num) {
        b[0] = (unsigned char) (num >> 56);
        b[1] = (unsigned char) (num >> 48);
        b[2] = (unsigned char) (num >> 40);
        b[3] = (unsigned char) (num >> 32);
        b[4] = (unsigned char) (num >> 24);
        b[5] = (unsigned char) (num >> 16);
        b[6] = (unsigned char) (num >> 8);
        b[7] = (unsigned char) (num);
    }

    int bytesToInt(const unsigned char *bytes) {
        lfloat buf;
        memcpy(buf.byte, bytes, 4);
        return buf.ivalue;
    }

    std::string floatToBinary(float f) {
        lfloat u;
        u.value = f;
        std::string str(32, '0');
        for (int i = 0; i < 32; i++) {
            str[31 - i] = (u.ivalue % 2) ? '1' : '0';
            u.ivalue >>= 1;
        }
        return str;
    }

    template<class T>
    void truncateArray(T data, size_t n, int byteLen, uchar *&binary) {
        lfloat bytes;
        int b;
        for (size_t i = 0; i < n; i++) {
            bytes.value = data[i];
            for (b = 4 - byteLen; b < 4; b++) {
                *binary++ = bytes.byte[b];
            }
        }
    }

    template<class T>
    void truncateArrayRecover(uchar *binary, size_t n, int byteLen, T *data) {
        lfloat bytes;
        bytes.ivalue = 0;
        int b;
        for (size_t i = 0; i < n; i++) {
            for (b = 4 - byteLen; b < 4; b++) {
                bytes.byte[b] = *binary++;
            }
            data[i] = bytes.value;
        }
    }

    std::vector<uchar> LeadingBitsEncode(float pre, float data) {
        lfloat lfBuf_pre;
        lfloat lfBuf_cur;

        lfBuf_pre.value = pre;
        lfBuf_cur.value = data;
        lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

        std::vector<uchar> bytes;
        int n = 0;
        if (lfBuf_pre.ivalue == 0) {
            n = 0;
        } else if (lfBuf_pre.ivalue >> 8 == 0) {
            n = 1;
        } else if (lfBuf_pre.ivalue >> 16 == 0) {
            n = 2;
        } else if (lfBuf_pre.ivalue >> 24 == 0) {
            n = 3;
        } else {
            n = 4;
        }

        for (int i = 0; i < n; i++) {
            bytes.push_back(lfBuf_cur.byte[i]);
        }
        return bytes;
    }

    float LeadingBitsDecode(float pre, std::vector<uchar> bytes) {
        lfloat lfBuf_pre;
        lfloat lfBuf_cur;

        lfBuf_pre.value = pre;
        lfBuf_cur = lfBuf_pre;

        for (int i = 0; i < bytes.size(); i++) {
            lfBuf_cur.byte[i] = bytes[i];
        }
        return lfBuf_cur.value;
    }

};
#endif //SZ3_BYTEUTIL_HPP
