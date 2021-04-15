#ifndef SZ3_FRONTEND_INTERFACE
#define SZ3_FRONTEND_INTERFACE

#include "def.hpp"

namespace SZ {


    namespace concepts {

        template<class T, uint N>
        class FrontendInterface {
        public:

            virtual ~FrontendInterface() = default;

            virtual std::vector<int> compress(T *data) = 0;

            virtual T *decompress(std::vector<int> &quant_inds) = 0;

            virtual void save(uchar *&c) = 0;

            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            virtual int get_radius() const = 0;

            virtual size_t get_num_elements() const = 0;

            virtual void print() = 0;

            virtual void clear() = 0;
        };

    }

}

#endif
