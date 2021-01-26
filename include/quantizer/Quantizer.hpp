#ifndef _SZ_QUANTIZER_HPP
#define _SZ_QUANTIZER_HPP


namespace SZ {
    namespace concepts {

        template<class T>
        class QuantizerInterface {
        public:

            virtual ~QuantizerInterface() = default;

            virtual void precompress_data() = 0;

            virtual void postcompress_data() = 0;

            virtual void predecompress_data() = 0;

            virtual void postdecompress_data() = 0;

            virtual void save(uchar *&c) const = 0;

            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            virtual int quantize(T data, T pred) = 0;

            virtual int quantize_and_overwrite(T &data, T pred) = 0;

            virtual T recover(T pred, int quant_index) = 0;

        };
    }
}

#endif
