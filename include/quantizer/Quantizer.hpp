#ifndef _SZ_QUANTIZER_HPP
#define _SZ_QUANTIZER_HPP

#include "utils/Concepts.hpp"

namespace SZ {
    namespace concepts {

        template<class T>
        class QuantizerInterface {
        public:

            virtual ~QuantizerInterface() = default;

            virtual void precompress_data() const = 0;

            virtual void postcompress_data() const = 0;

            virtual void predecompress_data() const = 0;

            virtual void postdecompress_data() const = 0;

            virtual void save(uchar *&c) const = 0;

            virtual void load(const uchar *&c, size_t &remaining_length) = 0;

            virtual int quantize(T data, T pred) = 0;

            virtual int quantize_and_overwrite(T &data, T pred) = 0;

            virtual T recover(T pred, int quant_index) = 0;

        };

        /**
         * Concept for the predictor class.
         *
         * Matches classes like the following:
         *
         * class my_quantizer {
         *  using iterator = std::multi_dimensional_range<T,N>::iterator
         *
         *  /// returns the prediction for the single element pointed to by the iterator
         *  int quantize(T);
         *
         *  /// pre-processing hook run at the beginning of compression
         *  void precompress_data();
         *
         *  /// post-processing hook run at the end of compression
         *  void postcompress_data();
         *
         *  /// pre-processing hook run before compressing each block
         *  /// post processing can be done either during postcompress_data or on the next call to precompress_block
         *  void precompress_block();
         *
         *  /// pre-processing hook run at the beginning of decompression
         *  void predecompress_data();
         *
         *  /// pre-processing hook run at the end of decompression
         *  void postdecompress_data();
         *
         *  /// pre-processing hook run before decompressing each block
         *  /// post processing can be done either during postcompress_data or on the next call to precompress_block
         *  void predecompress_block();
         *
         *  /// returns a string which represents the configuration of the preditor in a serialized form
         *  std::string save() const;
         *
         *  /// returns a predictor from a serialized form
         *  static my_quantizer load(const unsigned char*&, size_t& len);
         *
         * };
         */

        template<typename T, typename = void>
        struct is_quantizer : false_type {
        };
        template<typename T>
        struct is_quantizer<T, void_t<
                typename T::value_type,
                typename T::reference,
                decltype(std::declval<T>().quantize(
                        std::declval<typename T::value_type>(),
                        std::declval<typename T::value_type>()
                )),
                // decltype(std::declval<T>().quantize(
                //       std::declval<typename T::value_type>(),
                //       std::declval<typename T::value_type>(),
                //       std::declval<typename T::reference>())
                //     ),
                decltype(std::declval<T>().quantize_and_overwrite(
                        std::declval<typename T::reference>(),
                        std::declval<typename T::value_type>()
                )),
                decltype(std::declval<T>().recover(
                        std::declval<typename T::value_type>(),
                        std::declval<int>())
                ),
                decltype(std::declval<T>().save(
                        std::declval<unsigned char *&>())
                ),
                decltype(std::declval<T>().load(
                        std::declval<const unsigned char *&>(), std::declval<size_t &>())
                ),
                // decltype(T::load(std::declval<const unsigned char*&>(), std::declval<size_t&>())),
                decltype(std::declval<T>().precompress_data()),
                decltype(std::declval<T>().postcompress_data()),
                decltype(std::declval<T>().precompress_block()),
                decltype(std::declval<T>().predecompress_data()),
                decltype(std::declval<T>().postdecompress_data()),
                decltype(std::declval<T>().predecompress_block())//,
        >> : true_type {
            // static_assert(
            //     std::is_same<
            //       decltype(std::declval<T>().save()),
            //       std::string
            //     >::value, "save must return a string");
            // static_assert(
            //     std::is_same<
            //     decltype(std::declval<T>().quantize(
            //           std::declval<typename T::value_type>(),
            //           std::declval<typename T::value_type>(),
            //           std::declval<typename T::reference>())),
            //     int
            //     >::value, "quantize must return an int");
            static_assert(
                    std::is_same<
                            decltype(std::declval<T>().quantize_and_overwrite(
                                    std::declval<typename T::reference>(),
                                    std::declval<typename T::value_type>())),
                            int
                    >::value, "quantize must return an int");
            static_assert(
                    std::is_same<
                            decltype(std::declval<T>().quantize(
                                    std::declval<typename T::value_type>(),
                                    std::declval<typename T::value_type>())),
                            int
                    >::value, "quantize must return an int");
            static_assert(
                    std::is_same<
                            decltype(std::declval<T>().recover(
                                    std::declval<typename T::value_type>(),
                                    std::declval<int>())),
                            typename T::value_type
                    >::value, "recover must return an T::value_type");
        };

    }
}

#endif
