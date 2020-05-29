#ifndef _SZ_SZ_ZONE_COMPRESSOR_HPP
#define _SZ_SZ_ZONE_COMPRESSOR_HPP

#include "predictor/Predictor.hpp"
#include "quantizer/IntegerQuantizer.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "encoder/HuffmanEncoder.hpp"
#include "utils/Iterator.hpp"
#include "utils/Lossless.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/Config.hpp"
#include "def.hpp"
#include <cstring>

namespace SZ {
// type T
    template<class T, size_t N, class Predictor = concepts::VirtualPredictor<T, N>,
            class Quantizer = LinearQuantizer<T>, class Encoder = HuffmanEncoder<int> >
    class SZ_Zone_Compressor {
    public:
        // static_assert(concepts::is_predictor<Predictor>::value, "must implement the predictor interface");
        static_assert(concepts::is_quantizer<Quantizer>::value, "must implement the quatizer interface");

        SZ_Zone_Compressor(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer, Encoder encoder) :
                predictor(predictor), quantizer(quantizer), encoder(encoder), block_size(conf.block_size), stride(conf.stride),
                global_dimensions(conf.dims), num_elements_total(conf.num) {}

        // compress given the error bound
        std::vector<uchar *> compress(T *data, std::vector<size_t> &compressed_size, int num_zone) {
            std::array<size_t, N> zone_dimension = global_dimensions;
            zone_dimension[0] = ceil(zone_dimension[0] / num_zone);
            size_t num_elements_per_zone = 1;
            for (const auto &d:zone_dimension) {
                num_elements_per_zone *= d;
            }
            std::vector<int> quant_inds(num_elements_total);
            size_t quant_count = 0;
            std::vector<uchar *> compressed_data(num_zone + 1);
            std::vector<uchar *> compressed_data_pos(num_zone + 1);
            for (int z = 0; z < num_zone; z++) {
                compressed_data[z] = new uchar[2 * num_elements_per_zone * sizeof(T)];
                compressed_data_pos[z] = compressed_data[z];
            }
            compressed_data[num_zone] = new uchar[num_zone * num_elements_per_zone * sizeof(T)];
            compressed_data_pos[num_zone] = compressed_data[num_zone];

            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);

            for (int z = 0; z < num_zone; z++) {

                if (z == num_zone - 1) {
                    zone_dimension[0] = global_dimensions[0] - zone_dimension[0] * (num_zone - 1);
                }

                auto zone_data = data + z * num_elements_per_zone;
                auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(zone_data,
                                                                                             std::begin(zone_dimension),
                                                                                             std::end(zone_dimension), stride,
                                                                                             0);

                auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(zone_data,
                                                                                             std::begin(zone_dimension),
                                                                                             std::end(zone_dimension), 1, 0);
                std::array<size_t, N> intra_block_dims;
                predictor.precompress_data(inter_block_range->begin());
                quantizer.precompress_data();

                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    // std::cout << *block << " " << lp.predict(block) << std::endl;
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_current_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1 &&
                                               zone_dimension[i] - cur_index * stride < block_size) ?
                                              zone_dimension[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_current_index_vector());
                    T dec_data = 0;
                    predictor.precompress_block(intra_block_range);
                    predictor.precompress_block_commit();
                    quantizer.precompress_block();
//          	reg_count += predictor.get_sid();
                    {
                        auto intra_begin = intra_block_range->begin();
                        auto intra_end = intra_block_range->end();
                        for (auto element = intra_begin; element != intra_end; ++element) {
//                            quantizer.quantize_and_overwrite(*element, predictor.predict(element));
                            quant_inds[quant_count++] = quantizer.quantize_and_overwrite(*element,
                                                                                         predictor.predict(element));
                        }
                    }
                }
                predictor.postcompress_data(inter_block_range->begin());
                quantizer.postcompress_data();

                write(zone_dimension.data(), N, compressed_data_pos[z]);
                write(block_size, compressed_data_pos[z]);
                predictor.save(compressed_data_pos[z]);
                quantizer.save(compressed_data_pos[z]);

                predictor.clear();
                quantizer.clear();

            }

            clock_gettime(CLOCK_REALTIME, &end);
//            std::cout << "Predition & Quantization time = "
//                      << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
//                      << "s" << std::endl;

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos[num_zone]);
            for (int z = 0; z < num_zone; z++) {
                auto elements = (z == num_zone - 1) ? num_elements_total - num_elements_per_zone * (num_zone - 1)
                                                    : num_elements_per_zone;
                encoder.encode(&quant_inds[num_elements_per_zone * z], elements, compressed_data_pos[z]);
            }
            encoder.postprocess_encode();

            std::vector<uchar *> lossless_data(num_zone + 1);
            for (int z = 0; z < num_zone + 1; z++) {
                lossless_data[z] = lossless_compress(compressed_data[z],
                                                     compressed_data_pos[z] - compressed_data[z],
                                                     compressed_size[z]);
                delete[] compressed_data[z];
            }

            return lossless_data;
        }

        void decompress_encoder(uchar const *lossless_compressed_data, const size_t length) {
            auto compressed_data = lossless_decompress(lossless_compressed_data, length);
            uchar const *compressed_data_pos = compressed_data;
            size_t remaining_length = length;
            encoder.load(compressed_data_pos, remaining_length);
            delete[] compressed_data;
        }

        T *decompress_zone(uchar const *lossless_compressed_data, const size_t length, std::vector<T> &dec_data) {
            //exec decompress_encoder() first
            assert(encoder.isLoaded());

            auto compressed_data = lossless_decompress(lossless_compressed_data, length);
            uchar const *compressed_data_pos = compressed_data;
            size_t remaining_length = length;
            std::array<size_t, N> zone_dimensions;
            read(zone_dimensions.data(), N, compressed_data_pos, remaining_length);
            auto num_elements = 1;
            for (const auto &d : zone_dimensions) {
                num_elements *= d;
//                std::cout << d << " ";
            }
//            std::cout << std::endl;
            uint block_size = 0;
            read(block_size, compressed_data_pos, remaining_length);
            predictor.load(compressed_data_pos, remaining_length);
            // std::cout << "load predictor done\n";fflush(stdout);
            quantizer.load(compressed_data_pos, remaining_length);
            // std::cout << "load quantizer done\n";fflush(stdout);
//            encoder.load(compressed_data_pos, remaining_length);
            // size_t unpred_data_size = 0;
            // read(unpred_data_size, compressed_data_pos, remaining_length);
            // T const * unpred_data_pos = (T const *) compressed_data_pos;
            // compressed_data_pos += unpred_data_size * sizeof(T);
            auto quant_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();
            delete[] compressed_data;

            // std::cout << "load encoder done\n";fflush(stdout);
//    	std::cout << quant_inds[157684267] << std::endl;
            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::array<size_t, N> intra_block_dims;
//            auto dec_data = std::make_unique<T[]>(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.data(),
                                                                                         std::begin(zone_dimensions),
                                                                                         std::end(zone_dimensions), block_size,
                                                                                         0);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.data(),
                                                                                         std::begin(zone_dimensions),
                                                                                         std::end(zone_dimensions), 1, 0);

            predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();

//            std::cout << "start decompression" << std::endl;
            size_t reg_count = 0;
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_current_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1) ? zone_dimensions[i] - cur_index * block_size
                                                                      : block_size;
                    }
                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_current_index_vector());

                    predictor.predecompress_block(intra_block_range);
                    quantizer.predecompress_block();
//	          reg_count += predictor.get_sid();
                    // std::cout << "dimensions: " << intra_block_range->get_dimensions(0) << " " << intra_block_range->get_dimensions(1) << " " << intra_block_range->get_dimensions(2) << std::endl;
                    // std::cout << "index: " << block.get_current_index(0) << " " << block.get_current_index(1) << " " << block.get_current_index(2) << std::endl;
                    {
                        auto intra_begin = intra_block_range->begin();
                        auto intra_end = intra_block_range->end();
                        for (auto element = intra_begin; element != intra_end; ++element) {
                            *element = quantizer.recover(predictor.predict(element), *(quant_inds_pos++));
                        }
                    }
                }
            }
//            std::cout << "reg_count = " << reg_count << std::endl;
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
//            return dec_data.release();
            return 0;
        }


    private:
        Predictor predictor;
        Quantizer quantizer;
        Encoder encoder;
        uint block_size;
        uint stride;
        size_t num_elements_total;
        std::array<size_t, N> global_dimensions;
    };


    template<class T, uint N, class Predictor, class Quantizer, class Encoder>
    SZ_Zone_Compressor<T, N, Predictor, Quantizer, Encoder>
    make_sz_zone_compressor(const Config<T, N> &conf,
                            Predictor predictor,
                            Quantizer quantizer,
                            Encoder encoder
    ) {
        return SZ_Zone_Compressor<T, N, Predictor, Quantizer, Encoder>(conf, predictor, quantizer, encoder);
    }
}
#endif
