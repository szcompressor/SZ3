#ifndef SZ3_FRONTEND_INTERFACE
#define SZ3_FRONTEND_INTERFACE

#include "SZ3/def.hpp"
#include <vector>

namespace SZ3::concepts {


    /**
     * Frontend is the combination of predictor and quantizer
     * It handles original data  <-->  quantized error (usually in integer)
     *
     * This will be changed to DecompositionTypeII
     * DecompositionTypeI:
     *                 E.g: transformation, predictors with no quantizer fused
     *                 Example: regression, zfp, mgard, wavelet (cdf97 in sperr)
     * DecompositionTypeII:
     *                 predictors which internally calls quantizer APIs.
     *                 example: Lorenzo, interpolation
     *                 default quantizer need to be provided to reduce the initialization complexity for users
     * @tparam T original data
     * @tparam N data dimension
     */
    template<class T, uint N>
    class FrontendInterface {
    public:

        virtual ~FrontendInterface() = default;

        /**
         * predict the data and quantize the error
         * @param data original input
         * @return quantized prediction error
         */
        virtual std::vector<int> compress(T *data) = 0;

        /**
         * reverse of compress(), reconstruct the data
         * @param quant_inds quantized prediction error
         * @param dec_data place to write the reconstructed data
         * @return same value with dec_data
         */
        virtual T *decompress(std::vector<int> &quant_inds, T *dec_data) = 0;

        /**
         * serialize the frontend and store it to a buffer
         * @param c One large buffer is pre-allocated, and the start location of the serialized frontend in the buffer is indicated by c.
         *          After saving the frontend to the buffer, this function should change c to indicate the next empty location in the buffer
         */
        virtual void save(uchar *&c) = 0;

        /**
         * deserialize the frontend from a buffer
         * @param c start location of the frontend in the buffer
         * @param remaining_length the remaining length of the buffer
         */
        virtual void load(const uchar *&c, size_t &remaining_length) = 0;

        virtual size_t size_est() { return 0; };

        virtual int get_radius() { return 0; };

        virtual void print() {};

        virtual void clear() {};
    };


}

#endif
