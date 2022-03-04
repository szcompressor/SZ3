#ifndef SZ3_QOI_INTERFACE
#define SZ3_QOI_INTERFACE
/**
 * Interface for some specific quantities of interest (QoIs)
 * Created by Xin Liang on 12/06/2021.
 */

namespace SZ {


    namespace concepts {

        template<class T, uint N>
        class QoIInterface {
        public:

            virtual ~QoIInterface() = default;

            virtual T interpret_eb(T data) const = 0;

            virtual void update_tolerance(T data, T dec_data) = 0;

            virtual void print() = 0;

        };

    }

}

#endif
