#ifndef SZ3_QOI_INFO
#define SZ3_QOI_INFO

#include "QoI.hpp"
#include "XSquare.hpp"
#include "LogX.hpp"
#include "RegionalAverage.hpp"
#include "Isoline.hpp"
#include <vector>

namespace SZ {

    template<class T, SZ::uint N>
    std::shared_ptr<concepts::QoIInterface<T, N>> GetQOI(const Config &conf){
        switch(conf.qoi){
            case 1:
                return std::make_shared<SZ::QoI_X_Square<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 2:
                return std::make_shared<SZ::QoI_Log_X<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 3:
                return std::make_shared<SZ::QoI_RegionalAverage<T, N>>(conf.qoiEB, conf.absErrorBound);
            case 4:
            	std::vector<T> values;
            	for(int i=0; i<conf.isovalues.size(); i++){
            		values.push_back(conf.isovalues[i]);
            	}
                return std::make_shared<SZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound);
        }
        return NULL;
    }

}
#endif