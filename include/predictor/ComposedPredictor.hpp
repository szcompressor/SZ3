#ifndef _SZ_COMPOSED_PREDICTOR_HPP
#define _SZ_COMPOSED_PREDICTOR_HPP

#include <cassert>
#include "def.hpp"
#include "utils/Iterator.hpp"
#include "utils/Compat.hpp"
#include "predictor/Predictor.hpp"
#include <iostream>
#include <memory>

namespace SZ{
  template <class T, uint N>
  class VirtualPredictor{
  public:
    using Range = multi_dimensional_range<T, N>;
    using iterator = typename multi_dimensional_range<T, N>::iterator;
    virtual ~VirtualPredictor()=default;
    virtual void precompress_data(const iterator&) = 0;
    virtual void postcompress_data(const iterator&) = 0;
    virtual void predecompress_data(const iterator&) = 0;
    virtual void postdecompress_data(const iterator&) = 0;
    virtual void precompress_block(const std::shared_ptr<Range>&) = 0;
    virtual void predecompress_block(const std::shared_ptr<Range>&) = 0;
    virtual void save(uchar*& c) const = 0;
    virtual void load(const uchar*& c, size_t& remaining_length) = 0;
    virtual T predict(const iterator& iter) const noexcept = 0;
    virtual T estimate_error(const iterator& iter) const noexcept = 0;
    virtual void print()=0;
  };

  template <class T, uint N, class Base>
  class RealPredictor : public VirtualPredictor<T, N>, public Base{
  public:
    using Range = multi_dimensional_range<T, N>;
    using iterator = typename multi_dimensional_range<T, N>::iterator;
    template <class... Params>
    RealPredictor(Params&&... p){
    	Base(p...);
    }
    void precompress_data(const iterator& iter) {
      Base::precompress_data(iter);
    }
    void postcompress_data(const iterator& iter) {
      Base::postcompress_data(iter);
    }
    void predecompress_data(const iterator& iter) {
      Base::predecompress_data(iter);
    }
    void postdecompress_data(const iterator& iter) {
      Base::postdecompress_data(iter);
    }
    void precompress_block(const std::shared_ptr<Range>& range){
      Base::precompress_block(range);
    }
    void predecompress_block(const std::shared_ptr<Range>& range){
      Base::predecompress_block(range);
    }
    void save(uchar*& c) const{
      Base::save(c);
    }
    void load(const uchar*& c, size_t& remaining_length){
      Base::load(c, remaining_length);
    }
    inline T predict(const iterator& iter) const noexcept{
      return Base::predict(iter);
    }
	inline T estimate_error(const iterator& iter) const noexcept {
	  return Base::estimate_error(iter);
	}
    void print(){
    	Base::print();
    }
  };

template <class T, uint N>
class ComposedPredictor{
public:
    using Range = multi_dimensional_range<T, N>;
    using iterator = typename multi_dimensional_range<T, N>::iterator;
	void precompress_data(const iterator& iter) const noexcept{
		for(const auto& p:predictors){
			p->precompress_data(iter);
		}
	}
	void postcompress_data(const iterator& iter) const noexcept{
		for(const auto& p:predictors){
			p->postcompress_data(iter);
		}		
	}
	void predecompress_data(const iterator& iter) const noexcept{
		for(const auto& p:predictors){
			p->predecompress_data(iter);
		}				
	}
	void postdecompress_data(const iterator& iter) const noexcept{
		for(const auto& p:predictors){
			p->postdecompress_data(iter);
		}				
	}
	void precompress_block(const std::shared_ptr<Range>& range){
		for(const auto& p:predictors){
			p->precompress_block(range);
		}
		auto sample_range = range;
		// change sample range
		std::vector<double> err(predictors.size(), 0);
		for(auto iter = sample_range->begin(); iter != sample_range->end(); iter ++){
			for(int i=0; i<predictors.size(); i++){
				err[i] += predictors[i]->estimate_error(iter);
			}
		}
		sid = std::distance(err.begin(), std::min_element(err.begin(), err.end()));
		selection.push_back(sid);
		// std::cout << sid << std::endl;
	}
	void predecompress_block(const std::shared_ptr<Range>& range){
		for(const auto& p:predictors){
			p->predecompress_block(range);
		}
		sid = selection[current_index ++];
	}
	void save(uchar*& c) const{
		for(const auto& p:predictors){
			p->save(c);
		}
		// store selection
		*reinterpret_cast<size_t*>(c) = selection.size();
	    c += sizeof(size_t);
		memcpy(c, selection.data(), selection.size()*sizeof(int));
		c += selection.size()*sizeof(int);
		std::cout << "selection size: " << selection.size() << std::endl;
	}
	void load(const uchar*& c, size_t& remaining_length){
		for(const auto& p:predictors){
			p->load(c, remaining_length);
		}
		// load selection
		size_t selection_size = *reinterpret_cast<const size_t*>(c);
	    c += sizeof(size_t);
	    this->selection = std::vector<int>(reinterpret_cast<const int*>(c), reinterpret_cast<const int*>(c) + selection_size);
	    c += selection_size * sizeof(T);		
	}
    inline T predict(const iterator& iter) const noexcept{
      return predictors[sid]->predict(iter);
    }
	template <typename P1>
	void instantiate(P1 p1) {
    std::unique_ptr<VirtualPredictor<T,N>> p = compat::make_unique<RealPredictor<T, N, P1>>();
		uchar* buf_pos = buf;
		p1.save(buf_pos);
		const uchar* buf_pos2 = buf;
		size_t len = 512;
		p->load(buf_pos2, len);
		predictors.emplace_back(std::move(p));
	}
	template <typename P1>
	void unpack(P1 p1) {
		instantiate(p1);
	}
	template <typename P1, typename... Rest>
	void unpack(P1 p1, Rest... Rs) {
		instantiate<P1>(p1);
		unpack(Rs...);
	}
	template <class... Predictors>
	ComposedPredictor(Predictors&&... Ps){
		// static_assert(sizeof...(Predictors) == 0, "Number of predictors must be larger than 1");
		unpack(Ps...);
		for(const auto& p:predictors){
			p->print();
		}
	}
private:
	std::vector<std::unique_ptr<VirtualPredictor<T, N>>> predictors;
	std::vector<int> selection;
	int sid;							// selected index
	size_t current_index = 0;			// for decompression only
	uchar buf[512];
};

}


#endif
