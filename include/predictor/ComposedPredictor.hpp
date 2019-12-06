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
    virtual void print() const=0;
  };

  template <class T, uint N, class Base>
  class RealPredictor : public VirtualPredictor<T, N>, public Base{
  public:
    using Range = multi_dimensional_range<T, N>;
    using iterator = typename multi_dimensional_range<T, N>::iterator;
    RealPredictor(std::shared_ptr<Base> p){
    	base = p;
    }
    void precompress_data(const iterator& iter) {
      base->precompress_data(iter);
    }
    void postcompress_data(const iterator& iter) {
      base->postcompress_data(iter);
    }
    void predecompress_data(const iterator& iter) {
      base->predecompress_data(iter);
    }
    void postdecompress_data(const iterator& iter) {
      base->postdecompress_data(iter);
    }
    void precompress_block(const std::shared_ptr<Range>& range){
      base->precompress_block(range);
    }
    void predecompress_block(const std::shared_ptr<Range>& range){
      base->predecompress_block(range);
    }
    void save(uchar*& c) const{
      base->save(c);
    }
    void load(const uchar*& c, size_t& remaining_length){
      base->load(c, remaining_length);
    }
    inline T predict(const iterator& iter) const noexcept{
      return base->predict(iter);
    }
	inline T estimate_error(const iterator& iter) const noexcept {
	  return base->estimate_error(iter);
	}
    void print() const{
    	base->print();
    }
  private:
  	std::shared_ptr<Base> base;
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
    {
      auto sample_begin =sample_range->begin();
      auto sample_end = sample_range->end();
      for(auto iter = sample_begin; iter != sample_end; iter ++){
        for(int i=0; i<predictors.size(); i++){
          err[i] += predictors[i]->estimate_error(iter);
        }
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
    void print() const{
		for(const auto& p:predictors){
			p->print();
		}    	
    }
	template <typename P1>
	void instantiate(P1 p1) {
		predictors.push_back(std::move(p1));
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
		unpack(Ps...);
	}
private:
	std::vector<std::shared_ptr<VirtualPredictor<T, N>>> predictors;
	std::vector<int> selection;
	int sid;							// selected index
	size_t current_index = 0;			// for decompression only
	uchar buf[512];
};

}


#endif
