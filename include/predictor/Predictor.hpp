#ifndef _PREDICTOR_HPP
#define _PREDICTOR_HPP

namespace SZ{

// N-d lorenzo predictor
template <class T, int N>
class LorenzoPredictor{
protected:
	int id;	
public:
	Predictor();
	~Predictor();
	void predict() = 0;
};

}
#endif
