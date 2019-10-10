#ifndef _SZ_ENCODER_HPP
#define _SZ_ENCODER_HPP

namespace SZ{

// an encode that transform T1 data array to T2 data array
template <typename T1, typename T2>
class Encoder{
protected:
	int id;
public:
	Encoder();
	~Encoder();
	virtual void encode(const T1 * data, T2 * output) = 0;
	// virtual void serialize() = 0;
	// virtual void deserialize() = 0;
};

}
#endif
