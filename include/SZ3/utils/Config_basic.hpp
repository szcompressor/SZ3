#ifndef SZ_Config_Basic_HPP
#define SZ_Config_Basic_HPP

namespace SZ3{

    //enum ALGO { ALGO_LORENZO_REG, ALGO_INTERP_LORENZO, ALGO_INTERP, ALGO_NOPRED, ALGO_LOSSLESS };


	namespace concepts {

	class BasicConfig{
	 	public:
	 		virtual void loadcfg(const std::string &cfgpath) = 0;
	 		virtual static size_t size_est() = 0;
	 		virtual void save(unsigned char *&c) = 0;
	 		virtual void load(const unsigned char *&c) = 0;
	 		bool openmp = false;
	 		uint8_t cmprAlgo = 1;//ALGO_INTERP_LORENZO
	 		size_t num = 0;
	 		int QoZ = -1;
	}
}

#endif // SZ_Config_Basic_HPP
