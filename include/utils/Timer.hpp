//
// Created by Kai Zhao on 10/30/20.
//

#ifndef SZ_TIMER_HPP
#define SZ_TIMER_HPP

#include <ctime>
#include <string>
#include <iostream>

namespace SZ {
    class Timer {
    public:
        Timer() = default;

        void start() {
            clock_gettime(CLOCK_REALTIME, &begin);
        }

        void stop(const std::string &msg) {
            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << msg << " time = "
                      << (double) (end.tv_sec - begin.tv_sec) + (double) (end.tv_nsec - begin.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;
            fflush(stdout);
        }

    private:
        struct timespec begin, end;
    };
};


#endif //SZ_TIMER_HPP