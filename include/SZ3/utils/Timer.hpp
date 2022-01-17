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

        Timer(bool initstart) {
            if (initstart) {
                start();
            }
        }

        void start() {
            clock_gettime(CLOCK_REALTIME, &begin);
        }

        double stop() {
            clock_gettime(CLOCK_REALTIME, &end);
            return (double) (end.tv_sec - begin.tv_sec) + (double) (end.tv_nsec - begin.tv_nsec) / (double) 1000000000;
        }

        double stop(const std::string &msg) {
            double seconds = stop();
#if SZ3_DEBUG_TIMINGS
            std::cout << msg << " time = " << seconds << "s" << std::endl;
            fflush(stdout);
#endif
            return seconds;
        }

    private:
        struct timespec begin, end;
    };
};


#endif //SZ_TIMER_HPP
