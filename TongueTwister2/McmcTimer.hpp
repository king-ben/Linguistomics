#ifndef McmcTimer_hpp
#define McmcTimer_hpp

#include <chrono>
#include <cstddef>
typedef std::chrono::steady_clock::time_point TimePoint;



class McmcTimer {

    public:
                            McmcTimer(void);
                           ~McmcTimer(void) = default;
        void                end(void);
        void                start(void);
        void                printStatus(size_t currentCycle, size_t numCycles) const;
        
    private:
        TimePoint           startTime;
        TimePoint           endTime;
        bool                running;
};

#endif
