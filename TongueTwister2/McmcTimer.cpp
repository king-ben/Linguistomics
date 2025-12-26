#include <cstdio>
#include "McmcTimer.hpp"



McmcTimer::McmcTimer(void) : running(false) {

}

void McmcTimer::start(void) {

    startTime = std::chrono::steady_clock::now();
    running = true;
}

void McmcTimer::end(void) {

    endTime = std::chrono::steady_clock::now();
    running = false;
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    long hours   = elapsed / 3600;
    long minutes = (elapsed % 3600) / 60;
    long seconds = elapsed % 60;
    
    if (hours == 0)
        std::printf("   * Chain completed in %02ld:%02ld\n\n", minutes, seconds);
    else 
        std::printf("   * Chain completed in %02ld:%02ld:%02ld\n\n", hours, minutes, seconds);
}

void McmcTimer::printStatus(size_t currentCycle, size_t numCycles) const {

    if (running == false || currentCycle == 0)
        return;
    
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
    
    // estimate total time based on progress
    double fractionComplete = static_cast<double>(currentCycle) / static_cast<double>(numCycles);
    double estimatedTotalMs = static_cast<double>(elapsed) / fractionComplete;
    double remainingMs = estimatedTotalMs - static_cast<double>(elapsed);
    
    long remainingSec = static_cast<long>(remainingMs / 1000.0);
    if (remainingSec < 0)
        remainingSec = 0;
    
    long hours   = remainingSec / 3600;
    long minutes = (remainingSec % 3600) / 60;
    long seconds = remainingSec % 60;
    
    if (hours == 0)
        std::printf("%02ldm:%02lds remaining; ", minutes, seconds);
    else
        std::printf("%02ldh:%02ldm:%02lds remaining; ", hours, minutes, seconds);

    long totalSec = static_cast<long>(estimatedTotalMs / 1000.0);
    hours   = totalSec / 3600;
    minutes = (totalSec % 3600) / 60;
    seconds = totalSec % 60;

    if (hours == 0)
        std::printf("%02ldm:%02lds total", minutes, seconds);
    else
        std::printf("%02ldh:%02ldm:%02lds total", hours, minutes, seconds);
}
