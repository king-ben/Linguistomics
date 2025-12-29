#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <atomic>
#include <condition_variable>
#include <mutex> 
#include <string>
#include <thread>
#include <vector>
class Alignment;



class ThreadTask {

    public:
                                ThreadTask(void);
        virtual                ~ThreadTask(void);
        virtual void            run(void) = 0;
};

class ReadAlignmentTask : public ThreadTask {

    public: 
                                                ReadAlignmentTask(void) = delete;
                                                ReadAlignmentTask(std::string fn, size_t os, std::vector<std::vector<Alignment*>>* a, std::vector<std::string>* an);
        virtual                                ~ReadAlignmentTask(void);
        
        void                                    run(void) override {
                                                    readAlignment();
                                                }
        
        
    protected:
        void                                    readAlignment(void);
        
    private:
        static std::mutex                       coutMutex;  // shared across all instances
        size_t                                  offset;
        std::string                             fileName;
        std::vector<std::vector<Alignment*>>*   alignments;
        std::vector<std::string>*               alignmentNames;
};

class ThreadPool {

    public:
        explicit                ThreadPool(void);
                               ~ThreadPool(void);
        void                    pushTask(ThreadTask* task);
        void                    wait(void);
        int                     getThreadCount(void) const { return threadCount; }

    private:
        static constexpr size_t queueCapacity = 1024;           // cache-friendly circular buffer for tasks
        static constexpr size_t queueMask = queueCapacity - 1;  // note the size must be power of 2 for fast modulo via bitwise AND
        ThreadTask*             taskQueue[queueCapacity];
        size_t                  queueHead;                      // next slot to read from
        size_t                  queueTail;                      // next slot to write to
        size_t                  queueSize;                      // current number of tasks in queue
        int                     threadCount;
        std::atomic<size_t>     tasksInFlight;                  // tasks queued + tasks being executed
        std::atomic<bool>       running;
        std::mutex              mutex;                          // single mutex protects the queue; condition variables share it
        std::condition_variable taskAvailable;                  // signaled when queue non-empty
        std::condition_variable allComplete;                    // signaled when tasksInFlight hits 0
        std::thread*            threads;
        void                    worker(void);
        ThreadTask*             popTask(void);
};


#endif
