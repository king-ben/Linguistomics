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
                                                ReadAlignmentTask(std::string fn, size_t os, 
                                                                  std::vector<std::vector<Alignment*>>* a, 
                                                                  std::vector<std::string>* an,
                                                                  std::atomic<int>* counter);
        virtual                                ~ReadAlignmentTask(void);
        
        void                                    run(void) override { readAlignment(); }
        
    protected:
        void                                    readAlignment(void);
        
    private:
        size_t                                  offset;
        std::string                             fileName;
        std::vector<std::vector<Alignment*>>*   alignments;
        std::vector<std::string>*               alignmentNames;
        std::atomic<int>*                       completedCount;
};

class ThreadPool {

    public:
        explicit                ThreadPool(void);
                               ~ThreadPool(void);
        void                    pushTask(ThreadTask* task);
        void                    wait(void);
        int                     getThreadCount(void) const { return threadCount; }

    private:
        static constexpr size_t queueCapacity = 1024;
        static constexpr size_t queueMask = queueCapacity - 1;
        ThreadTask*             taskQueue[queueCapacity];
        size_t                  queueHead;
        size_t                  queueTail;
        size_t                  queueSize;
        int                     threadCount;
        std::atomic<size_t>     tasksInFlight;
        std::atomic<bool>       running;
        std::mutex              mutex;
        std::condition_variable taskAvailable;
        std::condition_variable allComplete;
        std::thread*            threads;
        void                    worker(void);
        ThreadTask*             popTask(void);
};


#endif
